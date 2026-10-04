"""Read what we need from remote GRIB2 files without downloading them whole.

- Records are located from the .idx inventory and fetched with HTTP Range.
  NOMADS accepts several ranges in one request (one hit per file); AWS does
  not, so there adjacent records are merged and fetched range by range.
- Point fields: value at the grid point nearest the site.
- Area fields: maximum over every grid point within each radius of the site
  (grid-point masks are computed once per grid and reused).
- Grid-relative Lambert u/v are rotated to earth-relative.
"""
from __future__ import annotations

import math
import os
import re
import tempfile
import threading
import time

import numpy as np

from .common import SESSION, log

try:
    import eccodes
    eccodes.codes_grib_multi_support_on()
except ImportError:
    eccodes = None

POINT = {
    "u10": r":UGRD:10 m above ground:",
    "v10": r":VGRD:10 m above ground:",
    "ceil": r":HGT:cloud ceiling:",
    "vis": r":VIS:surface:",
    "t2": r":TMP:2 m above ground:",
    "td2": r":DPT:2 m above ground:",
}
AREA = {
    "refd": r":REFD:1000 m above ground:",
    "refc": r":REFC:entire atmosphere",
    "ltng": r":LTNG:entire atmosphere",
}
_SKIP = re.compile(r"(\bmax\b|\bmin\b|\bave\b|\bacc\b|prob|%)", re.I)
NM = 1852.0


class Missing(Exception):
    """File or inventory not published (yet)."""


class _RateLimit:
    """Evenly spaced requests; NOMADS blocks IPs above ~120 hits/minute."""

    def __init__(self, per_minute: float):
        self.gap = 60.0 / per_minute
        self.next = 0.0
        self.lock = threading.Lock()

    def wait(self):
        with self.lock:
            now = time.monotonic()
            t = max(now, self.next)
            self.next = t + self.gap
        if t > now:
            time.sleep(t - now)


NOMADS_LIMIT = _RateLimit(90)


def _nomads(url: str) -> bool:
    return "nomads.ncep.noaa.gov" in url


def _get(url: str, **kw):
    if _nomads(url):
        NOMADS_LIMIT.wait()
    return SESSION.get(url, **kw)


def exists(url: str) -> bool:
    try:
        if _nomads(url):
            NOMADS_LIMIT.wait()
        return SESSION.head(url, timeout=15, allow_redirects=True).status_code == 200
    except Exception:
        return False


def read_idx(idx_url: str):
    r = _get(idx_url, timeout=30)
    if r.status_code in (403, 404):
        raise Missing(idx_url)
    r.raise_for_status()
    recs = []
    for line in r.text.splitlines():
        parts = line.split(":")
        if len(parts) < 3:
            continue
        try:
            start = int(parts[1])
        except ValueError:
            continue
        recs.append((parts[0], start, ":" + ":".join(parts[2:])))
    out = []
    for i, (no, start, desc) in enumerate(recs):
        end = next((recs[j][1] - 1 for j in range(i + 1, len(recs)) if recs[j][1] > start), None)
        out.append((no, start, end, desc))
    return out


def match_fields(inv, wanted: dict[str, str]) -> dict[str, tuple]:
    found = {}
    for name, pat in wanted.items():
        rx = re.compile(pat)
        for rec in inv:
            m = rx.search(rec[3])
            if m and not _SKIP.search(rec[3][m.end():]):
                found[name] = rec
                break
    return found


# ---------------------------------------------------------------- decoding
_MASKS: dict = {}
_MASK_LOCK = threading.Lock()


def _grid_key(gid):
    return tuple(eccodes.codes_get(gid, k) for k in (
        "gridType", "Ni", "Nj", "latitudeOfFirstGridPointInDegrees", "longitudeOfFirstGridPointInDegrees"))


def _masks(gid, sites, radius_nm):
    """{site_id: grid indices within radius_nm of the site} (computed once per grid)."""
    key = (_grid_key(gid), tuple((x["id"], x["lat"], x["lon"]) for x in sites), radius_nm)
    with _MASK_LOCK:
        if key in _MASKS:
            return _MASKS[key]
    lats = np.asarray(eccodes.codes_get_array(gid, "latitudes"))
    lons = np.asarray(eccodes.codes_get_array(gid, "longitudes"))
    out = {}
    for x in sites:
        p1, p2 = math.radians(x["lat"]), np.radians(lats)
        dl = np.radians(((lons - x["lon"] + 180) % 360) - 180)
        a = np.sin((p2 - p1) / 2) ** 2 + math.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
        d = 2 * 6371000.0 * np.arcsin(np.sqrt(np.clip(a, 0, 1))) / NM
        out[x["id"]] = np.nonzero(d <= radius_nm)[0]
    with _MASK_LOCK:
        _MASKS[key] = out
    return out


def _nearest(gid, lat, lon):
    for lo in (lon, lon % 360.0):
        try:
            r = eccodes.codes_grib_find_nearest(gid, lat, lo)[0]
            return float(r["value"] if isinstance(r, dict) else getattr(r, "value"))
        except Exception:
            continue
    raise RuntimeError("nearest-point lookup failed")


def _meta(gid) -> dict:
    meta = {"grid": eccodes.codes_get(gid, "gridType")}
    try:
        meta["rel"] = int(eccodes.codes_get(gid, "uvRelativeToGrid"))
    except Exception:
        meta["rel"] = 0
    if meta["grid"] == "lambert":
        meta["lov"] = float(eccodes.codes_get(gid, "LoVInDegrees"))
        meta["latin1"] = float(eccodes.codes_get(gid, "Latin1InDegrees"))
    return meta


def _decode(blob: bytes, jobs: list, sites, radius) -> list:
    """jobs[i] = 'point' | 'area' for the i-th message in blob."""
    out = []
    fd, path = tempfile.mkstemp(suffix=".grib2")
    try:
        with os.fdopen(fd, "wb") as f:
            f.write(blob)
        with open(path, "rb") as f:
            i = 0
            while True:
                gid = eccodes.codes_grib_new_from_file(f)
                if gid is None:
                    break
                try:
                    job = jobs[i] if i < len(jobs) else None
                    miss = eccodes.codes_get(gid, "missingValue")
                    if job == "point":
                        res = {}
                        for x in sites:
                            v = _nearest(gid, x["lat"], x["lon"])
                            res[x["id"]] = None if abs(v - miss) < 1e-6 or abs(v) > 1e10 else v
                        out.append((res, _meta(gid)))
                    elif job == "area":
                        vals = np.asarray(eccodes.codes_get_values(gid))
                        res = {}
                        for sid, idx in _masks(gid, sites, radius).items():
                            sub = vals[idx]
                            sub = sub[(np.abs(sub - miss) > 1e-6) & (np.abs(sub) < 1e10)]
                            res[sid] = float(sub.max()) if sub.size else None
                        out.append((res, {}))
                    else:
                        out.append((None, {}))
                    i += 1
                finally:
                    eccodes.codes_release(gid)
    finally:
        os.unlink(path)
    return out


def _multipart(resp) -> list[tuple[int, bytes]]:
    """Parse a multipart/byteranges body into [(start, bytes)]."""
    ctype = resp.headers.get("Content-Type", "")
    if "multipart/byteranges" not in ctype:
        cr = resp.headers.get("Content-Range", "")
        m = re.search(r"bytes (\d+)-", cr)
        return [(int(m.group(1)) if m else 0, resp.content)]
    boundary = re.search(r"boundary=\"?([^\";]+)\"?", ctype).group(1).encode()
    parts = []
    for chunk in resp.content.split(b"--" + boundary):
        if b"\r\n\r\n" not in chunk:
            continue
        head, body = chunk.split(b"\r\n\r\n", 1)
        m = re.search(rb"Content-Range:\s*bytes (\d+)-(\d+)", head, re.I)
        if not m:
            continue
        a, b = int(m.group(1)), int(m.group(2))
        parts.append((a, body[: b - a + 1]))
    return parts


def fetch(grib_url: str, sites: list, radius_nm: float) -> tuple[dict, dict]:
    """values[field] = {site_id: value}: nearest grid point for POINT fields, maximum within
    radius_nm for AREA fields. A field present in the file but undefined at a point gives None.
    Raises Missing."""
    if eccodes is None:
        raise RuntimeError("eccodes not installed")
    inv = read_idx(grib_url + ".idx")
    kinds = {**{k: "point" for k in POINT}, **{k: "area" for k in AREA}}
    found = sorted(match_fields(inv, {**POINT, **AREA}).items(), key=lambda kv: kv[1][1])
    return _fetch_found(grib_url, found, kinds, sites, radius_nm)


def fetch_select(grib_url: str, select, sites: list, radius_nm: float) -> dict:
    """select(inventory) -> [(name, record, "point"|"area")]."""
    if eccodes is None:
        raise RuntimeError("eccodes not installed")
    inv = read_idx(grib_url + ".idx")
    chosen = select(inv)
    found = sorted(((n, r) for n, r, _ in chosen), key=lambda kv: kv[1][1])
    vals, _ = _fetch_found(grib_url, found, {n: k for n, _, k in chosen}, sites, radius_nm)
    return vals


def _fetch_found(grib_url, found, kinds, sites, radius_nm):
    if not found:
        return {}, {}
    # contiguous runs of whole messages; a sub-message record stands alone
    groups: list[list] = []
    for name, rec in found:
        whole = "." not in rec[0]
        g = groups[-1] if groups else None
        if g and whole and g[-1][2] and g[-1][1][2] is not None and rec[1] == g[-1][1][2] + 1:
            g.append((name, rec, whole))
        else:
            groups.append([(name, rec, whole)])

    def rng(g):
        a, b = g[0][1][1], g[-1][1][2]
        return f"{a}-{b}" if b is not None else f"{a}-"

    blobs: dict[int, bytes] = {}
    if _nomads(grib_url) and len(groups) > 1 and all(g[-1][1][2] is not None for g in groups):
        r = _get(grib_url, headers={"Range": "bytes=" + ",".join(rng(g) for g in groups)}, timeout=90)
        if r.status_code in (403, 404, 416):
            raise Missing(grib_url)
        if r.status_code == 200:
            raise RuntimeError("server ignored byte ranges")
        r.raise_for_status()
        for a, body in _multipart(r):
            blobs[a] = body
    for g in groups:
        if g[0][1][1] in blobs:
            continue
        r = _get(grib_url, headers={"Range": "bytes=" + rng(g)}, timeout=90)
        if r.status_code in (403, 404, 416):
            raise Missing(grib_url)
        r.raise_for_status()
        blobs[g[0][1][1]] = r.content

    vals, meta = {}, {}
    for g in groups:
        blob = blobs.get(g[0][1][1])
        if blob is None:
            continue
        if len(g) == 1 and not g[0][2]:
            k = int(g[0][1][0].split(".")[1]) - 1
            jobs = [None] * k + [kinds[g[0][0]]]
            dec = _decode(blob, jobs, sites, radius_nm)
            pairs = [(g[0][0], dec[k] if k < len(dec) else (None, {}))]
        else:
            dec = _decode(blob, [kinds[n] for n, _, _ in g], sites, radius_nm)
            pairs = [(n, dec[i] if i < len(dec) else (None, {})) for i, (n, _, _) in enumerate(g)]
        for n, (v, m) in pairs:
            vals[n], meta[n] = v, m
    return vals, meta


def earth_relative(u, v, meta, lon):
    if not meta or not meta.get("rel"):
        return u, v
    if meta.get("grid") != "lambert":
        log.warning("grid-relative winds on %s grid not handled; dropping", meta.get("grid"))
        return None
    lov = meta["lov"] if meta["lov"] <= 180 else meta["lov"] - 360
    a = math.sin(math.radians(meta["latin1"])) * math.radians(lon - lov)
    c, s = math.cos(a), math.sin(a)
    return c * u + s * v, -s * u + c * v
