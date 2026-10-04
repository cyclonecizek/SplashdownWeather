"""Open-Meteo sources sampled at each site plus a 6-point ring at the precip radius.

Record keys "<var>@<site>": w = 10 m wind (kt) at the site; p = 1 h precip >= precip_in at
the site or any ring point (0/1); dd = 2 m dew point depression (C) at the site;
cs = low-ceiling stand-in (0/1) from dew point depression and low cloud cover, used only
where no source with direct ceiling guidance covers a window. Global ensembles have no lightning, ceiling or visibility
output, so they don't feed those rows.
"""
from __future__ import annotations

import json
import math
import os
import re
import time

from .common import SESSION, Context, SourceResult, log
from .marine import score as ml_score

PLEVS = (1000, 975, 950, 925, 850)
STD_HEIGHT_M = {1000: 110, 975: 330, 950: 560, 925: 780, 850: 1460}   # used if geopotential is absent

ENS_URL = "https://ensemble-api.open-meteo.com/v1/ensemble"
DET_URL = "https://api.open-meteo.com/v1/forecast"
_KEY = re.compile(r"^(wind_speed_10m|precipitation|temperature_2m|dew_point_2m|cloud_cover_low)(?:_member(\d+))?$")


def _points(ctx):
    r = float(ctx.cons["radius_nm"])
    pts = []
    for s in ctx.sites:
        pts.append((s["id"], s["lat"], s["lon"]))
        for k in range(6):
            b = math.radians(60 * k)
            pts.append((s["id"], round(s["lat"] + r / 60 * math.cos(b), 4),
                        round(s["lon"] + r / 60 * math.sin(b) / math.cos(math.radians(s["lat"])), 4)))
    return pts


def _fetch(url, model, ctx, pts, variables):
    days = max(1, math.ceil((ctx.t_end - ctx.now) / 86400) + 1)
    variables = list(variables)
    for _ in range(len(variables)):
        r = SESSION.get(url, timeout=180, params={
            "latitude": ",".join(str(p[1]) for p in pts), "longitude": ",".join(str(p[2]) for p in pts),
            "models": model, "hourly": ",".join(variables), "wind_speed_unit": "kn", "precipitation_unit": "inch",
            "timeformat": "unixtime", "timezone": "GMT", "past_days": 1, "forecast_days": min(days, 16)})
        bad = [v for v in variables if v in r.text] if r.status_code == 400 else []
        if not bad:
            break
        log.info("%s: dropping unsupported %s", model, bad[0])
        variables.remove(bad[0])
    r.raise_for_status()
    js = r.json()
    return js if isinstance(js, list) else [js]


def _members(locs, pts, ctx, only):
    c = ctx.cons
    out = {}
    by_site = {}
    for (sid, _, _), loc in zip(pts, locs):
        by_site.setdefault(sid, []).append(loc)
    for sid, site_locs in by_site.items():
        times = site_locs[0]["hourly"]["time"]
        per = {}
        for li, loc in enumerate(site_locs):
            for key, arr in loc["hourly"].items():
                m = _KEY.match(key)
                if m:
                    mid = f"m{int(m.group(2)):02d}" if m.group(2) else "m00"
                    per.setdefault(mid, {}).setdefault(m.group(1), [None] * len(site_locs))[li] = arr
        for mid, f in per.items():
            for i, t in enumerate(times):
                if not ctx.in_window(t):
                    continue
                rec = {}
                w = f.get("wind_speed_10m")
                if "w" in only and w and w[0] is not None and w[0][i] is not None:
                    rec[f"w@{sid}"] = round(w[0][i], 2)
                tt, td = f.get("temperature_2m"), f.get("dew_point_2m")
                if "d" in only and tt and td and tt[0] is not None and td[0] is not None \
                        and tt[0][i] is not None and td[0][i] is not None:
                    dd = max(0.0, tt[0][i] - td[0][i])
                    rec[f"dd@{sid}"] = round(dd, 2)
                    # Ceiling stand-in: cloud base of lifted surface air (about lcl_m_per_c metres per
                    # degree of dew point depression) below the ceiling limit, with low cloud cover at
                    # or above standin_low_cloud_pct when the model provides it.
                    if "cs" in only:
                        lc = f.get("cloud_cover_low")
                        low = lc[0][i] if lc and lc[0] is not None else None
                        base_m = dd * float(c.get("lcl_m_per_c", 125))
                        hit = base_m < c["ceiling_ft"] * 0.3048 and (low is None or low >= c.get("standin_low_cloud_pct", 70))
                        rec[f"cs@{sid}"] = int(hit)
                pr = f.get("precipitation")
                if "p" in only and pr:
                    vals = [a[i] for a in pr if a is not None and a[i] is not None]
                    if vals:
                        rec[f"p@{sid}"] = int(max(vals) >= c["precip_in"])
                if rec:
                    out.setdefault(mid, {}).setdefault(int(t), {}).update(rec)
    return out


def _profiles(url, model, ctx):
    """Marine layer records {mid: {t: {mlb@, mls@, ml@}}} from pressure-level temperatures at
    each site (Open-Meteo returns whichever levels the model has)."""
    pts = [(s["id"], s["lat"], s["lon"]) for s in ctx.sites]
    variables = ["temperature_2m", "dew_point_2m", "relative_humidity_1000hPa"]
    for p in PLEVS:
        variables += [f"temperature_{p}hPa", f"geopotential_height_{p}hPa"]
    locs = _fetch(url, model, ctx, pts, variables)
    key = re.compile(r"^(temperature_2m|dew_point_2m|relative_humidity_1000hPa|temperature_(\d+)hPa|geopotential_height_(\d+)hPa)(?:_member(\d+))?$")
    out = {}
    for (sid, _, _), loc in zip(pts, locs):
        h = loc["hourly"]
        per = {}
        for k, arr in h.items():
            m = key.match(k)
            if m:
                mid = f"m{int(m.group(4)):02d}" if m.group(4) else "m00"
                per.setdefault(mid, {})[m.group(1)] = arr
        for mid, f in per.items():
            for i, t in enumerate(h["time"]):
                if not ctx.in_window(t):
                    continue
                at = lambda name: (f.get(name) or [None] * (i + 1))[i]
                t2, td2 = at("temperature_2m"), at("dew_point_2m")
                levels = [(2.0, t2)] if t2 is not None else []
                for p in PLEVS:
                    tp = at(f"temperature_{p}hPa")
                    if tp is None:
                        continue
                    z = at(f"geopotential_height_{p}hPa")
                    z = z if z is not None else STD_HEIGHT_M[p]
                    if z > 10:
                        levels.append((float(z), tp))
                rec = ml_score(levels, ctx.cons, dd2m=(t2 - td2) if None not in (t2, td2) else None,
                               rh_low=at("relative_humidity_1000hPa"))
                if rec:
                    out.setdefault(mid, {}).setdefault(int(t), {}).update({f"{k}@{sid}": v for k, v in rec.items()})
    return out


def _source(scfg, ctx, url):
    only = set(scfg.get("only", ["w", "p", "d", "cs", "ml"]))
    variables = [v for k, v in (("w", "wind_speed_10m"), ("p", "precipitation")) if k in only]
    if "d" in only:
        variables += ["temperature_2m", "dew_point_2m"]
    if "cs" in only:
        variables += ["cloud_cover_low"]
    path = os.path.join(ctx.root, "cache", f"om_{scfg['id']}.json")
    try:
        with open(path) as f:
            old = json.load(f)
        if ctx.now - old["fetched"] < ctx.om_refresh_h * 3600:
            mem = {m: {int(t): r for t, r in s.items() if ctx.in_window(int(t))} for m, s in old["members"].items()}
            return SourceResult(mem, cycle=time.strftime("%H:%MZ", time.gmtime(old["fetched"])) + " (cached)",
                                note=f"{len(mem)} members via Open-Meteo", status="ok" if mem else "missing")
    except (OSError, ValueError, KeyError):
        pass
    pts = _points(ctx) if "p" in only else [(s["id"], s["lat"], s["lon"]) for s in ctx.sites]
    mem = _members(_fetch(url, scfg["model"], ctx, pts, variables), pts, ctx, only)
    if "ml" in only:
        try:
            for mid, series in _profiles(url, scfg["model"], ctx).items():
                for t, rec in series.items():
                    mem.setdefault(mid, {}).setdefault(t, {}).update(rec)
        except Exception as e:
            log.warning("%s: marine layer profiles failed: %s", scfg["id"], e)
    if mem:
        with open(path, "w") as f:
            json.dump({"fetched": ctx.now, "members": mem}, f, separators=(",", ":"))
    return SourceResult(mem, cycle=time.strftime("%H:%MZ", time.gmtime(ctx.now)),
                        note=f"{len(mem)} members via Open-Meteo", status="ok" if mem else "missing")


def openmeteo_ens(scfg, ctx):
    return _source(scfg, ctx, ENS_URL)


def openmeteo_det(scfg, ctx):
    return _source(scfg, ctx, DET_URL)
