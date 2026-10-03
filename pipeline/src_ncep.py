"""Hi-res NCEP guidance (HRRR time-lagged, HREF members, NAM 3 km) and NBM probabilities,
turned into per-site records for each member-hour.

Record keys "<var>@<site>":
    w     10 m wind speed (kt)
    p     precip within the radius: 1 km reflectivity (composite if absent) >= precip_dbz  (0/1)
    l     lightning within the radius: lightning field > ltng_threshold                    (0/1)
    c     ceiling below ceiling_ft at the site (no ceiling counts as 0)                     (0/1)
    v     visibility below vis_sm at the site                                               (0/1)
NBM records use probabilities (0-1) for p, l, c and v.
"""
from __future__ import annotations

import math
import re
import time
from concurrent.futures import ThreadPoolExecutor

from .common import Context, SourceResult, floor_hour, iso, log
from .grib import Missing, earth_relative, exists, fetch, fetch_select

MS_TO_KT = 1.943844
SM_M = 1609.344
FT_M = 0.3048


def _fmt(tmpl, cycle, fh):
    g = time.gmtime(cycle)
    return tmpl.format(ymd=time.strftime("%Y%m%d", g), hh=f"{g.tm_hour:02d}", fh=fh)


def _raw(url, cycle, ctx):
    cached = ctx.cache.get(url)
    if cached is not None:
        return cached
    try:
        vals, meta = fetch(url, ctx.sites, float(ctx.cons["radius_nm"]))
    except Missing:
        return None
    out = {}
    for site in ctx.sites:
        sid, o = site["id"], {}
        u, v = (vals.get("u10") or {}).get(sid), (vals.get("v10") or {}).get(sid)
        if u is not None and v is not None:
            uv = earth_relative(u, v, meta.get("u10"), site["lon"])
            if uv:
                o["w"] = round(math.hypot(*uv) * MS_TO_KT, 2)
        for k in ("refd", "refc", "ltng", "vis"):
            if k in vals and vals[k].get(sid) is not None:
                o[k] = vals[k][sid]
        if "ceil" in vals:                    # field present: undefined means no ceiling
            o["ceil"] = vals["ceil"].get(sid)
        out[sid] = o
    if any(out.values()):
        ctx.cache.put(url, cycle, out)
    return out


def record(raw, ctx):
    if not raw:
        return None
    c, rec = ctx.cons, {}
    for site in ctx.sites:
        sid, o = site["id"], raw.get(site["id"]) or {}
        if "w" in o:
            rec[f"w@{sid}"] = o["w"]
        refl = o.get("refd", o.get("refc"))
        if refl is not None:
            rec[f"p@{sid}"] = int(refl >= c["precip_dbz"])
        if o.get("ltng") is not None:
            rec[f"l@{sid}"] = int(o["ltng"] > c["ltng_threshold"])
        if "ceil" in o:
            h = o["ceil"]
            rec[f"c@{sid}"] = int(h is not None and 0 <= h < c["ceiling_ft"] * FT_M)
        if o.get("vis") is not None:
            rec[f"v@{sid}"] = int(o["vis"] < c["vis_sm"] * SM_M)
    return rec or None


def _pick(bases, tmpl, cycles):
    for base in bases:
        for c in cycles:
            if exists(f"{base}/{_fmt(tmpl, c, 1)}.idx"):
                return base, c
    return None


def _run(scfg, ctx, tasks, workers):
    def run(task):
        mid, url, c, valid = task
        try:
            return task, record(_raw(url, c, ctx), ctx)
        except Exception as e:
            log.warning("%s %s: %s", scfg["id"], url.rsplit("/", 1)[-1], e)
            return task, None
    members = {}
    with ThreadPoolExecutor(workers) as ex:
        for (mid, _, _, valid), rec in ex.map(run, tasks):
            if rec:
                members.setdefault(mid, {})[valid] = rec
    return members


def tle(scfg, ctx):
    recent = [floor_hour(ctx.now) - k * 3600 for k in range(10)]
    picked = _pick([scfg["base"]], scfg["file"], recent)
    if not picked:
        return SourceResult({}, status="missing", note="no recent cycle found")
    base, latest = picked
    cycles = [latest - k * 3600 for k in range(int(scfg.get("lag_cycles", 6)))]
    cycles += [c for c in (latest - k * 3600 for k in range(30))
               if time.gmtime(c).tm_hour % 6 == 0 and c not in cycles][: int(scfg.get("synoptic_extra", 2))]
    tasks = []
    for c in cycles:
        mx = scfg.get("long_fh", 48) if time.gmtime(c).tm_hour % 6 == 0 else scfg.get("short_fh", 18)
        for fh in range(0, mx + 1):
            if ctx.in_window(c + fh * 3600):
                tasks.append((time.strftime("%d/%HZ", time.gmtime(c)), f"{base}/{_fmt(scfg['file'], c, fh)}", c, c + fh * 3600))
    mem = _run(scfg, ctx, tasks, 16)
    return SourceResult(mem, cycle=iso(latest), note=f"{len(mem)}/{len(cycles)} cycles",
                        status="ok" if len(mem) == len(cycles) else ("partial" if mem else "missing"))


def multi_model(scfg, ctx):
    tasks, notes, missing = [], [], []
    for comp in scfg["components"]:
        hours = set(comp.get("cycles", [0, 6, 12, 18]))
        cands = [c for c in (floor_hour(ctx.now) - k * 3600 for k in range(60)) if time.gmtime(c).tm_hour in hours]
        picked = _pick([comp["base"]], comp["file"], cands)
        if not picked:
            missing.append(comp["id"])
            continue
        base, latest = picked
        use = [c for c in cands if c <= latest][: int(comp.get("lag", 2))]
        notes.append(f"{comp['id']} " + ", ".join(time.strftime("%HZ", time.gmtime(c)) for c in use))
        for c in use:
            for fh in range(0, int(comp.get("max_fh", 48)) + 1):
                if ctx.in_window(c + fh * 3600):
                    tasks.append((f"{comp['id']} {time.strftime('%d/%HZ', time.gmtime(c))}",
                                  f"{base}/{_fmt(comp['file'], c, fh)}", c, c + fh * 3600))
    mem = _run(scfg, ctx, tasks, int(scfg.get("workers", 4)))
    note = "; ".join(notes + ([f"missing: {', '.join(missing)}"] if missing else []))
    return SourceResult(mem, note=note, status="missing" if not mem else ("partial" if missing else "ok"))


# ---------------------------------------------------------------- NBM probabilities
_PERIOD = re.compile(r":(TSTM|APCP):surface:(\d+)-(\d+) hour")
_THR = re.compile(r"prob\s*([<>])\s*([\d.]+)", re.I)
_INST = re.compile(r":(CEIL|VIS):[^:]*:(\d+) hour fcst:", re.I)


def nbm_prob(scfg, ctx):
    """One member of probabilities (0-1). TSTM (1/3/6 h) and 1 h APCP > 0.254 mm: highest
    value within the radius. CEIL and VIS: at the site, using the NBM threshold closest
    to the limit. Period values are applied to every hour of the period."""
    c = ctx.cons
    cands = [x for x in (floor_hour(ctx.now) - k * 3600 for k in range(36))
             if time.gmtime(x).tm_hour in set(scfg.get("cycles", [0, 6, 12, 18]))]
    picked = _pick([scfg["base"]], scfg["file"], cands)
    if not picked:
        return SourceResult({}, status="missing", note="no recent cycle found")
    base, cycle = picked
    maxdur = int(scfg.get("max_period_h", 6))
    want = {"CEIL": c["ceiling_ft"] * FT_M, "VIS": c["vis_sm"] * SM_M}
    used = {}

    def select(inv):
        out, inst = [], {"CEIL": [], "VIS": []}
        for rec in inv:
            d = rec[3]
            th = _THR.search(d)
            if not th:
                continue
            m = _PERIOD.search(d)
            if m:
                var, a, b = m.group(1), int(m.group(2)), int(m.group(3))
                if var == "TSTM" and 0 < b - a <= maxdur and th.group(1) == ">":
                    out.append((f"l{a}-{b}", rec, "area"))
                elif var == "APCP" and b - a == 1 and th.group(1) == ">" and abs(float(th.group(2)) - 0.254) < 0.01:
                    out.append((f"p{a}-{b}", rec, "area"))
                continue
            m = _INST.search(d)
            if m and th.group(1) == "<":
                inst[m.group(1).upper()].append((float(th.group(2)), rec))
        for var, opts in inst.items():
            if opts:
                thr, rec = min(opts, key=lambda o: abs(o[0] - want[var]))
                used[var] = thr
                out.append((("c" if var == "CEIL" else "v") + f"<{thr:g}", rec, "point"))
        return out

    def run(fh):
        url = f"{base}/{_fmt(scfg['file'], cycle, fh)}"
        got = ctx.cache.get(url + "#nbm")
        if got is None:
            try:
                vals = fetch_select(url, select, ctx.sites, float(c["radius_nm"]))
            except Missing:
                return None
            except Exception as e:
                log.warning("%s f%03d: %s", scfg["id"], fh, e)
                return None
            got = {k: d for k, d in vals.items() if isinstance(d, dict)}
            ctx.cache.put(url + "#nbm", cycle, got)
        return fh, got

    fhs = [fh for fh in range(1, int(scfg.get("max_fh", 192)) + 1) if ctx.in_window(cycle + fh * 3600)]
    series, thr_note = {}, set()
    with ThreadPoolExecutor(16) as ex:
        for res in ex.map(run, fhs):
            if not res:
                continue
            fh, got = res
            for name, bysite in got.items():
                if name[0] in "cv" and "<" in name:
                    thr_note.add(f"{'ceiling' if name[0] == 'c' else 'visibility'} < {name[2:]} m")
                    hours = [fh]
                    key = name[0]
                else:
                    a, b = map(int, name[1:].split("-"))
                    hours = range(a + 1, b + 1)
                    key = name[0]
                for h in hours:
                    t = cycle + h * 3600
                    if not ctx.in_window(t):
                        continue
                    rec = series.setdefault(t, {})
                    for sid, v in bysite.items():
                        if v is not None:
                            k = f"{key}@{sid}"
                            rec[k] = max(rec.get(k, 0.0), min(1.0, max(0.0, v / 100.0)))
    series = {t: r for t, r in series.items() if r}
    note = "thunder and 1 h precip within the radius" + (f"; {', '.join(sorted(thr_note))} at the site" if thr_note else "")
    return SourceResult({"m00": series} if series else {}, cycle=iso(cycle), note=note,
                        status="ok" if series else "missing")
