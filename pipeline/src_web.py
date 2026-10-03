"""Open-Meteo sources sampled at each site plus a 6-point ring at the precip radius.

Record keys "<var>@<site>": w = 10 m wind (kt) at the site; p = 1 h precip >= precip_in at
the site or any ring point (0/1). Global ensembles have no lightning, ceiling or visibility
output, so they don't feed those rows.
"""
from __future__ import annotations

import json
import math
import os
import re
import time

from .common import SESSION, Context, SourceResult, log

ENS_URL = "https://ensemble-api.open-meteo.com/v1/ensemble"
DET_URL = "https://api.open-meteo.com/v1/forecast"
_KEY = re.compile(r"^(wind_speed_10m|precipitation)(?:_member(\d+))?$")


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
    r = SESSION.get(url, timeout=180, params={
        "latitude": ",".join(str(p[1]) for p in pts), "longitude": ",".join(str(p[2]) for p in pts),
        "models": model, "hourly": ",".join(variables), "wind_speed_unit": "kn", "precipitation_unit": "inch",
        "timeformat": "unixtime", "timezone": "GMT", "past_days": 1, "forecast_days": min(days, 16)})
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
                pr = f.get("precipitation")
                if "p" in only and pr:
                    vals = [a[i] for a in pr if a is not None and a[i] is not None]
                    if vals:
                        rec[f"p@{sid}"] = int(max(vals) >= c["precip_in"])
                if rec:
                    out.setdefault(mid, {}).setdefault(int(t), {}).update(rec)
    return out


def _source(scfg, ctx, url):
    only = set(scfg.get("only", ["w", "p"]))
    variables = [v for k, v in (("w", "wind_speed_10m"), ("p", "precipitation")) if k in only]
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
    if mem:
        with open(path, "w") as f:
            json.dump({"fetched": ctx.now, "members": mem}, f, separators=(",", ":"))
    return SourceResult(mem, cycle=time.strftime("%H:%MZ", time.gmtime(ctx.now)),
                        note=f"{len(mem)} members via Open-Meteo", status="ok" if mem else "missing")


def openmeteo_ens(scfg, ctx):
    return _source(scfg, ctx, ENS_URL)


def openmeteo_det(scfg, ctx):
    return _source(scfg, ctx, DET_URL)
