"""GFS MOS (MAV) and LAMP (LAV) ceiling and visibility at the nearest near-sea-level
station to each site, from the Iowa Environmental Mesonet MOS service.

Also dd@site, the 2 m dew point depression (C), from the MOS temperature and dew point.

MOS text guidance gives a best-category forecast, not threshold probabilities, so each
run counts as one member (the latest `runs` runs make a small time-lagged set):
    c@site = 1 when the ceiling category is below ceiling_ft   (cat 1 <200 ft, 2 200-400 ft)
    v@site = 1 when the visibility category is below vis_sm    (cat 1 <1/2 mi, 2 1/2-<1 mi)
For each site the first station in its candidate list with data is used.
"""
from __future__ import annotations

import math
import time
from datetime import datetime, timezone

from .common import SESSION, Context, SourceResult, log

API = "https://mesonet.agron.iastate.edu/api/1/mos.json"
# Category upper bounds (exclusive) used by MAV and LAMP text guidance
CIG_TOP_FT = {1: 200, 2: 500, 3: 1000, 4: 2000, 5: 3100, 6: 6600, 7: 12100, 8: 1e9}
VIS_TOP_SM = {1: 0.5, 2: 1.0, 3: 2.0, 4: 3.0, 5: 6.0, 6: 6.5, 7: 1e9}


def _ts(s):
    if s is None:
        return None
    if isinstance(s, (int, float)):
        return int(s / 1000) if s > 1e11 else int(s)
    s = str(s).replace("Z", "+00:00").replace(" ", "T")
    dt = datetime.fromisoformat(s)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return int(dt.timestamp())


def _num(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def _get(station, model, runtime):
    r = SESSION.get(API, timeout=30, params={"station": station, "model": model,
                                             "runtime": time.strftime("%Y-%m-%d %H:00Z", time.gmtime(runtime))})
    if r.status_code != 200:
        return []
    js = r.json()
    rows = js.get("data", js) if isinstance(js, dict) else js
    return rows if isinstance(rows, list) else []


def _dist_nm(a, b):
    p1, p2 = math.radians(a[0]), math.radians(b[0])
    dl = math.radians(b[1] - a[1])
    x = math.sin((p2 - p1) / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * 3440.065 * math.asin(math.sqrt(x))


def mos(scfg: dict, ctx: Context) -> SourceResult:
    model = scfg["model"]                      # GFS (MAV) or LAV (LAMP)
    c = ctx.cons
    cig_lim, vis_lim = float(c["ceiling_ft"]), float(c["vis_sm"])
    nruns = int(scfg.get("runs", 2))
    base = int(ctx.now // 21600 * 21600)
    runtimes = [base - k * 21600 for k in range(6)]   # 00/06/12/18Z runs, newest first
    members: dict[str, dict] = {}
    used = []
    for site in ctx.sites:
        cands = scfg["stations"].get(site["id"], [])
        for st in cands:
            st_id = st["id"]
            got = 0
            for rt in runtimes:
                rows = _get(st_id, model, rt)
                if not rows:
                    continue
                mid = time.strftime("%d/%HZ", time.gmtime(rt))
                for row in rows:
                    t = _ts(row.get("ftime_utc") or row.get("ftime"))   # ftime is station local time
                    if t is None or not ctx.in_window(t):
                        continue
                    rec = members.setdefault(mid, {}).setdefault(t, {})
                    cig, vis = _num(row.get("cig")), _num(row.get("vis"))
                    tmp, dpt = _num(row.get("tmp")), _num(row.get("dpt"))
                    if tmp is not None and dpt is not None and -60 < dpt <= tmp + 1 < 140:
                        rec[f"dd@{site['id']}"] = round(max(0.0, (tmp - dpt) * 5 / 9), 2)   # F to C
                    if cig in CIG_TOP_FT:
                        rec[f"c@{site['id']}"] = int(CIG_TOP_FT[cig] <= cig_lim)
                    if vis in VIS_TOP_SM:
                        rec[f"v@{site['id']}"] = int(VIS_TOP_SM[vis] <= vis_lim)
                got += 1
                if got >= nruns:
                    break
            if got:
                d = _dist_nm((site["lat"], site["lon"]), (st["lat"], st["lon"]))
                used.append(f"{site['id']}: {st_id} ({d:.0f} nm)")
                break
        else:
            used.append(f"{site['id']}: none of {', '.join(s['id'] for s in cands)} had data")
    members = {m: {t: r for t, r in s.items() if r} for m, s in members.items()}
    members = {m: s for m, s in members.items() if s}
    ok = all("none of" not in u for u in used)
    return SourceResult(members, note="; ".join(used), cycle=max(members) if members else "",
                        status=("ok" if ok else "partial") if members else "missing")
