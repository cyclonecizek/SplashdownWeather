"""Build docs/data/board.json for the splashdown board.

    python -m pipeline.run [--only hrrr,gefs]
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time

import yaml

from . import src_mos, src_ncep, src_web
from .common import Context, PointCache, SourceResult, floor_hour, iso, log

KINDS = {"tle": src_ncep.tle, "multi_model": src_ncep.multi_model, "nbm_prob": src_ncep.nbm_prob,
         "nbm_qmd": src_ncep.nbm_qmd, "mos": src_mos.mos, "ensprob": src_ncep.ensprob,
         "openmeteo_ens": src_web.openmeteo_ens, "openmeteo_det": src_web.openmeteo_det}
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def build(cfg, only=None):
    now = int(time.time())
    t_start = floor_hour(now) - int(cfg.get("hours_back", 12)) * 3600
    t_end = floor_hour(now) + (int(cfg["days"]) + 1) * 86400
    timeline = list(range(t_start, t_end + 1, 3600))
    idx = {t: i for i, t in enumerate(timeline)}
    cache = PointCache(os.path.join(ROOT, "cache", "points.json"))
    cache.prune(now - 4 * 86400)
    sites = cfg["sites"]
    ctx = Context(now=now, lat=sites[0]["lat"], lon=sites[0]["lon"], t_start=t_start, t_end=t_end, cache=cache,
                  cons=cfg["constraints"], root=ROOT, om_refresh_h=float(cfg.get("openmeteo_refresh_hours", 3)),
                  sites=sites)
    sources = []
    for scfg in cfg["sources"]:
        if (only and scfg["id"] not in only) or not scfg.get("enabled", True):
            continue
        t0 = time.time()
        try:
            res = KINDS[scfg["kind"]](scfg, ctx)
        except Exception as e:
            log.exception("%s failed", scfg["id"])
            res = SourceResult({}, status="error", note=f"{type(e).__name__}: {e}"[:200])
        members = []
        for mid, series in sorted(res.members.items()):
            v = {}
            for t, rec in series.items():
                i = idx.get(int(t))
                if i is None:
                    continue
                for k, x in rec.items():
                    v.setdefault(k, [None] * len(timeline))[i] = round(float(x), 2) if isinstance(x, float) else x
            if v:
                members.append({"id": mid, "v": v})
        el = round(time.time() - t0, 1)
        log.info("%-10s %-8s %3d members %6.1fs  %s", scfg["id"], res.status, len(members), el, res.note)
        sources.append({"id": scfg["id"], "label": scfg.get("label", scfg["id"]), "family": scfg.get("family", "global"),
                        "role": scfg.get("role", ""),
                        "weight": float(scfg.get("weight", 1)), "status": res.status, "cycle": res.cycle,
                        "note": res.note, "seconds": el, "members": members})
    cache.save()
    return {"generated": iso(now), "generated_unix": now, "display_tz": cfg["display_tz"], "days": int(cfg["days"]),
            "window_hours": int(cfg.get("window_hours", 3)), "min_window_coverage": cfg.get("min_window_coverage", 0.5),
            "sites": sites, "constraints": cfg["constraints"], "times": timeline, "sources": sources}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=os.path.join(ROOT, "config.yaml"))
    ap.add_argument("--only", default="")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    with open(a.config) as f:
        cfg = yaml.safe_load(f)
    data = build(cfg, set(filter(None, a.only.split(","))) or None)
    if not any(s["members"] for s in data["sources"]):
        log.error("no members from any source; keeping previous output")
        return 1
    out = os.path.join(ROOT, "docs", "data", "board.json")
    with open(out + ".tmp", "w") as f:
        json.dump(data, f, separators=(",", ":"))
    os.replace(out + ".tmp", out)
    log.info("wrote board.json (%.0f kB)", os.path.getsize(out) / 1024)
    return 0


if __name__ == "__main__":
    sys.exit(main())
