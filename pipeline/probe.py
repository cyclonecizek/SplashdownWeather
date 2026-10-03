"""Check each source: where it resolved and which fields exist.

    python -m pipeline.probe
"""
from __future__ import annotations

import os
import time

import yaml

from .common import SESSION, floor_hour, iso
from .grib import AREA, POINT, Missing, match_fields, read_idx
from .run import ROOT
from .src_ncep import _fmt


def probe(base, tmpl, cycles, show, fh=1):
    for c in cycles:
        url = f"{base}/{_fmt(tmpl, c, fh)}"
        try:
            inv = read_idx(url + ".idx")
        except Missing:
            continue
        except Exception as e:
            print(f"    {base}: {e}")
            return
        print(f"    OK {iso(c)}  {url}")
        show(inv)
        return
    print(f"    -- nothing found under {base}")


def show_fields(inv):
    found = match_fields(inv, {**POINT, **AREA})
    for k in (*POINT, *AREA):
        print(f"       {k:5s} {found[k][3] if k in found else '-- not in inventory (rows needing it skip this model)'}")


def show_nbm(inv):
    lines = [r[3] for r in inv if "prob" in r[3].lower() and any(f":{v}:" in r[3] for v in ("TSTM", "APCP", "CEIL", "VIS"))]
    print("\n".join(f"       {x}" for x in lines) or "       -- no TSTM/APCP/CEIL/VIS probability records")


def main():
    with open(os.path.join(ROOT, "config.yaml")) as f:
        cfg = yaml.safe_load(f)
    now = int(time.time())
    syn = lambda hours: [c for c in (floor_hour(now) - k * 3600 for k in range(36)) if time.gmtime(c).tm_hour in set(hours)]
    for s in cfg["sources"]:
        print(f"\n[{s['id']}] {s.get('label', '')} ({s['kind']})")
        if s["kind"] == "tle":
            probe(s["base"], s["file"], [floor_hour(now) - k * 3600 for k in range(10)], show_fields)
        elif s["kind"] == "multi_model":
            for comp in s["components"]:
                print(f"  {comp['id']}:")
                probe(comp["base"], comp["file"], syn(comp.get("cycles", [0, 6, 12, 18])), show_fields)
        elif s["kind"] == "nbm_prob":
            probe(s["base"], s["file"], syn(s.get("cycles", [0, 6, 12, 18])), show_nbm, fh=6)
        else:
            url = "https://ensemble-api.open-meteo.com/v1/ensemble" if s["kind"] == "openmeteo_ens" else "https://api.open-meteo.com/v1/forecast"
            site = cfg["sites"][0]
            r = SESSION.get(url, timeout=60, params={"latitude": site["lat"], "longitude": site["lon"], "models": s["model"],
                            "hourly": "wind_speed_10m,precipitation", "forecast_days": 2})
            print(f"    HTTP {r.status_code}" + ("" if r.ok else f": {r.text[:200]}"))
            if r.ok:
                h = r.json()["hourly"]
                for v in ("wind_speed_10m", "precipitation"):
                    print(f"    {v:16s} {sum(x is not None for x in h.get(v, []))} non-null hours")


if __name__ == "__main__":
    main()
