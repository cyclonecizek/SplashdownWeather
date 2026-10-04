"""Marine layer diagnostics from a temperature profile (situational awareness).

For each member-hour: find the lowest temperature inversion in the lowest `max_m` metres.
    mlb  inversion base height (m above sea level ~ above the surface offshore)
    mls  inversion strength (C): warming from base to the top of that inversion layer
    ml   1 when the classic fog / low stratus setup is present: base at or below
         ml_base_max_m, strength at least ml_strength_min_c, and moist air under it
         (2 m dew point depression at most ml_dd_max_c, or low-level RH at least ml_rh_min)
"""
from __future__ import annotations


def inversion(levels, max_m=2000.0, min_step_c=0.2):
    """levels: [(height_m, temp_c)]. Returns (base_m, strength_c) or (None, 0.0) if none."""
    lv = sorted((z, t) for z, t in levels if z is not None and t is not None and 0 <= z <= max_m)
    for k in range(len(lv) - 1):
        if lv[k + 1][1] - lv[k][1] > min_step_c:
            j = k + 1
            while j + 1 < len(lv) and lv[j + 1][1] >= lv[j][1]:
                j += 1
            return lv[k][0], lv[j][1] - lv[k][1]
    return None, 0.0


def score(levels, cons, dd2m=None, rh_low=None):
    """Record fragment {mlb, mls, ml} for one member-hour, or {} if the profile is too thin."""
    if sum(1 for z, t in levels if z is not None and t is not None) < 3:
        return {}
    base, strength = inversion(levels, float(cons.get("ml_max_m", 2000)))
    moist = None
    if dd2m is not None:
        moist = dd2m <= float(cons.get("ml_dd_max_c", 2.5))
    if rh_low is not None:
        moist = (moist or False) or rh_low >= float(cons.get("ml_rh_min", 90))
    out = {"mls": round(strength, 2)}
    if base is not None:
        out["mlb"] = round(base, 0)
    out["ml"] = int(base is not None and base <= float(cons.get("ml_base_max_m", 500))
                    and strength >= float(cons.get("ml_strength_min_c", 3)) and bool(moist))
    return out
