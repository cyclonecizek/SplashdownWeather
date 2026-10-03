# Splashdown Weather

Hourly-updated splashdown weather for three Southern California offshore sites
(OCN, LOS, SAND), built from HRRR, the HREF member models, NAM 3 km, NBM and the
ECMWF, AIFS, GEFS, ICON and GEM ensembles.

- **Wind (10 m):** each hour the ensemble mean and 99th percentile across all
  members decide Super GO (<= 7.25 kt and <= 10.0 kt), GO (<= 7.75 kt and
  <= 11.5 kt) or NO-GO. A window takes its worst hour.
- **Weather:** chance of precip within 10 nm, lightning within 10 nm (explicit
  lightning output only), ceiling below 500 ft and visibility below 1 sm at the
  site, plus the chance of any of them.
- **Page:** status for your window (or the next 3 hours) with a side-by-side site
  comparison, a 3-hour window grid out 7 days, the hourly wind plume and hourly
  weather chances, and source toggles.

Sites, limits and sources are in `config.yaml`.

## Set up

New repo, push with git, Pages from `/docs`, Actions workflow permissions set
to read and write, then run the workflow once and read the probe log. It lists
which fields each hi-res model has (wind, ceiling, visibility, reflectivity,
lightning) and which NBM probability thresholds were found.
