# Transect metadata for the Kerguelen2015 pyRAM transmission loss table

One row per transect in
[kerguelen2015TL](https://briansmiller.github.io/callDensity/reference/kerguelen2015TL.md),
giving the azimuth, source location, and the range at which that
transect's bathymetry (`cmems_mod_glo_phy_my_0.083deg_static`) runs onto
land. This is a different bathymetry product from whatever the AWR
manuscript's own AcTUP/WOA18-based
[`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
calls use, so it finds a different coastline and cannot be substituted
for that pipeline's own truncation vector, or vice versa.

## Usage

``` r
kerguelen2015TLmeta
```

## Format

A data.frame with 24 rows, one per transect:

- azimuths:

  Transect bearing, degrees, 0 to 345 in 15 degree steps.

- lonStart, latStart:

  Source location (WGS84).

- TruncationDistance_m:

  Range in metres at which this transect leaves the modelled bathymetry.
  1,000,000 where the full 1000 km transect stays in water.

- freq_hz, zr_m, call_type, site, season, year:

  Simulation parameters: 25 Hz, 20 m receiver depth, combined ABW/fin
  call type, Kerguelen2015, summer 2015.

- ssp_source, bathy_source:

  GLORYS12V1 and the static CMEMS bathymetry product used for the land
  mask.

## Source

`S:/manuscripts/2026-CopernicusPyRAM/data/output/metadata/TL_Kerguelen2015_combined_25hz_GLORYS12_summer_25Hz_metadata.csv`.
Bundled with
[kerguelen2015TL](https://briansmiller.github.io/callDensity/reference/kerguelen2015TL.md)
in `data/kerguelen2015pyram.rda`; rebuild both via
`data-raw/kerguelen2015pyram.R`.

## Details

Build the truncation vector for
[kerguelen2015TL](https://briansmiller.github.io/callDensity/reference/kerguelen2015TL.md)
by matching on azimuth, not row position:

    transectAz <- as.numeric(sub("^tl", "", names(kerguelen2015TL)[-1]))
    trunc      <- kerguelen2015TLmeta$TruncationDistance_m[
                    match(transectAz, kerguelen2015TLmeta$azimuths)]
