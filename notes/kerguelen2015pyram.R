# data-raw/kerguelen2015pyram.R
#
# One-off script to bundle the pyRAM TL table and its truncation-distance
# metadata for Kerguelen2015 (summer, 25 Hz, GLORYS12V1) into the package, so
# the noiseLevels vignette no longer depends on files living on Brian's S:
# drive. Run this once, from the package root, whenever the source CSVs
# change (e.g. CopernicusPyRAM is re-run for this site/season).
#
# This script is NOT run automatically -- it is not part of the package build,
# it just produces the .rda files under data/. See R/data.R for the
# documentation these two objects are meant to satisfy.

tlFile <- file.path("s:/manuscripts/2026-CopernicusPyRAM/data/output",
                    "TL_Kerguelen2015_combined_25hz_GLORYS12_summer_25Hz.csv")
metaFile <- file.path("s:/manuscripts/2026-CopernicusPyRAM/data/output/metadata",
                      "TL_Kerguelen2015_combined_25hz_GLORYS12_summer_25Hz_metadata.csv")

kerguelen2015TL     <- read.csv(tlFile)
kerguelen2015TLmeta <- read.csv(metaFile)

# Same sanity checks the vignette runs on load -- fail here, not silently
# downstream, if the source files ever change shape.
stopifnot(
  "TL table should have a range column plus 24 transects" =
    ncol(kerguelen2015TL) == 25,
  "TL table's first column should be named range_m" =
    names(kerguelen2015TL)[1] == "range_m",
  "Metadata should have one row per transect" =
    nrow(kerguelen2015TLmeta) == ncol(kerguelen2015TL) - 1,
  "Every transect column must have a matching azimuth in the metadata" =
    !anyNA(match(as.numeric(sub("^tl", "", names(kerguelen2015TL)[-1])),
                kerguelen2015TLmeta$azimuths))
)

# Bundled into one file rather than usethis::use_data()'s default of one
# object per .rda. That default (documented in R Packages, Wickham & Bryan)
# exists so a dataset's name always matches the file that defines it -- a real
# convention, not just tidiness for its own sake. Deviating from it here is a
# deliberate one-off for a single paired example, not a template: the TL/
# truncation split is an artefact of this dataset's own two source files, not
# a shape the package's functions require. pDetGivenNL and friends take a TL
# data.frame and a truncation vector as separate arguments regardless of how
# either is stored on disk, so nothing about the function interface is
# affected by bundling the two objects together here.
save(kerguelen2015TL, kerguelen2015TLmeta,
    file = "data/kerguelen2015pyram.rda")

# If data/kerguelen2015TL.rda and data/kerguelen2015TLmeta.rda exist from an
# earlier run of this script (two usethis::use_data() calls), delete them --
# otherwise the same two objects would be lazy-loaded from two places.
old <- c("data/kerguelen2015TL.rda", "data/kerguelen2015TLmeta.rda")
old <- old[file.exists(old)]
if (length(old) > 0) {
  message("Removing superseded file(s): ", paste(old, collapse = ", "))
  file.remove(old)
}
