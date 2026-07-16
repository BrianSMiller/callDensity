# notes/nlSweep.R
#
# Four-parameter sweep over detector 1 and detector 2, independently.
#
# Question: does nlFromDetections recover the true noise level, and does that
# recovery improve p_a and Dc, regardless of the location and scale of either
# detector.
#
# The previous sweep (test-nlSweep.R) moved both detectors together with a
# fixed 1 dB offset and never varied scale. That is the one configuration
# blind to the two things that could break the method: detector scale, and
# mismatch between the filter that selected the noise sample and the curve
# fitted to it. This sweep opens both gaps deliberately.
#
# The mismatch, stated plainly: capHist2snrInfo returns detector 1's
# detections. So detector 1 is the filter that skewed the noise sample. But
# detector 1 is assumed to be ground truth, so its misses are never seen and
# its curve cannot be fitted. The only curve available is detector 2's. Every
# NL correction here uses a det2-derived curve to undo a det1-shaped bias.
#
# False positives are off. With c = 0 the (1-c) term in the density equation
# is 1, so any error in Dc is an error in p_a, which is an error in NL.
# Nothing else can contribute.
#
# All cells share one call realisation via a fixed seed, so the design is
# paired: differences across the grid are the detectors, not the calls. It
# also means the whole grid rests on one draw. A second seed is a later check.
#
# Output: notes/nlSweep_results.csv, one row per (cell, NL variant).
# Resumable. Analysis lives elsewhere. This script only makes the table.
#
# Requires helper-fixtures.R to return truePa.
#
# B. Miller / Claude, 2026-07-16

library(callDensity)

# make_capture_history and make_tl_lookup live with the tests. Run this script
# from the package root.
source("tests/testthat/helper-fixtures.R")

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------

OUT_CSV   <- "notes/nlSweep_results.csv"
MODELTYPE <- "glm"
SEED      <- 42

# Grid. Start small, widen once the cost per cell is known.
grid <- expand.grid(
  det1location = c(0, 4),
  det1scale    = c(1, 4),
  det2location = c(0, 4),
  det2scale    = c(1, 4),
  KEEP.OUT.ATTRS = FALSE
)

# Full grid, for when the cost is known. Uncomment to use.
# grid <- expand.grid(
#   det1location = c(-4, 0, 4, 8),
#   det1scale    = c(0.5, 1, 2, 4),
#   det2location = c(-4, 0, 4, 8),
#   det2scale    = c(0.5, 1, 2, 4),
#   KEEP.OUT.ATTRS = FALSE
# )

# --------------------------------------------------------------------------
# One cell
# --------------------------------------------------------------------------

#' Run one cell of the sweep.
#'
#' Builds a capture history for the given detector pair, fits detector 2's
#' curve treating detector 1 as ground truth (observer ground truth), derives
#' four candidate NL distributions, and runs cde with each.
#'
#' @param det1location,det1scale Logistic parameters for detector 1.
#' @param det2location,det2scale Logistic parameters for detector 2.
#' @param TL        TL lookup table.
#' @param modelType Passed to fitDetFun and cde.
#' @param seed      Passed to make_capture_history. Fixed across cells.
#' @return Data frame, one row per NL variant.
nlSweepCell <- function(det1location, det1scale,
                        det2location, det2scale,
                        TL, modelType = MODELTYPE, seed = SEED) {

  h <- make_capture_history(det1location = det1location,
                            det1scale    = det1scale,
                            det2location = det2location,
                            det2scale    = det2scale,
                            fdr          = 0,
                            seed         = seed)

  # SNRinfo is detector 1's detections; detFun is detector 2's curve.
  SNRinfo <- capHist2snrInfo(h$capHistTab, "year")
  detFun  <- fitDetFun(SNRinfo, modelType = modelType)

  obsNlMean <- mean(SNRinfo$NoiseRL, na.rm = TRUE)

  nlVariants <- list(
    # Reference. Isolates NL error from everything else cde does.
    true  = h$NL,
    # Do nothing. The noise as measured at detector 1's detections, biased low
    # because quiet periods give high SNR which gives detections.
    naive = data.frame(mean = obsNlMean,
                       sd   = stats::sd(SNRinfo$NoiseRL, na.rm = TRUE),
                       sampleSize = sum(!is.na(SNRinfo$NoiseRL))),
    # The incumbent. Adds the SNR at which detFun reaches 0.5. That quantity
    # is a property of the detector. The bias is not.
    old   = suppressWarnings(nlFromSnrInfo(SNRinfo, detFun)),
    # The candidate. Inverts the pDetGivenNL weighting.
    new   = nlFromDetections(SNRinfo, detFun, h$SL, TL)
  )

  rows <- lapply(names(nlVariants), function(v) {
    nl <- nlVariants[[v]]
    r  <- cde(Nc         = h$Nc,
              capHistTab = h$capHistTab,
              snrDetFun  = detFun,
              SL         = h$SL,
              TL         = TL,
              A          = h$A,
              NL         = nl,
              modelType  = modelType)
    data.frame(
      det1location = det1location, det1scale = det1scale,
      det2location = det2location, det2scale = det2scale,
      variant      = v,
      nDet         = nrow(SNRinfo),
      obsNlMean    = obsNlMean,
      nlMean       = nl$mean,
      nlSd         = nl$sd,
      pa           = r$pa,
      Dc           = r$Dc,
      c            = r$c,
      trueNlMean   = h$NL$mean,
      trueNlSd     = h$NL$sd,
      truePa       = h$truePa,
      trueDc       = h$trueDc,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------

cellKey <- function(d) {
  sprintf("%g_%g_%g_%g", d$det1location, d$det1scale,
          d$det2location, d$det2scale)
}

runSweep <- function(grid, out_csv = OUT_CSV) {

  dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)

  TL <- make_tl_lookup(rangeStep = 1000)

  done <- character(0)
  if (file.exists(out_csv)) {
    prev <- utils::read.csv(out_csv, stringsAsFactors = FALSE)
    done <- unique(cellKey(prev))
    message(sprintf("Resuming: %d of %d cells already done.",
                    length(done), nrow(grid)))
  }

  for (i in seq_len(nrow(grid))) {
    g   <- grid[i, ]
    key <- cellKey(g)
    if (key %in% done) next

    message(sprintf("[%d/%d] det1(loc=%g, scale=%g)  det2(loc=%g, scale=%g)",
                    i, nrow(grid), g$det1location, g$det1scale,
                    g$det2location, g$det2scale))

    t0 <- Sys.time()
    res <- tryCatch(
      nlSweepCell(g$det1location, g$det1scale,
                  g$det2location, g$det2scale, TL = TL),
      error = function(e) {
        message(sprintf("  FAILED: %s", conditionMessage(e)))
        NULL
      }
    )
    if (is.null(res)) next

    message(sprintf("  %.1f s", as.numeric(difftime(Sys.time(), t0,
                                                    units = "secs"))))

    utils::write.table(res, out_csv, sep = ",", row.names = FALSE,
                       col.names = !file.exists(out_csv),
                       append = file.exists(out_csv))
  }

  invisible(utils::read.csv(out_csv, stringsAsFactors = FALSE))
}

# --------------------------------------------------------------------------
# Monte Carlo noise check. Run this before trusting the grid.
#
# cde estimates p_a by Monte Carlo integration. If that noise is comparable to
# the differences between NL variants, the grid measures the integrator rather
# than the method.
# --------------------------------------------------------------------------

mcNoiseCheck <- function(reps = 5, TL = make_tl_lookup(rangeStep = 1000)) {
  h <- make_capture_history(det1location = 1, det1scale = 2,
                            det2location = 2, det2scale = 4,
                            fdr = 0, seed = SEED)
  SNRinfo <- capHist2snrInfo(h$capHistTab, "year")
  detFun  <- fitDetFun(SNRinfo, modelType = MODELTYPE)

  pa <- vapply(seq_len(reps), function(i) {
    cde(Nc = h$Nc, capHistTab = h$capHistTab, snrDetFun = detFun,
        SL = h$SL, TL = TL, A = h$A, NL = h$NL, modelType = MODELTYPE)$pa
  }, numeric(1))

  message(sprintf("p_a over %d identical calls: mean %.5f, sd %.5f, range %.5f",
                  reps, mean(pa), stats::sd(pa), diff(range(pa))))
  pa
}

# --------------------------------------------------------------------------

if (sys.nframe() == 0L) {
  s <- runSweep(grid)
  print(head(s))
}
