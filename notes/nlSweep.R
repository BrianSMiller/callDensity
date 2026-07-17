# notes/nlSweep.R
#
# Four-parameter sweep over detector 1 and detector 2, independently, under two
# methods.
#
# Two questions:
#   1. Does nlFromDetections recover the true noise level regardless of the
#      location and scale of either detector.
#   2. Does adjudicated capture-recapture do better than observer ground truth.
#
# OG and CR are not two code paths in cde. They are two things you can put in
# detect_table1, which cde always treats as truth:
#
#   OG: detect_table1 is detector 1's raw detections. Detector 1 misses real
#       calls, falseDiscoveryRate counts those misses as detector 2's false
#       positives, and c inflates. The noise sample is filtered by detector 1's
#       curve, which is unfittable because detector 1 is assumed perfect. So
#       the NL correction has to use detector 2's curve to undo a
#       detector-1-shaped bias. That mismatch is structural.
#
#   CR: detect_table1 is the adjudicated verdict, exactly as mchToCR builds it.
#       A false positive is then a real false positive. The sample is calls
#       flagged by at least one observer, so the selecting filter is the union,
#       and the VGLM gives every component of it. The curve that did the
#       filtering is computable, so the correction can be right.
#
# All cells share one call realisation via a fixed seed, so the grid is paired:
# differences are the detectors, not the calls.
#
# FDR is zero. Under OG, c still comes out large, which is the point. A second
# sweep at FDR > 0 is the obvious next one, and tests whether CR recovers a
# non-zero c rather than merely reporting zero.
#
# Output: notes/nlSweep_results.csv, one row per (cell, method, NL variant).
# No resume. Analysis lives elsewhere. This script only makes the table.
#
# B. Miller / Claude, 2026-07-16

devtools::load_all()
source("tests/testthat/helper-fixtures.R")

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------

OUT_CSV <- "notes/nlSweep_results.csv"
SEED    <- 42
FDR     <- 0

grid <- expand.grid(
  det1location = c(0, 4),
  det1scale    = c(1, 4),
  det2location = c(0, 4),
  det2scale    = c(1, 4),
  KEEP.OUT.ATTRS = FALSE
)

# Full grid, once the cost per cell is known.
# grid <- expand.grid(
#   det1location = c(-4, 0, 4, 8),
#   det1scale    = c(0.5, 1, 2, 4),
#   det2location = c(-4, 0, 4, 8),
#   det2scale    = c(0.5, 1, 2, 4),
#   KEEP.OUT.ATTRS = FALSE
# )

# --------------------------------------------------------------------------
# Union detection function from a fitted VGLM
# --------------------------------------------------------------------------

#' Probability that at least one observer detects, as a plain function of SNR.
#'
#' This is the filter that selected an adjudicated sample: a call is in the
#' sample if anyone flagged it. VGAM computes it already, as type.fitted
#' "onempall0" (one minus the probability all observers scored zero).
#' predictDetFun.vglm uses it internally to unconditional-ise the per-observer
#' probability, but does not return it.
#'
#' Returned as a closure so it can be handed to nlFromDetections, which passes
#' it to pDetGivenNL, which has an is.function() branch.
#'
#' @param model A vglm fitted by fitDetFun(modelType = "vglm").
#' @param snrGrid SNR values at which to evaluate before interpolating.
#' @return function(snr) giving p_union. Flat extrapolation outside snrGrid.
pUnionFromVglm <- function(model, snrGrid = seq(-100, 100, by = 0.1)) {
  pu <- VGAM::predictvglm(model,
                          newdata     = data.frame(SNR = snrGrid),
                          type        = "response",
                          type.fitted = "onempall0")
  stats::approxfun(snrGrid, as.vector(pu), rule = 2)
}

# --------------------------------------------------------------------------
# Method preparation
#
# Each returns: ch (what cde sees), detFun (detector 2's curve, for p_a),
# selFun (the curve that selected the noise sample, for nlFromDetections),
# snrInfo (the noise sample), and modelType (for the cde result row).
# --------------------------------------------------------------------------

#' Observer ground truth. Detector 1 is asserted to be truth.
prepOG <- function(h, modelType = "glm") {
  ch      <- h$capHistTab
  snrInfo <- capHist2snrInfo(ch, "year")
  detFun  <- fitDetFun(snrInfo, modelType = modelType)

  # Detector 1 filtered the sample, but detector 1's curve cannot be fitted,
  # because its misses are never observed. Detector 2's curve is all there is.
  # This substitution is the structural limitation of OG, not a shortcut.
  list(ch = ch, snrInfo = snrInfo, detFun = detFun, selFun = detFun,
       modelType = modelType)
}

#' Adjudicated capture-recapture. The verdict is truth.
#'
#' Mirrors what mchToCR does with real data: subset to the rows anyone flagged,
#' then put the verdict in detect_table1. In the simulation the verdict is
#' groundTruth, which is what adjudication would recover.
prepCR <- function(h) {
  ch <- h$capHistTab

  # Adjudication only ever sees what someone flagged.
  ch <- ch[ch$detect_table1 == 1 | ch$detect_table2 == 1, ]

  # Keep both observers' raw detections; the VGLM needs them.
  ch$detect_observer1 <- as.integer(ch$detect_table1)
  ch$detect_observer2 <- as.integer(ch$detect_table2)

  # The verdict becomes truth, exactly as in mchToCR.
  ch$detect_table1 <- as.integer(ch$groundTruth1)

  # capHist2snrInfo drops the per-observer columns, so it cannot feed a VGLM.
  # Build the sample directly instead. Same filter, same quantities.
  g       <- ch[ch$detect_table1 == 1, ]
  snrInfo <- data.frame(
    detect_observer1 = g$detect_observer1,
    detect_observer2 = g$detect_observer2,
    Detected         = g$detect_observer2,
    CallRL           = as.numeric(g$signalRMSdB),
    NoiseRL          = as.numeric(g$noiseRMSdB),
    SNR              = as.numeric(g$signalRMSdB - g$noiseRMSdB),
    t                = g$t,
    month            = g$month,
    season           = g$season
  )
  snrInfo <- snrInfo[is.finite(snrInfo$SNR), ]

  # posbernoulli.t conditions on being seen by at least one observer, which is
  # exactly how the sample was drawn. A plain glm on this sample would flatter
  # detector 2, because a call both observers missed is not in it.
  detFun <- fitDetFun(snrInfo, modelType = "vglm",
                      yColNames     = c("detect_observer1", "detect_observer2"),
                      whichObserver = "detect_observer2")

  list(ch = ch, snrInfo = snrInfo, detFun = detFun,
       selFun = pUnionFromVglm(detFun), modelType = "vglm")
}

# --------------------------------------------------------------------------
# One cell, one method
# --------------------------------------------------------------------------

runVariants <- function(prep, h, TL, method, cell) {

  obsNlMean <- mean(prep$snrInfo$NoiseRL, na.rm = TRUE)

  nlVariants <- list(
    # Reference. Isolates the NL contribution from everything else.
    true  = h$NL,
    # Do nothing. Biased low: quiet gives high SNR gives detection.
    naive = data.frame(mean = obsNlMean,
                       sd   = stats::sd(prep$snrInfo$NoiseRL, na.rm = TRUE),
                       sampleSize = sum(!is.na(prep$snrInfo$NoiseRL))),
    # Incumbent. Adds the detector's 50% SNR, which is a property of the
    # detector, not of the bias.
    old   = nlFromSnrInfo(prep$snrInfo, prep$detFun),
    # Candidate. Inverts the pDetGivenNL weighting, using the curve that
    # actually selected the sample.
    new   = nlFromDetections(prep$snrInfo, prep$selFun, h$SL, TL)
  )

  rows <- lapply(names(nlVariants), function(v) {
    nl <- nlVariants[[v]]
    r  <- cde(Nc         = h$Nc,
              capHistTab = prep$ch,
              snrDetFun  = prep$detFun,
              SL         = h$SL,
              TL         = TL,
              A          = h$A,
              NL         = nl,
              modelType  = prep$modelType)
    data.frame(
      cell,
      method       = method,
      variant      = v,
      nDet         = nrow(prep$snrInfo),
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

#' Run one cell of the sweep under both methods.
nlSweepCell <- function(det1location, det1scale, det2location, det2scale,
                        TL, seed = SEED, fdr = FDR) {

  h <- make_capture_history(det1location = det1location,
                            det1scale    = det1scale,
                            det2location = det2location,
                            det2scale    = det2scale,
                            fdr          = fdr,
                            seed         = seed)

  cell <- data.frame(det1location = det1location, det1scale = det1scale,
                     det2location = det2location, det2scale = det2scale)

  rbind(
    runVariants(prepOG(h), h, TL, "OG", cell),
    runVariants(prepCR(h), h, TL, "CR", cell)
  )
}

# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------

runSweep <- function(grid, out_csv = OUT_CSV) {

  dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)

  # No resume. The cell key knows nothing about which code produced the row, so
  # a stale CSV silently survives exactly the changes that invalidate it.
  if (file.exists(out_csv)) file.remove(out_csv)

  TL <- make_tl_lookup(rangeStep = 1000)

  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]

    message(sprintf("[%d/%d] det1(loc=%g, scale=%g)  det2(loc=%g, scale=%g)",
                    i, nrow(grid), g$det1location, g$det1scale,
                    g$det2location, g$det2scale))

    t0  <- Sys.time()
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

if (sys.nframe() == 0L) {
  s <- runSweep(grid)
  print(utils::head(s))
}
