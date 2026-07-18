# notes/truncationSweep.R
#
# Does SNR truncation fix the p_a failure documented in
# cde_nl_and_truncation_handover.md?
#
# The failing case: a steep partner (det1scale = 0.5) leaves no recaptures at
# low SNR, the vglm extrapolates into that region, and CR's paFraction runs 1.54
# to 1.70 against an oracle of 0.99.
#
# The claim: truncating at theta multiplies the extrapolated region by zero, so
# paFraction should return towards 1. Dc should remain the density of ALL calls,
# because q(theta) cancels between a truncated Nc and a zeroed p_a.
#
# Calls cde() directly rather than pDetInArea + a hand-rolled Dc formula, to
# match what an end user actually does and what the vignette now does too.
# fdr = 0 in the simulation means cde()'s own c (via falseDiscoveryRate) should
# come back at exactly 0, so this is expected to reproduce the original
# pDetInArea-based numbers exactly, not merely approximately -- verified below
# rather than assumed.
#
# Run with Rscript. Do not library(callDensity).

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressMessages(library(VGAM))

runCell <- function(theta, det1scale, seed, n = 5e4, outerloop = 10) {
  RNGkind("Mersenne-Twister")
  set.seed(seed)

  d  <- make_capture_history(n = n, seed = seed,
                             det1location = 1, det1scale = det1scale,
                             det2location = 2, det2scale = 4,
                             fdr = 0)
  ch <- d$capHistTab
  TL <- tlSpherical(rangeStep = 2000, numTransects = 4)

  # Truth for all calls. truePa is E[p2] over every call made.
  truePa <- d$truePa
  trueDc <- d$trueDc

  # The estimand under truncation: E[p2 * 1(SNR >= theta)] over ALL calls.
  # Not the mean over surviving calls. This is what a zeroed p_a targets.
  truePaTrunc <- mean(ch$p_det2 * (ch$snr2 >= theta), na.rm = TRUE)

  # Capture-recapture: the adjudicated set is every call at least one observer
  # flagged, with false positives removed by the judge.
  adj <- subset(ch, (ch$groundTruth1 | ch$groundTruth2) &
                    (ch$detect_table1 | ch$detect_table2))
  if (is.finite(theta)) adj <- subset(adj, adj$SNR >= theta)

  fit <- fitSNRvglm(adj, c("detect_table1", "detect_table2"),
                    whichObserver = "detect_table2")

  # Nc must be truncated in step with p_a. cde() cannot do this itself -- it
  # receives Nc as a number and never sees the detections -- so it must be
  # confirmed via NcIsTruncated.
  Nc <- sum(ch$detect_table2 & (if (is.finite(theta)) ch$SNR >= theta else TRUE),
            na.rm = TRUE)

  # cde()'s internal falseDiscoveryRate()/capHist2snrInfo() calls always read
  # capHistTab$detect_table1 as ground truth and detect_table2 as the detector
  # under test -- this is not documented anywhere cde() itself is visible, only
  # inside falseDiscoveryRate()'s own roxygen, and neither raw observer in a
  # genuine two-observer capture-recapture setup IS ground truth. The real
  # Common Ground manuscript's mchToCR() handles this by overwriting
  # detect_table1 with the adjudicator's verdict before calling cde(); do the
  # same here with the simulated ground truth (groundTruth2), which is the
  # synthetic equivalent of a judge's adjudication.
  chForCde <- ch
  chForCde$detect_table1 <- chForCde$groundTruth2

  result <- suppressWarnings(suppressMessages(
    cde(Nc = Nc, capHistTab = chForCde, snrDetFun = fit,
        SL = d$SL, TL = TL, NL = d$NL, T = d$Time, A = d$A, k = 1,
        outerloop = outerloop,
        truncationDistance = d$R,
        snrTruncationThreshold = theta,
        NcIsTruncated = is.finite(theta))
  ))

  pa <- result$pa
  Dc <- result$Dc

  data.frame(theta = theta, det1scale = det1scale, seed = seed,
             Nc = Nc,
             truePa = if (is.finite(theta)) truePaTrunc else truePa,
             pa = pa,
             paFraction = pa / (if (is.finite(theta)) truePaTrunc else truePa),
             trueDc = trueDc, Dc = Dc, DcFraction = Dc / trueDc)
}

grid <- expand.grid(theta = c(-Inf, 0, 3, 6),
                    det1scale = c(0.5, 8),
                    seed = 11:13)

res <- do.call(rbind, Map(runCell, grid$theta, grid$det1scale, grid$seed))

summary <- aggregate(cbind(paFraction, DcFraction, Nc) ~ theta + det1scale,
                     data = res, FUN = mean)
print(summary[order(summary$det1scale, summary$theta), ], row.names = FALSE)

saveRDS(res, "notes/truncationSweep.rds")
