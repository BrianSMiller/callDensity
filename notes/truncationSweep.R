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

  # Nc must be truncated in step with p_a.
  Nc <- sum(ch$detect_table2 & (if (is.finite(theta)) ch$SNR >= theta else TRUE),
            na.rm = TRUE)

  pa <- pDetInArea(fit, SL = d$SL, TLlookup = TL, NL = d$NL,
                   output.resolution.m = 100,
                   outerloop = outerloop,
                   truncationDistance = d$R,
                   snrTruncationThreshold = theta)$overall

  Dc <- Nc / (1 * d$A * pa * d$Time)

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
