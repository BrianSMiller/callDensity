#' Probability of detection in the study area at a fixed noise level
#'
#' Answers the question that `pDetInArea` answers, but with the noise level
#' held still rather than drawn from a distribution: if the noise were always
#' exactly `nl`, what fraction of the calls in the study area would be
#' detected?
#'
#' Evaluated over a grid of noise levels this gives a curve, and that curve is
#' what makes the noise levels measured at detections biased. Quiet periods are
#' over-represented among detections by exactly the ratio the curve describes.
#' See `nlFromDetections` for the use, and the noiseLevels vignette for the
#' derivation.
#'
#' Transmission loss enters only as a lookup table, so this makes no assumption
#' about the propagation model. Under spherical spreading the curve falls by one
#' decade per 10 dB, but that is a consequence rather than an input.
#'
#' @param nl Vector of noise levels in dB.
#' @param detFun Either a detFun object from `fitDetFun`, or a plain function
#'   of SNR returning probability of detection.
#' @param SL List or data.frame containing the distribution of source levels,
#'   with elements named mean and sd.
#' @param TL Data.frame of transmission losses. The first column contains
#'   ranges in metres, the remaining columns contain TL in dB for each radial
#'   transect at that range. Same format as `cde` expects.
#' @param truncationDistance Scalar or vector of truncation distances in
#'   metres. If a vector, one value per transect. Cells beyond the truncation
#'   distance carry no weight, and the returned probability is relative to the
#'   truncated area.
#' @param nSLnodes Number of quadrature nodes used to average over the source
#'   level distribution. Default 41.
#' @param binWidth Width in dB of the transmission loss bins. The only
#'   approximation in this function. Default 0.25, which is conservative:
#'   results are typically stable to four significant figures at 1 dB.
#'
#' @returns Numeric vector of the same length as `nl`, each between 0 and 1.
#'
#' @importFrom stats dnorm approxfun
#' @export
pDetGivenNL <- function(nl, detFun, SL, TL,
                        truncationDistance = max(TL[[1]]),
                        nSLnodes = 41,
                        binWidth = 0.25) {

  TL <- as.data.frame(TL)
  stopifnot("TL needs a range column and at least one transect" = ncol(TL) >= 2)
  stopifnot("SL needs elements named mean and sd" =
              all(c("mean", "sd") %in% names(SL)))

  range_m <- TL[[1]]
  tlMat   <- as.matrix(TL[, -1, drop = FALSE])

  # Area weights. Each range step stands for an annulus, and the area of an
  # annulus grows with r, so a step at 200 km carries ten times the weight of a
  # step at 20 km. Cells beyond the truncation distance carry no weight.
  truncationDistance <- rep_len(truncationDistance, ncol(tlMat))
  wMat <- outer(range_m, truncationDistance, function(r, td) r * (r <= td))

  # Collapse the table by transmission loss. A cell at 300 km on one transect
  # and a cell at 280 km on another contribute identically if their TL is the
  # same, so all that matters is the total area weight at each TL. This turns
  # hundreds of thousands of cells into a few hundred bins. It is bookkeeping,
  # not an approximation of the physics; only binWidth costs anything.
  tlVec <- as.vector(tlMat)
  wVec  <- as.vector(wMat)
  ok    <- is.finite(tlVec) & wVec > 0
  tlVec <- tlVec[ok]
  wVec  <- wVec[ok]

  brks  <- seq(floor(min(tlVec)), ceiling(max(tlVec)) + binWidth, by = binWidth)
  idx   <- .bincode(tlVec, brks, include.lowest = TRUE)
  binW  <- as.vector(tapply(wVec,
                            factor(idx, levels = seq_len(length(brks) - 1)),
                            sum))
  binW[is.na(binW)] <- 0
  binTL <- brks[-length(brks)] + binWidth / 2

  keep   <- binW > 0
  binTL  <- binTL[keep]
  binW   <- binW[keep]
  totalW <- sum(binW)

  # Source levels. Calls are not all equally loud, so average over the SL
  # distribution using a weighted grid.
  slNodes <- seq(SL$mean - 4 * SL$sd, SL$mean + 4 * SL$sd, length.out = nSLnodes)
  slWt    <- dnorm(slNodes, SL$mean, SL$sd)
  slWt    <- slWt / sum(slWt)

  # The detection function depends on one number, SNR. Evaluate it once on a
  # grid and interpolate, rather than calling predict() repeatedly.
  pLookup <- if (is.function(detFun)) {
    detFun
  } else {
    snrGrid <- seq(-100, 100, by = 0.1)
    pGrid   <- predictDetFun(detFun, newdata = data.frame(SNR = snrGrid))$fit
    approxfun(snrGrid, pGrid, rule = 2)
  }

  vapply(nl, function(thisNL) {
    acc <- 0
    for (i in seq_along(slNodes)) {
      # sonar equation, once per TL bin rather than once per cell
      acc <- acc + slWt[i] * sum(pLookup(slNodes[i] - binTL - thisNL) * binW)
    }
    acc / totalW
  }, numeric(1))
}
