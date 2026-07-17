#' Predict the mean noise level you would measure at detections
#'
#' Given a candidate true noise distribution, predict the mean of the noise
#' levels that would be measured at detected calls. Those measurements are
#' biased low, because a quiet period gives a high SNR and a high SNR gives a
#' detection, so detections over-represent quiet periods.
#'
#' The weighting is `pDetGivenNL`: if the chance of detection at 80 dB is 2.5
#' times the chance at 84 dB, then 80 dB periods appear among the detections
#' 2.5 times more often than they occur.
#'
#' @param mu Candidate true mean noise level in dB.
#' @param sigma Standard deviation of the true noise distribution in dB.
#' @param detFun Detection function. See `pDetGivenNL`.
#' @param SL Source level distribution, with elements named mean and sd.
#' @param TL Transmission loss table. See `pDetGivenNL`.
#' @param truncationDistance Scalar or one value per transect, in metres.
#' @param nNodes Number of nodes used to integrate over the noise distribution.
#' @param ... Passed to `pDetGivenNL`, e.g. `binWidth`.
#'
#' @returns Scalar. Always less than or equal to `mu`.
#' @importFrom stats dnorm
#' @export
predictSampledNL <- function(mu, sigma, detFun, SL, TL,
                             truncationDistance = max(TL[[1]]),
                             nNodes = 121, ...) {
  # No variation in noise means no tilt, so nothing to correct.
  if (!is.finite(sigma) || sigma <= 0) return(mu)

  nl <- seq(mu - 6 * sigma, mu + 6 * sigma, length.out = nNodes)
  wt <- dnorm(nl, mu, sigma) *
    pDetGivenNL(nl, detFun, TL = TL, SL = SL,
                truncationDistance = truncationDistance, ...)
  if (sum(wt) <= 0) {
    stop("Nothing is detectable anywhere near this noise level. ",
         "Check that SL, TL and NL are in consistent units.")
  }
  sum(nl * wt) / sum(wt)
}


#' Estimate the noise level distribution from noise measured at detections
#'
#' Noise levels measured at detections are biased low, because detections
#' over-represent quiet periods. This undoes that bias.
#'
#' Replaces `nlFromSnrInfo`, which corrected the same bias by adding the SNR at
#' which the detection function reaches 0.5. That quantity is a property of the
#' detector. The bias is a property of the propagation geometry and the noise
#' variance. The two coincide only by chance. See the noiseLevels vignette.
#'
#' The standard deviation is taken from the measurements directly. Filtering by
#' detection shifts the mean but barely narrows the distribution, so the
#' measured standard deviation is close to the truth even though the measured
#' mean is not. This leaves one unknown, found by `stats::uniroot`.
#'
#' @param snrInfo Table of SNR information containing a column of noise level
#'   measurements in dB.
#' @param snrDetFun Detection function, as passed to `pDetInArea`.
#' @param SL Source level distribution, with elements named mean and sd.
#' @param TL Transmission loss table. First column ranges in metres, remaining
#'   columns TL in dB per radial transect.
#' @param truncationDistance Scalar or one value per transect, in metres.
#' @param nlColumn Name of the noise level column. Default "NoiseRL".
#' @param searchWidth Width in dB of the interval searched above the measured
#'   mean. The bias cannot be negative, so the search runs upwards only.
#' @param ... Passed to `pDetGivenNL`, e.g. `binWidth`.
#'
#' @returns Data.frame with one row and columns mean, sd and sampleSize, the
#'   same format as `nlFromSnrInfo` and `noiseLevelDistribution`.
#'
#' @importFrom stats sd uniroot
#' @export
nlFromDetections <- function(snrInfo, snrDetFun, SL, TL,
                             truncationDistance = max(TL[[1]]),
                             nlColumn = "NoiseRL",
                             searchWidth = 25, ...) {

  if (is.null(snrInfo[[nlColumn]])) {
    stop(sprintf("snrInfo has no column named '%s'", nlColumn))
  }
  nlObs   <- snrInfo[[nlColumn]]
  obsMean <- mean(nlObs, na.rm = TRUE)
  sigma   <- sd(nlObs, na.rm = TRUE)
  n       <- sum(!is.na(nlObs))

  if (!is.finite(obsMean)) stop("No usable noise level measurements.")

  # With no spread there is no tilt, so the measurements are unbiased.
  if (!is.finite(sigma) || sigma <= 0) {
    return(data.frame(mean = obsMean, sd = 0, sampleSize = n))
  }

  gap <- function(mu) {
    predictSampledNL(mu, sigma, snrDetFun, SL, TL, truncationDistance, ...) -
      obsMean
  }

  # The bias is never negative, so the true mean is at or above the measured
  # one. gap() is therefore <= 0 at the lower end and rises with mu.
  hi <- gap(obsMean + searchWidth)
  if (hi < 0) {
    stop(sprintf(paste("The noise bias appears to exceed searchWidth (%g dB).",
                       "Either the study area is far larger than the detection",
                       "range, or SL, TL and NL are not in consistent units.",
                       "Increase searchWidth if you believe the bias is",
                       "really this large."), searchWidth))
  }

  mu <- uniroot(gap, interval = c(obsMean, obsMean + searchWidth))$root

  data.frame(mean = mu, sd = sigma, sampleSize = n)
}
