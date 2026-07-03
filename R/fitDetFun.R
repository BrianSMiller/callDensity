# ============================================================
# Unified SNR detection model fitting
# ============================================================

#' Fit an SNR detection function
#'
#' @description
#' Fits a detection function relating probability of detection to SNR.
#'
#' Supported model types are:
#'
#' * glm
#' * gam
#' * scam
#' * vglm (capture-recapture)
#'
#' @param SNRinfo data.frame containing SNR and detection information.
#' @param modelType character.
#' @param numKnots spline basis dimension for GAM/SCAM.
#' @param yColNames observer columns for capture-recapture models.
#' @param whichObserver observer used for prediction.
#'
# ============================================================
# Unified SNR detection model fitting
# ============================================================

#' Fit an SNR detection function
#'
#' @description
#' Fits a detection function relating probability of detection to SNR.
#'
#' Supported model types are:
#'
#' * glm
#' * gam
#' * scam
#' * vglm (capture-recapture)
#'
#' @param SNRinfo data.frame containing SNR and detection information.
#' @param modelType character.
#' @param numKnots spline basis dimension for GAM/SCAM.
#' @param yColNames observer columns for capture-recapture models.
#' @param whichObserver observer used for prediction.
#'
#' @export
fitDetFun <- function(
    SNRinfo,
    modelType = c("gam", "glm", "scam", "vglm"),
    numKnots = 3,
    yColNames = c("detect_observer1", "detect_observer2"),
    whichObserver = NULL) {

  modelType <- match.arg(modelType)

  res <- switch(

    modelType,

    glm =

      stats::glm(
        Detected ~ SNR,
        data = SNRinfo,
        family = stats::binomial()
      ),

    gam =

      mgcv::gam(
        Detected ~ s(SNR, k = numKnots),
        data = SNRinfo,
        family = stats::binomial()
      ),

    scam =

      scam::scam(
        Detected ~ s(SNR, k = numKnots, bs = "mpi"),
        data = SNRinfo,
        family = stats::binomial()
      ),

    vglm = {

      fit <- VGAM::vglm(

        as.matrix(SNRinfo[, yColNames]) ~ SNR,

        VGAM::posbernoulli.t(parallel.t = TRUE ~ 0),

        data = SNRinfo
      )

      if (is.null(whichObserver))
        whichObserver <- tail(yColNames, 1)

      fit@extra$whichObserver <- whichObserver
      fit@extra$yColNames <- yColNames

      fit
    }

  )

  attr(res, "modelType") <- modelType

  if (!isS4(res)) {
    class(res) <- c("detFun", class(res))
  }

  res
}
