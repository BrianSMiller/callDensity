# ============================================================
# Unified SNR detection model fitting interface
# ============================================================

#' Fit SNR detection function (unified interface)
#'
#' @description
#' Unified wrapper for fitting SNR-based detection functions using:
#' - glm
#' - gam (mgcv)
#' - scam
#' - vglm (VGAM capture-recapture model)
#'
#' @param SNRinfo Data.frame containing:
#'   - SNR
#'   - Detected OR observer columns (for vglm)
#'
#' @param modelType One of:
#'   "glm", "gam", "scam", "vglm"
#'
#' @param yColNames Observer columns for vglm models
#'
#' @param whichObserver Observer used for marginal detection prediction
#'
#' @param numKnots Smoothing basis size for GAM/SCAM
#'
#' @return Fitted model object
#'
#' @export
fitDetFun <- function(SNRinfo,
                      modelType = c("glm", "gam", "scam", "vglm"),
                      yColNames = c("detect_observer1", "detect_observer2"),
                      whichObserver = NULL,
                      numKnots = 3) {

  modelType <- match.arg(modelType)

  # ------------------------------------------------------------
  # Single-observer models (binomial GLM/GAM/SCAM)
  # ------------------------------------------------------------

  if (modelType %in% c("glm", "gam", "scam")) {

    fml <- stats::as.formula("Detected ~ SNR")

    res <- switch(modelType,

                  glm = stats::glm(
                    fml,
                    data = SNRinfo,
                    family = stats::binomial(link = "logit")
                  ),

                  gam = mgcv::gam(
                    fml + mgcv::s(SNR, k = numKnots),
                    data = SNRinfo,
                    family = stats::binomial
                  ),

                  scam = scam::scam(
                    Detected ~ scam::s(SNR, k = numKnots, bs = "mpi"),
                    data = SNRinfo,
                    family = stats::binomial
                  )
    )

    attr(res, "modelType") <- modelType

    return(res)
  }

  # ------------------------------------------------------------
  # Capture-recapture model (VGAM vglm)
  # ------------------------------------------------------------

  if (modelType == "vglm") {

    res <- VGAM::vglm(
      as.matrix(SNRinfo[, yColNames]) ~ SNR,
      VGAM::posbernoulli.t(parallel.t = TRUE ~ 0),
      data = SNRinfo
    )

    attr(res, "modelType") <- "vglm"

    if (is.null(whichObserver)) {
      whichObserver <- yColNames[length(yColNames)]
    }

    res@extra$whichObserver <- whichObserver
    res@extra$yColNames <- yColNames

    return(res)
  }

  stop("Unsupported modelType")
}
