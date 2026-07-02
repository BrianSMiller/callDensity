# ============================================================
# Prediction engine for SNR detection functions
# ============================================================

#' Predict detection function from fitted SNR models
#'
#' @description
#' Unified prediction interface for SNR-based detection models.
#' Returns fitted detection probabilities over a sequence of SNR values
#' (or user-supplied newdata), optionally with 95% confidence intervals.
#'
#' Supports:
#' - glm, gam, scam (mgcv/scam-style models)
#' - VGAM vglm / vgam (posbernoulli.t models)
#'
#' @param model Fitted model from:
#'   fitSNRdetectionFunc(), fitSNRvglm(), or fitSNRvgam()
#'
#' @param newdata Data.frame with column `SNR`.
#'   If NULL, a regular grid is generated from the model frame.
#'
#' @param ci Logical. If TRUE, compute 95% confidence intervals.
#'
#' @param npoints If `newdata = NULL`, number of SNR grid points.
#'
#' @param nsim Number of simulations for VGAM confidence intervals.
#'
#' @return A data.frame with:
#'   - SNR
#'   - fit
#'   - lower (optional)
#'   - upper (optional)
#'
#' @export
predictDetFun <- function(model,
                          newdata = NULL,
                          ci = TRUE,
                          npoints = 300,
                          nsim = 500) {

  if (is.null(newdata)) {
    rng <- range(model.frame(model)$SNR, na.rm = TRUE)
    newdata <- data.frame(
      SNR = seq(rng[1], rng[2], length.out = npoints)
    )
  }

  cls <- class(model)

  # ------------------------------------------------------------
  # GLM / GAM / SCAM
  # ------------------------------------------------------------
  if (inherits(model, c("glm", "gam", "scam"))) {

    fit <- stats::predict(model,
                          newdata = newdata,
                          type = "response")

    out <- data.frame(
      SNR = newdata$SNR,
      fit = fit
    )

    if (ci) {
      pr <- stats::predict(model,
                           newdata = newdata,
                           type = "link",
                           se.fit = TRUE)

      crit <- 1.96

      link_fit <- pr$fit
      se <- pr$se.fit

      out$lower <- stats::plogis(link_fit - crit * se)
      out$upper <- stats::plogis(link_fit + crit * se)
    }

    return(out)
  }

  # ------------------------------------------------------------
  # VGAM (vglm / vgam)
  # ------------------------------------------------------------
  if (inherits(model, c("vglm", "vgam"))) {

    whichObserver <- model@extra$whichObserver

    # point prediction
    preds0 <- VGAM::predict(
      model,
      type = "response",
      newdata = newdata,
      type.fitted = "onempall0"
    )

    preds1 <- VGAM::predict(
      model,
      type = "response",
      newdata = newdata
    )

    ix <- which(colnames(preds1) == whichObserver)

    fit <- as.numeric(preds0 * preds1[, ix])

    out <- data.frame(
      SNR = newdata$SNR,
      fit = fit
    )

    # --------------------------------------------------------
    # CI via simulation (approximate, but stable)
    # --------------------------------------------------------
    if (ci) {

      beta_hat <- coef(model)
      V <- vcov(model)

      sims <- MASS::mvrnorm(nsim, mu = beta_hat, Sigma = V)

      sim_mat <- matrix(NA_real_, nrow = nsim, ncol = nrow(newdata))

      for (i in seq_len(nsim)) {

        tmp <- model
        tmp@coefficients <- sims[i, ]

        p0 <- VGAM::predict(
          tmp,
          type = "response",
          newdata = newdata,
          type.fitted = "onempall0"
        )

        p1 <- VGAM::predict(
          tmp,
          type = "response",
          newdata = newdata
        )

        sim_mat[i, ] <- p0 * p1[, ix]
      }

      out$lower <- apply(sim_mat, 2, stats::quantile, 0.025, na.rm = TRUE)
      out$upper <- apply(sim_mat, 2, stats::quantile, 0.975, na.rm = TRUE)
    }

    return(out)
  }

  stop("Unsupported model class: ", paste(cls, collapse = ", "))
}
