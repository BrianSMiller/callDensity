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
#' @export
predictDetFun <- function(model, ...) {
  UseMethod("predictDetFun")
}

#' @export
predictDetFun.glm <- function(model,
                              newdata = NULL,
                              ci = TRUE,
                              level = 0.95,
                              npoints = 300,
                              ...) {

  if (is.null(newdata)) {

    rng <- range(model.frame(model)$SNR, na.rm = TRUE)

    newdata <- data.frame(
      SNR = seq(rng[1], rng[2], length.out = npoints)
    )
  }

  fit <- stats::predict(
    model,
    newdata = newdata,
    type = "response"
  )

  out <- data.frame(
    SNR = newdata$SNR,
    fit = fit
  )

  if (ci) {

    pr <- stats::predict(
      model,
      newdata = newdata,
      type = "link",
      se.fit = TRUE
    )

    crit <- stats::qnorm((1 + level) / 2)

    out$lower <- stats::plogis(pr$fit - crit * pr$se.fit)
    out$upper <- stats::plogis(pr$fit + crit * pr$se.fit)

  }

  out
}

#' @export
predictDetFun.gam <- predictDetFun.glm

#' @export
predictDetFun.scam <- predictDetFun.glm

#' @export
predictDetFun.vglm <- function(model,
                               newdata = NULL,
                               ci = FALSE,
                               level = 0.95,
                               npoints = 300,
                               whichObserver = model@extra$whichObserver,
                               ...) {

  if (is.null(newdata)) {

    rng <- range(model@x[, "SNR"], na.rm = TRUE)

    newdata <- data.frame(
      SNR = seq(rng[1], rng[2], length.out = npoints)
    )
  }

  preds.any <- VGAM::predict(
    model,
    type = "response",
    newdata = newdata,
    type.fitted = "onempall0"
  )

  preds.obs <- VGAM::predict(
    model,
    type = "response",
    newdata = newdata
  )

  ix <- match(whichObserver, colnames(preds.obs))

  fit <- preds.any * preds.obs[, ix]

  out <- data.frame(
    SNR = newdata$SNR,
    fit = fit
  )

  if (ci) {

    warning(
      "Confidence intervals are not yet implemented for vglm models.",
      call. = FALSE
    )

    out$lower <- NA_real_
    out$upper <- NA_real_
  }

  out
}

#' @export
predictDetFun.default <- function(model, ...) {

  stop(
    "No predictDetFun() method for objects of class ",
    paste(class(model), collapse = "/"),
    call. = FALSE
  )
}

#' @export
predict.detFun <- function(object, ...) {
  # Strip the detFun class and re-dispatch to the underlying model's predict method
  class(object) <- class(object)[class(object) != "detFun"]
  predict(object, ...)
}

# ============================================================
# Unified multi-model prediction engine
# ============================================================

#' Predict detection functions for multiple models
#'
#' @description
#' Evaluates one or more fitted detection models on a shared SNR grid and
#' returns a long-format data frame suitable for ggplot or base plotting.
#'
#' @param models Named list of fitted models.
#' @param newdata Optional data.frame with SNR column. If NULL, a grid is created.
#' @param npoints Number of SNR grid points if newdata is NULL.
#' @param ci Logical; include confidence intervals where available.
#'
#' @return data.frame with columns: SNR, model, fit, lower, upper
#'
#' @export
predictDetFunList <- function(models,
                              newdata = NULL,
                              ci = TRUE,
                              npoints = 300,
                              ...) {

  if (is.null(newdata)) {

    rng <- range(unlist(
      lapply(models, function(m)
        model.frame(m)$SNR)
    ), na.rm = TRUE)

    newdata <- data.frame(
      SNR = seq(rng[1], rng[2], length.out = npoints)
    )
  }

  out <- do.call(

    rbind,

    lapply(names(models), function(nm) {

      p <- predictDetFun(
        models[[nm]],
        newdata = newdata,
        ci = ci,
        ...
      )

      p$model <- nm
      p

    })
  )

  rownames(out) <- NULL

  out
}
