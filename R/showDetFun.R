# ============================================================
# Plot SNR detection function with observed SNR distributions
# ============================================================

#' Plot SNR detection function with observed SNR distributions
#'
#' @description
#' Visualises a fitted SNR detection function together with the distribution
#' of detected and missed events. By default draws a rug along the y=1
#' (detected) and y=0 (missed) reference lines. Optionally overlays mirrored
#' density curves or histograms above and below those lines.
#'
#' @param model Fitted detection model from fitDetFun().
#'
#' @param SNRinfo Data frame used to fit the model. Must contain columns
#'   \code{SNR} and \code{Detected} (logical or 0/1).
#'
#' @param distribution One of \code{"none"} (default), \code{"density"}, or
#'   \code{"histogram"}. Controls whether a mirrored distribution is drawn
#'   above y=1 (detected) and below y=0 (missed). A rug is always drawn
#'   unless \code{rug = FALSE}.
#'
#' @param rug Logical. If \code{TRUE} (default), draws tick marks along the
#'   y=1 line for detected events and along y=0 for missed events.
#'   For large datasets, a random subsample of \code{rug.max} points is used.
#'
#' @param rug.max Maximum number of rug ticks per group. Default 500.
#'   Set to \code{Inf} to show all points (slow for n > 2000).
#'
#' @param mirror.height Maximum height of mirrored distributions as a
#'   fraction of the y-axis range. Distributions sit outside the [0,1] box.
#'   Default 0.20.
#'
#' @param rug.alpha Alpha for rug ticks. Default 0.4.
#'
#' @param npoints Number of SNR grid points for the fitted curve. Default 300.
#'
#' @param show.counts Logical. If \code{TRUE} (default), adds a subtitle
#'   reporting the number of detections out of total observations.
#'
#' @return A ggplot object.
#'
#' @export
showDetFun <- function(model, ...) {
  if (inherits(model, c("glm", "gam", "scam", "vglm", "vgam"))) {
    showDetFun.detFun(model, ...)
  } else {
    showDetFun.list(model, ...)
  }
}

#' @export
showDetFun.detFun <- function(model,
                               SNRinfo,
                               distribution  = c("none", "density", "histogram"),
                               rug           = TRUE,
                               rug.max       = 500,
                               mirror.height = 0.20,
                               rug.alpha     = 0.4,
                               npoints       = 300,
                               show.counts   = TRUE) {

  distribution <- match.arg(distribution)
  requireNamespace("ggplot2")

  # ------------------------------------------------------------------
  # Data preparation
  # ------------------------------------------------------------------

  obs          <- SNRinfo
  obs$Detected <- as.logical(obs$Detected)
  det_snr      <- obs$SNR[obs$Detected  & !is.na(obs$SNR)]
  miss_snr     <- obs$SNR[!obs$Detected & !is.na(obs$SNR)]
  n_det        <- length(det_snr)
  n_total      <- sum(!is.na(obs$SNR))

  # subsample rug points for large datasets
  subsample <- function(x, n) {
    if (is.finite(n) && length(x) > n) x[sample(length(x), n)] else x
  }
  rug_det  <- subsample(det_snr,  rug.max)
  rug_miss <- subsample(miss_snr, rug.max)

  # ------------------------------------------------------------------
  # Fitted curve
  # ------------------------------------------------------------------

  pred <- predictDetFun(
    model   = model,
    newdata = data.frame(SNR = seq(min(obs$SNR, na.rm = TRUE),
                                    max(obs$SNR, na.rm = TRUE),
                                    length.out = npoints)),
    ci = TRUE
  )

  # ------------------------------------------------------------------
  # Base plot: y-axis shows only 0 and 1; distributions sit outside
  # the box via coord_cartesian
  # ------------------------------------------------------------------

  p <- ggplot2::ggplot(pred, ggplot2::aes(SNR, fit)) +
    ggplot2::geom_hline(yintercept = c(0, 1),
                        linetype   = 2,
                        colour     = "grey60") +
    ggplot2::scale_y_continuous(
      breaks = c(0, 1),
      labels = c("0", "1"),
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(
      ylim = c(0 - mirror.height, 1 + mirror.height),
      clip = "off"
    ) +
    ggplot2::labs(x = "SNR (dB)", y = "Detection probability") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  # confidence band
  if (all(c("lower", "upper") %in% names(pred))) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower, ymax = upper),
        alpha       = 0.2,
        fill        = "grey60",
        inherit.aes = TRUE
      )
  }

  p <- p + ggplot2::geom_line(linewidth = 1.2, colour = "black")

  # ------------------------------------------------------------------
  # Rug: subsampled ticks along y=1 (detected) and y=0 (missed)
  # ------------------------------------------------------------------

  if (rug) {
    tick_h <- mirror.height * 0.3   # tick height in data units

    if (length(rug_det) > 0) {
      p <- p +
        ggplot2::geom_segment(
          data        = data.frame(SNR = rug_det),
          ggplot2::aes(x = SNR, xend = SNR,
                       y = 1, yend = 1 + tick_h),
          colour      = "steelblue",
          alpha       = rug.alpha,
          linewidth   = 0.3,
          inherit.aes = FALSE
        )
    }
    if (length(rug_miss) > 0) {
      p <- p +
        ggplot2::geom_segment(
          data        = data.frame(SNR = rug_miss),
          ggplot2::aes(x = SNR, xend = SNR,
                       y = 0, yend = 0 - tick_h),
          colour      = "firebrick",
          alpha       = rug.alpha,
          linewidth   = 0.3,
          inherit.aes = FALSE
        )
    }
  }

  # ------------------------------------------------------------------
  # Optional mirrored distributions
  # Detected sits above y=1, missed sits below y=0.
  # Both share a common normalisation so relative counts are honest.
  # ------------------------------------------------------------------

  if (distribution == "density" &&
      length(det_snr) > 1 && length(miss_snr) > 1) {

    dens_det  <- density(det_snr,  na.rm = TRUE)
    dens_miss <- density(miss_snr, na.rm = TRUE)

    global_max <- max(dens_det$y, dens_miss$y)
    scale_fun  <- function(y) y / global_max * mirror.height

    det_df  <- data.frame(SNR = dens_det$x,
                           ymin = 1,
                           ymax = 1 + scale_fun(dens_det$y))
    miss_df <- data.frame(SNR = dens_miss$x,
                           ymin = 0 - scale_fun(dens_miss$y),
                           ymax = 0)

    p <- p +
      ggplot2::geom_ribbon(
        data        = det_df,
        ggplot2::aes(x = SNR, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill        = "steelblue",
        alpha       = 0.35
      ) +
      ggplot2::geom_ribbon(
        data        = miss_df,
        ggplot2::aes(x = SNR, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill        = "firebrick",
        alpha       = 0.35
      )
  }

  if (distribution == "histogram") {

    breaks <- pretty(obs$SNR, n = 25)
    h_det  <- hist(det_snr,  breaks = breaks, plot = FALSE)
    h_miss <- hist(miss_snr, breaks = breaks, plot = FALSE)

    global_max <- max(h_det$counts, h_miss$counts)
    scale_fun  <- function(counts) counts / global_max * mirror.height

    bw <- diff(h_det$mids)[1]

    hist_df <- data.frame(
      SNR       = h_det$mids,
      det_hi    = 1 + scale_fun(h_det$counts),
      miss_lo   = 0 - scale_fun(h_miss$counts)
    )

    p <- p +
      # detected: bars from y=1 upward
      ggplot2::geom_rect(
        data        = hist_df,
        ggplot2::aes(xmin = SNR - bw * 0.45, xmax = SNR + bw * 0.45,
                     ymin = 1, ymax = det_hi),
        fill        = "steelblue",
        alpha       = 0.4,
        inherit.aes = FALSE
      ) +
      # missed: bars from y=0 downward
      ggplot2::geom_rect(
        data        = hist_df,
        ggplot2::aes(xmin = SNR - bw * 0.45, xmax = SNR + bw * 0.45,
                     ymin = miss_lo, ymax = 0),
        fill        = "firebrick",
        alpha       = 0.4,
        inherit.aes = FALSE
      )
  }

  # ------------------------------------------------------------------
  # Optional subtitle: detection count
  # ------------------------------------------------------------------

  if (show.counts) {
    p <- p + ggplot2::labs(
      subtitle = sprintf("%d / %d events detected", n_det, n_total)
    )
  }

  p
}

#' @export
showDetFun.list <- function(
    model,
    newdata   = NULL,
    ci        = FALSE,
    npoints   = 300,
    linewidth = 1.1,
    palette   = NULL,
    ...
) {

  pred <- predictDetFunList(
    model,
    newdata = newdata,
    ci      = ci,
    npoints = npoints,
    ...
  )

  p <- ggplot2::ggplot(
    pred,
    ggplot2::aes(SNR, fit, colour = model)
  ) +
    ggplot2::geom_line(linewidth = linewidth) +
    ggplot2::labs(
      x      = "SNR (dB)",
      y      = "Detection probability",
      colour = NULL
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_bw()

  if (ci && all(c("lower", "upper") %in% names(pred))) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower, ymax = upper, fill = model),
        alpha  = 0.15,
        colour = NA
      )
  }

  if (!is.null(palette)) {
    p <- p +
      ggplot2::scale_colour_manual(values = palette) +
      ggplot2::scale_fill_manual(values = palette)
  }

  p
}
