# ============================================================
# Plot SNR detection function with mirrored distributions
# ============================================================

#' Plot SNR detection function with observed SNR distributions
#'
#' @description
#' Visualises fitted SNR detection functions together with the distribution
#' of detected vs undetected events. Supports mirrored densities or histograms.
#'
#' @param model Fitted detection model from:
#'   fitSNRdetectionFunc(), fitSNRvglm()
#'
#' @param SNRinfo Data.frame used to fit the model. Must contain:
#'   - SNR
#'   - Detected (logical or 0/1)
#'
#' @param distribution One of "density", "histogram", "none"
#'
#' @param mirror.height Height scaling for mirrored distributions
#'
#' @param npoints Number of SNR grid points for prediction
#'
#' @return ggplot object
#'
#' @export
showDetFun <- function(model, ...) {
  UseMethod("showDetFun")
}

#' @export
showDetFun.detFun <- function(model,
                       SNRinfo,
                       distribution = c("density", "histogram", "none"),
                       mirror.height = 0.20,
                       npoints = 300) {

  distribution <- match.arg(distribution)
  requireNamespace("ggplot2")

  # ------------------------------------------------------------
  # Prediction
  # ------------------------------------------------------------

  pred <- predictDetFun(
    model = model,
    newdata = data.frame(
      SNR = seq(min(SNRinfo$SNR, na.rm = TRUE),
                max(SNRinfo$SNR, na.rm = TRUE),
                length.out = npoints)
    ),
    ci = TRUE
  )

  # ------------------------------------------------------------
  # Base plot
  # ------------------------------------------------------------

  p <- ggplot2::ggplot(pred, ggplot2::aes(SNR, fit)) +

    ggplot2::geom_line(linewidth = 1.2, colour = "black") +

    ggplot2::labs(
      x = "SNR (dB)",
      y = "Detection probability"
    ) +

    # FIXED AXIS: ONLY 0 and 1 visible scale
    ggplot2::scale_y_continuous(
      limits = c(-mirror.height, 1 + mirror.height),
      breaks = c(0, 1),
      labels = c("0", "1"),
      expand = c(0, 0)
    ) +

    ggplot2::theme_bw()

  # confidence band (still fine in this system)
  if (all(c("lower", "upper") %in% names(pred))) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower, ymax = upper),
        alpha = 0.2,
        fill = "grey60"
      )
  }

  # ------------------------------------------------------------
  # Observations
  # ------------------------------------------------------------

  obs <- SNRinfo
  obs$Detected <- as.logical(obs$Detected)

  # ------------------------------------------------------------
  # MIRRORED DENSITIES (FIXED ANCHORING)
  # ------------------------------------------------------------

  if (distribution == "density") {

    # ------------------------------------------------------------
    # GLOBAL NORMALISATION (key change)
    # ------------------------------------------------------------

    det_vals  <- obs$SNR[obs$Detected]
    miss_vals <- obs$SNR[!obs$Detected]

    dens_det <- density(det_vals, na.rm = TRUE)
    dens_miss <- density(miss_vals, na.rm = TRUE)

    global_max <- max(
      dens_det$y,
      dens_miss$y,
      na.rm = TRUE
    )

    scale_fun <- function(y) y / global_max * mirror.height

    # convert to data.frames for ggplot

    det_df <- data.frame(
      SNR = dens_det$x,
      y   = 1 + scale_fun(dens_det$y)
    )

    miss_df <- data.frame(
      SNR = dens_miss$x,
      y   = 0 - scale_fun(dens_miss$y)
    )

    p <- p +

      ggplot2::geom_area(
        data = det_df,
        ggplot2::aes(SNR, y),
        inherit.aes = FALSE,
        fill = "steelblue",
        alpha = 0.35
      ) +

      ggplot2::geom_area(
        data = miss_df,
        ggplot2::aes(SNR, y),
        inherit.aes = FALSE,
        fill = "firebrick",
        alpha = 0.35
      )
  }

  # ------------------------------------------------------------
  # HISTOGRAM VERSION (consistent anchoring)
  # ------------------------------------------------------------

  if (distribution == "histogram") {

    breaks <- pretty(obs$SNR, n = 25)

    h.det <- hist(obs$SNR[obs$Detected], breaks = breaks, plot = FALSE)
    h.miss <- hist(obs$SNR[!obs$Detected], breaks = breaks, plot = FALSE)

    hist.df <- data.frame(
      SNR = h.det$mids,
      det = h.det$density,
      miss = h.miss$density
    )

    hist.df$det  <- hist.df$det  / max(hist.df$det, na.rm = TRUE) * mirror.height
    hist.df$miss <- hist.df$miss / max(hist.df$miss, na.rm = TRUE) * mirror.height

    bw <- diff(hist.df$SNR)[1]

    p <- p +

      # DETECTED above 1
      ggplot2::geom_col(
        data = hist.df,
        ggplot2::aes(SNR, 1 + det),
        width = bw,
        fill = "steelblue",
        alpha = 0.4,
        inherit.aes = FALSE
      ) +

      # MISSED below 0
      ggplot2::geom_col(
        data = hist.df,
        ggplot2::aes(SNR, 0 - miss),
        width = bw,
        fill = "firebrick",
        alpha = 0.4,
        inherit.aes = FALSE
      )
  }

  # ------------------------------------------------------------
  # Reference lines (ONLY 0 and 1 axis structure)
  # ------------------------------------------------------------

  p +
    ggplot2::geom_hline(yintercept = c(0, 1),
                        linetype = 2,
                        colour = "grey60")
}

#' @export
showDetFun.list <- function(
    model,
    newdata = NULL,
    ci = FALSE,
    npoints = 300,
    linewidth = 1.1,
    palette = NULL,
    ...
) {

  pred <- predictDetFunList(
    model,
    newdata = newdata,
    ci = ci,
    npoints = npoints,
    ...
  )

  p <- ggplot2::ggplot(
    pred,
    ggplot2::aes(
      SNR,
      fit,
      colour = model
    )
  ) +
    ggplot2::geom_line(
      linewidth = linewidth
    ) +
    ggplot2::labs(
      x = "SNR (dB)",
      y = "Detection probability",
      colour = NULL
    ) +
    ggplot2::coord_cartesian(
      ylim = c(0, 1)
    ) +
    ggplot2::theme_bw()

  if (ci && all(c("lower", "upper") %in% names(pred))) {

    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = lower,
          ymax = upper,
          fill = model
        ),
        alpha = 0.15,
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
