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
showDetFun <- function(model,
                       SNRinfo,
                       distribution = c("density", "histogram", "none"),
                       mirror.height = 0.20,
                       npoints = 300) {
  
  distribution <- match.arg(distribution)
  
  requireNamespace("ggplot2")
  
  # ------------------------------------------------------------
  # Prediction via unified engine
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
  
  p <- ggplot2::ggplot(pred,
                       ggplot2::aes(SNR, fit)) +
    ggplot2::geom_line(linewidth = 1.2, colour = "black") +
    ggplot2::labs(
      x = "SNR (dB)",
      y = "Detection probability"
    ) +
    ggplot2::coord_cartesian(
      ylim = c(-mirror.height * 1.1, 1 + mirror.height * 1.1)
    ) +
    ggplot2::theme_bw()
  
  # Add confidence ribbon if available
  if (all(c("lower", "upper") %in% names(pred))) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower, ymax = upper),
        alpha = 0.2,
        fill = "grey60"
      )
  }
  
  # ------------------------------------------------------------
  # Prepare observed data
  # ------------------------------------------------------------
  
  obs <- SNRinfo
  obs$Detected <- as.logical(obs$Detected)
  
  # ------------------------------------------------------------
  # Mirrored distributions
  # ------------------------------------------------------------
  
  if (distribution == "density") {
    
    p <- p +
      
      ggplot2::geom_density(
        data = subset(obs, Detected),
        ggplot2::aes(x = SNR,
                     y = after_stat(scaled) * mirror.height + 1),
        inherit.aes = FALSE,
        fill = "steelblue",
        alpha = 0.35,
        colour = NA
      ) +
      
      ggplot2::geom_density(
        data = subset(obs, !Detected),
        ggplot2::aes(x = SNR,
                     y = -after_stat(scaled) * mirror.height),
        inherit.aes = FALSE,
        fill = "firebrick",
        alpha = 0.35,
        colour = NA
      )
  }
  
  if (distribution == "histogram") {
    
    breaks <- pretty(obs$SNR, n = 25)
    
    h.det <- hist(obs$SNR[obs$Detected],
                  breaks = breaks,
                  plot = FALSE)
    
    h.miss <- hist(obs$SNR[!obs$Detected],
                   breaks = breaks,
                   plot = FALSE)
    
    hist.df <- data.frame(
      SNR = h.det$mids,
      det = h.det$density,
      miss = h.miss$density
    )
    
    hist.df$det  <- hist.df$det  / max(hist.df$det,  na.rm = TRUE) * mirror.height
    hist.df$miss <- hist.df$miss / max(hist.df$miss, na.rm = TRUE) * mirror.height
    
    bw <- diff(hist.df$SNR)[1]
    
    p <- p +
      
      ggplot2::geom_col(
        data = hist.df,
        ggplot2::aes(SNR, 1 + det),
        width = bw,
        fill = "steelblue",
        alpha = 0.4,
        inherit.aes = FALSE
      ) +
      
      ggplot2::geom_col(
        data = hist.df,
        ggplot2::aes(SNR, -miss),
        width = bw,
        fill = "firebrick",
        alpha = 0.4,
        inherit.aes = FALSE
      )
  }
  
  # ------------------------------------------------------------
  # Reference lines
  # ------------------------------------------------------------
  
  p +
    ggplot2::geom_hline(yintercept = 0,
                        linetype = 2,
                        colour = "grey60") +
    ggplot2::geom_hline(yintercept = 1,
                        linetype = 2,
                        colour = "grey60")
}