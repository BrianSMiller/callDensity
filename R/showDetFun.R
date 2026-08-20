# ============================================================
# Plot SNR detection function with observed SNR distributions
# ============================================================

#' Plot SNR detection function with observed SNR distributions
#'
#' @description
#' Visualises one or more fitted SNR detection functions together with the
#' distribution of detected and missed events. Dispatches on the class of
#' \code{model}: a single fitted model (glm/gam/scam/vglm/vgam) draws one
#' curve (\code{showDetFun.detFun}); a named list of fitted models overlays
#' all of them, each in its own colour (\code{showDetFun.list}).
#'
#' By default draws a rug along the y=1 (detected) and y=0 (missed)
#' reference lines. Optionally overlays density curves or histograms of the
#' detected/missed SNR distribution above and below the main panel, via
#' \pkg{patchwork}.
#'
#' @param model A single fitted detection model from \code{fitDetFun()}, or
#'   a *named* list of them (names become the legend/colour labels).
#' @param ... Passed to \code{showDetFun.detFun()} or \code{showDetFun.list()}.
#'
#' @param SNRinfo Data frame used to fit the model(s). Must contain columns
#'   \code{SNR} and \code{Detected} (logical or 0/1). Optional: if omitted
#'   (\code{NULL}, the default), it's recovered directly from each fitted
#'   model, since that's what it was fit on anyway. For the list form,
#'   this can also be a named list of data frames, one per model (matched
#'   by name) -- use this only if you want to show *different* data than
#'   what a model was actually fit on. A single shared data.frame is used
#'   for every model otherwise.
#'
#' @param distribution One of \code{"none"} (default), \code{"density"}, or
#'   \code{"histogram"}. Adds detected/missed SNR distributions above/below
#'   the main panel. For the list form, one histogram/density per model is
#'   overlaid, coloured to match its curve.
#'
#' @param rug Logical. If \code{TRUE} (default), draws tick marks along the
#'   y=1 line for detected events and along y=0 for missed events.
#'   For large datasets, a random subsample of \code{rug.max} points is used.
#'
#' @param rug.max Maximum number of rug ticks per group. Default 500.
#'   Set to \code{Inf} to show all points (slow for n > 2000).
#'
#' @param rug.alpha Alpha for rug ticks. Default 0.4.
#'
#' @param panel.heights Relative heights of the top/main/bottom panels when
#'   \code{distribution != "none"}. Default \code{c(1, 4, 1)}.
#'
#' @param hist.alpha Fill alpha for histogram/density panels. Default 0.4
#'   for a single model; for a list of models, default \code{1/length(model)}
#'   so overlapping fills stay legible.
#'
#' @param palette Optional named vector of colours, one per model (list form
#'   only). Default \code{NULL} uses ggplot2's default discrete palette,
#'   applied identically across all three panels so a model's line,
#'   histogram, and rug always match -- even if that model has, say, zero
#'   missed events and would otherwise be silently dropped from just the
#'   bottom panel's own colour scale.
#'
#' @param npoints Number of SNR grid points for the fitted curve(s). Default 300.
#'
#' @param xlim Optional length-2 numeric vector restricting the SNR range
#'   shown. Applied consistently to every panel (and to the fitted curve's
#'   own prediction grid) at construction time, so the top/main/bottom
#'   panels can't drift out of sync the way they would if a scale were
#'   changed on just one panel after the fact -- see Details.
#'
#' @param show.counts Logical, default \code{TRUE}. For a single model, adds
#'   a subtitle reporting the number of detections out of total
#'   observations. For a list of models, a subtitle doesn't scale past one
#'   model, so instead each legend entry gets its own \code{(detections/events)}
#'   count appended, with the legend titled accordingly -- e.g. "GLM (62/800)"
#'   under a "Model (# detections / # events)" legend title.
#'
#' @return A ggplot object, or a patchwork object when \code{distribution != "none"}.
#'
#' @details
#' When \code{distribution != "none"}, the return value is a
#' \pkg{patchwork} object made of three independently-built ggplot panels
#' (top/main/bottom) sharing one x-range fixed at construction time --
#' patchwork doesn't keep panel axes "linked" the way, say, MATLAB's
#' \code{linkaxes} does; each panel is just a static plot glued to its
#' neighbours at render time. So rescaling only one panel after the fact
#' (e.g. \code{p[[2]] <- p[[2]] + xlim(...)}) only affects that one panel,
#' and the others silently keep their original range. Pass \code{xlim}
#' here instead, so every panel is built from the same range to begin with.
#'
#' Adding a plain ggplot layer (\code{geom_line()}, \code{theme()}, etc) to
#' that patchwork object with \code{+} applies it to whichever panel was
#' supplied \emph{last} to \code{patchwork::wrap_plots()} -- which is
#' deliberately the main panel here (supplied after top and bottom, then
#' placed back in the visual middle via a \code{design} string), since
#' that's almost always the panel people want to add curves/annotations to.
#' So \code{p + geom_line(...)} lands on the main panel directly -- no
#' indexing needed. To target the top or bottom panel specifically instead,
#' index it directly by its position in \code{wrap_plots()}'s argument list
#' (which is top, bottom, main -- \emph{not} their top-to-bottom visual
#' order): \code{p[[1]] <- p[[1]] + ...} for the top panel, or
#' \code{p[[2]] <- p[[2]] + ...} for the bottom one.
#'
#' \code{+} only ever reaches that one last-supplied panel, regardless of
#' what's being added -- there's no way to make some additions (like a
#' post-hoc \code{xlim()}) broadcast to all three panels while others stay
#' targeted at just the main one, using \code{+} for both. For anything
#' that should apply to every panel at once (an \code{xlim()} change
#' included), use patchwork's own broadcast operator \code{&} instead:
#' \code{p & ggplot2::xlim(-10, 10)}. That said, changing the range post hoc
#' this way re-clips already-built panels rather than rebuilding them, so
#' prefer passing \code{xlim} up front where possible.
#'
#' @export
showDetFun <- function(model, ...) {
  if (inherits(model, c("glm", "gam", "scam", "vglm", "vgam"))) {
    showDetFun.detFun(model, ...)
  } else {
    showDetFun.list(model, ...)
  }
}

# ------------------------------------------------------------------
# Internal helpers shared by both the single-model and list forms
# ------------------------------------------------------------------

#' Split an SNRinfo data.frame into detected/missed SNR vectors
#'
#' @param obs Data.frame with columns \code{SNR} and \code{Detected}.
#'
#' @return A list with elements \code{det} and \code{miss}.
#'
#' @keywords internal
detMissSNR <- function(obs) {
  detected <- as.logical(obs$Detected)
  list(
    det  = obs$SNR[detected  & !is.na(obs$SNR)],
    miss = obs$SNR[!detected & !is.na(obs$SNR)]
  )
}

#' Kernel density estimate pre-clipped to a shared x-range
#'
#' Passing \code{from}/\code{to} directly to \code{density()} keeps the
#' curve within the shared panel x-limits from the start, rather than
#' computing the usual wider KDE (which by default extends \code{cut}
#' bandwidths past the data) and letting ggplot2 silently drop and warn
#' about the clipped tails once \code{scale_x_continuous(limits=)} is
#' applied.
#'
#' @param x Numeric vector.
#' @param xlim Length-2 numeric vector, the shared panel x-range.
#'
#' @return A data.frame with columns \code{SNR} and \code{y}.
#'
#' @keywords internal
clippedDensity <- function(x, xlim) {
  d <- stats::density(x, na.rm = TRUE, from = xlim[1], to = xlim[2])
  data.frame(SNR = d$x, y = d$y)
}

#' Evenly spaced histogram breaks spanning exactly xlim
#'
#' Used instead of \code{pretty()} so bins never extend past the shared
#' panel x-limits -- \code{pretty()}'s round numbers routinely do, which
#' silently drops the outermost bin (and warns) once
#' \code{scale_x_continuous(limits=)} clips it.
#'
#' @param xlim Length-2 numeric vector, the shared panel x-range.
#' @param n Number of bins. Default 25.
#'
#' @return A numeric vector of \code{n + 1} break points.
#'
#' @keywords internal
evenBreaks <- function(xlim, n = 25) seq(xlim[1], xlim[2], length.out = n + 1)

#' Blank theme for the top/bottom distribution panels
#'
#' Keeps everything void (no ticks, text, or gridlines) except the y-axis
#' title, so "Detected"/"Missed" render as a normal rotated axis title
#' anchored to the same left-hand column as the main panel's own y-axis
#' title -- not off on the right edge, where its position would drift
#' with legend width (as it did with plot.tag positioning).
#'
#' @keywords internal
distTheme <- function() {
  ggplot2::theme_void() +
    ggplot2::theme(
      plot.margin  = ggplot2::margin(2, 5, 2, 5),
      axis.title.y = ggplot2::element_text(size = 9, face = "bold",
                                            colour = "grey30",
                                            angle = 90, vjust = 0.5)
    )
}

#' Subsample a vector of rug points for large datasets
#' @keywords internal
subsampleRug <- function(x, n) {
  if (is.finite(n) && length(x) > n) x[sample(length(x), n)] else x
}

#' Resolve SNRinfo into a named list matching modelNames
#'
#' Accepts a single shared data.frame (used for every model), a named list
#' of data.frames (one per model, for models fit to different underlying
#' data), or NULL (recover it from each model itself via
#' \code{extractSNRinfo()}) -- and always returns the list shape.
#'
#' @keywords internal
resolveSNRinfoList <- function(SNRinfo, model, modelNames) {
  if (is.null(SNRinfo)) {
    return(stats::setNames(lapply(model, extractSNRinfo), modelNames))
  }
  if (is.data.frame(SNRinfo)) {
    return(stats::setNames(rep(list(SNRinfo), length(modelNames)), modelNames))
  }
  if (is.list(SNRinfo)) {
    missing <- setdiff(modelNames, names(SNRinfo))
    if (length(missing) > 0) {
      stop("showDetFun: SNRinfo list is missing entries for: ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    return(SNRinfo[modelNames])
  }
  stop("showDetFun: SNRinfo must be a data.frame or a named list of ",
       "data.frames, one per model.", call. = FALSE)
}

#' Recover an SNRinfo-shaped data.frame directly from a fitted model
#'
#' Every model type \code{fitDetFun()} produces is fit against SNR and a
#' detected/missed outcome, so both are still there to pull back out --
#' this is what makes the separate \code{SNRinfo} argument optional.
#'
#' glm/gam/scam are fit via a plain formula against a data.frame
#' (\code{Detected ~ SNR}), so \code{model.frame()} already returns exactly
#' those two columns -- the same mechanism \code{predictDetFunList()}
#' already relies on for its default SNR grid.
#'
#' vglm's capture-recapture formula uses a multi-observer response matrix
#' rather than a single \code{Detected} column, so it needs its own path:
#' SNR comes from the \code{@@x} design matrix and Detected from the
#' \code{@@y} response matrix, picking out \code{whichObserver}'s column --
#' the same slot \code{predictDetFun.vglm()} already uses for its own SNR
#' range. Two list entries can share one underlying vglm fit and differ
#' only in \code{@@extra$whichObserver} (as in the CommonGround vignette,
#' comparing observer 1 vs observer 2 on the same model) -- this reads
#' whichObserver fresh from each model object, so that distinction is
#' preserved correctly.
#'
#' @param model A fitted detFun object.
#'
#' @return A data.frame with columns \code{SNR} and \code{Detected}.
#'
#' @keywords internal
extractSNRinfo <- function(model) {
  if (isS4(model)) {
    if (!all(c("x", "y") %in% methods::slotNames(model))) {
      stop("showDetFun: don't know how to recover SNRinfo from an object ",
           "of class ", paste(class(model), collapse = "/"), " -- pass ",
           "SNRinfo explicitly.", call. = FALSE)
    }
    whichObserver <- model@extra$whichObserver
    if (is.null(whichObserver)) {
      stop("showDetFun: model@extra$whichObserver is missing, so the ",
           "detected/missed outcome can't be recovered automatically -- ",
           "pass SNRinfo explicitly.", call. = FALSE)
    }
    return(data.frame(
      SNR      = model@x[, "SNR"],
      Detected = as.logical(model@y[, whichObserver])
    ))
  }

  mf <- stats::model.frame(model)
  if (!all(c("SNR", "Detected") %in% names(mf))) {
    stop("showDetFun: model.frame(model) doesn't contain both 'SNR' and ",
         "'Detected' columns -- pass SNRinfo explicitly.", call. = FALSE)
  }
  data.frame(SNR = mf$SNR, Detected = as.logical(mf$Detected))
}

# ------------------------------------------------------------------
# Single fitted model
# ------------------------------------------------------------------

#' @export
showDetFun.detFun <- function(model,
                               SNRinfo       = NULL,
                               distribution  = c("none", "density", "histogram"),
                               rug           = TRUE,
                               rug.max       = 500,
                               rug.alpha     = 0.4,
                               panel.heights = c(1, 4, 1),
                               hist.alpha    = 0.4,
                               npoints       = 300,
                               xlim          = NULL,
                               show.counts   = TRUE) {

  distribution <- match.arg(distribution)
  requireNamespace("ggplot2")
  requireNamespace("patchwork")

  if (is.null(SNRinfo)) SNRinfo <- extractSNRinfo(model)

  # ------------------------------------------------------------------
  # Data preparation
  # ------------------------------------------------------------------
  snr      <- detMissSNR(SNRinfo)
  n_det    <- length(snr$det)
  n_total  <- sum(!is.na(SNRinfo$SNR))
  x_limits <- if (!is.null(xlim)) xlim else range(SNRinfo$SNR, na.rm = TRUE)

  # ------------------------------------------------------------------
  # Fitted curve (main middle panel)
  # ------------------------------------------------------------------
  pred <- predictDetFun(
    model   = model,
    newdata = data.frame(SNR = seq(x_limits[1], x_limits[2], length.out = npoints)),
    ci = TRUE
  )

  p_main <- ggplot2::ggplot(pred, ggplot2::aes(SNR, fit)) +
    ggplot2::geom_hline(yintercept = c(0, 1), linetype = 2, colour = "grey60") +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_x_continuous(limits = x_limits, expand = c(0.02, 0)) +
    ggplot2::labs(x = "SNR (dB)", y = "Detection probability") +
    ggplot2::theme_bw()

  if (all(c("lower", "upper") %in% names(pred))) {
    p_main <- p_main +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                           alpha = 0.2, fill = "grey60", inherit.aes = TRUE)
  }
  p_main <- p_main + ggplot2::geom_line(linewidth = 1.2, colour = "black")

  # ------------------------------------------------------------------
  # Rug: subsampled ticks along y=1 (detected) and y=0 (missed)
  # ------------------------------------------------------------------
  if (rug) {
    rug_det  <- subsampleRug(snr$det,  rug.max)
    rug_miss <- subsampleRug(snr$miss, rug.max)

    if (length(rug_det) > 0) {
      p_main <- p_main + ggplot2::geom_rug(
        data = data.frame(SNR = rug_det), ggplot2::aes(x = SNR), sides = "t",
        colour = "steelblue", alpha = rug.alpha, inherit.aes = FALSE)
    }
    if (length(rug_miss) > 0) {
      p_main <- p_main + ggplot2::geom_rug(
        data = data.frame(SNR = rug_miss), ggplot2::aes(x = SNR), sides = "b",
        colour = "firebrick", alpha = rug.alpha, inherit.aes = FALSE)
    }
  }

  # ------------------------------------------------------------------
  # Optional top/bottom distribution panels
  # ------------------------------------------------------------------
  p_top    <- ggplot2::ggplot() + ggplot2::theme_void()
  p_bottom <- ggplot2::ggplot() + ggplot2::theme_void()

  if (distribution == "density" && length(snr$det) > 1 && length(snr$miss) > 1) {
    det_df  <- clippedDensity(snr$det,  x_limits)
    miss_df <- clippedDensity(snr$miss, x_limits)

    p_top <- ggplot2::ggplot(det_df, ggplot2::aes(x = SNR, y = y)) +
      ggplot2::geom_area(fill = "steelblue", alpha = hist.alpha) +
      ggplot2::scale_x_continuous(limits = x_limits, expand = c(0.02, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::labs(y = "Detected") +
      distTheme()

    p_bottom <- ggplot2::ggplot(miss_df, ggplot2::aes(x = SNR, y = y)) +
      ggplot2::geom_area(fill = "firebrick", alpha = hist.alpha) +
      ggplot2::scale_x_continuous(limits = x_limits, expand = c(0.02, 0)) +
      ggplot2::scale_y_continuous(trans = "reverse", expand = c(0, 0)) +
      ggplot2::labs(y = "Missed") +
      distTheme()
  }

  if (distribution == "histogram") {
    breaks <- evenBreaks(x_limits)

    p_top <- ggplot2::ggplot(data.frame(SNR = snr$det), ggplot2::aes(x = SNR)) +
      ggplot2::geom_histogram(breaks = breaks, fill = "steelblue", alpha = hist.alpha,
                              colour = "white", linewidth = 0.2) +
      ggplot2::scale_x_continuous(limits = x_limits, expand = c(0.02, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::labs(y = "Detected") +
      distTheme()

    p_bottom <- ggplot2::ggplot(data.frame(SNR = snr$miss), ggplot2::aes(x = SNR)) +
      ggplot2::geom_histogram(breaks = breaks, fill = "firebrick", alpha = hist.alpha,
                              colour = "white", linewidth = 0.2) +
      ggplot2::scale_x_continuous(limits = x_limits, expand = c(0.02, 0)) +
      ggplot2::scale_y_continuous(trans = "reverse", expand = c(0, 0)) +
      ggplot2::labs(y = "Missed") +
      distTheme()
  }

  # ------------------------------------------------------------------
  # Subtitle: detection count
  # ------------------------------------------------------------------
  if (show.counts) {
    lbl <- sprintf("%d / %d events detected", n_det, n_total)
    if (distribution != "none") {
      p_top <- p_top + ggplot2::labs(subtitle = lbl)
    } else {
      p_main <- p_main + ggplot2::labs(subtitle = lbl)
    }
  }

  # ------------------------------------------------------------------
  # Panel composition
  # ------------------------------------------------------------------
  if (distribution == "none") return(p_main)

  # p_main is supplied LAST (with p_bottom before it) so that a later
  # `p + geom_line(...)`-style addition -- patchwork's "+" always targets
  # whichever panel was supplied last -- lands on the main panel, which is
  # what people almost always want to tweak. The design string places
  # panels visually as top/main/bottom regardless of that supply order:
  # tags are assigned in supply order (p_top=A, p_bottom=B, p_main=C), and
  # "A\nC\nB" reads those tags out top-to-bottom as top/main/bottom.
  # `heights` follows this visual row order, not supply order, so
  # panel.heights still means (top, main, bottom) as documented.
  patchwork::wrap_plots(p_top, p_bottom, p_main, design = "A\nC\nB",
                        heights = panel.heights) &
    ggplot2::theme(plot.margin = ggplot2::margin(2, 5, 2, 5))
}

# ------------------------------------------------------------------
# Named list of fitted models
# ------------------------------------------------------------------

#' @export
showDetFun.list <- function(
    model,
    SNRinfo       = NULL,
    newdata       = NULL,
    distribution  = c("none", "density", "histogram"),
    ci            = FALSE,
    rug           = TRUE,
    rug.max       = 500,
    rug.alpha     = 0.4,
    panel.heights = c(1, 4, 1),
    hist.alpha    = NULL,
    npoints       = 300,
    linewidth     = 1.1,
    palette       = NULL,
    xlim          = NULL,
    show.counts   = TRUE,
    ...
) {
  distribution <- match.arg(distribution)
  requireNamespace("ggplot2")

  if (is.null(names(model)) || any(!nzchar(names(model)))) {
    stop("showDetFun: 'model' must be a *named* list when passing multiple ",
         "detection functions -- the names become the legend/colour ",
         "labels.", call. = FALSE)
  }
  modelNames <- names(model)

  if (is.null(newdata) && !is.null(xlim)) {
    newdata <- data.frame(SNR = seq(xlim[1], xlim[2], length.out = npoints))
  }

  # Fixed palette applied identically to every panel (line, ribbon,
  # histogram/density fill, rug), so a model's colour can't drift between
  # panels -- e.g. a model with zero missed events would otherwise be
  # silently dropped from just the bottom panel's own default colour scale.
  if (is.null(palette)) {
    palette <- stats::setNames(scales::hue_pal()(length(modelNames)), modelNames)
  } else if (is.null(names(palette))) {
    palette <- stats::setNames(palette, modelNames)
  }

  if (is.null(hist.alpha)) hist.alpha <- 1 / length(modelNames)

  # ------------------------------------------------------------------
  # Fitted curves (main middle panel)
  # ------------------------------------------------------------------
  pred <- predictDetFunList(model, newdata = newdata, ci = ci, npoints = npoints, ...)
  pred$model <- factor(pred$model, levels = modelNames)

  # One shared x-range for every panel, computed once, so they can't
  # diverge -- an explicit xlim wins, then whatever newdata spanned
  # (user-supplied or built from xlim above), falling back to the fitted
  # curves' own range.
  x_limits <- if (!is.null(xlim)) {
    xlim
  } else if (!is.null(newdata)) {
    range(newdata$SNR)
  } else {
    range(pred$SNR, na.rm = TRUE)
  }

  p_main <- ggplot2::ggplot(pred, ggplot2::aes(SNR, fit, colour = model)) +
    ggplot2::geom_line(linewidth = linewidth) +
    ggplot2::scale_x_continuous(limits = x_limits, expand = c(0.02, 0)) +
    ggplot2::labs(x = "SNR (dB)", y = "Detection probability") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_bw()

  if (ci && all(c("lower", "upper") %in% names(pred))) {
    p_main <- p_main +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = model),
                           alpha = 0.15, colour = NA) +
      ggplot2::scale_fill_manual(values = palette, limits = modelNames, guide = "none")
  }

  # SNR/Detected data (needed for rug, the distribution panels, and now
  # also the legend counts) is only worth resolving -- and only worth an
  # extractSNRinfo() call per model, when SNRinfo wasn't supplied -- if
  # something actually asked for it.
  needsSNR <- show.counts || distribution != "none" || rug

  legendLabels <- stats::setNames(modelNames, modelNames)
  legendName   <- NULL
  snrList        <- NULL
  detMissByModel <- NULL

  if (needsSNR) {
    # --------------------------------------------------------------
    # Per-model detected/missed SNR
    # --------------------------------------------------------------
    snrList <- resolveSNRinfoList(SNRinfo, model, modelNames)
    detMissByModel <- stats::setNames(
      lapply(modelNames, function(nm) detMissSNR(snrList[[nm]])),
      modelNames
    )

    if (show.counts) {
      counts <- vapply(detMissByModel, function(s) length(s$det), integer(1))
      totals <- vapply(snrList[modelNames], function(d) sum(!is.na(d$SNR)), integer(1))
      legendLabels <- stats::setNames(sprintf("%s (%d/%d)", modelNames, counts, totals),
                                      modelNames)
      legendName <- "Model (# detections / # events)"
    }
  }

  p_main <- p_main +
    ggplot2::scale_colour_manual(values = palette, limits = modelNames,
                                 labels = legendLabels, name = legendName)

  if (!needsSNR) return(p_main)

  # ------------------------------------------------------------------
  # Rug: one coloured tick set per model, along y=1 (top) / y=0 (bottom)
  # ------------------------------------------------------------------
  if (rug) {
    for (nm in modelNames) {
      rug_det  <- subsampleRug(detMissByModel[[nm]]$det,  rug.max)
      rug_miss <- subsampleRug(detMissByModel[[nm]]$miss, rug.max)

      if (length(rug_det) > 0) {
        p_main <- p_main + ggplot2::geom_rug(
          data = data.frame(SNR = rug_det), ggplot2::aes(x = SNR), sides = "t",
          colour = palette[[nm]], alpha = rug.alpha, inherit.aes = FALSE)
      }
      if (length(rug_miss) > 0) {
        p_main <- p_main + ggplot2::geom_rug(
          data = data.frame(SNR = rug_miss), ggplot2::aes(x = SNR), sides = "b",
          colour = palette[[nm]], alpha = rug.alpha, inherit.aes = FALSE)
      }
    }
  }

  # ------------------------------------------------------------------
  # Optional top/bottom distribution panels: one histogram/density per
  # model, overlaid (position = "identity") and coloured by the same
  # palette as the curves, at hist.alpha so overlaps stay legible.
  # ------------------------------------------------------------------
  p_top    <- ggplot2::ggplot() + ggplot2::theme_void()
  p_bottom <- ggplot2::ggplot() + ggplot2::theme_void()

  stackByModel <- function(extract, minLen = 0) {
    rows <- lapply(modelNames, function(nm) {
      x <- extract(detMissByModel[[nm]])
      if (length(x) <= minLen) return(NULL)
      data.frame(SNR = x, model = nm)
    })
    out <- do.call(rbind, rows)
    if (!is.null(out)) out$model <- factor(out$model, levels = modelNames)
    out
  }

  if (distribution == "density") {
    detDF  <- stackByModel(function(s) s$det,  minLen = 1)
    missDF <- stackByModel(function(s) s$miss, minLen = 1)

    densify <- function(df) {
      if (is.null(df)) return(NULL)
      do.call(rbind, lapply(split(df, df$model), function(g) {
        cbind(clippedDensity(g$SNR, x_limits), model = g$model[1])
      }))
    }
    detDF  <- densify(detDF)
    missDF <- densify(missDF)

    if (!is.null(detDF)) {
      p_top <- ggplot2::ggplot(detDF, ggplot2::aes(x = SNR, y = y, fill = model)) +
        ggplot2::geom_area(position = "identity", alpha = hist.alpha, colour = NA) +
        ggplot2::scale_fill_manual(values = palette, limits = modelNames, guide = "none") +
        ggplot2::scale_x_continuous(limits = x_limits, expand = c(0.02, 0)) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::labs(y = "Detected") +
        distTheme()
    }
    if (!is.null(missDF)) {
      p_bottom <- ggplot2::ggplot(missDF, ggplot2::aes(x = SNR, y = y, fill = model)) +
        ggplot2::geom_area(position = "identity", alpha = hist.alpha, colour = NA) +
        ggplot2::scale_fill_manual(values = palette, limits = modelNames, guide = "none") +
        ggplot2::scale_x_continuous(limits = x_limits, expand = c(0.02, 0)) +
        ggplot2::scale_y_continuous(trans = "reverse", expand = c(0, 0)) +
        ggplot2::labs(y = "Missed") +
        distTheme()
    }
  }

  if (distribution == "histogram") {
    breaks <- evenBreaks(x_limits)
    detDF  <- stackByModel(function(s) s$det)
    missDF <- stackByModel(function(s) s$miss)

    if (!is.null(detDF)) {
      p_top <- ggplot2::ggplot(detDF, ggplot2::aes(x = SNR, fill = model)) +
        ggplot2::geom_histogram(breaks = breaks, position = "identity",
                                alpha = hist.alpha, colour = "white", linewidth = 0.2) +
        ggplot2::scale_fill_manual(values = palette, limits = modelNames, guide = "none") +
        ggplot2::scale_x_continuous(limits = x_limits, expand = c(0.02, 0)) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::labs(y = "Detected") +
        distTheme()
    }
    if (!is.null(missDF)) {
      p_bottom <- ggplot2::ggplot(missDF, ggplot2::aes(x = SNR, fill = model)) +
        ggplot2::geom_histogram(breaks = breaks, position = "identity",
                                alpha = hist.alpha, colour = "white", linewidth = 0.2) +
        ggplot2::scale_fill_manual(values = palette, limits = modelNames, guide = "none") +
        ggplot2::scale_x_continuous(limits = x_limits, expand = c(0.02, 0)) +
        ggplot2::scale_y_continuous(trans = "reverse", expand = c(0, 0)) +
        ggplot2::labs(y = "Missed") +
        distTheme()
    }
  }

  if (distribution == "none") return(p_main)

  # See showDetFun.detFun() above for why p_main is supplied last here.
  patchwork::wrap_plots(p_top, p_bottom, p_main, design = "A\nC\nB",
                        heights = panel.heights) &
    ggplot2::theme(plot.margin = ggplot2::margin(2, 5, 2, 5))
}
