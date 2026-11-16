# colour.R -- palettes and the coloured summary tables.
#
# Colour is off everywhere by default. When it is switched on (jamovi: the
# "Colour" checkbox; R: `colored = TRUE`, `highlight = TRUE`), fit indices,
# reliabilities and standardized loadings are placed on a ten-step ramp from
# "worst" to "best", as the 2024 cfawiz() did in the console. The position on
# the ramp comes from psw_color_index() (R/print.R): fixed bands per kind of
# statistic, e.g. CFI below .70 is step 1 and .99 or more is step 10. The
# bands are a reading aid, not cutoffs with any authority.
#
# jamovi tables can only mark a cell as "negative" (red), so the ramp is shown
# in a separate Html result ("Colour Summary") instead of in the tables.

# Ten colours each, worst first. The default has to survive the common forms
# of colour blindness; "redblue" is the original cfawiz() ramp.
psw_palettes <- list(
  purpleorange = c("#7F3B08", "#B35806", "#E08214", "#FDB863", "#FEE0B6",
                   "#D8DAEB", "#B2ABD2", "#8073AC", "#542788", "#2D004B"),
  cividis = c("#FFEA46", "#E4CF5B", "#C4B56C", "#A69D75", "#8A8779",
              "#707173", "#575C6D", "#39486B", "#00336F", "#00204D"),
  redblue = c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
              "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"),
  grey = c("#F2F2F2", "#E0E0E0", "#CCCCCC", "#B3B3B3", "#999999",
           "#808080", "#666666", "#4D4D4D", "#333333", "#1A1A1A"))

psw_ramp <- function(palette = NULL) {
  if (is.null(palette) || !palette %in% names(psw_palettes)) palette <- "purpleorange"
  psw_palettes[[palette]]
}

# Colour of one value of one kind ("cfi", "rmsea", "loading", ...); NA if the kind has no ramp.
psw_value_colour <- function(kind, value, palette = NULL) {
  idx <- if (is.na(value)) NA_integer_ else psw_color_index(kind, value)
  if (is.na(idx)) NA_character_ else psw_ramp(palette)[max(1L, min(idx, 10L))]
}

# Black or white text, whichever reads better on a background colour.
psw_text_on <- function(bg) {
  rgb <- grDevices::col2rgb(bg) / 255
  if (sum(c(0.299, 0.587, 0.114) * rgb) > 0.55) "#1A1A1A" else "#FFFFFF"
}

#' One coloured table as HTML (internal)
#' @param title Table title.
#' @param values Numeric matrix; row and column names are the labels.
#' @param kinds Character matrix of the same shape (or a vector recycled by
#'   row): the kind of statistic in each cell, see `psw_color_index()`.
#' @keywords internal
#' @noRd
psw_heat_table <- function(title, values, kinds, palette = NULL) {
  kinds <- matrix(kinds, nrow(values), ncol(values))
  cell <- "padding:3px 10px;text-align:right;"
  head <- paste0("<tr><th></th>", paste0("<th style='", cell, "font-weight:normal'>",
                                         psw_jmv_escape(colnames(values)), "</th>", collapse = ""), "</tr>")
  body <- vapply(seq_len(nrow(values)), function(i) {
    tds <- vapply(seq_len(ncol(values)), function(j) {
      v <- values[i, j]
      if (is.na(v)) return(sprintf("<td style='%s'></td>", cell))
      bg <- psw_value_colour(kinds[i, j], v, palette)
      style <- if (is.na(bg)) cell else sprintf("%sbackground:%s;color:%s;", cell, bg, psw_text_on(bg))
      sprintf("<td style='%s'>%s</td>", style, sub("^(-?)0\\.", "\\1.", sprintf("%.3f", v)))
    }, "")
    sprintf("<tr><td style='padding:3px 10px 3px 0'>%s</td>%s</tr>", psw_jmv_escape(rownames(values)[i]),
            paste(tds, collapse = ""))
  }, "")
  sprintf("<p style='margin:14px 0 4px 0'><b>%s</b></p><table style='border-collapse:collapse'>%s%s</table>",
          psw_jmv_escape(title), head, paste(body, collapse = ""))
}

psw_heat_legend <- function(palette = NULL) {
  sw <- paste0("<span style='display:inline-block;width:18px;height:12px;background:", psw_ramp(palette),
               "'></span>", collapse = "")
  paste0("<p style='margin-top:12px;font-size:smaller'>worse ", sw, " better. Ten fixed bands per statistic ",
         "(CFI and TLI from .70 to .99, RMSEA from .30 to .01, SRMR from .35 to .03, loadings and AVE from .1 to .9, ",
         "\u03C9 from .50 to .95). The bands are a reading aid, not cutoffs.</p>")
}

# Coloured summary of a set of models that share factors and items: one column
# per model (a single CFA, or one per estimator).
#   fit: data frame with one row per model (cfi, tli, rmsea, srmr); labels: its column names in the output
#   loadings / reliability: long data frames with a column `by` naming the model
psw_heat_summary <- function(fit, labels, loadings, reliability, by, palette = NULL) {
  wide <- function(long, rows, value) {
    key <- do.call(paste, c(long[rows], sep = ": "))
    m <- matrix(NA_real_, length(unique(key)), length(labels), dimnames = list(unique(key), labels))
    m[cbind(match(key, unique(key)), match(long[[by]], labels))] <- long[[value]]
    m
  }
  idx <- c(cfi = "CFI", tli = "TLI", rmsea = "RMSEA", srmr = "SRMR")
  fm <- t(as.matrix(fit[names(idx)]))
  dimnames(fm) <- list(unname(idx), labels)
  out <- psw_heat_table("Fit", fm, names(idx), palette)
  if (!is.null(loadings) && nrow(loadings))
    out <- paste0(out, psw_heat_table("Standardized loadings", wide(loadings, c("factor", "item"), "std"), "loading", palette))
  if (!is.null(reliability) && nrow(reliability))
    out <- paste0(out,
                  psw_heat_table("Composite reliability (\u03C9)", wide(reliability, "factor", "omega"), "omega", palette),
                  psw_heat_table("AVE", wide(reliability, "factor", "ave"), "ave", palette))
  paste0(out, psw_heat_legend(palette))
}
