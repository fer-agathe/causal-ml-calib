library(ggplot2)
library(tikzDevice)

#' Theme for ggplot2
#'
#' @param ... arguments passed to the theme function
#' @export
#' @importFrom ggplot2 element_rect element_text element_blank element_line unit
#'   rel
theme_paper <- function (...) {
  ggthemes::theme_base() +
    theme(
      legend.background = element_rect(
        fill = "transparent", linetype="solid", colour ="black"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.key = element_blank(),
      panel.background = element_rect(fill = NA),
      strip.background = element_rect(fill = NA, color = NA)
    )
}


#' From plot created with {tikzDevice}, create a standalone latex document
#' and compile it with pdflatex to save the plot as pdf
#'
#' @param filename Name of the tex file (WITHOUT THE EXTENSION) that contains
#'  the tikzpicture.
#' @param path_to_latex Path to LaTeX engine (Defaults to
#'   `/Library/TeX/texbin/`).
#' @param interpreter By default, use pdflatex (`pdflatex`).
#' @param path Path to the destination folder.
#' @param keep_tex should the tex file (only the one from the standalone doc)
#'  be kept after compilation? Defaults to `FALSE`.
#' @param verbose A logical value indicating whether diagnostic messages are
#'   printed when measuring dimensions of strings. Defaults to `FALSE`.
#' @param ignore.stdout A logical (not NA) indicating whether messages written
#'   to ‘stdout’  should be ignored. Defaults to `TRUE`.
#' @param crop If `TRUE` (default to `FALSE`), the PDF is cropped using pdfcrop.
#'
plot_to_pdf <- function(filename,
                        path_to_latex = "/Library/TeX/texbin/",
                        interpreter = "pdflatex",
                        path = "./",
                        keep_tex = FALSE,
                        verbose = FALSE,
                        ignore.stdout = TRUE,
                        crop = FALSE) {
  content <- paste0(
    "\\documentclass{standalone}
      \\usepackage{amsmath,amssymb,amsthm,mathtools,graphicx}
      \\usepackage{array,dcolumn}
      \\usepackage{tikz}
      \\usetikzlibrary{arrows.meta, positioning, calc}
      %\\usepackage{dsfont}
      %\\usepackage{fontspec}
      \\renewcommand{\\familydefault}{\\rmdefault}
      %\\usepackage{natbib}
      \\usepackage{microtype}
      %\\usepackage{newtxtext,newtxmath}
      %\\usepackage{times,mathpazo}
      \\usepackage{pgfplots}
      \\usetikzlibrary{pgfplots.groupplots}
      \\usepackage{xcolor}
      \\begin{document}

      \\input{",
    path, filename,
    ".tex}

      \\end{document}"
  )

  # The file which will import the graph in tex format
  fileConn <- file(paste0(path, filename, "_tmp.tex"))
  writeLines(content, fileConn)
  close(fileConn)

  # Process tex file to get the PDF
  system(
    paste0(
      path_to_latex,
      interpreter, " -shell-escape -synctex=1 -interaction=nonstopmode  ",
      path,
      filename, "_tmp.tex"),
    ignore.stdout = TRUE
  )
  if (crop == TRUE) {
    system(
      paste0(
        "pdfcrop ", filename, "_tmp.pdf ", filename, "_tmp.pdf"
      )
    )
  }
  if(!path %in%  c(".", "./", "/"))
    system(paste0("mv ", filename, "_tmp.pdf ", path))
  system(paste0("rm ", filename, "_tmp.aux"))
  system(paste0("rm ", filename, "_tmp.log"))
  system(paste0("rm ", filename, "_tmp.synctex.gz"))
  if (!keep_tex) {
    system(paste0("rm ", path, filename, "_tmp.tex"))
  }
  system(paste0("mv ", path, filename, "_tmp.pdf ", path, filename, ".pdf"))
}

#' Save a ggplot2 plot as PDF, using LaTeX tikz
#'
#' @param plot A ggplot2 object.
#' @param path_to_latex Path to LaTeX engine (Defaults to
#'   `/Library/TeX/texbin/`).
#' @param interpreter By default, use pdflatex (`pdflatex`).
#' @param path Path to the destination folder.
#' @param filename File name (without the extension).
#' @param keep_tex should the tex file be kept after compilation? Defaults to
#'   `FALSE`.
#' @param width Width in inches (default to 15).
#' @param height Height in inches (default to 15).
#' @param verbose A logical value indicating whether diagnostic messages are
#'   printed when measuring dimensions of strings. Defaults to `FALSE`.
#' @param ignore.stdout A logical (not NA) indicating whether messages written
#'   to ‘stdout’  should be ignored. Defaults to `TRUE`.
#' @param crop If `TRUE` (default to `FALSE`), the PDF is cropped using pdfcrop.
#'
#' @importFrom tikzDevice tikz
#' @importFrom grDevices dev.off
#' @export
#' @md
#'
ggplot2_to_pdf <- function(plot,
                           path_to_latex = "/Library/TeX/texbin/",
                           interpreter = "pdflatex",
                           path = "./",
                           filename,
                           keep_tex = FALSE,
                           width = 15,
                           height = 15,
                           verbose = FALSE,
                           ignore.stdout = TRUE,
                           crop = FALSE) {

  content <- paste0(
    "\\documentclass{standalone}
\\usepackage[utf8]{inputenc}
\\usepackage[T1]{fontenc}
\\usepackage{amsmath,amssymb,amsthm,mathtools,graphicx}
\\usepackage{array,dcolumn}
%\\usepackage{dsfont}
%\\usepackage{fontspec}
%\\setmainfont{Noto Sans}
\\usepackage{tikz}
\\usetikzlibrary{arrows.meta,positioning,calc}
\\usepackage{nicefrac}
\\usepackage{microtype}
\\usepackage{pgfplots}
\\usetikzlibrary{pgfplots.groupplots}
\\usepackage{xcolor}
\\renewcommand{\\familydefault}{\\rmdefault}

\\begin{document}

\\input{", path, filename, "_content.tex}

\\end{document}"
  )

  ## Write wrapper TeX file
  fileConn <- file(paste0(path, filename, ".tex"))
  writeLines(content, fileConn)
  close(fileConn)

  ## Temporarily set the preamble used by tikzDevice
  old_decl <- getOption("tikzDocumentDeclaration")

  options(
    tikzDocumentDeclaration =
      "\\documentclass[11pt]{article}
\\usepackage[utf8]{inputenc}
\\usepackage[T1]{fontenc}
\\usepackage{tikz}
\\usepackage{amsmath,amssymb,amsthm,mathtools}
\\newcommand{\\E}{\\mathbb{E}}
\\newcommand{\\Cov}{\\mathrm{Cov}}
\\newcommand{\\one}{\\mathbf{1}}"
  )

  on.exit(
    options(tikzDocumentDeclaration = old_decl),
    add = TRUE
  )

  ## Export plot to TikZ
  tikz(
    file = paste0(path, filename, "_content.tex"),
    width = width,
    height = height,
    verbose = verbose
  )

  print(plot)
  dev.off()

  ## Move raster image produced by ggplot (if any)
  name_scale <- paste0(filename, "_content_ras1.png")
  scale_exists <- file.exists(name_scale)

  if (scale_exists && !path %in% c(".", "./", "/")) {
    system(paste0("mv ", name_scale, " ", path))
  }

  ## Compile PDF
  system(
    paste0(
      path_to_latex,
      interpreter,
      " -shell-escape -synctex=1 -interaction=nonstopmode ",
      path,
      filename,
      ".tex"
    ),
    ignore.stdout = ignore.stdout
  )

  if (crop) {
    system(paste0("pdfcrop ", path, filename, ".pdf ", path, filename, ".pdf"))
  }

  if (!path %in% c(".", "./", "/")) {
    system(paste0("mv ", filename, ".pdf ", path))
  }

  ## Cleanup
  system(paste0("rm ", filename, ".aux"))
  system(paste0("rm ", filename, ".log"))
  system(paste0("rm ", filename, ".synctex.gz"))

  if (!keep_tex) {
    system(paste0("rm ", path, filename, ".tex"))
    system(paste0("rm ", path, filename, "_content.tex"))
  }

  if (scale_exists) {
    system(paste0("rm ", path, "/", name_scale))
  }
}
