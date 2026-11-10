# Re-run the precomputed vignettes. They need CmdStan, so CRAN-style builds and
# CI only render the stored output (vignettes/*.Rmd); the sources with live
# code are vignettes/*.Rmd.orig. Install this version of the package first,
# then from the package root:  Rscript data-raw/vignettes.R [name ...]
# Figures are written next to the vignettes as <name>-<chunk>-<n>.png.
args <- commandArgs(TRUE)
setwd("vignettes")
src <- if (length(args)) paste0(args, ".Rmd.orig") else list.files(pattern = "\\.Rmd\\.orig$")
for (f in src) {
  t0 <- Sys.time()
  knitr::knit(f, sub("\\.orig$", "", f), envir = new.env())
  message(f, ": ", format(round(Sys.time() - t0, 1)))
}
