# Builds data/mep_io.rda from the public MEP_IOcurve repository
# (https://github.com/ShayOfir/MEP_IOcurve, data/cleaned_MEP_data.csv;
# Wonderman Bar Sela et al. 2026, PLOS ONE 21(2):e0343725).
# Usage: Rscript data-raw/mep_io.R <path to cleaned_MEP_data.csv>
# Same exclusions and intensity scaling as load_and_preprocess_MEP_data() there.
emg <- read.csv(commandArgs(TRUE)[1])
d <- subset(emg, toUse == 1 & patient_code != 303 &
                 !(patient_code == 301 & coil == "fig8" & side == "rt"))
vmax <- function(coil) ifelse(coil == "fig8", 1670, 1700)  # stimulator output (V)
mep_io <- data.frame(
  # anonymous study codes 301..311 recoded to S01..S10 (sorted order kept)
  subj      = factor(sprintf("S%02d", match(d$patient_code, sort(unique(d$patient_code))))),
  coil      = factor(d$coil, levels = c("fig8", "h7", "rf")),
  side      = factor(d$side, levels = c("lt", "rt"), labels = c("left", "right")),
  intensity = vmax(d$coil)^2 * d$machine_setting / 1700^2,   # % of the RF device's output
  y         = d$PeakToPeak                                    # MEP peak-to-peak (mV)
)
usethis::use_data(mep_io, overwrite = TRUE, compress = "xz")
