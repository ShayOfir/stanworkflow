#' Motor evoked potential input-output curves
#'
#' MEP amplitudes recorded from the first dorsal interosseous muscle of 10
#' healthy adults stimulated with three TMS coils over both hemispheres, at
#' several intensities each (about five trials per intensity). Used in
#' `vignette("mep")`. Subject codes are anonymous; no demographic data are
#' included.
#'
#' @format A data frame with 1881 rows and 5 columns:
#' \describe{
#'   \item{subj}{subject, factor `S01`-`S10` (the study codes, recoded)}
#'   \item{coil}{coil, factor: `fig8` (Figure-of-8), `h7` (H7), `rf`
#'     (multidirectional rotational-field coil array)}
#'   \item{side}{stimulated hemisphere, factor `left`/`right`. Subject `S01`
#'     has no Figure-of-8 data for the right hemisphere.}
#'   \item{intensity}{stimulation intensity on a common scale, \% of maximal
#'     stimulator output (MSO): the machine setting for `h7` and `rf`; for
#'     `fig8`, whose stimulator has a lower maximal voltage (1670 vs 1700 V),
#'     the setting times (1670/1700)^2, as in the original analysis}
#'   \item{y}{MEP peak-to-peak amplitude (mV)}
#' }
#' @source Wonderman Bar Sela O, Ofir-Geva S, Pell GS, Roth Y, Friedman J,
#'   Muhana A, et al. (2026). MEPs elicited by multidirectional
#'   rotational-field TMS show marked differences compared to unidirectional
#'   Figure-of-8 and H7 coils. PLOS ONE 21(2): e0343725.
#'   \doi{10.1371/journal.pone.0343725}. Data and original model:
#'   <https://github.com/ShayOfir/MEP_IOcurve> (BSD-3-Clause), file
#'   `data/cleaned_MEP_data.csv`, with the exclusions of that analysis; built
#'   by `data-raw/mep_io.R`.
"mep_io"
