# Shared console palette for MAIHDA print methods.
#
# Built entirely on cli, so every colour degrades to plain text wherever ANSI is
# unsupported (knitr/vignettes, R CMD check, testthat, NO_COLOR) -- the styling adds
# nothing to captured/plain output. The hues match the plot palette so the console
# and the figures agree, and they are deliberately VALENCE-NEUTRAL: colour encodes
# emphasis/decidedness (look-here vs. a neutral conclusion vs. de-emphasis), never
# good-vs-bad. In particular nothing is coloured green/red, so a null ("negligible")
# result is not dressed up as a success and a finding is not dressed up as an alarm.
#
#   accent    #D55E00  the plot highlight colour: headline values, "relevant",
#                      flagged counts -- i.e. "look here"
#   secondary #0072B2  the plot data colour: a neutral firm conclusion ("negligible")
#   muted     grey     notes, captions, footers, fit diagnostics, undecided states
#   warn      yellow   genuine caveats only (singular/boundary fits)
#   bold               section titles / labels
maihda_palette <- function() {
  list(
    accent    = cli::make_ansi_style("#D55E00"),
    secondary = cli::make_ansi_style("#0072B2"),
    muted     = cli::col_silver,
    warn      = cli::col_yellow,
    bold      = cli::style_bold
  )
}
