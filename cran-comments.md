## Submission

This version fixes the test failures CRAN reported on the r-devel flavours on 2026-11-15 (to be corrected before 2026-12-04).

From R-devel r90878, `model.matrix.default()` raises the error for a factor with a single level itself, where `contrasts<-` raised it before. `compare_maihda_groups()` recognised that failure by the call that raised it, so under R-devel its warning lost the sentence naming the single-level stratum dimension, and three expectations in `tests/testthat/test-audit-2026-10-12-compare-failure-wording.R` failed. The package now recognises both forms. No estimate was affected.

## Test environments

* local: Windows 11 x64, R 4.5.1 (`R CMD check --as-cran`)
* GitHub Actions: ubuntu-latest (R-devel r90926, R-release, R-oldrel-1), windows-latest (R-release), macos-latest (R-release)
* win-builder: R-devel, R-release r90926 

## R CMD check results

0 errors | 0 warnings | 1 note

* checking for future file timestamps ... NOTE
  unable to verify current time

  This note appears only in the local check, where the time server could not be reached; it is not about the package.

The tests that failed on CRAN's r-devel flavours pass on R-devel r90926 (GitHub Actions, ubuntu-latest).

## Reverse dependencies

I checked the one reverse dependency, guideR 0.12.0, with `R CMD check` against this version: Status OK, no new problems.

## Changes that may affect users

Besides the fix, this version adds intervals for the variance components to `summary()` and `tidy()` and changes which lme4 fit is kept, and which of lme4's convergence warnings are passed on, for fits that lme4's convergence checks flag. All of these are listed in NEWS.md, and none is used by guideR.
