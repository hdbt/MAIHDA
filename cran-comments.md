## Test environments

* local: Windows 11 x64, R 4.5.1 (`R CMD check --as-cran`)
* win-builder: R-devel, R-release
* mac-builder: R-release
* GitHub Actions: ubuntu-latest (R-devel, R-release, R-oldrel-1), windows-latest (R-release), macos-latest (R-release)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking for future file timestamps ... NOTE
  unable to verify current time

  This note appears only in the local check, where the time server could not be reached; it is not about the package.

## Reverse dependencies

I checked the one reverse dependency, guideR 0.12.0, with `R CMD check` against this version: Status OK, no new problems.

## Changes that may affect users

This release removes functions deprecated in earlier releases (`calculate_pvc()` since 0.2.1; `plot_comparison()`, `plot_group_comparison()` and `decomposition = "cross-classified"` since 0.1.11). `"cross-classified"` now stops with an error naming its replacements, and reading the old `pvc` element of a PCV result warns.

The Gaussian PCV now uses each model's fitted (REML) between-stratum variance by default; `estimation = "ML"` restores the previous behaviour. All of these are listed in NEWS.md, and none is used by guideR.
