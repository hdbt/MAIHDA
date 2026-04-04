## Test environments

* local: R 4.3.3
* GitHub Actions (ubuntu-latest): R-devel, R-release, R-oldrel-1
* GitHub Actions (windows-latest): R-release
* GitHub Actions (macos-latest): R-release
* win-builder: R-devel, R-release

## R CMD check results

0 errors | 0 warnings | 0 notes

## Resubmission Notes

This is a resubmission. In this version I have addressed the feedback from Konstanze Lauseker:

* Added formatted literature references (with DOIs) to the Description field of the DESCRIPTION file.
* Added missing `\value` tags to exported print methods (`print.maihda_model.Rd`, `print.maihda_strata.Rd`, `print.maihda_summary.Rd`, `print.pvc_result.Rd`).
* Replaced `\dontrun{}` with `\donttest{}` in examples, and completely unwrapped fast-running examples such as `make_strata()`.

## Downstream dependencies

There are currently no downstream dependencies for this package.
