# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do your testing?
# - If you are developing a package, run your tests using covr::report(),
#   devtools::test(), or by pressing Ctrl/Cmd + Shift + T.
# - If you are testing a script, run your tests using testthat::test_file()

library(testthat)
library(MAIHDA)

test_check("MAIHDA")
