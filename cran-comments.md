This version fixes several bugs and adds minor functionality. This version
also works around a pre-check FAIL due to the 'effects' package 4.2-3 that
was submitted to CRAN the same day---and therefore installed on the Windows
check server---dropping an export of Effects.default.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* Local M2 macOS install, R 4.5.1
* Intel macOS (on github-actions), R-release
* Ubuntu 24.04.2 (on github-actions), R-release
* Windows (on github-actions), R-release
* Windows (winbuilder), R-devel
* Windows (winbuilder), R-release
* Windows (winbuilder), R-oldrelease
