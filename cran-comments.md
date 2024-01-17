This version fixes some minor bugs and is meant to trigger a rebuild of the
sdmTMB binaries on CRAN with the latest version of Matrix, which had breaking
ABI changes that affected TMB.
 
## R CMD check results

0 errors | 0 warnings | 1 notes

checking installed package size ... NOTE
    installed size is  6.1Mb
    sub-directories of 1Mb or more:
      libs   4.2Mb
      
* This is due to compiled code; we cannot reduce the size further.

## Test environments

* Local M2 macOS install, R 4.3.2
* Intel macOS (on github-actions), R-release
* Ubuntu 22.04.3 (on github-actions), R-release
* Ubuntu 22.04.3 (on github-actions), R-devel
* Ubuntu 22.04.3 (on github-actions), R-oldrel
* Windows (on github-actions), R 4.3.2
* Windows (winbuilder), R-devel

With sanitizer checks:
 
* Ubuntu 22.04.3 (on github-actions), R-devel with valgrind

