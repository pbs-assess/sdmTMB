This version fixes several issues that have arisen on CRAN as per an email
from Brian Ripley on 2023-10-07 and a note about rgdal on 2023-10-18.
 
* Removes dependencies on 'rgdal'
* Fixes issues caused by updates to 'sf', 'INLA', and 'glmmTMB'

This version no longer depends on 'INLA'. It also contains several bug fixes
and additions to functionality.
  
## R CMD check results

0 errors | 0 warnings | 1 notes

checking installed package size ... NOTE
    installed size is  6.0Mb
    sub-directories of 1Mb or more:
      libs   4.1Mb
      
* This is due to compiled code; we cannot reduce the size further.

## Test environments

* local macOS install, R 4.3.1
* Windows (on github-actions), R 4.3.1
* Ubuntu 22.04.3 (on github-actions), R-devel
* Windows (winbuilder), R-devel
* Windows (winbuilder), R-release

With sanitizer checks:
 
* Ubuntu 22.04.3 (on github-actions), R-devel with valgrind
* Ubuntu 22.04.3 (on github-actions), R-devel with clang-asan
