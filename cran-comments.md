## Resubmission

This is a resubmission to fix 'Additional issues' on CRAN checks as requested.

We have removed calls from the rstan package, which appear to have caused
valgrind/gcc-ASAN/clang-ASAN errors on CRAN checks. We have tested the new
package with valgrind on two systems and clang-ASAN on the latest Fedora
attempting to match the CRAN setup and our checks pass.

Compared to the previously published version 0.2.1:

* We removed 2 Latin-1 strings in package data.
* We added 'future' and 'lme4' to Suggests to avoid 'Rd cross-references'
  'Undeclared packages' NOTE.
* We removed selected vignettes and reduced the size of an included dataset.
  The majority of the size is due to compiled code.
  
## R CMD check results

0 errors | 0 warnings | 2 notes

Suggests or Enhances not in mainstream repositories:
  INLA
Availability using Additional_repositories specification:
  INLA   yes   https://inla.r-inla-download.org/R/stable
...
Package unavailable to check Rd xrefs: 'INLA'
  
* This is correct and intended; INLA is not on CRAN.

checking installed package size ... NOTE
    installed size is  5.5Mb
    sub-directories of 1Mb or more:
      libs   3.7Mb
      
* This is due to compiled code; we cannot reduce the size further.

## Test environments

* local macOS install, R 4.2.2
* Windows (on github-actions), R 4.2.2
* Ubuntu 20.04.4 (on github-actions), R-devel
* Windows (winbuilder), R devel
* Windows Server 2022 (R-hub), R-devel, 64 bit

With sanitizer checks:
 
* Fedora release 38 (on docker), R-devel with clang ASAN
* Ubuntu 20.04.4 (on github-actions), R-devel with valgrind
* Ubuntu 20.04.4 (on docker), R-devel with valgrind level 2 instrumentation
