## Resubmission

This is a resubmission to fix 'Additional issues' on CRAN checks. This version:

* Corrects a memory issue identified by valgrind/gcc-ASAN/clan-ASAN tests
  
* Removes 2 Latin-1 strings in package data

* Adds 'future' and 'lme4' to Suggests to avoid 'Rd cross-references'
  'Undeclared packages' NOTE
  
> Result: NOTE
>   installed size is 8.9Mb
>   sub-directories of 1Mb or more:
>   data 1.7Mb
>   libs 5.4Mb

* We removed selected vignettes and reduced the size of an included dataset.
  The majority of the size is due to compiled code.

> Result: NOTE
>     Package suggested but not available for checking: ‘INLA’

* This is as intended; INLA is not on CRAN.
    
## R CMD check results

0 errors | 0 warnings | 2 notes

Possibly misspelled words in DESCRIPTION:
  GLMMs (3:46, 49:71)
  SDMs (53:43)
  SPDE (3:35, 50:75)
  Spatiotemporal (3:20)
  al (53:66)
  et (53:63)
  spatiotemporal (49:37)

* These are not misspelled.

Suggests or Enhances not in mainstream repositories:
  INLA
Availability using Additional_repositories specification:
  INLA   yes   https://inla.r-inla-download.org/R/stable
...
Package unavailable to check Rd xrefs: 'INLA'
  
* This is correct and intended; INLA is not on CRAN.

## Test environments

* local macOS install, R 4.2.2
* Windows (on github-actions), R 4.2.2
* Ubuntu 20.04.4 (on github-actions), R devel
* Windows (winbuilder), R devel
* Windows (winbuilder), R release
* Windows Server 2022 (R-hub), R-devel, 64 bit
