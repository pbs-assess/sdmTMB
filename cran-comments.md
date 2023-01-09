## Resubmission

This is a resubmission. In this version, we have fixed the following:

* Added \value to .Rd files for ggplot2_installed.Rd, inla_installed.Rd,
  plot_smooth.Rd
  
* Removed an internal function (R/make_enum.R) that wrote by default to 
  the user's home filespace.

> The LICENSE file is only needed if you have
> additional restrictions to the license which
> you have not? In that case omit the file and its
> reference in the DESCRIPTION file.

* We have moved the contents of LICENSE to inst/COPYRIGHTS, removed the
  reference to LICENSE in DESCRIPTION, and added 'Copyright: inst/COPYRIGHTS'
  in DESCRIPTION.

## R CMD check results

0 errors | 0 warnings | 3 notes

* This is a new release.

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
