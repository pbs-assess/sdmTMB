.onAttach <- function(libname, pkgname) {
  if (utils::packageVersion("Matrix") >= package_version('1.3.0')) {
    packageStartupMessage(
      "There is currently a bug involving versions of the Matrix package >= 1.3-0 ",
      "that results in sdmTMB models eventually failing to fit with errors about ",
      "`par.random : replacement has length zero`. ",
      "If you encounter this, you will want to install an older version of the ",
      "Matrix package (1.2-18 is suggested). ",
      "Download the version from here: `https://cran.r-project.org/src/contrib/Archive/Matrix/` ",
      "and, from the necessary working directory, install with ",
      "`install.packages('Matrix_1.2-18.tar.gz', type = 'source', repos = NULL)`."
    )
  }
}
