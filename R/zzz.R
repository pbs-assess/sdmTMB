.onAttach <- function(libname, pkgname) {
  if (utils::packageVersion("Matrix") >= package_version('1.3.0')) {
    packageStartupMessage(
      "There is currently a bug involving versions of the Matrix package >= 1.3-0 \n",
      "and the version of TMB on CRAN that results in some TMB (and sdmTMB) models \n",
      "eventually failing to fit with errors about: \n",
      "`par.random : replacement has length zero`. \n",
      "If you encounter this, you will want to install the latest versions from source:\n",
      "\n",
      "install.packages('Matrix', type = 'source')\n",
      "remotes::install_github('kaskr/adcomp/TMB')\n",
      "remotes::install_github('pbs-assess/sdmTMB')"
    )
  }
}
