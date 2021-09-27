.onUnload <- function(libpath) {
  library.dynam.unload("sdmTMB", libpath)
}
