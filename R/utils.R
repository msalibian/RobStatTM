#' @import fit.models

.onLoad <- function(libname, pkgname) {
  ##--------------- begin {fit.models} -----------------
  requireNamespace("fit.models")
  FM.add.class <- fit.models::fmclass.add.class
  FM.register  <- fit.models::fmclass.register

  FM.add.class("lmfm", "lmrobM", warn = F)
  FM.add.class("lmfm", "lmrobdetMM", warn = F)
  FM.add.class("lmfm", "lmrobdetDCML", warn = F)

  FM.add.class("covfm", "covClassic", warn = F)
  FM.add.class("covfm", "covRob", warn = F)
  
  FM.register(fmclass = "pcompfm",
              classes = c("prcomp", "prcompRob"),
              validation.function = NULL)
  ##--------------- end {fit.models} -------------------

  invisible()
}

#.onUnload <- function (libpath) {
#  library.dynam.unload("RobStatTM", libpath)
#}
