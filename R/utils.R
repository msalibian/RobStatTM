#' @import fit.models

.onAttach <- function(libname, pkgname) {
  fit.models::fmclass.add.class("lmfm", "lmrobM", warn = F)
  fit.models::fmclass.add.class("lmfm", "lmrobdetMM", warn = F)
  fit.models::fmclass.add.class("lmfm", "lmrobdetDCML", warn = F)
  
  fit.models::fmclass.register(fmclass = "covfm",
                   classes = c("covClassic", "covRob"),
                   validation.function = NULL)
                                 
  fit.models::fmclass.register(fmclass = "pcompfm",
                   classes = c("prcomp", "prcompRob"),
                   validation.function = NULL)

  invisible()
}
