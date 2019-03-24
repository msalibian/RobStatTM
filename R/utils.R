.onAttach <- function(libname, pkgname) {
  fmclass.add.class("lmfm", "lmrobM", warn = F)
  fmclass.add.class("lmfm", "lmrobdetMM", warn = F)
  fmclass.add.class("lmfm", "lmrobdetDCML", warn = F)
  
  fmclass.register(fmclass = "covfm",
                   classes = c("covClassic", "covRob"),
                   validation.function = NULL)
                                 
  fmclass.register(fmclass = "pcompfm",
                   classes = c("prcomp", "pcaRobS"),
                   validation.function = NULL)

  invisible()
}