# Version 1.0.3:
- The Shiny User Interface capability has been removed, and correspondingly
the Shiny Interface to RobStatTM Vignette has also been removed. The
reason for doing so is that the Shiny capability resulted in too many package dependencies.
We anticipate that Greg Brownson, the creator of the RobSTatTM Shiny UI, will make 
that capability available in an independent package complement to RobStatTM. 
- The package fit.models has been removed as a dependency, i.e., it is no longer listed as Depends.
The reason for this is that including fit.models in RobStatTM created too many dependencies.  The
fit.model package stand-alone functionality works very well with RobStatTM.  Correspondingly, the
use of fit.models has been removed from the "Vignette for Command Line Use of RobStatTM.pdf"
document, and is now provided in the separate vignette  "fit.models using RobStatTM.pdf".
- In `src/lmrob.c` add USE_FC_LEN_T and use FCONE when calling BLAS and LAPACK Fortran functions 
(https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Fortran-character-strings)
- Re-ordered the elements in the object returned by `covRob` and `Multirobu`, renamed the 
argument `cor` to `corr`, and now the correlation matrix (if the argument `corr = TRUE`) is returned in 
the element `cor`. 
- Renamed the returned entry `weigths` to `wts` and the order of other entries
in the object returned by `covClassic`
- Removed comment on the help page of `lmrobdetMM.RFPE` referring to it being for internal use. 
This was not correct, the function can be used directly.
- `lmrobdetMM.RFPE` now includes the argument `bothVals`. If set to `TRUE` the function returns
two numbers, that when added together yield the estimated robust prediction errors. 
- Objects returned by `covClassic` and `covRob` now include an element `call` with the
matched called to the corresponding function. 
- The function `lmrobLinTest` has been renamed to `lmrobdetLinTest`
- Fixed a bug producing undesired behaviour when an exact fit (more than half of the 
data lying perfectly on a line) was detected. 
- The help page for lmrobdetMM was revised to describe all entries in the 
returned object. 


# Version 1.0.1:
- Fix an issue with exact fits if the M-scale estimate is (close to) 0.
- Includes help pages for all datasets.

