# RobStatTM

This repository contains a development version of the companion package to the upcoming 2nd edition of
the book [Robust Statistics: Theory and Methods](https://www.wiley.com/go/maronna/robust), by Ricardo Maronna, Doug Martin, Victor Yohai and Matias Salibian-Barrera.

* The latest official version of the package is available on [CRAN](https://cran.r-project.org/package=RobStatTM). You should probably use [that version](https://cran.r-project.org/package=RobStatTM).
* To install the "development" version on GitHub use 
```
devtools::install_github("msalibian/RobStatTM")
```
* The scripts reproducing the examples and figures in the book can be found in the folder `inst/scripts`.

#### Bug reports

A bug is a reproducible problem that is caused by the code in the package.
Good bug reports are extremely helpful. Please follow the guidelines below to submit your
bug report.

##### Guidelines for bug reports:

We use [GitHub issues](https://guides.github.com/features/issues/) to track and solve potential bugs in our package. When submitting your bug report please:

1. Use our package's [GitHub issue search](https://github.com/msalibian/RobStatTM/issues) to check
whether your issue / bug has already been reported.

2. Check if your issue / bug has been fixed by trying to reproduce it using the latest version of the package.

3. Isolate the problem by creating a **minimal reproducible example** (see below)

4. Create an [issue](https://guides.github.com/features/issues/) for this repository. Refer to [this page](https://help.github.com/en/articles/creating-an-issue) for instructions on how to create a GitHub issue.

A good bug report should not require others to contact you to find more information. Please
try to be as detailed as possible in your report. What is your environment? What steps will
reproduce the issue? What outcome did you expect and what outcome did you get?

###### Example:

> A short and descriptive bug report title
>
> A summary of the issue and the OS environment in which it occurs.
> Include the steps required to reproduce the bug.
>
> 1. This is the first step
> 2. This is the second step
> 3. Further steps, etc.
>
> Any other information you want to share that is relevant to the issue being
> reported. This might include the lines of code that you have identified as
> causing the bug, and potential solutions.


##### Minimal reproducible examples

(This section is adapted from [Rob Hyndman's notes on minimal reproducible examples](https://robjhyndman.com/hyndsight/minimal-reproducible-examples/)).

A Minimal reproducible example (MRE) is intended to reproduce an error using the smallest amount of
code possible. To check that your MRE code is reproducible, try running it in a fresh R
session before you submit the issue. Using minimal reproducible examples
saves package developers time in wading through messy code that is not
relevant to the apparent bug.

A MRE should consist of a single R script file that can be run without error in a fresh R
session, and should contain the following three sections:


  * Packages to be loaded.
  * The shortest amount of code that reproduces the problem.
  * The output of `sessionInfo()` as a comment.

Please remove anything that is not necessary to reproduce the problem.

Try to use one of the built-in datasets if possible. If you need to include
some data, then use `dput()` so
the data can be included as part of the same text file. In
most cases, you do not need to include all of
your data, just a small subset that will allow the problem to be reproduced.

If you randomly generate some data, use `set.seed(somenumber)`.

Please spend time adding comments so we can understand your code quickly.
