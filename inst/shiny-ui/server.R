# server.R
# 
# Author: Gregory Brownson
#
# Description: Backend portion of the Shiny 
#              interface for loading data and
#              selecting location and scale
#              parameters

options(warn = -1L, scipen = -1L)

pkgs <- c("DT", "fit.models", "ggplot2", "grid", "gridExtra", "gtools",
          "PerformanceAnalytics", "robust", "robustbase", "shiny", "xts")

missing.packages <- pkgs[sapply(pkgs, function(p) { length(find.package(p, quiet = T)) == 0 })]
if (length(missing.packages) > 0) {
  cat(paste("The following packages are missing:", missing.packages, ".\n"))
  cat("Installing missing packages!")
  install.packages(missing.packages)
}

lapply(pkgs, library, character.only = TRUE, quietly = T)

library(RobStatTM, quietly = T)

# Set theme for ggplot
thm <- theme_bw() +
       theme(plot.title = element_text(hjust = 0.5),
             axis.text = element_text(size = 12),
             axis.title = element_text(size = 14, face = "bold"))

theme_set(thm)

# Add models to respective class and create new fit.models class when necessary

if (is.null(fit.models:::e$fmreg$lmfm)) {
  fmclass.register("lmfm", c("lmrobM", "lmrobdetMM", "lmrobdetDCML"))
} else {
  fmclass.add.class("lmfm", "lmrobM", warn = F)
  fmclass.add.class("lmfm", "lmrobdetMM", warn = F)
  fmclass.add.class("lmfm", "lmrobdetDCML", warn = F)
}

if (is.null(fit.models:::e$fmreg$covfm)) {
  fmclass.register("covfm", c("covClassic", "covRob"))
} else {
  fmclass.add.class("covfm", "covClassic", warn = F)
  fmclass.add.class("covfm", "covRob", warn = F)
}

if (is.null(fit.models:::e$fmreg$pcompfm)) {
  fmclass.register("pcompfm", c("prcomp", "prcompRob"))
} else {
  fmclass.add.class("pcompfm", "prcomp", warn = F)
  fmclass.add.class("pcompfm", "prcompRob", warn = F)
}

# Functions which interface with rest of package and makes things easier for fit.models
locScaleClassic <- function(x) {
  z <- list()
  
  z$disper <- sd(x)
  z$mu <- mean(x)
  z$std.mu <- sd(x) / sqrt(length(x))
  
  z
}

lm <- function(form, ...) {
  fit <- stats::lm(form, ...)
  fit$call <- form
  return(fit)
}

lmrobM <- function(form, ...) {
  fit <- RobStatTM::lmrobM(form, ...)
  fit$call <- form
  return(fit)
}

lmrobdetDCML <- function(form, ...) {
  fit <- RobStatTM::lmrobdetDCML(form, ...)
  fit$call <- form
  return(fit)
}

lmrobdetMM <- function(form, ...) {
  fit <- RobStatTM::lmrobdetMM(form, ...)
  fit$call <- form
  return(fit)
}

covClassic <- function(data, data.name, corr = F, ...) {
  z <- RobStatTM::covClassic(data, ...)
  z$corr <- corr
  z$call <- call("covClassic", data = as.name(data.name))
  
  if(corr) {
    z$sdev <- sqrt(diag(z$cov))
    z$cov <- z$cov / (z$sdev %o% z$sdev)
  }
  
  return(z)
}

covRob <- function(data, data.name, corr = F, ...) {
  z <- RobStatTM::covRob(data, ...)
  z$corr <- corr
  z$call <- call("covRob", data = as.name(data.name))
  class(z) <- "covRob"
  
  if(corr) {
    z$sdev <- sqrt(diag(z$cov))
    z$cov <- z$cov / (z$sdev %o% z$sdev)
  }
  
  return(z)
}

covRobMM <- function(data, data.name, corr = F, ...) {
  z <- RobStatTM::covRob(data, type = "MM", ...)
  z$corr <- corr
  z$call <- call("covRobMM", data = as.name(data.name))
  class(z) <- "covRob"
  
  if(corr) {
    z$sdev <- sqrt(diag(z$cov))
    z$cov <- z$cov / (z$sdev %o% z$sdev)
  }
  
  return(z)
}

covRobRocke <- function(data, data.name, corr = F, ...) {
  z <- RobStatTM::covRob(data, type = "Rocke", ...)
  z$corr <- corr
  z$call <- call("covRobRocke", data = as.name(data.name))
  class(z) <- "covRob"
  
  if (corr) {
    z$sdev <- sqrt(diag(z$cov))
    z$cov <- z$cov / (z$sdev %o% z$sdev)
  }
  
  return(z)
}

prcompRob <- function(data, ncomp = ncol(data), desprop = 1.0, ...) {
  X = pcaRobS(data, ncomp = ncomp, desprop = desprop)
  
  # X.scaled <- scale(data, center = center, scale = scale)
  
  n <- ncol(X$eigvec)
  
  pc.index <- sapply(1L:n, function(i) { paste("PC", i) })
  
  z <- list()
  
  z$sdev     <- sqrt(diag(var(X$repre)))
  z$rotation <- X$eigvec
  colnames(z$rotation) <- pc.index
  
  z$center   <- X$mu
  z$scale    <- FALSE
  z$x        <- X$repre
  
  class(z) <- "prcompRob"
  
  z
}

# Custom summary and html output functions for methods above
htmlText <- function(x, ...) {
  UseMethod("htmlText", x)
}

htmlCoefTable <- function(coefs, lines, digits = max(3, getOption("digits") - 3)) {
    temp <- sapply(colnames(coefs), function(name) {  paste0("<th scope=\"col\" align=\"center\">", name, "</th>") })
    
    lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
    
    temp <- sapply(1:nrow(coefs),
                   function (i, coef) {
                     line <- paste0("<tr><th scope=\"row\">", row.names(coef)[i], "</th>")
                     temp <- sapply(coef[i, 1:3], function(val) { paste0("<td align=\"right\">", format(signif(val, digits = digits)), "</td>") })
                     paste0(line, paste(temp, collapse = ''), "<td align=\"right\">", format(signif(coef[i, 4], digits = digits)),
                            " ", stars.pval(coef[i, 4]),"</td></tr>")
                   },
                   coefs)
    
    lines <- paste0(lines, paste(temp, collapse = ''), "<tfoot><td colspan=\"5\">Signif. codes: ",
                    attr(stars.pval(0), 'legend'), "</td></tfoot></table></div>")
    
    return(lines)
}

htmlText.summary.lm <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  # Get residuals
  resid <- x$residuals
  df <- x$df
  resid.df <- df[2L]
  
  lines <- ""
  
  # Display residual statistics
  if (resid.df > 5L) {
  	nam <- c("Min", "1Q", "Median", "3Q", "Max")
  	resid.quantile <- if (length(dim(resid)) == 2L) {
              	        structure(apply(t(resid), 1L, quantile),
              		      dimnames = list(nam, dimnames(resid)[[2L]]))
                      }	else  {
                          zz <- zapsmall(quantile(resid), digits + 1L)
                          structure(zz, names = nam)
                      }
  	lines <- paste0(lines, "<div style=\"overflow:auto;\"><table cellpadding:\"10\">",
                    "<caption><font color=\"#000000\"><b>Residual Statistics:</b></font></caption>",
                    "<tr>")
  	
  	temp <- sapply(nam, function(name) {  paste0("<th scope=\"col\" align=\"center\">", name, "</th>") })
  	
  	lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
  	
  	temp <- sapply(resid.quantile, function(val) { paste0("<td align=\"right\">", format(signif(val, digits = digits)), "</td>") })
  	
  	lines <- paste0(lines, "<tr>", paste(temp, collapse = ''), "</tr>")
  } else if (rdf > 0L) {
    lines <- paste0(lines, "<div style=\"overflow:auto;\"><table cellpadding:\"10\">",
                    "<caption><font color=\"#000000\"><b>Residual Statistics:</b></font></caption>",
                    "<tr>")
    
    temp <- sapply(1L:length(resid), function(i) {  paste0("<th scope=\"col\" align=\"center\">", paste("Res.", i), "</th>") })
  	
  	lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
    
    temp <- sapply(resid, function(val) { paste0("<td align=\"right\">", format(signif(val, digits = digits)), "</td>") })
  	
  	lines <- paste0(lines, "<tr>", paste(temp, collapse = ''), "</tr></table></div>")
  } else { # resid.df == 0 : perfect fit!
    lines <- paste0(lines, "<div><font color=\"#000000\"><b>", "ALL ", df[1L],
                    " residuals are 0: no residual degrees of freedom!", "</b></font></div>" )
  }
  
  # Display coefficient statistics
  if (length(x$aliased) == 0L) {
    lines <- paste0(lines, "<div><font color=\"#000000\"><b>No Coefficients</b></font></div>" )
  } else {
    if (nsingular <- df[3L] - df[1L]) {
      lines <- paste0(lines, "<div><font color=\"#000000\"><b>", "Coefficients: (", nsingular,
                      " not defined because of singularities)", "</b></font></div>")
    } else {
      lines <- paste0(lines, "<div style=\"overflow:auto;\"><table cellpadding:\"10\">",
                      "<caption><font color=\"#000000\"><b>Coefficients:</b></font></caption>",
                      "<tr><td></td>")
    }
    
    coefs <- x$coefficients
    
    if(!is.null(aliased <- x$aliased) && any(aliased)) {
      cn <- names(aliased)
      coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn, colnames(coefs)))
      coefs[!aliased, ] <- x$coefficients
    }
    
    lines <- htmlCoefTable(coefs, lines, digits = digits)
  }
  
  # Display regression statistics
  lines <- paste0(lines, "<div><p><font color=\"#000000\">", "<br>Residual standard error: ", format(signif(x$sigma, digits)),
                  " on ", resid.df, " degrees of freedom.", "<br>" )
  
  if(nzchar(mess <- naprint(x$na.action))) {
    lines <- paste0(lines, "(", mess, ")<br>")
  }
  
  if (!is.null(x$fstatistic)) {
    
    lines <- paste0(lines, "<div><p><font color=\"#000000\">", "Multiple R-squared:", formatC(x$r.squared, digits = digits),
                    ",  Adjusted R-squared: ", formatC(x$adj.r.squared, digits = digits), "<br>",
                    "F-statistic: ", formatC(x$fstatistic[1L], digits = digits), " on ", x$fstatistic[2L], " and ",
              	    x$fstatistic[3L], " DF,  p-value: ",
              	    format.pval(pf(x$fstatistic[1L], x$fstatistic[2L],
                                   x$fstatistic[3L], lower.tail = FALSE),
                                   digits = digits))
  }
  
  lines <- paste0(lines, "</font></p></div>")
  
#   correl <- x$correlation
#   
#   if (!is.null(correl)) {
#   	p <- NCOL(correl)
#   	if (p > 1L) {
# 	    cat("\nCorrelation of Coefficients:\n")
# 	    if(is.logical(symbolic.cor) && symbolic.cor) {# NULL < 1.7.0 objects
# 		    print(symnum(correl, abbr.colnames = NULL))
#       } else {
#         correl <- format(round(correl, 2), nsmall = 2, digits = digits)
#         correl[!lower.tri(correl)] <- ""
#         print(correl[-1, -p, drop=FALSE], quote = FALSE)
#       }
# 	  }
#   }
#   
#   cat("\n")#- not in S
  return(paste0("<div style=\"border: 1px solid #ccc;
                                  border-radius: 6px;
                                  padding: 0px 5px;
                                  margin: 5px -10px;
                                  background-color: #f5f5f5;\">",
                lines,
                "</div>"))
}

htmlText.summary.lmrobdetMM <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  lines <- ""
  
  control <- x$control
  # Residual statistics
  resid <- x$residuals
  df <- x$df
  rdf <- df[2L]
  if (!is.null(x$weights) && diff(range(x$weights))) {
    lines <- paste0(lines, "<div style=\"overflow:auto;\"><table cellpadding:\"10\">",
                    "<caption><font color=\"#000000\"><b>Weighted Residual Statistics:</b></font></caption>",
                    "<tr>")
  } else {
    lines <- paste0(lines, "<div style=\"overflow:auto;\"><table cellpadding:\"10\">",
                    "<caption><font color=\"#000000\"><b>Residual Statistics:</b></font></caption>",
                    "<tr>")
  }

  if (rdf > 5L) {
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
  	resid.quantile <- if (length(dim(resid)) == 2L) {
              	        structure(apply(t(resid), 1L, quantile),
              		      dimnames = list(nam, dimnames(resid)[[2L]]))
                      }	else  {
                          zz <- zapsmall(quantile(resid), digits + 1L)
                          structure(zz, names = nam)
                      }
  	
  	
  	temp <- sapply(nam, function(name) {  paste0("<th scope=\"col\" align=\"center\">", name, "</th>") })
  	
  	lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
  	
  	temp <- sapply(resid.quantile, function(val) { paste0("<td align=\"right\">", format(signif(val, digits = digits)), "</td>") })
  	
  	lines <- paste0(lines, "<tr>", paste(temp, collapse = ''), "</tr>")
  } else if (rdf > 0L) {
    lines <- paste0(lines, "<div style=\"overflow:auto;\"><table cellpadding:\"10\">",
                    "<caption><font color=\"#000000\"><b>Residual Statistics:</b></font></caption>",
                    "<tr>")
    
    temp <- sapply(1L:length(resid), function(i) {  paste0("<th scope=\"col\" align=\"center\">", paste("Res.", i), "</th>") })
  	
  	lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
    
    temp <- sapply(resid, function(val) { paste0("<td align=\"right\">", format(signif(val, digits = digits)), "</td>") })
  	
  	lines <- paste0(lines, "<tr>", paste(temp, collapse = ''), "</tr></table></div>")
  } else { # resid.df == 0 : perfect fit!
    lines <- paste0(lines, "<div><font color=\"#000000\"><b>", "ALL ", df[1L],
                    " residuals are 0: no residual degrees of freedom!", "</b></font></div>" )
  }
  
  # Coefficient statistics
  if( length(x$aliased) ) {
    if( !(x$converged) ) {
      if (x$scale == 0) {
        lines <- paste0(lines, "<div style=\"overflow:auto;\"><p>Exact fit detected!<br></p>",
                        "<table cellpadding:\"10\"><caption><font color=\"#000000\"><b>",
                        "Coefficients:</b></font></caption><tr><td></td>")
      } else {
        if (control$method == "S") {
          lines <- paste0(lines, "<div style=\"overflow:auto;\"><p>Algorithm did not converge!<br></p>",
                          "<table cellpadding:\"10\"><caption><font color=\"#000000\"><b>",
                          "Coefficients of the *initial* S-estimator:</b></font></caption><tr><td></td>")
        } else {
          lines <- paste0(lines, "<div style=\"overflow:auto;\"><p>Algorithm did not converge!<br></p>",
                          "<table cellpadding:\"10\"><caption><font color=\"#000000\"><b>",
                          "Coefficients of the ", control$method, "-estimator:</b></font></caption><tr><td></td>")
        }
      }
      
      lines <- htmlCoefTable(coefs, lines, digits = digits)
    } else {
      if (nsingular <- df[3L] - df[1L]) {
        lines <- paste0(lines, "<div><font color=\"#000000\"><b>", "Coefficients: (", nsingular,
                      " not defined because of singularities)", "</b></font></div>")
      } else {
        lines <- paste0(lines, "<div style=\"overflow:auto;\"><table cellpadding:\"10\">",
                      "<caption><font color=\"#000000\"><b>Coefficients:</b></font></caption>",
                      "<tr><td></td>")
      }
      
      coefs <- x$coefficients
      
      if(!is.null(aliased <- x$aliased) && any(aliased)) {
        cn <- names(aliased)
        coefs <- matrix(NA, length(aliased), 4, dimnames=list(cn, colnames(coefs)))
        coefs[!aliased, ] <- x$coefficients
      }
      
      lines <- htmlCoefTable(coefs, lines, digits = digits)

      lines <- paste0(lines, "<div><p><font color=\"#000000\">", "<br>Robust residual standard error: ", format(signif(x$scale, digits)),"<br>")
      
      if (!is.null(x$r.squared) && x$df[1] != attr(x$terms, "intercept")) {
        lines <- paste0(lines, "Multiple R-squared:  ", formatC(x$r.squared, digits = digits),
                        ",  Adjusted R-squared: ", formatC(x$adj.r.squared, digits = digits), "<br>", "Convergence in ", x$iter,
                        " IRWLS iterations")
      }
    }
  } else {
    lines <- paste0(lines, "<div><font color=\"#000000\"><b>No Coefficients</b></font></div>")
  }
  
  return(paste0("<div style=\"border: 1px solid #ccc;
                                  border-radius: 6px;
                                  padding: 0px 5px;
                                  margin: 5px -10px;
                                  background-color: #f5f5f5;\">",
                lines,
                "</div>"))
}

htmlText.summary.lmfm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  lines <- ""
  n.models <- length(x)
  mod.names <- names(x)
  fancy.names <- format(paste(mod.names, ":", sep = ""), justify = "right")
  
  resid.q <- t(sapply(x, function(u) quantile(residuals(u), na.rm = TRUE)))
  colnames(resid.q) <- c("Min", "1Q", "Median", "3Q", "Max")

  lines <- paste0(lines, "<div style=\"overflow-x:auto;\"><table cellpadding:\"10\">",
                  "<caption><font color=\"#000000\"><b>Residual Statistics:</b></font></caption><tr><td><b>Method</b></td>")
  
  temp <- sapply(colnames(resid.q), function(name) { paste0("<th scope=\"col\" align=\"center\">", name, "</th>") })
  
  lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
  
  temp <- sapply(resid.q[1, ],
                 function(val) {
                   paste0("<td align=\"right\"><font color=\"#FF0000\">", format(signif(val, digits = digits)), "</font></td>")
                 })
  
  lines <- paste0(lines, "<tr><td align=\"left\"><font color=\"#FF0000\">", mod.names[1], "</font></td>", paste(temp, collapse = ''), '</tr>')
  
  temp <- sapply(resid.q[2, ],
                 function(val) {
                   paste0("<td align=\"right\"><font color=\"#000000\">", format(signif(val, digits = digits)), "</font></td>")
                 })
  
  lines <- paste0(lines, "<tr><td align=\"left\"><font color=\"#000000\">", mod.names[2], "</font></td>", paste(temp, collapse = ''), "</tr></table></div>")

  coefs <- lapply(x, coef)
  variables <- lapply(coefs, rownames)
  
  var.common <- intersect(variables[[1]], variables[[2]])
  var.uncommon <- lapply(variables, setdiff, var.common)
  var.names <- c(var.common, unlist(var.uncommon))
  model.var.names <- lapply(var.uncommon, function(x, y) { c(y, x) }, var.common)
  n.coefs <- length(var.names)
  coefmat <- matrix(NA, n.coefs * n.models, 4)
  r.names <- rep(fancy.names, n.coefs)
  rownames(coefmat) <- format(r.names, justify = "right")
  names <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  colnames(coefmat) <- names

  for(i in 1:n.models) {
    mc1 <- coefs[[i]][which(model.var.names[[i]] %in% var.common), ]
    mc2 <- coefs[[i]][which(model.var.names[[i]] %in% var.uncommon[[i]]), ]
    coefmat[n.models * 0:(n.coefs - 1) + i, ][which(var.names %in% model.var.names[[i]]), ] <- matrix(rbind(mc1, mc2), ncol = 4)
  }

  lines <- paste0(lines, "<div style=\"overflow:auto;\"><table cellpadding:\"10\">",
                  "<caption><font color=\"#000000\"><b>Coefficients:</b></font></caption>",
                  "<tr><td></td><th scope=\"col\" align=\"center\">Method</th>")
  
  temp <- sapply(names, function(name) { paste0("<th scope=\"col\" align=\"center\">", name, "</th>") })
  
  lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
  
  vars <- var.names
  
  for (i in 1:n.coefs) {
    ind <- 2 * i - 1
    
    lines <- paste0(lines, "<tr><th scope=\"row\" rowspan=\"2\">", vars[i], "</th>",
                    "<td align=\"left\"><font color=\"#FF0000\">", mod.names[1L], "</font></td>")
    temp <- sapply(coefmat[ind, 1L:3L],
                   function(val) {
                     paste0("<td align=\"right\"><font color=\"#FF0000\">", format(signif(val, digits = digits)), "</font></td>")
                   })
    
    lines <- paste0(lines, paste(temp, collapse = ''), "<td align=\"right\"><font color=\"#FF0000\">",
                    format(signif(coefmat[ind, 4L], digits = digits, ...)), " ", stars.pval(coefmat[ind, 4L]),
                    "</font></td></tr><tr><td align=\"left\">", mod.names[2L], "</td>")
    
    
    temp <- sapply(coefmat[ind + 1L, 1L:3L],
                   function(val) {
                     paste0("<td align=\"right\"><font color=\"#000000\">", format(signif(val, digits = digits, ...)), "</font></td>")
                   })
    
    lines <- paste0(lines, paste(temp, collapse = ''), "<td align=\"right\"><font color=\"#000000\">",
                    format(signif(coefmat[ind + 1, 4L], digits = digits, ...)),
                    " ", stars.pval(coefmat[ind + 1, 4L]), "</font></td></tr>")
  }
  
  lines <- paste0(lines, "<tfoot><td colspan=\"5\">Signif. codes: ",
                  attr(stars.pval(0), 'legend'), "</td></tfoot></table></div>")

  lines <- paste0(lines, "<div><p><font color=\"#000000\">", "<br><b>Residual Scale Estimates:</b><br>")
  for(i in 1:n.models) {
    lines <- paste0(lines, "\t", fancy.names[i], " ", format(x[[i]]$sigma, digits = digits, ...), " on ",
                    x[[i]]$df[2L], " degrees of freedom<br>")
  }
  
  rsq <- sapply(x, function(u) u$r.squared)
  has.rsq <- which(!sapply(rsq, function(u) is.null(u) || is.na(u)))
  if(length(has.rsq)) {
    lines <- paste0(lines, "<br><b>Multiple R-squared:</b>")
    for(i in 1:n.models) {
      lines <- paste0(lines, "<br>  ", fancy.names[i], " ", format(rsq[i], digits = digits))
    }
  }
  
  lines <- paste0(lines, "</font></p></div>")

  # Add correlations in a later release?
  
  # correlations <- lapply(x, function(u) u$correlation)
  # if(all(!sapply(correlations, is.null))) {    
  #   if(any(sapply(correlations, NCOL) > 1)) {
  #     cat("Correlations:\n")
  #     for(i in 1:n.models) {
  #       if((p <- NCOL(correlations[[i]])) > 1) {        
  #         correl <- format(round(correlations[[i]], 2), nsmall = 2,
  #                          digits = digits, ...)
  #         correl[!lower.tri(correl)] <- ""
  #         cat(mod.names[i], ":\n", sep = "")
  #         print(correl[-1, -p, drop = FALSE], quote = FALSE, ...)
  #       }
  #       cat("\n")
  #     }
  #   }
  # }
  
  return(paste0("<div style=\"border: 1px solid #ccc;
                                  border-radius: 6px;
                                  padding: 0px 5px;
                                  margin: 5px -10px;
                                  background-color: #f5f5f5;\">",
                lines,
                "</div>"))
}

htmlText.summary.covClassic <- function(x, digits = max(3L, getOption("digits") - 3L),
                                        print.distance = FALSE, ...) {
  lines <- paste0("<div style=\"overflow:auto;\"><table cellpadding:\"10\">",
                  "<caption><font color=\"#000000\"><b>Classical Estimate of Covariance:</b></font></caption>",
                  "<tr><td></td>")
  
  temp <- sapply(colnames(x$cov), function(name) { paste0("<th scope=\"col\" align=\"center\">", name, "</th>") })
  
  lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
  
  for (i in 1:nrow(x$cov)) {
    lines <- paste0(lines, "<tr><th scope=\"row\">", row.names(x$cov)[i], "</th>")
    temp <- sapply(x$cov[i, ], function(val) { paste0("<td align=\"right\">", format(val, scientific = T, digits = 4), "</td>") })
    lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
  }
  
  lines <- paste0(lines, "</table></div>")
  
  lines <- paste0(lines, "<div style=\"overflow-x:auto;\"><table cellpadding:\"10\">",
                  "<caption><font color=\"#000000\"><b>Classical Estimate of Location:</b></font></caption><tr>")
  
  temp <- sapply(names(x$center), function(name) { paste0("<th scope=\"col\" align=\"center\">", name, "</th>") })
  
  lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
  
  temp <- sapply(x$center, function(val) { paste0("<td align=\"right\">", format(val, scientific = T, digits = 4), "</td>") })
  
  lines <- paste0(lines, "<tr>", paste(temp, collapse = ''), "</tr></table></div>")
  
  lines <- paste0(lines, "<div style=\"overflow-x:auto;\"><table>",
                  "<caption><font color=\"#000000\"><b>Eigenvalues:</b></font></caption><tr>")
  
  temp <- sapply(names(x$evals), function(name) { paste0("<th scope=\"col\ align=\"center\">", name, "</th>") })
  
  lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
  
  temp <- sapply(x$evals, function(val) { 
                   paste0("<td align=\"right\">", format(val, scientific = T, digits = 4), "</td>")
                 })
  
  lines <- paste0(lines, "<tr>", paste(temp, collapse = ''), "</tr></table></div>")
  
  return(paste0("<div style=\"border: 1px solid #ccc;
                                  border-radius: 6px;
                                  padding: 0px 5px;
                                  margin: 5px -10px;
                                  background-color: #f5f5f5;\">",
                lines,
                "</div>"))
}

htmlText.summary.covRob <- function(x, digits = max(3L, getOption("digits") - 3L),
                                      print.distance = FALSE, ...) {
  lines <- paste0("<div style=\"overflow:auto;\"><table cellpadding:\"10\">",
                  "<caption><font color=\"#000000\"><b>Robust Estimate of Covariance:</b></font></caption>",
                  "<tr><td></td>")
  
  temp <- sapply(colnames(x$cov), function(name) { paste0("<th scope=\"col\" align=\"center\">", name, "</th>") })
  
  lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
  
  for (i in 1:nrow(x$cov)) {
    lines <- paste0(lines, "<tr><th scope=\"row\">", row.names(x$cov)[i], "</th>")
    temp <- sapply(x$cov[i, ], function(val) { paste0("<td align=\"right\">", format(val, scientific = T, digits = 4), "</td>") })
    lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
  }
  
  lines <- paste0(lines, "</table></div>")
  
  lines <- paste0(lines, "<div style=\"overflow-x:auto;\"><table cellpadding:\"10\">",
                  "<caption><font color=\"#000000\"><b>Robust Estimate of Location:</b></font></caption><tr>")
  
  temp <- sapply(names(x$center), function(name) { paste0("<th scope=\"col\" align=\"center\">", name, "</th>") })
  
  lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
  
  temp <- sapply(x$center, function(val) { paste0("<td align=\"right\">", format(val, scientific = T, digits = 4), "</td>") })
  
  lines <- paste0(lines, "<tr>", paste(temp, collapse = ''), "</tr></table></div>")
  
  lines <- paste0(lines, "<div style=\"overflow-x:auto;\"><table>",
                  "<caption><font color=\"#000000\"><b>Eigenvalues:</b></font></caption><tr>")
  
  temp <- sapply(names(x$evals), function(name) { paste0("<th scope=\"col\ align=\"center\">", name, "</th>") })
  
  lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
  
  temp <- sapply(x$evals, function(val) { 
                   paste0("<td align=\"right\">", format(val, scientific = T, digits = 4), "</td>")
                 })
  
  lines <- paste0(lines, "<tr>", paste(temp, collapse = ''), "</tr></table></div>")
  
  return(paste0("<div style=\"border: 1px solid #ccc;
                                  border-radius: 6px;
                                  padding: 0px 5px;
                                  margin: 5px -10px;
                                  background-color: #f5f5f5;\">",
                lines,
                "</div>"))
}

htmlText.summary.covfm <- function(x, digits = max(3L, getOption("digits") - 3L),
                                     print.distance = FALSE, ...) {
  mod.names <- names(x)
  
  lines <- paste0("<div style=\"overflow:auto;\"><table cellpadding:\"10\">",
                  "<caption><font color=\"#000000\"><b>Comparison of Covariance/Correlation Estimates:</b></font></caption>",
                  "<tr><td></td><th scope=\"col\" align=\"center\">Method</th>")
  
  temp <- sapply(colnames(x[[1]]$cov), function(name) { paste0("<th scope=\"col\" align=\"center\">", name, "</th>") })
  
  lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
  
  
  for (i in 1:nrow(x[[1]]$cov)) {
    lines <- paste0(lines, "<tr><th scope=\"row\" rowspan=\"2\">", row.names(x[[1]]$cov)[i], "</th>",
                    "<td align=\"left\"><font color=\"#FF0000\">", mod.names[1], "</font></td>")
    temp <- sapply(x[[1]]$cov[i, ],
                   function(val) {
                     paste0("<td align=\"right\"><font color=\"#FF0000\">", format(signif(val, digits = digits)), "</font></td>")
                   })
    lines <- paste0(lines, paste(temp, collapse = ''), "</tr><tr><td align=\"left\">", mod.names[2], "</td>")
    
    temp <- sapply(x[[2]]$cov[i, ],
                   function(val) {
                     paste0("<td align=\"right\"><font color=\"#000000\">", format(signif(val, digits = digits)), "</font></td>")
                   })
    
    lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
  }
  
  lines <- paste0(lines, "</table></div>")
  
  lines <- paste0(lines, "<div style=\"overflow-x:auto;\"><table cellpadding:\"10\">",
                  "<caption><font color=\"#000000\"><b>Robust Estimate of Location:</b></font></caption><tr><td><b>Method</b></td>")
  
  temp <- sapply(names(x[[1]]$center), function(name) { paste0("<th scope=\"col\" align=\"center\">", name, "</th>") })
  
  lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
  
  temp <- sapply(x[[1]]$center,
                 function(val) {
                   paste0("<td align=\"right\"><font color=\"#FF0000\">", format(signif(val, digits = digits)), "</font></td>")
                 })
  
  lines <- paste0(lines, "<tr><td align=\"left\"><font color=\"#FF0000\">", mod.names[1], "</font></td>", paste(temp, collapse = ''), '</tr>')
  
  temp <- sapply(x[[2]]$center,
                 function(val) {
                   paste0("<td align=\"right\"><font color=\"#000000\">", format(signif(val, digits = digits)), "</font></td>")
                 })
  
  lines <- paste0(lines, "<tr><td align=\"left\"><font color=\"#000000\">", mod.names[2], "</font></td>", paste(temp, collapse = ''), "</tr></table></div>")
  
  lines <- paste0(lines, "<div style=\"overflow-x:auto;\"><table>",
                  "<caption><font color=\"#000000\"><b>Eigenvalues:</b></font></caption><tr><td><b>Method</b></td>")
  
  temp <- sapply(names(x[[1]]$evals), function(name) { paste0("<th scope=\"col\ align=\"center\">", name, "</th>") })
  
  lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
  
  temp <- sapply(x[[1]]$evals,
                 function(val) { 
                   paste0("<td align=\"right\"><font color=\"#FF0000\">", format(signif(val, digits = digits)), "</font></td>")
                 })
  
  lines <- paste0(lines, "<tr><td align=\"left\"><font color=\"#FF0000\">", mod.names[1], "</font></td>", paste(temp, collapse = ''), "</tr>")
  
  temp <- sapply(x[[2]]$evals,
                 function(val) { 
                   paste0("<td align=\"right\"><font color=\"#000000\">", format(signif(val, digits = digits)), "</font></td>")
                 })
  
  lines <- paste0(lines, "<tr><td align=\"left\"><font color=\"#000000\">", mod.names[2],
                  "</font></td>", paste(temp, collapse = ''), "</tr></table></div>")

  return(paste0("<div style=\"border: 1px solid #ccc;
                                  border-radius: 6px;
                                  padding: 0px 5px;
                                  margin: 5px -10px;
                                  background-color: #f5f5f5;\">",
                lines,
                "</div>"))
}

# Custom stuff for robust pca

summary.prcompRob <- function(object, ...) {
  chkDots(...)
  
  vars <- object$sdev^2
  vars <- vars / sum(vars)
  
  importance <- rbind("Standard deviation" = object$sdev,
                      "Proportion of Variance" = round(vars, 5L),
                      "Cumulative Proportion" = round(cumsum(vars), 5L))
  
  colnames(importance) <- colnames(object$rotation)
  
  object$importance <- importance
  
  class(object) <- "summary.prcompRob"
  object
}

htmlText.summary.prcompRob <- function(x,
                                      digits = max(3L, getOption("digits") - 3L),
                                      ...) {
  lines <- paste0("<div style=\"overflow:auto;\"><table cellpadding:\"10\">",
                  "<caption><font color=\"#000000\"><b>Importance of Components:</b></font></caption>",
                  "<tr><td></td>")
  temp <- sapply(colnames(x$importance), function(name) { paste0("<th scope=\"col\" align=\"center\">", name, "</th>") })
  
  lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
  
  for (i in 1:nrow(x$importance)) {
    lines <- paste0(lines, "<tr><th scope=\"row\">", row.names(x$importance)[i], "</th>")
    temp <- sapply(x$importance[i, ], function(val) { paste0("<td align=\"right\">", signif(val, 4), "</td>") })
    lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
  }
  
  paste0(lines, "</table></div>")
  
  return(paste0("<div style=\"border: 1px solid #ccc;
                                  border-radius: 6px;
                                  padding: 0px 5px;
                                  margin: 5px -10px;
                                  background-color: #f5f5f5;\">",
                lines,
                "</div>"))
}

htmlText.summary.prcomp <- function(x,
                                    digits = max(3L, getOption("digits") - 3L),
                                    ...) {
  lines <- paste0("<div style=\"overflow:auto;\"><table cellpadding:\"10\">",
                  "<caption><font color=\"#000000\"><b>Importance of Components:</b></font></caption>",
                  "<tr><td></td>")
  temp  <- sapply(colnames(x$importance), function(name) { paste0("<th scope=\"col\" align=\"center\">", name, "</th>") })
  
  lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
  
  for (i in 1L:nrow(x$importance)) {
    lines <- paste0(lines, "<tr><th scope=\"row\">", row.names(x$importance)[i], "</th>")
    temp  <- sapply(x$importance[i, ], function(val) { paste0("<td align=\"right\">", signif(val, 4), "</td>") })
    lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
  }
  
  lines <- paste0(lines, "</table></div>")
  
  return(paste0("<div style=\"border: 1px solid #ccc;
                                  border-radius: 6px;
                                  padding: 0px 5px;
                                  margin: 5px -10px;
                                  background-color: #f5f5f5;\">",
                lines,
                "</div>"))
}

summary.pcompfm <- function(object, ...) {
  object$Classic <- summary(object$Classic)
  object$Robust  <- summary(object$Robust)
  oldClass(object) <- "summary.pcompfm"
  object
}

# Summary method comparing
htmlText.summary.pcompfm <- function(x,
                                     digits = max(3L, getOption("digits") - 3L),
                                     ...) {
  mod.names <- names(x)
  
  lines <- paste0("<div style=\"overflow:auto;\"><table cellpadding:\"10\">",
                  "<caption><font color=\"#000000\"><b>Importance of Components:</b></font></caption>",
                  "<tr><td></td><td><b>Method</b></td>")
  
  temp  <- sapply(colnames(x[[1L]]$importance), function(name) { paste0("<th scope=\"col\" align=\"center\">", name, "</th>") })
  
  lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
  
  for (i in 1L:nrow(x[[1L]]$importance)) {
    lines <- paste0(lines, "<tr><th scope=\"row\" rowspan=\"2\">", row.names(x[[1L]]$importance)[i], "</th>",
                    "<td align=\"left\"><font color=\"#FF0000\">", mod.names[1L], "</font></td>")
    temp  <- sapply(x[[1L]]$importance[i, ],
                    function(val) {
                      paste0("<td align=\"right\"><font color=\"#FF0000\">", signif(val, 4L), "</font></td>")
                    })
    
    lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
    lines <- paste0(lines, "<td align=\"left\"><font color=\"#000000\">", mod.names[2L], "</font></td>")
    
    temp  <- sapply(x[[2L]]$importance[i, ],
                    function(val) {
                      paste0("<td align=\"right\"><font color=\"#000000\">", signif(val, 4L), "</font></td>")
                    })
    
    lines <- paste0(lines, paste(temp, collapse = ''), "</tr>")
  }
  
  lines <- paste0(lines, "</table></div>")
  
  return(paste0("<div style=\"border: 1px solid #ccc;
                                  border-radius: 6px;
                                  padding: 0px 5px;
                                  margin: 5px -10px;
                                  background-color: #f5f5f5;\">",
                lines,
                "</div>"))
}

# Back-end implementation of Shiny server
shinyServer(function(input, output) {
  
  # Need this to store reactive objects
  values <- reactiveValues()
  
  values$linRegress.active <- F
  values$covariance.active <- F
  values$pca.active        <- F
  
  values$linRegress.plots.active <- F
  values$covariance.plots.active <- F
  values$pca.plots.active        <- F
  
  
####################  
## Data Selection ##
####################
  
  # Display correct options to obtain data given
  # method selected by user
  output$select.data <- renderUI({
    # If no package is chosen, then do nothing
    if(is.null(input$library)) {
      return("")
    }
    
    # List of datasets
    lst <- data(package = input$library)
    
    item.lst <- lst$results[, "Item"][grepl(" \\(", lst$results[, "Item"]) == FALSE]
    
    # Create selection for datasets
    selectInput("dataset",
                label   = "Select Dataset",
                choices = item.lst)
  })
  
  observeEvent(input$display.table, {
    output$data.panel <- renderUI({
      fluidPage(
        # Create Data table
        fluidRow(
          wellPanel(DT::dataTableOutput("contents.table"))
        ),
        
        fluidRow(
          column(2, offset = 5,
            disabled(
              
              actionLink("data.info.link", label = "More Info")
            )
          )
        )
      )
    })
    
    output$contents.table <- DT::renderDataTable({
      contents_table()
    })
  })
  
  # On-click, load the data and return the data frame
  contents_table <- eventReactive(input$display.table, {
    if (input$source == "upload") {
      req(input$file)
      
      # Read data from file
      values$dat <- tryCatch({
                      read.csv(input$file$datapath,
                               header = input$header,
                               sep    = input$sep,
                               quote  = input$quote,
                               numerals = "allow.loss")
                    },
                    error = function(cond) {
                      message("There appears to be an error. Try checking the format of the CSV file.")
                      message("Here's the original error message:")
                      message(cond)
                      
                      return(NA)
                    })
      
      class(values$dat[[1]])
      
      if (input$data.ts == TRUE) {
        values$dat[, 1] <- as.Date(as.character(values$dat[, 1]))
        
        values$dat <- xts(values$dat[, -1], values$dat[, 1])
      }
      
      # If there are no headers, give data headers
      if (input$header == FALSE) {
        colnames(dat)[-1] <- paste0('X', 1L:(ncol(dat) - 1L))
      }
    } else {
      # If no dataset exists, return nothing
      if (is.null(input$dataset)) {
        return()
      }
      
      # Obtain specified dataset from specified package
      data(list = input$dataset, package = input$library)
      
      values$dat <- get(input$dataset)
    }
    
    if (is.vector(values$dat)) {
      data <- data.frame(values$dat)
      
      colnames(data) <- input$dataset
      
      values$dat <- data
      
      values$dat.numeric.variables <- colnames(values$dat)
    } else {
      # Get variable names
      if (all(class(values$dat) == "zoo")) {
        values$dat <- xts(values$dat)
      }
      
      values$dat.variables <- colnames(values$dat)
      
      num.index <- sapply(values$dat, is.numeric)
      
      values$dat.numeric.variables <- values$dat.variables[num.index]
      
      values$dat.numeric <- as.data.frame(values$dat[values$dat.numeric.variables])
      
      data <- values$dat
      
      data[, num.index] <- round(data[, num.index], 3L)
    }
    
    return (as.data.frame(data))
  })
  
  
####################
## Location/Scale ##
####################
  
    # Render variable input list
  output$locScale.select.variable <- renderUI({
    # If there is no data, do nothing
    if (is.null(dim(values$dat))) {
      return("")
    }
    
    # Render select input for variables
    selectInput("locScale.variable", "Variable",
                choices = values$dat.numeric.variables)
  })
  
  output$locScale.eff.options <- renderUI({
    if (any(input$locScale.psi == c("modopt", "optimal"))) {
      radioButtons("locScale.eff", "Asymptotic Efficiency",
                            choices  = c("0.85" = 0.85, "0.9" = 0.9, "0.95" = 0.95, "0.99" = 0.99),
                            selected = 0.99)
    } else {
      radioButtons("locScale.eff", "Asymptotic Efficiency",
                            choices  = c("0.85" = 0.85, "0.9" = 0.9, "0.95" = 0.95),
                            selected = 0.9)
    }
  })
  
  # On-click, find the estimators and create string object of results
  contents_estimators <- eventReactive(input$locScale.display, {
    if (is.null(dim(values$dat))) {
        return(paste0("<font color=\"#FF0000\"><b>", "ERROR: No data loaded!", "</b><font>"))
    }
    
    data <- as.numeric(values$dat[, input$locScale.variable])
    
    if (input$locScale.method == 'robust') {
      est <- locScaleM(x     = data,
                       psi   = input$locScale.psi,
                       eff   = as.numeric(input$locScale.eff))
      line1 <- paste0("<font color=\"#000000\">Location <strong>(SE)</strong>: ", round(est$mu, 4), " (<strong>", round(est$std.mu, 4), ")</strong><br>")
      line2 <- paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Scale: ", round(est$disper, 4),"</font>")

      return(paste0("<div style=\"border: 1px solid #ccc;
                                  border-radius: 6px;
                                  padding: 0px 5px;
                                  margin: 5px -10px;
                                  background-color: #f5f5f5;\">",
                    line1, line2,
                    "</div>"))
      
    } else if (input$locScale.method == 'classic') {
      est <- locScaleClassic(data)
      
      line1 <- paste0("<font color=\"#000000\">Location <strong>(SE)</strong>: ", round(est$mu, 4), " (<strong>", round(est$std.mu, 4), ")</strong><br>")
      line2 <- paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Scale: ", round(est$disper, 4),"</font>")

      return(paste0("<div style=\"border: 1px solid #ccc;
                                  border-radius: 6px;
                                  padding: 0px 5px;
                                  margin: 5px -10px;
                                  background-color: #f5f5f5;\">",
                      line1, line2,
                    "</div>"))
    } else {
      est1 <- locScaleClassic(data)
      
      est2 <- locScaleM(x     = data,
                        psi   = input$locScale.psi,
                        eff   = as.numeric(input$locScale.eff))
      
      line1 <- paste0("<font color=\"#000000\">Comparison of Location <strong>(SE)</strong>:<br>Classical&emsp;",
                      signif(est1$mu, 3L), " (<strong>", signif(est1$std.mu, 3L), ")</strong><br>Robust&nbsp;&nbsp;&nbsp;&emsp;",
                      signif(est2$mu, 3L), " (<strong>", signif(est2$std.mu, 3L), ")</strong><br><br>")
      line2 <- paste0("Comparison of Scale:<br> Classical&emsp;", signif(est1$disper, 3L), "<br> Robust&nbsp;&nbsp;&nbsp;&emsp;",
                      signif(est2$disper, 3L), "</font>")

      return(paste0("<div style=\"border: 1px solid #ccc;
                                  border-radius: 6px;
                                  padding: 0px 5px;
                                  margin: 5px -10px;
                                  background-color: #f5f5f5;\">",
                    line1, line2,
                    "</div>"))
    }
    
    # Get values for location and scale using 'locscaleM' function from
    # 'RobStatTM' package
    
    
    # Store results in string objects
    
    
  })
  
  # Display output
  output$locScale.Results <- renderText({
    contents_estimators()
  })
  
#######################
## Linear Regression ##
#######################

  ## Running Regression ##
  
  values$linRegress.methods <- c("LS", "M", "MM", "DCML")
  values$linRegress.functions  <- c("lm", "lmrobM", "lmrobdetMM", "lmrobdetDCML")
  
  output$linRegress.options <- renderUI({
    if (input$linRegress.second.method) {
      tabPanel("",
        fluidRow(
          column(6L,
            selectInput("linRegress.fit.option", "Method",
                        choices   = c("LS", "M", "MM", "DCML"),
                        selected  = "MM"),
            
            # List of dependent variables, must be selected
            uiOutput("linRegress.select.dependent"),
            
            # List of predictors to choose from
            uiOutput("linRegress.select.independent"),
            
            # String representing regression formula of form Y ~ X_0 + ... + X_n
            uiOutput("linRegress.formula"),
            
            # Interface to robust options
            uiOutput("linRegress.robust.control")
          ),
          
          column(6L,
            selectInput("linRegress.fit.option2", "Method",
                        choices   = c("LS", "M", "MM", "DCML"),
                        selected  = "MM"),
            
            # List of dependent variables, must be selected
            uiOutput("linRegress.select.dependent2"),
            
            # List of predictors to choose from
            uiOutput("linRegress.select.independent2"),
            
            # String representing regression formula of form Y ~ X_0 + ... + X_n
            uiOutput("linRegress.formula2"),
            
            # Interface to robust options
            uiOutput("linRegress.robust.control2")
          )
        )
      )
    } else {
      tabPanel("",
        selectInput("linRegress.fit.option", "Method",
                    choices   = c("LS", "M", "MM", "DCML"),
                    selected  = "MM"),
        
        # List of dependent variables, must be selected
        uiOutput("linRegress.select.dependent"),
        
        # List of predictors to choose from
        uiOutput("linRegress.select.independent"),
        
        # String representing regression formula of form Y ~ X_0 + ... + X_n
        uiOutput("linRegress.formula"),
        
        # Interface to robust options
        uiOutput("linRegress.robust.control")
      )
    }
  })

  output$linRegress.select.dependent <- renderUI({
    # If there is no data, do nothing
    if (is.null(dim(values$dat))) {
      return("")
    }
    
    # Render select input for variables
    selectInput("linRegress.dependent", "Dependent",
                choices = values$dat.numeric.variables)
  })

  output$linRegress.select.independent <- renderUI({
    # If no dependent variable is selected, do nothing
    if (is.null(input$linRegress.dependent)) {
      return("")
    }
    
    ind.vars <- setdiff(values$dat.numeric.variables, input$linRegress.dependent)
    
    # Render select input for variables
    selectInput("linRegress.independent", "Independent",
                choices  = ind.vars,
                selected = ind.vars[1],
                multiple = TRUE)
  })
  
  output$linRegress.formula <- renderUI({
    if (is.null(input$linRegress.dependent)) {
      return("")
    }
    
    formula.str <- paste0(input$linRegress.dependent, " ~ ")
    
    ind.vars <- input$linRegress.independent
    
    n <- length(ind.vars)
    
    if (n > 1L) {
      for (i in 1L:(n - 1L)) {
        formula.str <- paste0(formula.str, ind.vars[i], " + ")
      }
    }
    
    formula.str <- paste0(formula.str, ind.vars[n])
    
    textInput("linRegress.formula.text", "Formula",
              formula.str)
  })
  
  output$linRegress.select.dependent2 <- renderUI({
    # If there is no data, do nothing
    if (is.null(dim(values$dat))) {
      return("")
    }
    
    # Render select input for variables
    selectInput("linRegress.dependent2", "Dependent",
                choices = values$dat.numeric.variables)
  })

  output$linRegress.select.independent2 <- renderUI({
    # If no dependent variable is selected, do nothing
    if (is.null(input$linRegress.dependent2)) {
      return("")
    }
    
    ind.vars <- setdiff(values$dat.variables, input$linRegress.dependent2)
    
    # Render select input for variables
    selectInput("linRegress.independent2", "Independent",
                choices  = ind.vars,
                selected = ind.vars[1L],
                multiple = TRUE)
  })
  
  output$linRegress.formula2 <- renderUI({
    if (is.null(input$linRegress.dependent)) {
      return("")
    }
    
    formula.str <- paste(input$linRegress.dependent2, " ~ ")
    
    ind.vars <- input$linRegress.independent2
    
    n <- length(ind.vars)
    
    if (n > 1) {
      for (i in 1:(n - 1)) {
        formula.str <- paste(formula.str, ind.vars[i], " + ")
      }
    }
    
    formula.str <- paste(formula.str, ind.vars[n])
    
    textInput("linRegress.formula.text2", "Formula",
              formula.str)
  })
  
  output$linRegress.robust.control <- renderUI({
    if (input$linRegress.fit.option != "LS") {
      if (input$linRegress.second.method == TRUE) {
        req(input$linRegress.fit.option2)
        
        if(input$linRegress.fit.option2 != "LS") {
          tabPanel("",
            tags$hr(),
            
            h4("Robust Controls 1"),
            
            selectInput("linRegress.family", "Family",
                        choices = c("Bisquare" = "bisquare",
                                    "Opt."      = "optimal",
                                    "Mod. Opt." = "modopt"),
                        selected = "optimal"),
            
            numericInput("linRegress.eff", "Efficiency", value = 0.99, min = 0.80, max = 0.99, step = 0.01)
          )
        } else {
          tabPanel("",
            tags$hr(),
            
            h4("Robust Controls"),
            
            selectInput("linRegress.family", "Family",
                        choices = c("Bi-square" = "bisquare",
                                    "Opt."      = "optimal",
                                    "Mod. Opt." = "modopt"),
                        selected = "optimal"),
            
            numericInput("linRegress.eff", "Efficiency", value = 0.99, min = 0.80, max = 0.99, step = 0.01)
          )
        }
      } else {
        tabPanel("",
          tags$hr(),
          
          h4("Robust Controls"),
          
          selectInput("linRegress.family", "Family",
                      choices = c("Bi-square" = "bisquare",
                                  "Opt."      = "optimal",
                                  "Mod. Opt." = "modopt"),
                      selected = "optimal"),
          
          numericInput("linRegress.eff", "Efficiency", value = 0.99, min = 0.80, max = 0.99, step = 0.01)
        )
      }
    }
  })
  
  output$linRegress.robust.control2 <- renderUI({
    if (input$linRegress.fit.option2 != "LS") {
      tabPanel("",
        tags$hr(),
        
        h4("Robust Controls 2"),
        
        selectInput("linRegress.family2", "Family",
                    choices = c("Bi-square" = "bisquare",
                                "Opt."      = "optimal",
                                "Mod. Opt." = "modopt"),
                    selected = "optimal"),
        
        numericInput("linRegress.eff2", "Efficiency", value = 0.99, min = 0.80, max = 0.99, step = 0.01)
      )
    }
  })
  
  observeEvent(input$linRegress.display, {
    if (is.null(values$dat)) {
      output$linRegress.results <- renderText({
        return(paste0("<font color=\"#FF0000\"><b>", "ERROR: No data loaded!", "</b><font>"))
      })
    } else if (any(is.null(input$linRegress.independent) || is.null(input$linRegress.dependent))) {
      output$linRegress.results <- renderText({
        return(paste0("<font color=\"#FF0000\"><b>", "ERROR: Missing independent variables! Please add at least one.", "</b><font>"))
      })
    } else if (!is.numeric(values$dat[, input$linRegress.dependent])) {
      output$linRegress.results <- renderText({
        return(paste0("<font color=\"#FF0000\"><b>", invalid_response(), "</b><font>"))
      })
    }
      
    if (input$linRegress.second.method) {
      
      methods <- c(input$linRegress.fit.option, input$linRegress.fit.option2)
      
      index <- match(methods, values$linRegress.methods)
    
      n <- length(index)
    
      model <- sapply(index,
                      function(i, m) {
                        m[i]
                      },
                      values$linRegress.functions)
    } else {
      model <- values$linRegress.functions[match(input$linRegress.fit.option, values$linRegress.methods)]
    }
    
    fit <- vector(mode = "list", length = length(index))
    if (model[1] == "lm") {
      fit[[1]] <- do.call(model[1], list(as.formula(input$linRegress.formula.text), data = values$dat))
    } else {
      control <- lmrobdet.control(efficiency = input$linRegress.eff,
                                  family = input$linRegress.family,
                                  compute.rd = T)
    
      fit[[1]] <- do.call(model[1], list(as.formula(input$linRegress.formula.text),
                                         data    = values$dat,
                                         control = control))
    }
    
    fit[[1]]$call <- call(model[1], as.formula(input$linRegress.formula.text))
    
    if (input$linRegress.second.method) {
      if (model[2] == "lm") {
        fit[[2]] <- do.call(model[2], list(as.formula(input$linRegress.formula.text2), data = values$dat))
      } else {
        control <- lmrobdet.control(efficiency = input$linRegress.eff2,
                                    family = input$linRegress.family2,
                                    compute.rd = T)
      
        fit[[2]] <- do.call(model[2], list(as.formula(input$linRegress.formula.text2),
                                           data    = values$dat,
                                           control = control))
      }
        
      fit[[2]]$call <- call(model[2], as.formula(input$linRegress.formula.text2))
      
      
      if (model[2] == "lm" && model[1] != "lm") {
        model <- model[2:1]
        
        methods <- methods[2:1]
        
        fm <- fit.models(fit[[2]], fit[[1]])
      } else {
        fm <- fit.models(fit[[1]], fit[[2]])
      }
    } else {
      fm <- fit[[1]]
    }

    values$linRegress.active <- T
    
    if (input$linRegress.second.method) {
      if (input$linRegress.fit.option == input$linRegress.fit.option2) {
        values$linRegress.models <- c(paste(input$linRegress.fit.option, "1"), paste(input$linRegress.fit.option[1], "2"))
      } else {
        values$linRegress.models <- methods
      }
      names(fm) <- values$linRegress.models
      
      values$linRegress.fit <- fm[[1]]
      
      values$linRegress.fit2 <- fm[[2]]
      
      values$linRegress.num.fits <- 2L
    } else {
      values$linRegress.fit <- fm
      
      values$linRegress.models <- input$linRegress.fit.option
      values$linRegress.num.fits <- 1L
    }
    
    values$linRegress.fm <- fm
    
      
    output$linRegress.results <- renderText({
      return(htmlText(summary(values$linRegress.fm)))
    })
  })
  
  invalid_response <- eventReactive(input$linRegress.display, {
    return(paste("ERROR: Response variable is of",
                 class(values$dat[, input$linRegress.dependent]),
                 "type. Please select a response variable with numeric values"))
  })
  
  observeEvent(input$linRegress.display.plots, {
    values$linRegress.plots.active <- T
    
    plots <- vector(mode = "list")
    
    i <- 0L
    
    fit <- values$linRegress.fit
    
    fit.vals <- fitted(fit)
    
    # Residual v. fitted values
    if (input$linRegress.residual.fit == T) {
      i <- i + 1L
      
      dat <- data.frame(X = fit.vals, Y = fit$residuals)
      
      sigma <- 1L
      
      if (any(class(fit) == "lm")) {
        sigma <- sd(fit$residuals)
      } else {
        sigma <- fit$scale
      }
      
      title.name <- ifelse(any(class(fit) == "lm"), values$linRegress.models[1], paste("Robust", values$linRegress.models[1]))
      
      plots[[i]] <- ggplot(data = dat, aes(x = X, y = Y)) +
                      ggtitle(title.name) +
                      xlab("Fitted Values") +
                      theme(aspect.ratio = 3/4) +
                      ylab("Residuals") +
                      geom_point(color = "dodgerblue2", shape = 18, size = 2.5) +
                      geom_hline(yintercept = c(-2.5 * sigma, 0, 2.5 * sigma),
                                 linetype = 2L)
      
      if (input$include.rugplot == T) {
        plots[[i]] <- plots[[i]] + geom_rug()
      }
    }
    
    # Response v. fitted values
    if (input$linRegress.response.fit == T) {
      i <- i + 1L
      
      response <- fit$model[, 1]

      dat <- data.frame(X = fit.vals, Y = response)
      
      title.name <- ifelse(any(class(fit) == "lm"), values$linRegress.models[1], paste("Robust", values$linRegress.models[1]))

      plots[[i]] <- ggplot(data = dat, aes(x = X, y = Y)) +
                      ggtitle(title.name) +
                      theme(aspect.ratio = 3/4) +
                      xlab("Fitted Values") +
                      ylab("Response") +
                      geom_point(color = "dodgerblue2", shape = 18, size = 2.5)

      if (input$include.rugplot == T) {
        plots[[i]] <- plots[[i]] + geom_rug()
      }
    }

    # QQ Plot
    if (input$linRegress.qq == T) {
      i <- i + 1L

      dat <- data.frame(Res = sort(fit$residuals))
      
      title.name <- ifelse(any(class(fit) == "lm"), values$linRegress.models[1], paste("Robust", values$linRegress.models[1]))

      # Calculate slope and intercept for qqline
      
      y     <- quantile(fit$residuals, c(0.25, 0.75), type = 5)
      x     <- qnorm(c(0.25, 0.75))
      slope <- diff(y) / diff(x)
      int   <- y[1] - slope * x[1]
      
      if (input$linRegress.qq.env == T) {
        confidence.level <- 0.95
        n <- length(fit$residuals)
        P <- ppoints(n)
        z <- qnorm(P)
        
        zz <- qnorm(1 - (1 - confidence.level) / 2)
        SE <- (slope / dnorm(z)) * sqrt(P * (1 - P) / n)
        fit.vals <- int + slope * z
        dat$z <- z
        dat$lower <- fit.vals - zz * SE
        dat$upper <- fit.vals + zz * SE
      }

      # Normal QQ plot
      
      if (input$linRegress.qq.env == T) {
        plots[[i]] <- ggplot(data = dat, aes(x = z, y = Res)) +
                        ggtitle(title.name) +
                        xlab("Normal Quantiles") +
                        ylab("Ordered Residuals") +
                        theme(aspect.ratio = 3/4) +
                        geom_point(color = "dodgerblue2", shape = 18, size = 2.5) +
                        geom_abline(slope = slope, intercept = int) +
                        geom_ribbon(aes(ymin = lower, ymax = upper),
                                    alpha = 0.2, color = "dodgerblue2", fill = "dodgerblue2")
      } else {
        plots[[i]] <- ggplot(data = dat, aes(sample = Res)) +
                        ggtitle(title.name) +
                        xlab("Normal Quantiles") +
                        ylab("Ordered Residuals") +
                        theme(aspect.ratio = 3/4) +
                        geom_qq() +
                        geom_abline(slope = slope, intercept = int)
      }
    }
    
    # Standardized residuals vs. robust distances
    if (input$linRegress.resid.dist == T) {
      i <- i + 1
      
      title.name <- ifelse(any(class(fit) == "lm"), values$linRegress.models[1], paste("Robust", values$linRegress.models[1]))
      
      if (any(class(fit) == "lm")) {
        st.residuals <- rstandard(fit)
        
        tr <- terms(fit)
        
        mat <- model.matrix(fit)
        
        if(attr(tr, "intercept") == 1L) {
          mat <- mat[, -1, drop = F]
        }
        
        MD <- sqrt(mahalanobis(x      = mat,
                               center = colMeans(mat),
                               cov    = var(mat)))
        
        chi <- sqrt(qchisq(p = 1 - 0.025, df = fit$rank))
              
        dat <- data.frame(X = MD, Y = st.residuals)

        plots[[i]] <- ggplot(data = dat, aes(x = X, y = Y)) +
                         ggtitle(title.name) +
                         xlab("Distances") +
                         ylab("Standardized Residuals") +
                         theme(aspect.ratio = 3/4) +
                         geom_point(color = "dodgerblue2", shape = 18, size = 2.5) +
                         geom_hline(yintercept = c(-2.5, 0, 2.5),
                                    linetype = 2) +
                         geom_vline(xintercept = chi,
                                    linetype = 2)
      } else {
        st.residuals <- fit$residuals / fit$scale
        
        MD <- fit$MD
        
        chi <- sqrt(qchisq(p = 1 - 0.025, df = fit$rank))
              
        dat <- data.frame(X = MD, Y = st.residuals)

        plots[[i]] <- ggplot(data = dat, aes(x = X, y = Y)) +
                         ggtitle(title.name) +
                         xlab("Robust Distances") +
                         ylab("Robustly Standardized Residuals") +
                         theme(aspect.ratio = 3/4) +
                         geom_point(color = "dodgerblue2", shape = 18, size = 2.5) +
                         geom_hline(yintercept = c(-2.5, 0, 2.5),
                                    linetype = 2) +
                         geom_vline(xintercept = chi,
                                    linetype = 2)
      }
    }
    
    # Estimated residual density
    if (input$linRegress.residual.density == T) {
      i <- i + 1L
      
      dat <- data.frame(Res = fit$residuals)
      
      title.name <- ifelse(any(class(fit) == "lm"), values$linRegress.models[1], paste("Robust", values$linRegress.models[1]))
      
      plots[[i]] <- ggplot(data = dat) +
                      ggtitle(title.name) +
                      xlab("Residuals") +
                      ylab("Density") +
                      theme(aspect.ratio = 3/4) +
                      geom_histogram(aes(x = Res, y = ..density..),
                                     fill  = 'white',
                                     color = 'black',
                                     bins  = 35) +
                      geom_density(aes(x = Res),
                                   color = 'dodgerblue2',
                                   fill  = 'dodgerblue2',
                                   alpha = 0.1)
      
      if (input$include.rugplot == T) {
        plots[[i]] <- plots[[i]] + geom_rug()
      }
    }
    
    # Standardized residuals vs. index values
    if (input$linRegress.resid.index == T) {
      i <- i + 1
      
      title.name <- ifelse(any(class(fit) == "lm"), values$linRegress.models[1], paste("Robust", values$linRegress.models[1]))
      
      if (any(class(values$dat) == "zoo")) {
        indx <- index(values$dat)
      } else {
        indx <- 1:nrow(values$dat)
      }
      
      if (any(class(fit) == "lm")) {
        st.residuals <- rstandard(fit)
        
        dat <- data.frame(X = indx, Y = st.residuals)
      
        plots[[i]] <- ggplot(data = dat, aes(x = X, y = Y)) +
                        ggtitle(title.name) +
                        xlab("Index") +
                        ylab("Standardized Residuals") +
                        theme(aspect.ratio = 3/4) +
                        geom_point(color = "dodgerblue2", shape = 18, size = 2.5) +
                        geom_line()
      } else {
        st.residuals <- fit$residuals / fit$scale
        
        dat <- data.frame(X = indx, Y = st.residuals)
      
        plots[[i]] <- ggplot(data = dat, aes(x = X, y = Y)) +
                        ggtitle(title.name) +
                        xlab("Index") +
                        ylab("Robustly Standardized Residuals") +
                        theme(aspect.ratio = 3/4) +
                        geom_point(color = "dodgerblue2", shape = 18, size = 2.5) +
                        geom_line()
      }
    }
    
    ## 2 Regression Case
    
    if (values$linRegress.num.fits == 2) {
      j <- 0
      
      fit2 <- values$linRegress.fit2
      
      fit.vals <- fitted(fit2)
      
      # Residual v. fitted values
      if (input$linRegress.residual.fit == T) {
        j <- j + 1
      
        dat <- data.frame(X = fit.vals, Y = fit2$residuals)
        
        title.name <- ifelse(any(class(fit2) == "lm"), values$linRegress.models[2], paste("Robust", values$linRegress.models[2]))
        
        sigma <- 1
        
        if (any(class(fit2) == "lm")) {
          sigma <- sd(fit2$residuals)
        } else {
          sigma <- fit2$scale
        }
        
        plt <- ggplot(data = dat, aes(x = X, y = Y)) +
                 ggtitle(title.name) +
                 xlab("Fitted Values") +
                 ylab("Residuals") +
                 theme(aspect.ratio = 3/4) +
                 geom_point(color = "dodgerblue2", shape = 18, size = 2.5) +
                 geom_hline(yintercept = c(-2.5 * sigma, 0, 2.5 * sigma),
                            linetype = 2)
        
        if (input$include.rugplot == T) {
          plt <- plt + geom_rug()
        }
        
        x.min <- min(layer_scales(plots[[j]])$x$range$range, layer_scales(plt)$x$range$range)
        x.max <- max(layer_scales(plots[[j]])$x$range$range, layer_scales(plt)$x$range$range)

        y.min <- min(layer_scales(plots[[j]])$y$range$range, layer_scales(plt)$y$range$range)
        y.max <- max(layer_scales(plots[[j]])$y$range$range, layer_scales(plt)$y$range$range)
        
        plots[[j]] <- plots[[j]] + scale_y_continuous(limits = c(y.min, y.max)) +
                        scale_x_continuous(limits = c(x.min, x.max))
        plt <- plt + scale_y_continuous(limits = c(y.min, y.max)) +
                 scale_x_continuous(limits = c(x.min, x.max))
        
        plots[[j]] <- cbind(ggplotGrob(plots[[j]]), ggplotGrob(plt), size = "last")
      }
      
      # Response v. fitted values
      if (input$linRegress.response.fit == T) {
        j <- j + 1
  
        dat <- data.frame(X = fit.vals, Y = fit2$model[, 1])
        
        title.name <- ifelse(any(class(fit2) == "lm"), values$linRegress.models[2], paste("Robust", values$linRegress.models[2]))
  
        plt <- ggplot(data = dat, aes(x = X, y = Y)) +
                 ggtitle(title.name) +
                 xlab("Fitted Values") +
                 ylab("Response") +
                 theme(aspect.ratio = 3/4) +
                 geom_point(color = "dodgerblue2", shape = 18, size = 2.5)
  
        if (input$include.rugplot == T) {
          plt <- plt + geom_rug()
        }
        
        x.min <- min(layer_scales(plots[[j]])$x$range$range, layer_scales(plt)$x$range$range)
        x.max <- max(layer_scales(plots[[j]])$x$range$range, layer_scales(plt)$x$range$range)

        y.min <- min(layer_scales(plots[[j]])$y$range$range, layer_scales(plt)$y$range$range)
        y.max <- max(layer_scales(plots[[j]])$y$range$range, layer_scales(plt)$y$range$range)
        
        plots[[j]]  <- plots[[j]] + scale_y_continuous(limits = c(y.min, y.max),
                                                       expand = expand_scale(mult = c(0.05, 0.05))) +
                         scale_x_continuous(limits = c(x.min, x.max))
        
        plt <- plt + scale_y_continuous(limits = c(y.min, y.max),
                                        expand = expand_scale(mult = c(0.05, 0.05))) +
                 scale_x_continuous(limits = c(x.min, x.max))
        
        plots[[j]] <- cbind(ggplotGrob(plots[[j]]), ggplotGrob(plt), size = "last")
      }
  
      # QQ Plot
      if (input$linRegress.qq == T) {
        j <- j + 1
  
        dat <- data.frame(Res = sort(fit2$residuals))
        
        title.name <- ifelse(any(class(fit2) == "lm"), values$linRegress.models[2], paste("Robust", values$linRegress.models[2]))
  
        # Calculate slope and intercept for qqline 
        y     <- quantile(fit2$residuals, c(0.25, 0.75), type=5)
        x     <- qnorm(c(0.25, 0.75))
        slope <- diff(y) / diff(x)
        int   <- y[1] - slope * x[1]
        
        if (input$linRegress.qq.env == T) {
          confidence.level <- 0.95
          n <- length(fit2$residuals)
          P <- ppoints(n)
          z <- qnorm(P)
          
          zz <- qnorm(1 - (1 - confidence.level) / 2)
          SE <- (slope / dnorm(z)) * sqrt(P * (1 - P) / n)
          fit.vals <- int + slope * z
          dat$z <- z
          dat$lower <- fit.vals - zz * SE
          dat$upper <- fit.vals + zz * SE
          
         plt <- ggplot(data = dat, aes(x = z, y = Res)) +
                  ggtitle(title.name) +
                  xlab("Normal Quantiles") +
                  ylab("Ordered Residuals") +
                  theme(aspect.ratio = 3/4) +
                  geom_point(color = "dodgerblue2", shape = 18, size = 2.5) +
                  geom_abline(slope = slope, intercept = int) +
                  geom_ribbon(aes(ymin = lower, ymax = upper),
                              alpha = 0.2, color = "dodgerblue2", fill = "dodgerblue2")
        } else {
          plt <- ggplot(data = dat, aes(x = z, y = Res)) +
                   ggtitle(title.name) +
                   xlab("Normal Quantiles") +
                   ylab("Ordered Residuals") +
                   theme(aspect.ratio = 3/4) +
                   geom_qq() +
                   geom_abline(slope = slope, intercept = int)
        }
        
        x.min <- min(layer_scales(plots[[j]])$x$range$range, layer_scales(plt)$x$range$range)
        x.max <- max(layer_scales(plots[[j]])$x$range$range, layer_scales(plt)$x$range$range)

        y.min <- min(layer_scales(plots[[j]])$y$range$range, layer_scales(plt)$y$range$range)
        y.max <- max(layer_scales(plots[[j]])$y$range$range, layer_scales(plt)$y$range$range)
        
        plots[[j]]  <- plots[[j]] + scale_y_continuous(limits = c(y.min, y.max),
                                                       expand = expand_scale(mult = c(0.05, 0.05))) +
                        scale_x_continuous(limits = c(x.min, x.max))
        
        plt <- plt + scale_y_continuous(limits = c(y.min, y.max),
                                        expand = expand_scale(mult = c(0.05, 0.05))) +
                 scale_x_continuous(limits = c(x.min, x.max))
        
        plots[[j]] <- cbind(ggplotGrob(plots[[j]]), ggplotGrob(plt), size = "last")
      }
      
      # Standardized residuals vs. robust distances
      if (input$linRegress.resid.dist == T) {
        j <- j + 1
        
        title.name <- ifelse(any(class(fit2) == "lm"), values$linRegress.models[2], paste("Robust", values$linRegress.models[2]))
        
        if (any(class(fit2) == "lm")) {
          st.residuals <- rstandard(fit2)
          
          tr <- terms(fit2)
        
          mat <- model.matrix(fit2)
          
          if(attr(tr, "intercept") == 1) {
            mat <- mat[, -1, drop = F]
          }
          
          MD <- sqrt(mahalanobis(x      = mat,
                                 center = colMeans(mat),
                                 cov    = var(mat)))
          
          chi <- sqrt(qchisq(p = 1 - 0.025, df = fit2$rank))
                
          dat <- data.frame(X = MD, Y = st.residuals)
  
          plt <- ggplot(data = dat, aes(x = X, y = Y)) +
                   ggtitle(title.name) +
                   xlab("Distances") +
                   ylab("Standardized Residuals") +
                   theme(aspect.ratio = 3/4) +
                   geom_point(color = "dodgerblue2", shape = 18, size = 2.5) +
                   geom_hline(yintercept = c(-2.5, 0, 2.5),
                              linetype = 2) + 
                   geom_vline(xintercept = chi,
                              linetype = 2)
        } else {
          st.residuals <- fit2$residuals / fit2$scale
          
          MD <- fit2$MD
          
          chi <- sqrt(qchisq(p = 1 - 0.025, df = fit2$rank))
                
          dat <- data.frame(X = MD, Y = st.residuals)
  
          plt <- ggplot(data = dat, aes(x = X, y = Y)) +
                   ggtitle(title.name) +
                   xlab("Robust Distances") +
                   ylab("Robustly Standardized Residuals") +
                   theme(aspect.ratio = 3/4) +
                   geom_point(color = "dodgerblue2", shape = 18, size = 2.5) +
                   geom_hline(yintercept = c(-2.5, 0, 2.5),
                              linetype = 2) + 
                   geom_vline(xintercept = chi,
                              linetype = 2)
        }
        
        x.min <- min(layer_scales(plots[[j]])$x$range$range, layer_scales(plt)$x$range$range)
        x.max <- max(layer_scales(plots[[j]])$x$range$range, layer_scales(plt)$x$range$range)

        y.min <- min(layer_scales(plots[[j]])$y$range$range, layer_scales(plt)$y$range$range)
        y.max <- max(layer_scales(plots[[j]])$y$range$range, layer_scales(plt)$y$range$range)
        
        plots[[j]] <- plots[[j]] + scale_y_continuous(limits = c(y.min, y.max),
                                                       expand = expand_scale(mult = c(0.05, 0.05))) +
                        scale_x_continuous(limits = c(x.min, x.max))
        
        plt <- plt + scale_y_continuous(limits = c(y.min, y.max),
                                        expand = expand_scale(mult = c(0.05, 0.05))) +
                 scale_x_continuous(limits = c(x.min, x.max))
        
        plots[[j]] <- cbind(ggplotGrob(plots[[j]]), ggplotGrob(plt), size = "last")
      }
      
      # Estimated residual density
      if (input$linRegress.residual.density == T) {
        j <- j + 1
        
        dat <- data.frame(Res = fit2$residuals)
        
        title.name <- ifelse(any(class(fit2) == "lm"), values$linRegress.models[2], paste("Robust", values$linRegress.models[2]))
        
        plt <- ggplot(data = dat) +
                 ggtitle(title.name) +
                 xlab("Residuals") +
                 ylab("Density") +
                 theme(aspect.ratio = 3/4) +
                 geom_histogram(aes(x = Res, y = ..density..),
                                fill  = 'white',
                                color = 'black',
                                bins  = 35) +
                 geom_density(aes(x = Res),
                              color = 'dodgerblue2',
                              fill  = 'dodgerblue2',
                              alpha = 0.1)
        
        if (input$include.rugplot == T) {
          plt <- plt + geom_rug()
        }
        
        x.min <- min(layer_scales(plots[[j]])$x$range$range, layer_scales(plt)$x$range$range)
        x.max <- max(layer_scales(plots[[j]])$x$range$range, layer_scales(plt)$x$range$range)

        y.min <- min(layer_scales(plots[[j]])$y$range$range, layer_scales(plt)$y$range$range)
        y.max <- max(layer_scales(plots[[j]])$y$range$range, layer_scales(plt)$y$range$range)
        
        plots[[j]]  <- plots[[j]] + scale_y_continuous(limits = c(y.min, y.max),
                                                       expand = expand_scale(mult = c(0.01, 0.05))) +
                        scale_x_continuous(limits = c(x.min, x.max))
        
        plt <- plt + scale_y_continuous(limits = c(y.min, y.max),
                                        expand = expand_scale(mult = c(0.01, 0.05))) +
                 scale_x_continuous(limits = c(x.min, x.max))
        
        plots[[j]] <- cbind(ggplotGrob(plots[[j]]), ggplotGrob(plt), size = "last")
      }
      
      # Standardized residuals vs. index values
      if (input$linRegress.resid.index == T) {
        j <- j + 1
        
        title.name <- ifelse(any(class(fit2) == "lm"), values$linRegress.models[2], paste("Robust", values$linRegress.models[2]))
        
        indx <- NULL
        
        if (any(class(values$dat) == "zoo")) {
          indx <- index(values$dat)
        } else {
          indx <- 1:nrow(values$dat)
        }
        
        if (any(class(fit2) == "lm")) {
          st.residuals <- rstandard(fit2)
          
          dat <- data.frame(X = indx, Y = st.residuals)
        
          plt <- ggplot(data = dat, aes(x = X, y = Y)) + 
                   ggtitle(title.name) +
                   xlab("Index") +
                   ylab("Standardized Residuals") +
                   theme(aspect.ratio = 3/4) +
                   geom_point(color = "dodgerblue2", shape = 18, size = 2.5) +
                   geom_line()
        } else {
          st.residuals <- fit2$residuals / fit2$scale
          
          dat <- data.frame(X = indx, Y = st.residuals)
        
          plt <- ggplot(data = dat, aes(x = X, y = Y)) + 
                   ggtitle(title.name) +
                   xlab("Index") +
                   ylab("Robustly Standardized Residuals") +
                   theme(aspect.ratio = 3/4) +
                   geom_point(color = "dodgerblue2", shape = 18, size = 2.5) +
                   geom_line()
        }
        
        # x.min <- min(layer_scales(plots[[j]])$x$range$range, layer_scales(plt)$x$range$range)
        # x.max <- max(layer_scales(plots[[j]])$x$range$range, layer_scales(plt)$x$range$range)

        y.min <- min(layer_scales(plots[[j]])$y$range$range, layer_scales(plt)$y$range$range)
        y.max <- max(layer_scales(plots[[j]])$y$range$range, layer_scales(plt)$y$range$range)
        
        plots[[j]]  <- plots[[j]] + scale_y_continuous(limits = c(y.min, y.max),
                                                       expand = expand_scale(mult = c(0.05, 0.05)))
                        # scale_x_continuous(limits = c(x.min, x.max),
                        #                               expand = c(0, 0))
        
        plt <- plt + scale_y_continuous(limits = c(y.min, y.max),
                                        expand = expand_scale(mult = c(0.05, 0.05)))
                        # scale_x_continuous(limits = c(x.min, x.max),
                        #                               expand = c(0, 0))
        
        plots[[j]] <- cbind(ggplotGrob(plots[[j]]), ggplotGrob(plt), size = "last")
      }
      
      if (input$linRegress.overlaid.scatter == T) {
        j <- j + 1
        mat <- as.data.frame(fit$model)
        
        names <- colnames(mat)
        
        LR.names <- ifelse(values$linRegress.models == "LS", values$linRegress.models, paste("Robust", values$linRegress.models))
        
        colnames(mat) <- c("Y", "X")
        
        int1 <- fit$coefficients[1]
        int2 <- fit2$coefficients[1]
        
        slope1 <- fit$coefficients[2]
        slope2 <- fit2$coefficients[2]
        
        line.dat <- data.frame(slope = c(slope1, slope2), int = c(int1, int2), method = LR.names)
        
        title.name <- input$linRegress.formula.text
        
        plt <- ggplot() +
                 ggtitle(title.name) +
                 xlab(names[2]) +
                 ylab(names[1]) +
                 theme(aspect.ratio = 3/4) +
                 geom_point(data = mat, aes(x = X, y = Y), color = "dodgerblue2", shape = 18, size = 2.5) +
                 geom_abline(data = line.dat, aes(slope = slope, intercept = int, linetype = method, color = method), size = 1, show.legend = T) +
                 scale_color_manual(values = c("red", "black")) +
                 scale_linetype_manual(values = c("dashed", "solid"))
        
        plots[[j]] <- ggplotGrob(plt)
      }
      
      if (j == 0) {
        output$linRegress.plot.ui <- renderUI({
          fluidPage(
            fluidRow(verbatimTextOutput("no.selection"))
          )
        })
      } else {
        values$linRegress.plots        <- plots
        values$linRegress.num.plots    <- j
        values$linRegress.active.index <- 1
        values$linRegress.active.plot  <- plots[[1]]
        
        if (j > 1) {
          output$linRegress.plot.ui <- renderUI({
            fluidPage(
              wellPanel(
                plotOutput("linRegress.plot.output")
              ),
              
              fluidRow(
                column(1,
                       offset = 10,
                       actionButton("linRegress.next.plot",
                                    "",
                                    icon = icon("angle-right", "fa-2x"))
                )
              )
            )
          })
        } else {
          output$linRegress.plot.ui <- renderUI({
            fluidPage(
              wellPanel(
                plotOutput("linRegress.plot.output")
              )
            )
          })
        }
        
        output$linRegress.plot.output <- renderPlot({
          grid.draw(values$linRegress.active.plot)
        })
      }
    } else {
      if (i == 0) {
        output$linRegress.plot.ui <- renderUI({
          fluidPage(verbatimTextOutput("no.selection"))
        })
      } else {
        values$linRegress.plots        <- plots
        values$linRegress.num.plots    <- i
        values$linRegress.active.index <- 1
        values$linRegress.active.plot  <- plots[[1]]
        
        if (i > 1) {
          output$linRegress.plot.ui <- renderUI({
            fluidPage(
              wellPanel(
                plotOutput("linRegress.plot.output")
              ),
              
              fluidRow(
                column(1,
                       offset = 10,
                       actionButton("linRegress.next.plot",
                                    "",
                                    icon = icon("angle-right", "fa-2x"))
                )
              )
            )
          })
        } else {
          output$linRegress.plot.ui <- renderUI({
            fluidPage(
              wellPanel(
                plotOutput("linRegress.plot.output")
              )
            )
          })
        }
        
        output$linRegress.plot.output <- renderPlot({
          values$linRegress.active.plot
        })
      }
    }
  })
  
  output$no.selection <- renderPrint({
    cat("No plots selected.")
  })
  
  # On button press, show next plot(s)
  observeEvent(input$linRegress.next.plot, {
    if (values$linRegress.num.plots > 0) {
      if (values$linRegress.num.fits == 2) {
        values$linRegress.active.index <- values$linRegress.active.index %% values$linRegress.num.plots + 1
        values$linRegress.active.plot  <- values$linRegress.plots[[values$linRegress.active.index]]
        
        output$linRegress.plot.ui <- renderUI({
          if (values$linRegress.active.index == values$linRegress.num.plots) {
            fluidPage(
              wellPanel(
                plotOutput("linRegress.plot.output")
              ),
              
              fluidRow(
                column(1,
                       offset = 1,
                       actionButton("linRegress.prev.plot",
                                    "",
                                    icon = icon("angle-left", "fa-2x"))
                )
              )
            )
          } else {
            fluidPage(
              wellPanel(
                plotOutput("linRegress.plot.output")
              ),
              
              fluidRow(
                column(1,
                       offset = 1,
                       actionButton("linRegress.prev.plot",
                                    "",
                                    icon = icon("angle-left", "fa-2x"))
                ),
                
                column(1,
                       offset = 8,
                       actionButton("linRegress.next.plot",
                                    "",
                                    icon = icon("angle-right", "fa-2x"))
                )
              )
            )
          }
        })
        
        output$linRegress.plot.output <- renderPlot({
          grid.draw(values$linRegress.active.plot)
        })
      } else {
        
        values$linRegress.active.index <- values$linRegress.active.index %% values$linRegress.num.plots + 1
      
        values$linRegress.active.plot <- values$linRegress.plots[[values$linRegress.active.index]]
        
        output$linRegress.plot.ui <- renderUI({
          if (values$linRegress.active.index == values$linRegress.num.plots) {
            fluidPage(
              wellPanel(
                plotOutput("linRegress.plot.output")
              ),
            
              fluidRow(
                column(1,
                       offset = 1,
                       actionButton("linRegress.prev.plot",
                                    "",
                                    icon = icon("angle-left", "fa-2x"))
                )
              )
            )
          } else {
            fluidPage(
              wellPanel(
                plotOutput("linRegress.plot.output")
              ),
            
              fluidRow(
                column(1,
                       offset = 1,
                       actionButton("linRegress.prev.plot",
                                    "",
                                    icon = icon("angle-left", "fa-2x"))
                ),
                
                column(1,
                       offset = 8,
                       actionButton("linRegress.next.plot",
                                    "",
                                    icon = icon("angle-right", "fa-2x"))
                )
              )
            )
          }
        })
      
        output$linRegress.plot.output <- renderPlot({
          values$linRegress.active.plot
        })
      }
    }
  })
  
  # On button press, move to previous plot(s)
  observeEvent(input$linRegress.prev.plot, {
    if (values$linRegress.num.plots > 0) {
      if (values$linRegress.num.fits == 2) {
        values$linRegress.active.index <- values$linRegress.num.plots + (values$linRegress.active.index - 1) %% (-values$linRegress.num.plots)
        values$linRegress.active.plot  <- values$linRegress.plots[[values$linRegress.active.index]]
        
        output$linRegress.plot.ui <- renderUI({
          if (values$linRegress.active.index == 1) {
            fluidPage(
              wellPanel(
                plotOutput("linRegress.plot.output")
              ),
              
              fluidRow(
                column(1,
                       offset = 10,
                       actionButton("linRegress.next.plot",
                                    "",
                                    icon = icon("angle-right", "fa-2x"))
                )
              )
            )
          } else {
            fluidPage(
              wellPanel(
                plotOutput("linRegress.plot.output")
              ),
              
              fluidRow(
                column(1,
                       offset = 1,
                       actionButton("linRegress.prev.plot",
                                    "",
                                    icon = icon("angle-left", "fa-2x"))
                ),
                
                column(1,
                       offset = 8,
                       actionButton("linRegress.next.plot",
                                    "",
                                    icon = icon("angle-right", "fa-2x"))
                )
              )
            )
          }
        })
        
        output$linRegress.plot.output <- renderPlot({
          grid.draw(values$linRegress.active.plot)
        })
      } else {
        values$linRegress.active.index <- values$linRegress.num.plots + (values$linRegress.active.index - 1) %% (-values$linRegress.num.plots)
      
        values$linRegress.active.plot <- values$linRegress.plots[[values$linRegress.active.index]]
        
        output$linRegress.plot.ui <- renderUI({
          if (values$linRegress.active.index == 1) {
            fluidPage(
              wellPanel(
                plotOutput("linRegress.plot.output")
              ),
            
              fluidRow(
                column(1,
                       offset = 10,
                       actionButton("linRegress.next.plot",
                                    "",
                                    icon = icon("angle-right", "fa-2x"))
                )
              )
            )
          } else {
            fluidPage(
              wellPanel(
                plotOutput("linRegress.plot.output")
              ),
            
              fluidRow(
                column(1,
                       offset = 1,
                       actionButton("linRegress.prev.plot",
                                    "",
                                    icon = icon("angle-left", "fa-2x"))
                ),
                
                column(1,
                       offset = 8,
                       actionButton("linRegress.next.plot",
                                    "",
                                    icon = icon("angle-right", "fa-2x"))
                )
              )
            )
          }
        })
      
        output$linRegress.plot.output <- renderPlot({
          values$linRegress.active.plot
        })
      }
    }
  })
  
################
## Covariance ##
################
  
  output$covariance.select.variables <- renderUI({
    if (is.null(dim(values$dat))) {
      return("")
    }
    
    selectInput("covariance.variables", "Variables",
                choices  = values$dat.numeric.variables,
                selected = values$dat.numeric.variables,
                multiple = TRUE)
  })
  
  observeEvent(input$covariance.display, {
    if (is.null(values$dat)) {
      output$covariance.results <- renderText({
        return(paste0("<font color=\"#FF0000\"><b>", "ERROR: No data loaded!", "</b><font>"))
      })
    } else if (is.null(input$covariance.variables)) {
      output$covariance.results <- renderText({
        return(paste0("<font color=\"#FF0000\"><b>", "ERROR: Missing variables! Please input at least 2 variables.", "</b><font>"))
      })
    } else if (length(input$covariance.variables) == 1) {
      output$covariance.results <- renderText({
        return(paste0("<font color=\"#FF0000\"><b>", "ERROR: Insufficient number of variables! Please add more.", "</b><font>"))
      })
    } else {
      values$covariance.active <- TRUE
    
      corr <- FALSE
      
      if (input$covariance.type == "corr") {
        corr <- TRUE
      }
      
      data <- values$dat[, input$covariance.variables]
      
      if (input$covariance.method == "classic") {
        values$covariance.num <- 1
        
        
        values$covariance.fit <- covClassic(data,
                                            data.name = input$dataset,
                                            corr = corr)
        
      } else if (input$covariance.method == "rob") {
        values$covariance.num <- 1
        
        if (input$covariance.estimator == "MM") {
          values$covariance.fit <- covRobMM(data,
                                                 data.name = input$dataset,
                                                 corr = corr)
        } else if (input$covariance.estimator == "Rocke") {
          values$covariance.fit <- covRobRocke(data,
                                               data.name = input$dataset,
                                               corr = corr)
        } else {
          values$covariance.fit <- covRob(data,
                                          data.name = input$dataset,
                                          corr = corr)
        }
      } else {
        values$covariance.num <- 2
        
        if (input$covariance.estimator == "MM") {
          values$covariance.fit <- fit.models(c(Classic = "covClassic", Robust = "covRobMM"),
                                              data = data,
                                              data.name = input$dataset,
                                              corr = corr)
        } else if (input$covariance.estimator == "Rocke") {
          values$covariance.fit <- fit.models(c(Classic = "covClassic", Robust = "covRobRocke"),
                                              data = data,
                                              data.name = input$dataset,
                                              corr = corr)
        } else {
          values$covariance.fit <- fit.models(c(Classic = "covClassic", Robust = "covRob"),
                                              data = data,
                                              data.name = input$dataset,
                                              corr = corr)
        }
      }
      
      # Print summary method for covariance results
      output$covariance.results <- renderText({
        return(htmlText(summary(values$covariance.fit)))
      })
    }
  })
  
  observeEvent(input$covariance.display.plots, {
    values$covariance.plots.active <- TRUE
    
    i <- 0
    
    data <- values$dat[, input$covariance.variables]
    
    plots <- vector(mode = "list")
    
    if (values$covariance.num == 2) {
      fm <- values$covariance.fit
      
      corr <- fm[[1]]$corr
      
      if (input$covariance.eigen == T) {
        i <- i + 1
        
        eigen.vals <- if (corr) {
                        c(eigen(fm[[1]]$cov * (fm[[1]]$sdev %o% fm[[1]]$sdev))$values,
                          eigen(fm[[2]]$cov * (fm[[2]]$sdev %o% fm[[2]]$sdev))$values)
                      } else {
                        c(eigen(fm[[1]]$cov)$values, eigen(fm[[2]]$cov)$values)
                      }
        
        n <- length(eigen.vals) / 2
        
        indx <- rep(as.integer(1:n), 2)
        
        Method <- c(rep("Classic", n), rep("Robust", n))
        
        dat <- data.frame(indx, eigen.vals, Method)
        
        plt <- ggplot(data = dat) +
                 ggtitle("Eigenvalues") +
                 ylab("Eigenvalue") +
                 theme(aspect.ratio = 3/4) +
                 geom_line(aes(x = indx, y = eigen.vals, color = Method, linetype = Method), size = 1) +
                 geom_point(aes(x = indx, y = eigen.vals, color = Method), shape = 16, size = 2.5) +
                 scale_color_manual(values = c("red", "black")) +
                 scale_linetype_manual(values = c("dashed", "solid")) +
                 scale_x_continuous(name = "Factor Number", breaks = min(dat$indx):max(dat$indx))
       
        plots[[i]] <- ggplotGrob(plt)
      }
      
      if (input$covariance.mahalanobis == T) {
        i <- i + 1
        
        md <- if (corr) {
                lapply(values$covariance.fit,
                       function(x, mat) {
                         sqrt(mahalanobis(mat, x$center, x$cov * (x$sdev %o% x$sdev)))
                       },
                       data)
              } else {
                lapply(values$covariance.fit,
                       function(x, mat) {
                         sqrt(mahalanobis(mat, x$center, x$cov))
                       },
                       data)
              }
        
        N <- sapply(md, length)
        p <- sapply(values$covariance.fit, function(x) nrow(x$cov))
        
        thresh <- sqrt(qchisq(0.95, df = p))
        
        indx <- lapply(N, function(n) 1:n)
        
        dat = data.frame(indx = indx[[1]], MD = md[[1]])
        
        p1 <- ggplot(data = dat, aes(x = indx, y = MD)) +
                ggtitle("Classic") +
                xlab("Index") +
                ylab("Mahalanobis Distance") +
                geom_point(color = "dodgerblue2", shape = 18, size = 2.5) +
                geom_hline(yintercept = thresh[1], linetype = 2)
        
        dat = data.frame(indx = indx[[2]], MD = md[[2]])
        
        p2 <- ggplot(data = dat, aes(x = indx, y = MD)) +
                ggtitle("Robust") +
                xlab("Index") +
                ylab("Mahalanobis Distance") +
                geom_point(color = "dodgerblue2", shape = 18, size = 2.5) +
                geom_hline(yintercept = thresh[2], linetype = 2)
        
        y.min <- min(layer_scales(p1)$y$range$range, layer_scales(p2)$y$range$range)
        y.max <- max(layer_scales(p1)$y$range$range, layer_scales(p2)$y$range$range)
        
        p1 <- p1 + scale_y_continuous(limits = c(y.min, y.max),
                                      expand = expand_scale(mult = c(0.05, 0.05)))
        
        p2 <- p2 + scale_y_continuous(limits = c(y.min, y.max),
                                      expand = expand_scale(mult = c(0.05, 0.05)))
        
        plots[[i]] <- cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last")
      }
      
      if (input$covariance.ellipses.matrix == T) {
        i <- i + 1
        
        # Function to create ellipses taken from 'robust' package
        ellipse <- function(loc, A) {
          detA <- A[1, 1] * A[2, 2] - A[1, 2]^2
          dist <- sqrt(qchisq(0.95, 2))
          ylimit <- sqrt(A[2, 2]) * dist
          y <- seq(-ylimit, ylimit, 0.01 * ylimit)
          sqrt.discr <- detA/A[2,2]^2 * (A[2,2] * dist^2 - y^2)
          sqrt.discr[c(1, length(sqrt.discr))] <- 0.0
          sqrt.discr <- sqrt(sqrt.discr)
          b <- loc[1] + A[1, 2] / A[2, 2] * y
          x1 <- b - sqrt.discr
          x2 <- b + sqrt.discr
          y <- loc[2] + y
          rbind(cbind(x1, y), cbind(rev(x2), rev(y)))
        }
        
        fit  <- fm[[1]]
        
        fit2 <- fm[[2]]
        
        grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
        
        text.size <- ifelse(ncol(data) > 4, 4, 10)
        
        # data.frame with xy coordinates
        all <- lapply(1:nrow(grid),
                      function(i, data) {
                        xcol <- grid[i, "x"]
                        ycol <- grid[i, "y"]
                        data.frame(xvar = names(data)[ycol], yvar = names(data)[xcol],
                                   i = grid[i, "x"], j = grid[i, "y"],
                                   x = as.vector(data[, xcol]), y = as.vector(data[, ycol]))
                      },
                      data)
        
        all <- do.call("rbind", all)
        
        all.upper <- all
        
        all.upper[[6]][all.upper[[3]] <  all.upper[[4]]] <- "NA"
        all.upper[[5]][all.upper[[3]] <= all.upper[[4]]] <- "NA"
        
        all.upper$x <- suppressWarnings(as.numeric(as.character(all.upper$x)))
        all.upper$y <- suppressWarnings(as.numeric(as.character(all.upper$y)))
        
        all.upper$xvar <- factor(all.upper$xvar, levels = names(data))
        all.upper$yvar <- factor(all.upper$yvar, levels = names(data))
        
        xy.mapping <- aes_string(x = "x", y = "y")
        class(xy.mapping) <- "uneval"
        
        plt <- ggplot(all.upper, mapping = aes_string(x = "x", y = "y")) +
                 theme(axis.text.x  = element_blank(), axis.text.y  = element_blank(),
                       axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
                       axis.title.x = element_blank(), axis.title.y = element_blank(),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.spacing = unit(0, "lines")) +
                 facet_grid(xvar ~ yvar, scales = "free") +
                 theme(aspect.ratio = 1)
        
        if (fit$corr) {
          X <- fit$cov
        } else {
          s <- sqrt(diag(fit$cov))
          X <- fit$cov / (s %o% s) 
        }
        
        if (fit2$corr) {
          X2 <- fit2$cov
        } else {
          s <- sqrt(diag(fit2$cov))
          X2 <- fit2$cov / (s %o% s)
        }
        
        ellipse.pts <- lapply(1:nrow(grid),
                       function(i) {
                         row <- grid[i, "x"]
                         col <- grid[i, "y"]
                         if (row > col) {
                           cov <- X[row, col]
                           
                           pts <- c(seq(0.0, 2*pi, length.out = 181), NA)
                           xs <- cos(pts + acos(cov) / 2)
                           ys <- cos(pts - acos(cov) / 2)
                           
                           cov2 <- X2[row, col]
                           
                           xs2 <- cos(pts + acos(cov2) / 2)
                           ys2 <- cos(pts - acos(cov2) / 2)
                           
                           rbind(data.frame(xvar = names(data)[col], yvar = names(data)[row],
                                            i = row, j = col,
                                            x = xs,  y = ys,
                                            type = "classic"),
                                 data.frame(xvar = names(data)[col], yvar = names(data)[row],
                                            i = row, j = col,
                                            x = xs2,  y = ys2,
                                            type = "robust"))
                         } else {
                           rbind(data.frame(xvar = names(data)[col], yvar = names(data)[row],
                                            i = row, j = col,
                                            x = as.numeric(NA), y = as.numeric(NA),
                                            type = "classic"),
                                 data.frame(xvar = names(data)[col], yvar = names(data)[row],
                                            i = row, j = col,
                                            x = as.numeric(NA), y = as.numeric(NA),
                                            type = "robust"))
                         }
                       })
        
        ellipse.pts <- do.call("rbind", ellipse.pts)
        
        ellipse.pts$x <- suppressWarnings(as.numeric(as.character(ellipse.pts$x)))
        ellipse.pts$y <- suppressWarnings(as.numeric(as.character(ellipse.pts$y)))
        
        ellipse.pts$xvar <- factor(ellipse.pts$xvar, levels = names(data))
        ellipse.pts$yvar <- factor(ellipse.pts$yvar, levels = names(data))
        
        plt <- plt + geom_polygon(mapping = aes_string(x = "x", y = "y", color = "type", linetype = "type"),
                                  data    = ellipse.pts, na.rm = T, fill = NA, show.legend = T)
        
        # Calculate each variable's kernel density 
        densities <- with(all, simplify = FALSE, 
                          tapply(x, yvar, stats::density, na.rm = TRUE))
        densities <- lapply(densities, function(x) data.frame(x = x$x, y = x$y))
        names(densities) <- levels(all$xvar)
        
        # Melt densities
        densities <- reshape2::melt(densities, id.vars = c("x", "y"))
        names(densities) <- c("x", "y", "xvar")
        densities$xvar <- factor(densities$xvar, levels = levels(all$xvar))
        densities$yvar <- factor(densities$xvar, levels = levels(all$xvar))
        
        # Scale each variable's density estimate to match range of the variable
        densities$y.scaled <- NA
        densities$x.scaled <- NA
        for(v in levels(all$xvar)) {
          var.vals <- densities$x[densities$xvar == v]
          var.dens <- densities$y[densities$xvar == v]
          
          scaled.vals <- scales::rescale(var.vals, to = range(-1:1, na.rm = TRUE))
          scaled.dens <- scales::rescale(var.dens, to = range(-1:1, na.rm = TRUE))
          
          densities$x.scaled[densities$xvar == v] <- scaled.vals
          densities$y.scaled[densities$xvar == v] <- scaled.dens
        }
        
        # Add density line to plot
        plt <- plt + geom_line(data = densities, aes(x = x.scaled, y = y.scaled), color = "red3")
        
        corr.text <- lapply(1:nrow(grid),
                            function(i) {
                              row <- grid[i, "x"]
                              col <- grid[i, "y"]
                              if (row < col) {
                                corr  <- X[row, col]
                                corr2 <- X2[row, col]
                               
                                rbind(data.frame(xvar = names(data)[col], yvar = names(data)[row],
                                                 lab = as.character(formatC(round( corr, 2 ), format='f', digits=2)),
                                                 i = row, j = col,
                                                 x = 0.0,  y = 0.25,
                                                 type = "classic"),
                                      data.frame(xvar = names(data)[col], yvar = names(data)[row],
                                                 lab = as.character(formatC(round( corr2, 2 ), format='f', digits=2)),
                                                 i = row, j = col,
                                                 x = 0.0,  y = -0.25,
                                                 type = "robust"))
                              } else {
                                rbind(data.frame(xvar = names(data)[col], yvar = names(data)[row],
                                                 lab = as.character(NA),
                                                 i = row, j = col,
                                                 x = as.numeric(NA), y = as.numeric(NA),
                                                 type = "classic"),
                                      data.frame(xvar = names(data)[col], yvar = names(data)[row],
                                                 lab = as.character(NA),
                                                 i = row, j = col,
                                                 x = as.numeric(NA), y = as.numeric(NA),
                                                 type = "robust"))
                              }
                            })
        
        corr.text <- do.call("rbind", corr.text)
        
        corr.text$x <- suppressWarnings(as.numeric(as.character(corr.text$x)))
        corr.text$y <- suppressWarnings(as.numeric(as.character(corr.text$y)))
        
        corr.text$xvar <- factor(corr.text$xvar, levels = names(data))
        corr.text$yvar <- factor(corr.text$yvar, levels = names(data))
        
        plt <- plt + geom_text(data = corr.text, aes_string(x = "x", y = "y", label = "lab", color = "type"),
                               size = text.size, na.rm = T, show.legend = F, check_overlap = TRUE) +
                 scale_color_manual(values = c("red", "black")) +
                 scale_linetype_manual(values = c("dashed", "solid"))
        
        plots[[i]] <- ggplotGrob(plt)
      }
      
      if (input$covariance.chi.qqplot == T) {
        i <- i + 1
        
        md <- if(corr) {
                lapply(fm,
                       function(x, mat) {
                         sqrt(mahalanobis(mat, x$center, x$cov * (x$sdev %o% x$sdev)))
                       },
                       data)
              } else {
                lapply(fm,
                       function(x, mat) {
                         sqrt(mahalanobis(mat, x$center, x$cov))
                       },
                       data)
              }
        
        chisq.points <- sqrt(qchisq(ppoints(nrow(data)), ncol(data)))
        
        dat <- data.frame(x = chisq.points, y = sort(md[[1]]))
        
        p1 <- ggplot(data = dat, aes(x = x, y = y)) +
                ggtitle("Classic") +
                xlab("Chi-Squared Quantiles") +
                ylab("Sorted Distances") +
                geom_point(color = "dodgerblue2", shape = 18, size = 2.5) +
                geom_abline(slope = 1, intercept = 0, linetype = 2)
        
        dat <- data.frame(x = chisq.points, y = sort(md[[2]]))
        
        p2 <- ggplot(data = dat, aes(x = x, y = y)) +
                ggtitle("Robust") +
                xlab("Chi-Squared Quantiles") +
                ylab("Sorted Distances") +
                geom_point(color = "dodgerblue2", shape = 18, size = 2.5) +
                geom_abline(slope = 1, intercept = 0, linetype = 2)
        
        y.min <- min(layer_scales(p1)$y$range$range, layer_scales(p2)$y$range$range)
        y.max <- max(layer_scales(p1)$y$range$range, layer_scales(p2)$y$range$range)
        
        p1 <- p1 + scale_y_continuous(limits = c(y.min, y.max),
                                      expand = expand_scale(mult = c(0.05, 0.05)))
        
        p2 <- p2 + scale_y_continuous(limits = c(y.min, y.max),
                                      expand = expand_scale(mult = c(0.05, 0.05)))
        
        plots[[i]] <- cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last")
      }
      
      if (input$covariance.dist.dist == T) {
        i <- i + 1
        
        md <- lapply(fm,
                     function(x, mat) {
                       sqrt(mahalanobis(mat, x$center, x$cov))
                     },
                     data)
        
        dat <- data.frame(x = md[[1]], y = md[[2]])
        
        N <- sapply(md, length)
        p <- sapply(fm, function(x) nrow(x$cov))
        
        thresh <- sqrt(qchisq(0.95, df = p))
        
        min.val <- min(dat)
        max.val <- max(dat)
        
        plt <- ggplot(data = dat, aes(x = x, y = y)) +
                 ggtitle("Distance-Distance Plot") +
                 xlab("Classic Distance") +
                 ylab("Robust Distance") +
                 xlim(c(min.val, max.val)) +
                 ylim(c(min.val, max.val)) +
                 theme(aspect.ratio = 1) +
                 geom_point(color = "dodgerblue2", shape = 18, size = 2.5) +
                 geom_hline(yintercept = thresh[1], linetype = 2) +
                 geom_vline(xintercept = thresh[1], linetype = 2) +
                 geom_abline(aes(slope = 1, intercept = 0), linetype = 2)
        
        plots[[i]] <- ggplotGrob(plt)
      }
    } else {
      fit <- values$covariance.fit
      
      # Scree plot for covariance
      if (input$covariance.eigen == T) {
        i <- i + 1
        
        eigen.vals <- if (fit$corr == TRUE) {
                        eigen(fit$cov * (fit$sdev %o% fit$sdev))$values
                      } else {
                        eigen(fit$cov)$values
                      }
        
        
        n <- length(eigen.vals)
        
        indx <- 1:n
        
        dat <- data.frame(indx, eigen.vals)
        
        plt <- ggplot(data = dat) +
                 ggtitle("Eigenvalues") +
                 ylab("Eigenvalue") +
                 theme(aspect.ratio = 3/4) +
                 geom_line(aes(x = indx, y = eigen.vals), color = "dodgerblue2", size = 1) +
                 geom_point(aes(x = indx, y = eigen.vals), color = "dodgerblue2", shape = 5, size = 2.5) +
                 scale_x_continuous(name = "Factor Number", breaks = min(dat$indx):max(dat$indx))
       
        plots[[i]] <- ggplotGrob(plt)
      }
      
      if (input$covariance.mahalanobis == T) {
        i <- i + 1
        
        md <- if (fit$corr == TRUE) {
                sqrt(mahalanobis(data, fit$center, fit$cov * (fit$sdev %o% fit$sdev)))
              } else {
                sqrt(mahalanobis(data, fit$center, fit$cov))
              }
        
        n <- length(md)
        p <- nrow(fit$cov)
        
        thresh <- sqrt(qchisq(0.95, df = p))
        
        indx <- 1:n
        
        if (class(fit) == "covClassic") {
          title.name = "Classic Distance"
        } else {
          title.name = "Robust Distance"
        }
        
        dat = data.frame(indx = indx, MD = md)
        
        plt <- ggplot(data = dat, aes(x = indx, y = MD)) +
                 ggtitle(title.name) +
                 theme(aspect.ratio = 3/4) +
                 xlab("Index") +
                 ylab("Mahalanobis Distance") +
                 geom_point(color = "dodgerblue2", shape = 18, size = 2.5) +
                 geom_hline(yintercept = thresh[1], linetype = 2)
        
        plots[[i]] <- ggplotGrob(plt)
      }
      
      if (input$covariance.ellipses.matrix == T) {
        i <- i + 1
        
        ellipse <- function(loc, A) {
          detA <- A[1, 1] * A[2, 2] - A[1, 2]^2
          dist <- sqrt(qchisq(0.95, 2))
          ylimit <- sqrt(A[2, 2]) * dist
          y <- seq(-ylimit, ylimit, 0.01 * ylimit)
          sqrt.discr <- detA/A[2,2]^2 * (A[2,2] * dist^2 - y^2)
          sqrt.discr[c(1, length(sqrt.discr))] <- 0.0
          sqrt.discr <- sqrt(sqrt.discr)
          b <- loc[1] + A[1, 2] / A[2, 2] * y
          x1 <- b - sqrt.discr
          x2 <- b + sqrt.discr
          y <- loc[2] + y
          rbind(cbind(x1, y), cbind(rev(x2), rev(y)))
        }
        
        grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
        
        text.size <- ifelse(ncol(data) > 2, 4, 10)
        
        # data.frame with xy coordinates
        all <- lapply(1:nrow(grid),
                      function(i) {
                        xcol <- grid[i, "x"]
                        ycol <- grid[i, "y"]
                        data.frame(xvar = names(data)[ycol], yvar = names(data)[xcol],
                                   i = grid[i, "x"], j = grid[i, "y"],
                                   x = data[, xcol], y = data[, ycol])
                      })
        
        all <- do.call("rbind", all)
        
        all.upper <- all
        
        all.upper[[6]][all.upper[[3]] <  all.upper[[4]]] <- "NA"
        all.upper[[5]][all.upper[[3]] <= all.upper[[4]]] <- "NA"
        
        all.upper$x <- suppressWarnings(as.numeric(as.character(all.upper$x)))
        all.upper$y <- suppressWarnings(as.numeric(as.character(all.upper$y)))
        
        all.upper$xvar <- factor(all.upper$xvar, levels = names(data))
        all.upper$yvar <- factor(all.upper$yvar, levels = names(data))
        
        xy.mapping <- aes_string(x = "x", y = "y")
        class(xy.mapping) <- "uneval"
        
        plt <- ggplot(all.upper, mapping = aes_string(x = "x", y = "y")) +
                 theme(axis.text.x  = element_blank(), axis.text.y  = element_blank(),
                       axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
                       axis.title.x = element_blank(), axis.title.y = element_blank(),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.spacing = unit(0, "lines")) +
                 facet_grid(xvar ~ yvar, scales = "free") +
                 theme(aspect.ratio = 1)
        
        if (fit$corr == F) {
          s <- sqrt(diag(fit$cov))
          X <- fit$cov / (s %o% s) 
        } else {
          X <- fit$cov
        }
        
        ellipse.pts <- lapply(1:nrow(grid),
                              function(i) {
                                row <- grid[i, "x"]
                                col <- grid[i, "y"]
                                if (row > col) {
                                  cov <- X[row, col]
                                  
                                  pts <- c(seq(0.0, 2*pi, length.out = 181), NA)
                                  xs <- cos(pts + acos(cov) / 2)
                                  ys <- cos(pts - acos(cov) / 2)
                                  
                                  data.frame(xvar = names(data)[col], yvar = names(data)[row],
                                             i = row, j = col,
                                             x = xs,  y = ys)
                                } else {
                                  data.frame(xvar = names(data)[col], yvar = names(data)[row],
                                             i = row, j = col,
                                             x = as.numeric(NA), y = as.numeric(NA))
                                }
                              })
        
        ellipse.pts <- do.call("rbind", ellipse.pts)
        
        ellipse.pts$x <- suppressWarnings(as.numeric(as.character(ellipse.pts$x)))
        ellipse.pts$y <- suppressWarnings(as.numeric(as.character(ellipse.pts$y)))
        
        ellipse.pts$xvar <- factor(ellipse.pts$xvar, levels = names(data))
        ellipse.pts$yvar <- factor(ellipse.pts$yvar, levels = names(data))
        
        plt <- plt + geom_polygon(mapping = aes_string(x = "x", y = "y"), color = "black",
                                  data    = ellipse.pts, na.rm = T, fill = NA, show.legend = T)
        
        # Calculate each variable's kernel density 
        densities <- with(all, simplify = FALSE, 
                          tapply(x, yvar, stats::density, na.rm = TRUE))
        densities <- lapply(densities, function(x) data.frame(x = x$x, y = x$y))
        names(densities) <- levels(all$xvar)
        
        # Melt densities
        densities <- reshape2::melt(densities, id.vars = c("x", "y"))
        names(densities) <- c("x", "y", "xvar")
        densities$xvar <- factor(densities$xvar, levels = levels(all$xvar))
        densities$yvar <- factor(densities$xvar, levels = levels(all$xvar))
        
        # Scale each variable's density estimate to match range of the variable
        densities$y.scaled <- NA
        densities$x.scaled <- NA
        for(v in levels(all$xvar)) {
          var.vals <- densities$x[densities$xvar == v]
          var.dens <- densities$y[densities$xvar == v]
          
          scaled.vals <- scales::rescale(var.vals, to = range(-1:1, na.rm = TRUE))
          scaled.dens <- scales::rescale(var.dens, to = range(-1:1, na.rm = TRUE))
          
          densities$x.scaled[densities$xvar == v] <- scaled.vals
          densities$y.scaled[densities$xvar == v] <- scaled.dens
        }
        
        # Add density line to plot
        plt <- plt + geom_line(data = densities, aes(x = x.scaled, y = y.scaled), color = "red3")
        
        corr.text <- lapply(1:nrow(grid),
                            function(i) {
                              row <- grid[i, "x"]
                              col <- grid[i, "y"]
                              if (row < col) {
                                corr  <- X[row, col]
                               
                                data.frame(xvar = names(data)[col], yvar = names(data)[row],
                                           lab = as.character(formatC(round( corr, 2 ), format='f', digits=2)),
                                           i = row, j = col,
                                           x = 0.0,  y = -0.25)
                              } else {
                                data.frame(xvar = names(data)[col], yvar = names(data)[row],
                                           lab = as.character(NA),
                                           i = row, j = col,
                                           x = as.numeric(NA), y = as.numeric(NA))
                              }
                            })
        
        corr.text <- do.call("rbind", corr.text)
        
        corr.text$x <- suppressWarnings(as.numeric(as.character(corr.text$x)))
        corr.text$y <- suppressWarnings(as.numeric(as.character(corr.text$y)))
        
        corr.text$xvar <- factor(corr.text$xvar, levels = names(data))
        corr.text$yvar <- factor(corr.text$yvar, levels = names(data))
        
        plt <- plt + geom_text(data = corr.text, aes_string(x = "x", y = "y", label = "lab"),
                               size = text.size, na.rm = TRUE, check_overlap = TRUE)
        
        plots[[i]] <- ggplotGrob(plt)
      }
      
      if (input$covariance.chi.qqplot == T) {
        i <- i + 1
        
        md <- if (fit$corr == T) {
                sqrt(mahalanobis(data, fit$center, fit$cov * (fit$sdev %o% fit$sdev)))
              } else {
                sqrt(mahalanobis(data, fit$center, fit$cov))
              }
        
        chisq.points <- sqrt(qchisq(ppoints(nrow(data)), ncol(data)))
        
        dat <- data.frame(x = chisq.points, y = sort(md))
        
        plt <- ggplot(data = dat, aes(x = x, y = y)) +
                 ggtitle("Classic") +
                 xlab("Chi-Squared Quantiles") +
                 theme(aspect.ratio = 3/4) +
                 ylab("Sorted Distances") +
                 geom_point(color = "dodgerblue2", shape = 18, size = 2.5) +
                 geom_abline(slope = 1, intercept = 0, linetype = 2)
        
        plots[[i]] <- ggplotGrob(plt)
      }
    }
    
    
    if (i == 0) {
      output$covariance.plot.ui <- renderUI({
        fluidPage(
          fluidRow(verbatimTextOutput("no.selection"))
        )
      })
    } else {
      values$covariance.plots        <- plots
      values$covariance.num.plots    <- i
      values$covariance.active.index <- 1
      values$covariance.active.plot  <- plots[[1]]
      
      if (i > 1) {
        output$covariance.plot.ui <- renderUI({
          fluidPage(
            wellPanel(
              plotOutput("covariance.plot.output")
            ),
            
            fluidRow(
              column(1,
                     offset = 10,
                     actionButton("covariance.next.plot",
                                  "",
                                  icon = icon("angle-right", "fa-2x"))
              )
            )
          )
        })
      } else {
        output$covariance.plot.ui <- renderUI({
          fluidPage(
            wellPanel(
              plotOutput("covariance.plot.output")
            )
          )
        })
      }
      
      output$covariance.plot.output <- renderPlot({
        grid.draw(values$covariance.active.plot)
      })
    }
  })
  
  observeEvent(input$covariance.next.plot, {
    if (values$covariance.num.plots > 0) {
      values$covariance.active.index <- values$covariance.active.index %% values$covariance.num.plots + 1
      values$covariance.active.plot  <- values$covariance.plots[[values$covariance.active.index]]
      
      output$covariance.plot.ui <- renderUI({
        if (values$covariance.active.index == values$covariance.num.plots) {
          fluidPage(
            wellPanel(
              plotOutput("covariance.plot.output")
            ),
            
            fluidRow(
              column(1,
                     offset = 1,
                     actionButton("covariance.prev.plot",
                                  "",
                                  icon = icon("angle-left", "fa-2x"))
              )
            )
          )
        } else {
          fluidPage(
            wellPanel(
              plotOutput("covariance.plot.output")
            ),
            
            fluidRow(
              column(1,
                     offset = 1,
                     actionButton("covariance.prev.plot",
                                  "",
                                  icon = icon("angle-left", "fa-2x"))
              ),
              
              column(1,
                     offset = 8,
                     actionButton("covariance.next.plot",
                                  "",
                                  icon = icon("angle-right", "fa-2x"))
              )
            )
          )
        }
      })
      
      output$covariance.plot.output <- renderPlot({
        grid.draw(values$covariance.active.plot)
      })
    }
  })
  
  # On button press, move to previous plot(s)
  observeEvent(input$covariance.prev.plot, {
    if (values$covariance.num.plots > 0) {
      values$covariance.active.index <- values$covariance.num.plots + (values$covariance.active.index - 1) %% (-values$covariance.num.plots)
      values$covariance.active.plot  <- values$covariance.plots[[values$covariance.active.index]]
      
      output$covariance.plot.ui <- renderUI({
        if (values$covariance.active.index == 1) {
          fluidPage(
            wellPanel(
              plotOutput("covariance.plot.output")
            ),
            
            fluidRow(
              column(1,
                     offset = 10,
                     actionButton("covariance.next.plot",
                                  "",
                                  icon = icon("angle-right", "fa-2x"))
              )
            )
          )
        } else {
          fluidPage(
            wellPanel(
              plotOutput("covariance.plot.output")
            ),
            
            fluidRow(
              column(1,
                     offset = 1,
                     actionButton("covariance.prev.plot",
                                  "",
                                  icon = icon("angle-left", "fa-2x"))
              ),
              
              column(1,
                     offset = 8,
                     actionButton("covariance.next.plot",
                                  "",
                                  icon = icon("angle-right", "fa-2x"))
              )
            )
          )
        }
      })
      
      output$covariance.plot.output <- renderPlot({
        grid.draw(values$covariance.active.plot)
      })
    }
  })
  
#########################
## PCA ##
#########################
  
  output$pca.select.variables <- renderUI({
    if (is.null(dim(values$dat))) {
      return("")
    }
    
    selectInput("pca.variables", "Variables",
                choices  = values$dat.numeric.variables,
                selected = values$dat.numeric.variables,
                multiple = TRUE)
  })
  
  observeEvent(input$pca.display, {
    if (is.null(values$dat)) {
      output$pca.results <- renderText({
        return(paste0("<font color=\"#FF0000\"><b>", "ERROR: No data loaded!", "</b><font>"))
      })
    } else if (is.null(input$pca.variables)) {
      output$pca.results <- renderText({
        return(paste0("<font color=\"#FF0000\"><b>", "ERROR: Missing variables! Please input at least 2 variables.", "</b><font>"))
      })
    } else if (length(input$pca.variables) <= 2) {
      output$pca.results <- renderText({
        return(paste0("<font color=\"#FF0000\"><b>", "ERROR: Insufficient number of variables! Please add more.", "</b><font>"))
      })
    } else {
      
      values$pca.active <- TRUE
      data <- values$dat.numeric[, input$pca.variables]
      
      if (input$pca.standardize) {
        data <- scale(data, apply(data, 2, median), apply(data, 2, mad))
      }
      
      if (input$pca.method == "classic") {
        values$pca.num <- 1
        
        values$pca.fit <- prcomp(data)
      } else if (input$pca.method == "rob") {
        values$pca.num <- 1
        
        values$pca.fit <- prcompRob(data)
      } else {
        values$pca.num <- 2
        
        values$pca.fit <- fit.models(Classic = prcomp(data),
                                     Robust  = prcompRob(data,
                                     est = input$pca.estimator))
      }
      
      # Print summary method for pca results
      output$pca.results <- renderText({
        
        return(htmlText(summary(values$pca.fit)))
      })
    }
  })
  
  ## PCA - Plotting
  # observeEvent(input$pca.display.plots, {
  #   values$pca.plots.active <- TRUE
  #   
  #   i <- 0
  #   
  #   n <- min(5, ncol(data))
  #   
  #   plots <- vector(mode = "list")
  #   
  #   if (values$pca.num == 2) {
  #     fm <- values$pca.fit
  #     if (input$pca.scree == T) {
  #        i <- i + 1
  #       
  #       vars <- lapply(fm, function(x) { x$sdev^2 })
  #       vars <- sapply(vars, function(v) { v / sum(v) })
  #       
  #       indx <- sapply(1:n, function(i) { paste("Comp.", i) })
  #       
  #       indx <- rep(indx, 2)
  #       
  #       Method <- c(rep("Classic", n), rep("Robust", n))
  #       
  #       dat <- data.frame(indx, eigen.vals, Method)
  #       
  #       plt <- ggplot(data = dat) +
  #                ggtitle("Scree Plot") +
  #                ylab("Variance") +
  #                geom_line(aes(x = indx, y = eigen.vals, color = Method, linetype = Method), size = 1) +
  #                geom_point(aes(x = indx, y = eigen.vals, color = Method), shape = 16, size = 2.5) +
  #                scale_color_manual(values = c("red", "black")) +
  #                scale_linetype_manual(values = c("dashed", "solid"))
  #      
  #       plots[[i]] <- ggplotGrob(plt)
  #     }
  #   } else {
  #     fit <- values$pca.fit
  #     if (input$pca.scatter == T) {
  #       
  #     }
  #     
  #     if (input$pca.loadings == T) {
  #       dat <- fit$x
  #       
  #       grid <- expand.grid(x = 1:n, y = 1:n)
  #       
  #       text.size <- ifelse(ncol(data) > 2, 4, 10)
  #       
  #       # data.frame with xy coordinates
  #       all <- lapply(1:n,
  #                     function(i) {
  #                       xcol <- grid[i, "x"]
  #                       ycol <- grid[i, "y"]
  #                       data.frame(xvar = names(data)[ycol], yvar = names(data)[xcol],
  #                                  i = grid[i, "x"], j = grid[i, "y"],
  #                                  x = data[, xcol], y = data[, ycol])
  #                     })
  #       
  #       all <- do.call("rbind", all)
  #       
  #       all.upper <- all
  #       
  #       all.upper[[6]][all.upper[[3]] <  all.upper[[4]]] <- "NA"
  #       all.upper[[5]][all.upper[[3]] <= all.upper[[4]]] <- "NA"
  #       
  #       all.upper$x <- suppressWarnings(as.numeric(as.character(all.upper$x)))
  #       all.upper$y <- suppressWarnings(as.numeric(as.character(all.upper$y)))
  #       
  #       all.upper$xvar <- factor(all.upper$xvar, levels = names(data))
  #       all.upper$yvar <- factor(all.upper$yvar, levels = names(data))
  #     }
  #     
  #     if (input$pca.scree == T) {
  #       i <- i + 1
  #       
  #       vars <- fit$sdev^2
  #       vars <- vars / sum(vars)
  #       
  #       indx <- sapply(1:n, function(i) { paste("Comp.", i) })
  #       
  #       dat <- data.frame(indx, vars)
  #       
  #       plt <- ggplot(data = dat) +
  #                ggtitle("Scree Plot") +
  #                ylab("Variance") +
  #                geom_line(aes(x = indx, y = vars), color = "dodgerblue2", size = 1) +
  #                geom_point(aes(x = indx, y = vars), color = "dodgerblue2", shape = 5, size = 2.5)
  #               # scale_x_continuous(name = "Factor Number", breaks = min(dat$indx):max(dat$indx))
  #      
  #       plots[[i]] <- ggplotGrob(plt)
  #     }
  #   }
  # })
  
  # Reset all windows once new data set is loaded
  observeEvent(input$display.table, {
    if (values$linRegress.active) {
      output$linRegress.results.ui <- renderUI({ invisible() })
      
      values$linRegress.active <- FALSE
      
      if (values$linRegress.plots.active) {
        output$linRegress.plot.ui <- renderUI({ invisible() })
        
        values$linRegress.plots.active <- FALSE
      }
    }
    
    if (values$covariance.active) {
      output$covariance.results.ui <- renderUI({ invisible() })
      
      values$covariance.active <- FALSE
      
      if (values$covariance.plots.active) {
        output$covariance.plot.ui <- renderUI({ invisible() })
        
        values$covariance.plots.active <- FALSE
      }
    }
    
    if (values$pca.active) {
      output$pca.results.ui <- renderUI({ invisible() })
      
      values$pca.active <- FALSE
      
      if (values$pca.plots.active) {
        output$pca.plot.ui <- renderUI({ invisible() })
        
        values$pca.plots.active <- FALSE
      }
    }
    
    updateTabsetPanel(session  = getDefaultReactiveDomain(), "linear.tabs",
                      selected = "linear.model")
    
    updateTabsetPanel(session  = getDefaultReactiveDomain(), "covariance.tabs",
                      selected = "covariance.est")
    
    updateTabsetPanel(session  = getDefaultReactiveDomain(), "pca.tabs",
                      selected = "pca.est")
  })
  
  # Reset plotting window for linear regression
  observeEvent(input$linRegress.display, {
    if (values$linRegress.plots.active) {
      output$linRegress.plot.ui <- renderUI({ invisible() })
      
      values$linRegress.plots.active <- FALSE
    }
  })
  
  observeEvent(input$covariance.display, {
    if (values$covariance.plots.active) {
      output$covariance.plot.ui <- renderUI({ invisible() })
      
      values$covariance.plots.active <- FALSE
    }
  })
  
  observeEvent(input$pca.display, {
    if (values$pca.plots.active) {
      output$pca.plot.ui <- renderUI({ invisible() })
      
      values$pca.plots.active <- FALSE
    }
  })
  
  output$about.text <- renderText({"
    <div>
      <font color=\"#000000\"> Text for this section to be added on a later date. 
        Please see the <a href=\"https://github.com/GregoryBrownson/RobStatTM-GUI/blob/master/vignette/RobStatTM%20GUI_vignette.pdf\">vignette</a>
        for more information.</font>
    </div>
  "})
  
  output$help.text <- renderText({"
    <div>
      <font color=\"#000000\"> Please contact <a href=\"mailto:gregory.brownson@gmail.com\">Greg Brownson</a>
        for any questions about the application. Any comments or suggestions to improve on the
        <strong>RobStatTM</strong> shiny interface are encouraged as well!</font>
    </div>
  "})
})
