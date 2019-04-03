# autism
# Example 6.7
# Tables 6.8,6.9

library(RobStatTM)

# Must install robustvarComp
# Must install nlme
# Must install WWWGbook

data(autism, package='WWGbook')
autism <- autism[complete.cases(autism),]
completi <- table(autism$childid)==5
completi <- names(completi[completi])
indici <- as.vector(unlist(sapply(completi, function(x) which(autism$childid==x))))
ind <- rep(FALSE, nrow(autism))
ind[indici] <- TRUE
autism <- subset(autism, subset=ind) ## complete cases 41
attach(autism)
sicdegp.f <- factor(sicdegp)
age.f <- factor(age)
age.2 <- age - 2
sicdegp2 <- sicdegp
sicdegp2[sicdegp == 3] <- 0
sicdegp2[sicdegp == 2] <- 2
sicdegp2[sicdegp == 1] <- 1
sicdegp2.f <- factor(sicdegp2)
autism.updated <- subset(data.frame(autism, sicdegp2.f, age.2), !is.na(vsae))
autism.grouped <- nlme::groupedData(vsae ~ age.2 | childid, data=autism.updated, order.groups = FALSE)
p <- 5
n <- 41
z1 <- rep(1, p)
z2 <- c(0, 1, 3, 7, 11)
z3 <- z2^2
K <- list()
K[[1]] <- tcrossprod(z1,z1)
K[[2]] <- tcrossprod(z2,z2)
K[[3]] <- tcrossprod(z3,z3)
K[[4]] <- tcrossprod(z1,z2) + tcrossprod(z2,z1)
K[[5]] <- tcrossprod(z1,z3) + tcrossprod(z3,z1)
K[[6]] <- tcrossprod(z3,z2) + tcrossprod(z2,z3)
names(K) <- c("Int", "age", "age2", "Int:age", "Int:age2", "age:age2")

groups <- cbind(rep(1:p, each=n), rep((1:n), p))

## Composite Tau
tmp.ctrl <- robustvarComp::varComprob.control(lower=c(0.01,0.01,0.01,-Inf,-Inf,-Inf))
AutismCompositeTau <- robustvarComp::varComprob(vsae ~ age.2 + I(age.2^2)
                                 + sicdegp2.f + age.2:sicdegp2.f + I(age.2^2):sicdegp2.f,
                                 groups = groups, data = autism.grouped, varcov = K,
                                 control=tmp.ctrl)
# The Warning messages produced by the above function indicates a small change to the code
# that is recommended, but the results are not effected by this Warning.
# That code change will be made in a future release of RobStatTM
summary(AutismCompositeTau)

## Classic S
tmp2.ctrl <- robustvarComp::varComprob.control(method="S", psi="rocke", cov.init="covOGK", 
                                               lower=c(0.01,0.01,0.01,-Inf,-Inf,-Inf))
AutismS <- robustvarComp::varComprob(vsae ~ age.2 + I(age.2^2)
                      + sicdegp2.f + age.2:sicdegp2.f + I(age.2^2):sicdegp2.f,
                      groups = groups, data = autism.grouped, varcov = K,
                      control=tmp2.ctrl)
# The Warning message produced by the above function is a notification, and not a
# true "warning".  This will be fixed in a future release of RobStatTM
summary(AutismS)















