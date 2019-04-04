# oats.R
# EXAMPLE 4.2
# Figures 4.2, 4.4

library(RobStatTM)
data(oats)

## LS regression
oats1LS <- lm(response1 ~ variety+block, data=oats)
oats1LS_var <- lm(response1 ~ block, data=oats)
oats1LS_block <- lm(response1 ~ variety, data=oats)
oats2LS <- lm(response2 ~ variety+block, data=oats)
oats2LS_var <- lm(response2 ~ block, data=oats)
oats2LS_block <- lm(response2 ~ variety, data=oats)
sigma1LS <- summary(oats1LS)$sigma
sigma2LS <- summary(oats2LS)$sigma

## Classical ANOVA  tests
anov1_var  <- anova(oats1LS, oats1LS_var)
anov1_block <- anova(oats1LS, oats1LS_block)
anov2_var <- anova(oats2LS, oats2LS_var)
anov2_block <- anova(oats2LS, oats2LS_block)

## M regressions
cont <- lmrobdet.control(bb = 0.5, efficiency = 0.85, family = "bisquare")
oats1M <- lmrobM(response1 ~ variety+block, control=cont, data=oats)
oats1M_var <- lmrobM(response1 ~ block, control=cont, data=oats)
oats1M_block <- lmrobM(response1 ~ variety, control=cont, data=oats)
oats2M <- lmrobM(response2 ~ variety+block, control=cont, data=oats)
oats2M_var <- lmrobM(response2 ~ block, control=cont, data=oats)
oats2M_block <- lmrobM(response2 ~ variety, control=cont, data=oats)
sM2 <- oats2M$scale

## Robust ANOVA  tests
anov1M_var <- rob.linear.test(oats1M, oats1M_var)
anov1M_block <- rob.linear.test(oats1M, oats1M_block)
anov2M_var <- rob.linear.test(oats2M, oats2M_var)
anov2M_block <- rob.linear.test(oats2M, oats2M_block)

plot(oats2LS, which=2)
abline(h=0, lty=2)


tmp <- qqnorm(resid(oats2M)/sM2, ylab="Standardized residuals", pch=19, col='gray30')
qqline(resid(oats2M)/sM2)
abline(h=c(-2.5, 0, 2.5), lty=2)
w <- c(24,36,1,35,20)
text(tmp$x[w] + .1, tmp$y[w] + .1, w)

 
A=matrix(0,2,4); A[1,]= c(anov1_var[2,6], anov1M_var$F.pvalue, anov1_block[2,6], anov1M_block$F.pvalue)
 
A[2,]= c(anov2_var[2,6], anov2M_var$F.pvalue, anov2_block[2,6], anov2M_block$F.pvalue)
rownames(A)=c("Original", "Altered")
colnames(A)=c("F(rows)", "Robust(rows)", "F(cols)", "Robust(cols)")

"classical and robust p-values of ANOVA tests"
A
 
# Comment:  Due to changes in the codes that were made after the book's printing,
# not all p-values coincide with those in the example in the book. 






