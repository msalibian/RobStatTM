
# robust ratio likelihood test for linear assumptions
# input
#object1 lmRob object corresponding to  the complete model
#object2 lmRob object  corresponding to the the null hypothesis model 
#output
# $.test value of the F-statistics
# $.Fpvalue pvalue using the F distribution
# $chispvalue pvalue using the chisquare distribution



rob.linear.test=function(object1, object2)
{
	p = length(object1$coeff)
	q = length(object1$coeff) - length(object2$coeff)
	n = length(object1$resid)
	 con=object1$control
                   cc=con$tuning.psi
                  s=object1$scale
	a = sum(Mchi(object2$resid/s, cc,  psi="bisquare" ))
	b = sum(Mchi(object1$resid/s, cc, psi="bisquare"))
	c = sum(Mchi(object1$resid/s, cc, psi="bisquare",2))
	d = sum(Mchi(object1$resid/s, cc, psi="bisquare",1)^2)
	test = (2 * (a - b) * c)/d
	chisq.pvalue =1- pchisq(test, q)
	F.pvalue = 1-pf(test/q, q, n - p)
	df=c(q,n-p)
	output = list(test, chisq.pvalue, F.pvalue,df)
	names(output)=c("test","chisq.pvalue","f.pvalue","df")
	output
}
 