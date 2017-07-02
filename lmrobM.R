lmrobM=function (formula, tuning.psi=3.44)
{
 outL=rq(formula, tau=.5)
resL=sort(abs(outL$resid))
p=length(outL$coeff)
n=length(outL$resid)
s=resL[(n+p)/2]/.6745
initial=list(coefficients=outL$coeff,scale=s)
outM=lmrob (formula, tuning.psi=tuning.psi ,init=initial)
outM
}


