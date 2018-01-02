  TR2=function(R2, family, cc)
# given de clasic R^2 (R2) compute the robust (RR2)
#cc is the tuning constant
{ a=Erho (family,cc,1)
 b=Erho( family,cc,sqrt(1-R2))				
 RR2=(b-a)/(b*(1-a))
 RR2}

Erho =function(family,cc,z)	
# compute E(rho(u/z,family,cc))
{
if(family=="bisquare")
{
dd=cc*z
a0=2*pnorm(dd)-1
a2=-2*dd*dnorm(dd)+a0
a4=-2*(dd^3)*dnorm(dd)+3*a2
a6=-2*(dd^5)*dnorm(dd)+5*a4
ee=(a6/dd^6)+(3*a2/dd^2)-(3*a4/dd^4)+1-a0
} else 
 {ee=2*(integrate(hh, 0, cc[3]*z,family=family,cc=cc,z=z)$value+1-pnorm(cc[3]*z))
}
ee
}
  	 	
 





ff=function(x,y,family,cc)
{uu=TR2(x,family,cc)-y
uu}

hh=function(v,family, cc,z)
{uu=rho(v/z, family=family,cc,z)*dnorm(v)
uu}


INVTR2=function(RR2,family,cc)
# given the robust R^2 (RR2) compute the classic (R2)
#as the inverse of TR2 
#cc is the tuning constant
{aa=TR2(.99999,family,cc)
bb=TR2(.00001,family,cc)
 
if(RR2>.99)
{R2=1} 
if (RR2<bb) 
{R2=0} 
if((RR2<=.99)&(RR2>=bb))
{R2=uniroot(ff,c( bb/2, aa+((1-aa)/2)),y=RR2,family=family,cc=cc)$root}
R2=as.numeric(R2)
R2}

 
	
 


		




