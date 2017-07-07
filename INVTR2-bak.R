  TR2=function(R2,cc)
# given de clasic R^2 (R2) compute the robust (RR2)
#cc is the tuning constant
{ a=Erhobic(cc)
 b=Erhobic( cc*sqrt(1-R2))
 RR2=(b-a)/b
 RR2}

Erhobic=function(cc)
# compute E(rho(u,cc)), rho is the bisquare function
{a0=2*pnorm(cc)-1
a2=-2*cc*dnorm(cc)+a0
a4=-2*(cc^3)*dnorm(cc)+3*a2
a6=-2*(cc^5)*dnorm(cc)+5*a4
dd=(a6/cc^6)+(3*a2/cc^2)-(3*a4/cc^4)+1-a0
dd}


ff=function(x,y,cc)
{uu=TR2(x,cc)-y
uu}




INVTR2=function(RR2,cc)
# given de robust R^2 (RR2) compute the classic (R2)
#as the inverse of TR2 
#cc is the tuning constant
{
aa=TR2(.99999999,cc)
bb=TR2(.00000001,cc)
 
if(RR2>aa)
{R2=1} 
if (RR2<bb) 
{R2=0} 
if((RR2<=aa)&(RR2>=bb))
{R2=uniroot(ff,c( .000000001,.999999999),y=RR2,cc=cc)$root}
R2}


INVTR2(0.7367034,3.44)

 







