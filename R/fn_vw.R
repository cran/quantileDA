fn_vw=function(theta,data,k,cl,qq,lambda)
{VV=0
p=ncol(data)
for (i in 1:k) if (sum(cl==i)>0) {
  nn=sum(cl==i)
  xx=data[cl==i,,drop=FALSE]
  a=lambda*rowSums((theta+((1-2*theta)*(xx<t(matrix(qq[i,],p,nn)))))*abs(xx-t(matrix(qq[i,],p,nn))))
  VV=VV+sum(a)
}
numobs=length(cl)
VV=VV-numobs*log(lambda*theta*(1-theta))
return(VV)

}
