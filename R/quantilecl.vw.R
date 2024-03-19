quantilecl.vw<-function(train,test,cl,theta=NULL,cl.test=NULL){
  
  eps=1e-3
  it.max=50
  lambda=1
  qq.ratio=5
  k=length(unique(cl))
  if(is.character(cl)) cl<-as.numeric(as.factor(cl))
	if(is.character(cl.test)) cl.test<-as.numeric(as.factor(cl.test))
	p<-ncol(train)
	if(length(lambda)==1){lambda<-rep(lambda,p)}

	
	### init partition e quantiles
	
  	qq.best=qq=matrix(0,k,p)
  	VV=VV.test=Inf
	  VV.temp=NULL
	  n.train<-dim(train)[1]
	  n.test<-dim(test)[1]
  	QQ=QQ0=array(0,c(n.train,k,p))
	  QQ.test=QQ0.test=array(0,c(n.test,k,p))
	  me.train_seq=me.test_seq=me.test=NULL
	
	#######################################################################
	## initialization
	#######################################################################
	
  	## 1) equispaced interval
  	
	  if (!is.null(theta)) theta=c(matrix(theta,p,1))
	  
  	if (is.null(theta)) {
  	          theta=init.theta(train=train,cl=cl,k=k,p=p)
	            } 
	
  	### clustering step
  	ratio=5
  	h=0
  	nk=table(cl)
  	
  	      
    ### starting estimation 
  	
	 while ((ratio>eps) & (h<it.max)) {
	   
	   h=h+1
	   
	    ### b) compute the new barycenters
	    for (j in 1:p){
 	    for (i in 1:k) {
	    if (nk[i]>0) qq[i,j]=stats::quantile(train[cl==i,j],theta[j])
	    QQ0[,i,j]=(theta[j]+(1-2*theta[j])*(train[,j] < matrix(qq[i,j],n.train)))*abs(train[,j]-matrix(qq[i,j],n.train))
	    QQ0.test[,i,j]=(theta[j]+(1-2*theta[j])*(test[,j] < matrix(qq[i,j],n.test)))*abs(test[,j]-matrix(qq[i,j],n.test))
	    }
	    }
	    
	    ### d) estimate lambda
	   select<-function(QQ0,cl) return(QQ0[cbind(seq_along(cl),cl)])
	   den=colMeans(apply(QQ0,3,select,cl))
	   den=ifelse(den==0,eps,den)
	   lambda=1/den

	for (j in 1:p){
 	    for (i in 1:k) {
	    QQ[,i,j]=lambda[j]*QQ0[,i,j]-log(lambda[j]*theta[j]*(1-theta[j]))
	    QQ.test[,i,j]=lambda[j]*QQ0.test[,i,j]-log(lambda[j]*theta[j]*(1-theta[j]))
		}
	}
	   clpost.train=apply(apply(QQ,c(1,2),sum),1,which.min)
	   clpost.test=apply(apply(QQ.test,c(1,2),sum),1,which.min)

	   ### c) compute theta
	   for (j in 1:p) theta[j]=stats::optim(theta[j],fn_vw,method="L-BFGS-B",lower=0.0001,upper=0.9999,data=train[,j,drop=FALSE],k=k,cl=cl,qq=qq[,j,drop=FALSE],lambda=lambda[j])$par
	   
	   ### e) compute z
	   conta=conta.test=0
	   for (j in 1:p){
	     conta=conta+sum(QQ[cbind(seq_along(cl),cl,j)])
	     conta.test=conta.test+sum(QQ.test[cbind(seq_along(clpost.test),clpost.test,j)])
	   }
	   me.train<-misc(cl,clpost.train)
	   if (!is.null(cl.test)) me.test<-misc(cl.test,clpost.test)
  	   me.train_seq=c(me.train_seq,me.train)
  	   me.test_seq=c(me.test_seq,me.test)
  	 	
  	   VV[h]=conta
		  VV.test[h]=conta.test
  		ratio=abs(VV[h-1]-VV[h])/VV[h-1]
		if (h<qq.ratio) ratio=2*eps
	 }
  	
	names(theta)<-names(lambda)<-colnames(qq)<-colnames(train)
	out<-list(Vseq=VV, theta.choice=theta,me.train=me.train,me.test=me.test, 
            train=train,test=test,cl.train=clpost.train,cl.test=clpost.test,lambda=lambda)
	
	class(out)<-"quantileDA"
	return(out)
}





