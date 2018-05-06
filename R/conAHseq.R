#seq

conAHseq=function(Z,ld,n,k=6,alpha=0.05){
  pp=JointSum(Z,rep(1,length(Z)),N=n,XX=ld,YY0=diag(1),adj_Y=0)
  bbb=pp$beta
  ccc=pp$cov
  q=length(bbb)
  t0=t(bbb)%*%ginv(ccc)%*%bbb
  ppp=c(pchisq(t0[1,1],df=q,lower.tail=FALSE))
  
  if(ppp>alpha){
    return(list(pvalue=ppp,k=0))
  }
  
  for(s in 1:min(k-1,q-1)){
    J=combn(1:q,s)
    ppdie=c()
    for(jj in 1:ncol(J)){
      j=J[,jj]
      bbb2=bbb[-j]
      ccc2=ccc[-j,]
      if(is.matrix(ccc2)==TRUE){
        ccc2=ccc2[,-j]
      }else{
        ccc2=ccc2[-j]
      }
      t1=t(bbb2)%*%ginv(ccc2)%*%bbb2
      #ppdie=c(ppdie,pchisq(t1[1,1],df=s,lower.tail=FALSE))
      ppdie=c(ppdie,pchisq(t1[1,1],df=q-s,lower.tail=FALSE))
    }
    
    ppp=c(ppp,max(ppdie))
    if(max(ppdie)>alpha){
      return(list(pvalue=ppp,k=s))
    }
  }
  return(list(pvalue=ppp,k=s+1))
}
