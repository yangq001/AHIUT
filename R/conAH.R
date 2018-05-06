#Testing AH with IUT

conAH=function(Z,ld,n){
  pp=JointSum(Z,rep(1,length(Z)),N=n,XX=ld,YY0=diag(1),adj_Y=0)
  bbb=pp$beta
  ccc=pp$cov
  q=length(bbb)
  t0=t(bbb)%*%ginv(ccc)%*%bbb
  ppp=c(pchisq(t0[1,1],df=q,lower.tail=FALSE))
  for(j in 1:q){
    bbb2=bbb[-j]
    ccc2=ccc[-j,]
    if(is.matrix(ccc2)==TRUE){
      ccc2=ccc2[,-j]
    }else{
      ccc2=ccc2[-j]
    }
    t1=t(bbb2)%*%ginv(ccc2)%*%bbb2
    ppp=c(ppp,pchisq(t1[1,1],df=q-1,lower.tail=FALSE))
  }
  return(max(ppp))
}
