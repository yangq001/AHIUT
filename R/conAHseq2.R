#faster seq

conAHseq2=function(Z,ld,n,k=6,alpha=0.05){
  pp=JointSum(Z,rep(1,length(Z)),N=n,XX=ld,YY0=diag(1),adj_Y=0)
  bbb=pp$beta
  ccc=pp$cov
  q=length(bbb)
  t0=t(bbb)%*%ginv(ccc)%*%bbb
  ppp=c(pchisq(t0[1,1],df=q,lower.tail=FALSE))
  ok=c()
  
  if(ppp>alpha){
    return(list(snp=c(),pvalue=ppp,k=0))
  }
  
  for(s in 1:min(k-1,q-1)){
    if(length(ok)>0){
      J=(1:q)[-ok]
    }else{
      J=1:q
    }
    ppdie=ppsig=c()
    for(jj in J){
      j=c(ok,jj)
      bbb2=bbb[-j]
      ccc2=ccc[-j,]
      if(is.matrix(ccc2)==TRUE){
        ccc2=ccc2[,-j]
      }else{
        ccc2=ccc2[-j]
      }
      t1=t(bbb2)%*%ginv(ccc2)%*%bbb2
      
      bbb2=bbb[j]
      ccc2=ccc[j,]
      if(is.matrix(ccc2)==TRUE){
        ccc2=ccc2[,j]
      }else{
        ccc2=ccc2[j]
      }
      t2=t(bbb2)%*%ginv(ccc2)%*%bbb2
      
      ppdie=c(ppdie,pchisq(t1[1,1],df=q-s,lower.tail=FALSE))
      ppsig=c(ppsig,pchisq(t2[1,1],df=s,lower.tail=FALSE))
    }
    wow=which(ppsig==min(ppsig))[1]
    ok=c(ok,J[wow])
    ppp=c(ppp,max(ppdie))
    
    if(max(ppdie)>alpha){
      #J=(1:q)[-ok]
      #ppsig=c()
      #for(jj in J){
      #  bbb2=bbb[j]
      #  ccc2=ccc[j,]
      #  if(is.matrix(ccc2)==TRUE){
      #    ccc2=ccc2[,j]
      #  }else{
      #    ccc2=ccc2[j]
      #  }
      #  t2=t(bbb2)%*%ginv(ccc2)%*%bbb2
      #  ppsig=c(ppsig,pchisq(t2[1,1],df=s,lower.tail=FALSE))
      #}
      #wow=which(ppsig==min(ppsig))[1]
      #ok=c(ok,J[wow])
      return(list(snp=ok,pvalue=ppp,k=s))
    }
  }
  
  J=(1:q)[-ok]
  ppsig=c()
  for(jj in J){
    bbb2=bbb[j]
    ccc2=ccc[j,]
    if(is.matrix(ccc2)==TRUE){
      ccc2=ccc2[,j]
    }else{
      ccc2=ccc2[j]
    }
    t2=t(bbb2)%*%ginv(ccc2)%*%bbb2
    ppsig=c(ppsig,pchisq(t2[1,1],df=s,lower.tail=FALSE))
  }
  wow=which(ppsig==min(ppsig))[1]
  ok=c(ok,J[wow])
  return(list(snp=ok,pvalue=ppp,k=s+1))
}

