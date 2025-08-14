reg.link = function(x, coef.fix, logit=F, effect.fix=NULL) 
{ # coef.fix=alpha(for p) or gamma(for mu); coef.random=beta(for p) or nu(for mu)
  # logit=F (for mu) or T(for p) 
  if(is.null(effect.fix)){effect.fix=(x%*%coef.fix)[,1]}
  effect.total = effect.fix
  if(logit==T){return(effect.total-log1pexp(effect.total))} else {
    return(effect.total)
  }
}

loglik = function(z, rd, tn, x, alpha, gamma, phi, effect.fix.z=NULL, effect.fix.rd=NULL)
{
  if(is.null(effect.fix.z)){effect.fix.z=(x%*%alpha)[,1]}
  if(is.null(effect.fix.rd)){effect.fix.rd=(x%*%gamma)[,1]}
  p.log.vec=reg.link(x, coef.fix=alpha, logit=T, effect.fix=effect.fix.z) 
  mu.log.vec=reg.link(x, coef.fix=gamma, logit=F, effect.fix=effect.fix.rd) 
  mu.vec=exp(mu.log.vec)
  
  ll.gamma=ll.gamma.tn.bar=numeric(length(z))#=ll.gamma.tn=numeric(length(z))
  ll.gamma[z==1]=dgamma(rd[z==1],shape=1/phi,scale=mu.vec[z==1]*phi,log=TRUE)
  ll.gamma.tn.bar[z==0]=pgamma(tn[z==0],shape=1/phi,scale=mu.vec[z==0]*phi,lower.tail=FALSE,log.p=TRUE)
  p.bar.log.vec=log1mexp(-p.log.vec)
  
  ll.obs=z*(p.log.vec+ll.gamma)+(1-z)*(p.bar.log.vec+log1pexp(p.log.vec-p.bar.log.vec+ll.gamma.tn.bar))

  ll=sum(ll.obs)
  list(ll.obs=ll.obs, ll=ll)
}

rd.posterior.sim = function(z, rd, tn, x, gamma, phi, seed.eval=1001)
{
  set.seed(seed.eval)
  z.length=length(z)
  mu.log.vec=reg.link(x, coef.fix=gamma, logit=F) 
  mu.vec=exp(mu.log.vec)
  ll.gamma.tn.bar=pgamma(tn,shape=1/phi,scale=mu.vec*phi,lower.tail=FALSE,log.p=TRUE)
  rd.sim=numeric(z.length)
  rd.sim[(z==1)]=rd[(z==1)]
  rd.quantile.log=ll.gamma.tn.bar+log(1-runif(z.length))
  rd.sim[(z==0)]=qgamma(rd.quantile.log[(z==0)],shape=1/phi,scale=mu.vec[(z==0)]*phi,
                        lower.tail=FALSE,log.p=TRUE)
  return(rd.sim)
}

recur.z = function(z, rd, tn, x, alpha, gamma, phi, effect.fix.z=NULL, effect.fix.rd=NULL)
{
  if(is.null(effect.fix.z)){effect.fix.z=(x%*%alpha)[,1]}
  if(is.null(effect.fix.rd)){effect.fix.rd=(x%*%gamma)[,1]}
  p.log.vec=reg.link(x, coef.fix=alpha, logit=T, effect.fix=effect.fix.z) 
  mu.log.vec=reg.link(x, coef.fix=gamma, logit=F, effect.fix=effect.fix.rd) 
  p.vec=exp(p.log.vec)
  mu.vec=exp(mu.log.vec)
  gamma.tn.bar=pgamma(tn,shape=1/phi,scale=mu.vec*phi,lower.tail=FALSE,log.p=FALSE)
  
  z.e0=p.vec*gamma.tn.bar/(p.vec*gamma.tn.bar+(1-p.vec))
  z.e=rep(1,length(z))
  z.e[z==0]=z.e0[z==0]
  
  return(z.e)
}

recur.logit = function(z, rd, tn, x, alpha, gamma, phi,
                       iter.eval=20, sig.pen=NULL) ##alpha and beta
{
  alpha.length=length(alpha)
  param=alpha
  param.new=param
  param.old=param-Inf
  x.eval=x
  iter=0
  while ((iter<=iter.eval)&(sum((param.old-param.new)^2)>10^(-8)))
  {
    param.old=param.new
    alpha.new=param.new
    p.vec=exp(reg.link(x, coef.fix=alpha.new, logit=T))
    dQ = apply(sweep(x.eval,1,z-p.vec,FUN="*"),2,sum) -if(!is.null(sig.pen)){c(0, param.new[-1]/sig.pen^2)} else{0}
    dQ2 = -crossprod(sweep(x.eval,1,p.vec*(1-p.vec),FUN="*"),x.eval) -if(!is.null(sig.pen)){diag(c(0,1/sig.pen^2))} else{0}
    param.new=param.new-crossprod(dQ,solve(dQ2))[1,]
    iter = iter+1
  }
  alpha.new=param.new
  list(alpha=alpha.new, iter=iter)
}

recur.gamma = function(z, rd, tn, x, alpha, gamma, phi,
                       iter.eval=20, sig.pen=NULL) ##gamma and nu
{
  gamma.length=length(gamma)
  param=gamma
  param.new=param
  param.old=param-Inf
  x.eval=x
  iter=0
  while ((iter<=iter.eval)&(sum((param.old-param.new)^2)>10^(-8)))
  {
    param.old=param.new
    gamma.new=param.new
    mu.log.vec=reg.link(x, coef.fix=gamma.new, logit=F) 
    mu.vec=exp(mu.log.vec)
    dQ = apply(sweep(x.eval,1,z*(-1+rd/mu.vec),FUN="*"),2,sum)-if(!is.null(sig.pen)){c(0, param.new[-1]/sig.pen^2)} else{0}
    dQ2 = -crossprod(sweep(x.eval,1,z*rd/mu.vec,FUN="*"),x.eval)-if(!is.null(sig.pen)){diag(c(0,1/sig.pen^2))} else{0}
    param.new=param.new-crossprod(dQ,solve(dQ2))[1,]
    iter = iter+1
  }
  gamma.new=param.new
  list(gamma=gamma.new, iter=iter)
}

Q.phi = function(phi, term.1, term.2)
{
  return(-(term.1*(-phi*lgamma(1/phi)-log(phi))+term.2)/phi)
}

recur.phi = function(z, rd, tn, x, alpha, gamma, phi)
{
  mu.log.vec=reg.link(x, coef.fix=gamma, logit=F) 
  mu.vec=exp(mu.log.vec)
  term.1=sum(z)
  term.2=sum(z*(log(rd)-log(mu.vec)-rd/mu.vec))
  phi.opt=optimize(f=Q.phi, interval=c(0.5*phi,2*phi), term.1=term.1, term.2=term.2)$minimum
  return(phi.opt)
}


ecm = function(z, rd, tn, x, alpha, gamma, phi,
               iter.max=100,
               sig.pen=NULL, seed.num=999, print=TRUE)
{
  n=length(z)
  iter=0
  param.list=list()
  while(iter<iter.max)
  {
    #loglik.em.old = loglik.em
    alpha.old=alpha
    gamma.old=gamma
    phi.old=phi
    ## E-step
    rd.sim=rd.posterior.sim(z, rd, tn, x, gamma, phi, seed.eval=(seed.num*(iter+1)))

    
    z.e= recur.z(z, rd, tn, x, alpha, gamma, phi)
    
    ## M-step
    update.logit=recur.logit(z.e, rd.sim, tn, x, 
                             alpha, gamma, phi, iter.eval=10, sig.pen)
    alpha=update.logit$alpha
    update.gamma=recur.gamma(z.e, rd.sim, tn, x, 
                             alpha, gamma, phi, iter.eval=10, sig.pen)
    gamma=update.gamma$gamma
    update.phi=recur.phi(z.e, rd.sim, tn, x, alpha, gamma, phi)
    phi=update.phi
    
    iter = iter + 1
    
    if(print) 
    {
      cat("===== Iteration", iter, "; ECM step =====", "\n")
      cat("alpha:", "\n")
      print(alpha)
      cat("gamma:", "\n")
      print(gamma)
      cat("phi:", c(phi), sep="\t", "\n")
    }
    param.list[[iter]] = list(alpha=alpha, gamma=gamma, phi=phi)
  }
  return(param.list)
}

analysis.fit = function(z, z0, rd, tn, x, 
                        alpha, gamma, phi,
                        z.valid, z0.valid, rd.valid, tn.valid, x.valid,
                        z.test, z0.test, rd.test, tn.test, x.test,
                        seed.eval=1001) ## unreported claim estimation
{
  set.seed(seed.eval)
  p.log.vec=reg.link(x, coef.fix=alpha, logit=T) 
  mu.log.vec=reg.link(x, coef.fix=gamma, logit=F) 
  p.vec=exp(p.log.vec)
  mu.vec=exp(mu.log.vec)
  gamma.tn.bar=pgamma(tn,shape=1/phi,scale=mu.vec*phi,lower.tail=FALSE,log.p=FALSE)
  gamma.tn=pgamma(tn,shape=1/phi,scale=mu.vec*phi,lower.tail=TRUE,log.p=FALSE)
  
  z.e0=p.vec*gamma.tn.bar/(p.vec*gamma.tn.bar+(1-p.vec))
  ur.e.vec=numeric(length(z))
  ur.e.vec[z==0]=z.e0[z==0]
  ur.e.ind=ur.e.vec
  ur.e.agg=sum(ur.e.ind)
  
  p.obs=pmax(p.vec*gamma.tn,1e-100)
  p.ur=pmax(ur.e.vec,1e-100)
  p.com=pmax(p.vec,1e-100)
  
  ## Validation
  p.log.vec.valid=reg.link(x.valid, coef.fix=alpha, logit=T) 
  mu.log.vec.valid=reg.link(x.valid, coef.fix=gamma, logit=F) 
  p.vec.valid=exp(p.log.vec.valid)
  mu.vec.valid=exp(mu.log.vec.valid)
  gamma.tn.bar.valid=pgamma(tn.valid,shape=1/phi,scale=mu.vec.valid*phi,lower.tail=FALSE,log.p=FALSE)
  gamma.tn.valid=pgamma(tn.valid,shape=1/phi,scale=mu.vec.valid*phi,lower.tail=TRUE,log.p=FALSE)
  
  z.e0.valid=p.vec.valid*gamma.tn.bar.valid/(p.vec.valid*gamma.tn.bar.valid+(1-p.vec.valid))
  ur.e.vec.valid=numeric(length(z.valid))
  ur.e.vec.valid[z.valid==0]=z.e0.valid[z.valid==0]
  ur.e.ind.valid=ur.e.vec.valid
  ur.e.agg.valid=sum(ur.e.ind.valid)
  
  p.obs.valid=pmax(p.vec.valid*gamma.tn.valid,1e-100)
  p.ur.valid=pmax(ur.e.vec.valid,1e-100)
  p.com.valid=pmax(p.vec.valid,1e-100)
  
  ## Testing
  p.log.vec.test=reg.link(x.test, coef.fix=alpha, logit=T) 
  p.vec.test=exp(p.log.vec.test)
  p.com.test=pmax(p.vec.test,1e-100)
  
  analysis.list=list(ur.e.ind=ur.e.ind, ur.e.agg=ur.e.agg,
                     p.obs=p.obs, p.ur=p.ur, p.com=p.com,
                     ur.e.ind.valid=ur.e.ind.valid, ur.e.agg.valid=ur.e.agg.valid,
                     p.obs.valid=p.obs.valid, p.ur.valid=p.ur.valid, p.com.valid=p.com.valid,
                     p.com.test=p.com.test)
  return(analysis.list)
}