reg.link = function(x, coef.fix, coef.random, 
                    match.b, match.s, match.p, ww.b, ww.s, ww.p, 
                    logit=F, effect.fix=NULL) 
{ # coef.fix=alpha(for p) or gamma(for mu); coef.random=beta(for p) or nu(for mu)
  # logit=F (for mu) or T(for p) 
  if(is.null(effect.fix)){effect.fix=(x%*%coef.fix)[,1]}
  effect.random.b = (coef.random[1]*ww.b)[match.b]
  effect.random.s = (coef.random[2]*ww.s)[match.s]
  effect.random.p = (coef.random[3]*ww.p)[match.p]
  effect.random = effect.random.b + effect.random.s + effect.random.p
  effect.total = effect.fix + effect.random
  if(logit==T){return(effect.total-log1pexp(effect.total))} else {
    return(effect.total)
  }
}

loglik = function(z, rd, tn, x, alpha, beta, gamma, nu, phi, rho,
                  match.b, match.s, match.p, ww.b, ww.s, ww.p, ww.length.multi,
                  effect.fix.z=NULL, effect.fix.rd=NULL)
{
  if(is.null(effect.fix.z)){effect.fix.z=(x%*%alpha)[,1]}
  if(is.null(effect.fix.rd)){effect.fix.rd=(x%*%gamma)[,1]}
  p.log.vec=pmin(reg.link(x, coef.fix=alpha, coef.random=beta, 
                     match.b, match.s, match.p, ww.b, ww.s, ww.p, 
                     logit=T, effect.fix=effect.fix.z),-1e-100)
  mu.log.vec=reg.link(x, coef.fix=gamma, coef.random=nu, 
                      match.b, match.s, match.p, ww.b, ww.s, ww.p, 
                      logit=F, effect.fix=effect.fix.rd) 
  mu.vec=exp(mu.log.vec)
  
  ll.gamma=ll.gamma.tn.bar=numeric(length(z))
  ll.gamma[z==1]=dgamma(rd[z==1],shape=1/phi,scale=mu.vec[z==1]*phi,log=TRUE)
  ll.gamma.tn.bar[z==0]=pgamma(tn[z==0],shape=1/phi,scale=mu.vec[z==0]*phi,lower.tail=FALSE,log.p=TRUE)
  p.bar.log.vec=log1mexp(-p.log.vec)
  
  ll.obs=z*(p.log.vec+ll.gamma)+(1-z)*(p.bar.log.vec+log1pexp(p.log.vec-p.bar.log.vec+ll.gamma.tn.bar))
  ll.lat.multi=dmvnorm(cbind(ww.b[1:ww.length.multi],ww.s[1:ww.length.multi]),mean=rep(0,2),
                       sigma=array(c(1,rho,rho,1),dim=c(2,2)),log=T)
  ll.lat.uni.b=dnorm(ww.b[-(1:ww.length.multi)],mean=0,sd=1,log=T)
  ll.lat.uni.s=dnorm(ww.s[-(1:ww.length.multi)],mean=0,sd=1,log=T)
  ll.lat.p=dnorm(ww.p,mean=0,sd=1,log=T)
  ll=sum(ll.obs)+sum(ll.lat.multi)+sum(ll.lat.uni.b)+sum(ll.lat.uni.s)+sum(ll.lat.p)
  list(ll.obs=ll.obs, ll.lat.multi=ll.lat.multi, ll.lat.uni.b=ll.lat.uni.b,
       ll.lat.uni.s=ll.lat.uni.s, ll.lat.p=ll.lat.p, ll=ll)
}


ww.posterior.sim = function(z, rd, tn, x, alpha, beta, gamma, nu, phi, rho,
                            t.mat, match.b, match.s, match.p, ww.b, ww.s, ww.p, 
                            unique.b, unique.s, unique.p, ww.length.multi,
                            MH.sd.ww, num.sim=1000, seed.eval=1001)
{
  set.seed(seed.eval)
  effect.fix.z=(x%*%alpha)[,1]
  effect.fix.rd=(x%*%gamma)[,1]
  ww.length.b=length(ww.b)
  ww.length.s=length(ww.s)
  ww.length.p=length(ww.p)
  ww.MH.matrix.b=array(0,dim=c(ww.length.b,num.sim))
  ww.MH.matrix.s=array(0,dim=c(ww.length.s,num.sim))
  ww.MH.matrix.p=array(0,dim=c(ww.length.p,num.sim))
  ww.MH.matrix.b[,1]=ww.b
  ww.MH.matrix.s[,1]=ww.s
  ww.MH.matrix.p[,1]=ww.p
  
  loglik.list=loglik(z, rd, tn, x, alpha, beta, gamma, nu, phi, rho,
                     match.b, match.s, match.p, ww.b, ww.s, ww.p, ww.length.multi,
                     effect.fix.z, effect.fix.rd)
  loglik.b.dataframe = data.frame(ll=loglik.list$ll.obs, idx=t.mat[,1])
  loglik.s.dataframe = data.frame(ll=loglik.list$ll.obs, idx=t.mat[,2])
  loglik.p.dataframe = data.frame(ll=loglik.list$ll.obs, idx=t.mat[,3])
  ll.obs.agg.b  =setDT(loglik.b.dataframe)[,.(sum(ll)),by=.(idx)]
  ll.obs.agg.s  =setDT(loglik.s.dataframe)[,.(sum(ll)),by=.(idx)]
  ll.obs.agg.p  =setDT(loglik.p.dataframe)[,.(sum(ll)),by=.(idx)]
  match.ll.b = match(unique.b,ll.obs.agg.b$idx)
  match.ll.s = match(unique.s,ll.obs.agg.s$idx)
  match.ll.p = match(unique.p,ll.obs.agg.p$idx)
  
  for (s in 1:(num.sim-1))
  {
    ww.old.b=ww.MH.matrix.b[,s+1]=ww.MH.matrix.b[,s]
    ww.old.s=ww.MH.matrix.s[,s+1]=ww.MH.matrix.s[,s]
    ww.old.p=ww.MH.matrix.p[,s+1]=ww.MH.matrix.p[,s]
    ## buyer
    ww.trial.b=ww.old.b+rnorm(ww.length.b,mean=0,sd=MH.sd.ww)
    loglik.list.old.b=loglik(z, rd, tn, x, alpha, beta, gamma, nu, phi, rho,
                             match.b, match.s, match.p, ww.old.b, ww.old.s, ww.old.p, ww.length.multi,
                             effect.fix.z, effect.fix.rd)
    loglik.list.trial.b=loglik(z, rd, tn, x, alpha, beta, gamma, nu, phi, rho,
                               match.b, match.s, match.p, ww.trial.b, ww.old.s, ww.old.p, ww.length.multi,
                               effect.fix.z, effect.fix.rd)
    loglik.old.b.dataframe = data.frame(ll=loglik.list.old.b$ll.obs, idx=t.mat[,1])
    loglik.trial.b.dataframe = data.frame(ll=loglik.list.trial.b$ll.obs, idx=t.mat[,1])
    #ll.obs.agg.old.b  =aggregate(loglik.list.old.b$ll.obs,by=list(t.mat[,1]),FUN=sum)
    ll.obs.agg.old.b  =setDT(loglik.old.b.dataframe)[,.(sum(ll)),by=.(idx)]
    loglik.ind.old.b  =ll.obs.agg.old.b$V1[match.ll.b]+
      c(loglik.list.old.b$ll.lat.multi, loglik.list.old.b$ll.lat.uni.b)
    #ll.obs.agg.trial.b=aggregate(loglik.list.trial.b$ll.obs,by=list(t.mat[,1]),FUN=sum)
    ll.obs.agg.trial.b  =setDT(loglik.trial.b.dataframe)[,.(sum(ll)),by=.(idx)]
    loglik.ind.trial.b=ll.obs.agg.trial.b$V1[match.ll.b]+
      c(loglik.list.trial.b$ll.lat.multi, loglik.list.trial.b$ll.lat.uni.b)
    accept.idx=(runif(ww.length.b)<exp(loglik.ind.trial.b-loglik.ind.old.b))
    ww.MH.matrix.b[accept.idx,s+1]=ww.trial.b[accept.idx]
    ww.old.b=ww.MH.matrix.b[,s+1]
    ## seller
    ww.trial.s=ww.old.s+rnorm(ww.length.s,mean=0,sd=MH.sd.ww)
    loglik.list.old.s=loglik(z, rd, tn, x, alpha, beta, gamma, nu, phi, rho,
                             match.b, match.s, match.p, ww.old.b, ww.old.s, ww.old.p, ww.length.multi,
                             effect.fix.z, effect.fix.rd)
    loglik.list.trial.s=loglik(z, rd, tn, x, alpha, beta, gamma, nu, phi, rho,
                               match.b, match.s, match.p, ww.old.b, ww.trial.s, ww.old.p, ww.length.multi,
                               effect.fix.z, effect.fix.rd)
    loglik.old.s.dataframe = data.frame(ll=loglik.list.old.s$ll.obs, idx=t.mat[,2])
    loglik.trial.s.dataframe = data.frame(ll=loglik.list.trial.s$ll.obs, idx=t.mat[,2])
    #ll.obs.agg.old.s  =aggregate(loglik.list.old.s$ll.obs,by=list(t.mat[,2]),FUN=sum)
    ll.obs.agg.old.s  =setDT(loglik.old.s.dataframe)[,.(sum(ll)),by=.(idx)]
    loglik.ind.old.s  =ll.obs.agg.old.s$V1[match.ll.s]+
      c(loglik.list.old.s$ll.lat.multi, loglik.list.old.s$ll.lat.uni.s)
    #ll.obs.agg.trial.s=aggregate(loglik.list.trial.s$ll.obs,by=list(t.mat[,2]),FUN=sum)
    ll.obs.agg.trial.s=setDT(loglik.trial.s.dataframe)[,.(sum(ll)),by=.(idx)]
    loglik.ind.trial.s=ll.obs.agg.trial.s$V1[match.ll.s]+
      c(loglik.list.trial.s$ll.lat.multi, loglik.list.trial.s$ll.lat.uni.s)
    accept.idx=(runif(ww.length.s)<exp(loglik.ind.trial.s-loglik.ind.old.s))
    ww.MH.matrix.s[accept.idx,s+1]=ww.trial.s[accept.idx]
    ww.old.s=ww.MH.matrix.s[,s+1]
    ## policy
    ww.trial.p=ww.old.p+rnorm(ww.length.p,mean=0,sd=MH.sd.ww)
    loglik.list.old.p=loglik(z, rd, tn, x, alpha, beta, gamma, nu, phi, rho,
                             match.b, match.s, match.p, ww.old.b, ww.old.s, ww.old.p, ww.length.multi,
                             effect.fix.z, effect.fix.rd)
    loglik.list.trial.p=loglik(z, rd, tn, x, alpha, beta, gamma, nu, phi, rho,
                               match.b, match.s, match.p, ww.old.b, ww.old.s, ww.trial.p, ww.length.multi,
                               effect.fix.z, effect.fix.rd)
    loglik.old.p.dataframe = data.frame(ll=loglik.list.old.p$ll.obs, idx=t.mat[,3])
    loglik.trial.p.dataframe = data.frame(ll=loglik.list.trial.p$ll.obs, idx=t.mat[,3])
    #ll.obs.agg.old.p  =aggregate(loglik.list.old.p$ll.obs,by=list(t.mat[,3]),FUN=sum)
    ll.obs.agg.old.p  =setDT(loglik.old.p.dataframe)[,.(sum(ll)),by=.(idx)]
    loglik.ind.old.p  =ll.obs.agg.old.p$V1[match.ll.p]+
      loglik.list.old.p$ll.lat.p
    #ll.obs.agg.trial.p=aggregate(loglik.list.trial.p$ll.obs,by=list(t.mat[,3]),FUN=sum)
    ll.obs.agg.trial.p=setDT(loglik.trial.p.dataframe)[,.(sum(ll)),by=.(idx)]
    loglik.ind.trial.p=ll.obs.agg.trial.p$V1[match.ll.p]+
      loglik.list.trial.p$ll.lat.p
    accept.idx=(runif(ww.length.p)<exp(loglik.ind.trial.p-loglik.ind.old.p))
    ww.MH.matrix.p[accept.idx,s+1]=ww.trial.p[accept.idx]
    ww.old.p=ww.MH.matrix.p[,s+1]
  }
  accept.prob.b=apply(ww.MH.matrix.b,1,FUN=function(x){length(unique(x))})/num.sim
  accept.prob.s=apply(ww.MH.matrix.s,1,FUN=function(x){length(unique(x))})/num.sim
  accept.prob.p=apply(ww.MH.matrix.p,1,FUN=function(x){length(unique(x))})/num.sim
  list(ww.MH.matrix.b=ww.MH.matrix.b, ww.MH.matrix.s=ww.MH.matrix.s, ww.MH.matrix.p=ww.MH.matrix.p, 
       accept.prob.b=accept.prob.b, accept.prob.s=accept.prob.s, accept.prob.p=accept.prob.p)
}

rd.posterior.sim = function(z, rd, tn, x, gamma, nu, phi, rho, 
                            match.b, match.s, match.p, ww.b, ww.s, ww.p, 
                            seed.eval=1001)
{
  set.seed(seed.eval)
  z.length=length(z)
  mu.log.vec=reg.link(x, coef.fix=gamma, coef.random=nu, 
                      match.b, match.s, match.p, ww.b, ww.s, ww.p, logit=F) 
  mu.vec=exp(mu.log.vec)
  ll.gamma.tn.bar=pgamma(tn,shape=1/phi,scale=mu.vec*phi,lower.tail=FALSE,log.p=TRUE)
  rd.sim=numeric(z.length)
  rd.sim[(z==1)]=rd[(z==1)]
  rd.quantile.log=ll.gamma.tn.bar+log(1-runif(z.length))
  rd.sim[(z==0)]=qgamma(rd.quantile.log[(z==0)],shape=1/phi,scale=mu.vec[(z==0)]*phi,
                        lower.tail=FALSE,log.p=TRUE)
  return(rd.sim)
}

recur.z = function(z, rd, tn, x, alpha, beta, gamma, nu, phi, rho,
                   match.b, match.s, match.p, ww.b, ww.s, ww.p,
                   effect.fix.z=NULL, effect.fix.rd=NULL)
{
  if(is.null(effect.fix.z)){effect.fix.z=(x%*%alpha)[,1]}
  if(is.null(effect.fix.rd)){effect.fix.rd=(x%*%gamma)[,1]}
  p.log.vec=reg.link(x, coef.fix=alpha, coef.random=beta, 
                     match.b, match.s, match.p, ww.b, ww.s, ww.p, 
                     logit=T, effect.fix=effect.fix.z) 
  mu.log.vec=reg.link(x, coef.fix=gamma, coef.random=nu, 
                      match.b, match.s, match.p, ww.b, ww.s, ww.p,
                      logit=F, effect.fix=effect.fix.rd) 
  p.vec=exp(p.log.vec)
  mu.vec=exp(mu.log.vec)
  gamma.tn.bar=pgamma(tn,shape=1/phi,scale=mu.vec*phi,lower.tail=FALSE,log.p=FALSE)
  
  z.e0=p.vec*gamma.tn.bar/(p.vec*gamma.tn.bar+(1-p.vec))
  z.e=rep(1,length(z))
  z.e[z==0]=z.e0[z==0]
  
  return(z.e)
}

recur.logit = function(z, rd, tn, x, alpha, beta, gamma, nu, phi, rho,
                       match.b, match.s, match.p, ww.b, ww.s, ww.p,
                       iter.eval=20, sig.pen=NULL) ##alpha and beta
{
  alpha.length=length(alpha)
  param=c(alpha,beta)
  param.new=param
  param.old=param-Inf
  w.b=ww.b[match.b]
  w.s=ww.s[match.s]
  w.p=ww.p[match.p]
  x.eval=cbind(x, w.b, w.s, w.p)
  iter=0
  while ((iter<=iter.eval)&(sum((param.old-param.new)^2)>10^(-8)))
  {
    param.old=param.new
    alpha.new=param.new[1:alpha.length]
    beta.new=param.new[-(1:alpha.length)]
    p.vec=exp(reg.link(x, coef.fix=alpha.new, coef.random=beta.new,
                       match.b, match.s, match.p, ww.b, ww.s, ww.p, logit=T))
    dQ = apply(sweep(x.eval,1,z-p.vec,FUN="*"),2,sum) -if(!is.null(sig.pen)){c(0, param.new[-1]/sig.pen^2)} else{0}
    dQ2 = -crossprod(sweep(x.eval,1,p.vec*(1-p.vec),FUN="*"),x.eval) -if(!is.null(sig.pen)){diag(c(0,1/sig.pen^2))} else{0}
    param.new=param.new-crossprod(dQ,solve(dQ2))
    iter = iter+1
  }
  alpha.new=param.new[1:alpha.length]
  beta.new=param.new[-(1:alpha.length)]
  list(alpha=alpha.new, beta=beta.new, iter=iter)
}

recur.gamma = function(z, rd, tn, x, alpha, beta, gamma, nu, phi, rho,
                       match.b, match.s, match.p, ww.b, ww.s, ww.p,
                       iter.eval=20, sig.pen=NULL) ##gamma and nu
{
  gamma.length=length(gamma)
  param=c(gamma,nu)
  param.new=param
  param.old=param-Inf
  w.b=ww.b[match.b]
  w.s=ww.s[match.s]
  w.p=ww.p[match.p]
  x.eval=cbind(x, w.b, w.s, w.p)
  iter=0
  while ((iter<=iter.eval)&(sum((param.old-param.new)^2)>10^(-8)))
  {
    param.old=param.new
    gamma.new=param.new[1:gamma.length]
    nu.new=param.new[-(1:gamma.length)]
    mu.log.vec=reg.link(x, coef.fix=gamma.new, coef.random=nu.new, 
                        match.b, match.s, match.p, ww.b, ww.s, ww.p, logit=F) 
    mu.vec=exp(mu.log.vec)
    dQ = apply(sweep(x.eval,1,z*(-1+rd/mu.vec),FUN="*"),2,sum)-if(!is.null(sig.pen)){c(0, param.new[-1]/sig.pen^2)} else{0}
    dQ2 = -crossprod(sweep(x.eval,1,z*rd/mu.vec,FUN="*"),x.eval)-if(!is.null(sig.pen)){diag(c(0,1/sig.pen^2))} else{0}
    param.new=param.new-crossprod(dQ,solve(dQ2))
    iter = iter+1
  }
  gamma.new=param.new[1:gamma.length]
  nu.new=param.new[-(1:gamma.length)]
  list(gamma=gamma.new, nu=nu.new, iter=iter)
}

Q.phi = function(phi, term.1, term.2)
{
  return(-(term.1*(-phi*lgamma(1/phi)-log(phi))+term.2)/phi)
}

recur.phi = function(z, rd, tn, x, alpha, beta, gamma, nu, phi, rho,
                     match.b, match.s, match.p, ww.b, ww.s, ww.p)
{
  mu.log.vec=reg.link(x, coef.fix=gamma, coef.random=nu, 
                      match.b, match.s, match.p, ww.b, ww.s, ww.p, logit=F) 
  mu.vec=exp(mu.log.vec)
  term.1=sum(z)
  term.2=sum(z*(log(rd)-log(mu.vec)-rd/mu.vec))
  phi.opt=optimize(f=Q.phi, interval=c(0.5*phi,2*phi), term.1=term.1, term.2=term.2)$minimum
  return(phi.opt)
}

Q.rho = function(rho, MS.bb, MS.ss, MS.bs) {
  return((-rho^3+(1-MS.bb-MS.ss)*rho+MS.bs*(rho^2+1))^2)
}

recur.rho = function(rho, MS.bb, MS.ss, MS.bs){
  rho.opt=optimize(f=Q.rho, interval=c(max(-0.99,rho-0.05),min(0.99,rho+0.05)),
                   MS.bb=MS.bb, MS.ss=MS.ss, MS.bs=MS.bs)$minimum
  return(rho.opt)
}


ecm = function(z, rd, tn, x, alpha, beta, gamma, nu, phi, rho,
               t.mat, ww.b, ww.s, ww.p, 
               MH.sd.ww=1, iter.max=100, num.sim=200, num.burn=20, num.sim.e=10,
               sig.pen=NULL, seed.num=999, print=TRUE)
{
  n=length(z)
  unique.b=unique(t.mat[,1])
  unique.s=unique(t.mat[,2])
  unique.multi=unique.b[unique.b%in%unique.s]
  unique.uni.b=unique.b[!(unique.b%in%unique.s)]
  unique.uni.s=unique.s[!(unique.s%in%unique.b)]
  unique.b=c(unique.multi,unique.uni.b)
  unique.s=c(unique.multi,unique.uni.s)
  unique.p=unique(t.mat[,3])
  match.b=match(t.mat[,1],unique.b)
  match.s=match(t.mat[,2],unique.s)
  match.p=match(t.mat[,3],unique.p)
  ww.length.multi=length(unique.multi)
  iter=0
  param.list=list()
  while(iter<iter.max)
  {
    #loglik.em.old = loglik.em
    alpha.old=alpha
    beta.old=beta
    gamma.old=gamma
    nu.old=nu
    ww.b.old=ww.b
    ww.s.old=ww.s
    ww.p.old=ww.p
    phi.old=phi
    rho.old=rho
    ## E-step
    ww.list = ww.posterior.sim(z, rd, tn, x, alpha, beta, gamma, nu, phi, rho,
                               t.mat, match.b, match.s, match.p, ww.b, ww.s, ww.p, 
                               unique.b, unique.s, unique.p, ww.length.multi, 
                               MH.sd.ww, num.sim=num.sim, seed.eval=(seed.num*iter))
    ww.matrix.b = ww.list$ww.MH.matrix.b
    ww.matrix.s = ww.list$ww.MH.matrix.s
    ww.matrix.p = ww.list$ww.MH.matrix.p
    
    z.trial=rep(z,num.sim.e)
    tn.trial=rep(tn,num.sim.e)
    x.trial=do.call("rbind", replicate(num.sim.e, x, simplify = FALSE))
    rd.trial=t.mat.trial=ww.b.trial=ww.s.trial=ww.p.trial=c()
    t.mat.trial=array("a",dim=c(num.sim.e*n,3))
    for (s in 1:num.sim.e)
    {
      for (u in 1:3) {t.mat.trial[((s-1)*n+1):(s*n),u]=paste(t.mat[,u], "&", s)}
      ww.sel.idx=num.burn+ceiling((num.sim-num.burn)*s/num.sim.e)
      ww.b.append=ww.matrix.b[,ww.sel.idx]
      ww.s.append=ww.matrix.s[,ww.sel.idx]
      ww.p.append=ww.matrix.p[,ww.sel.idx]
      ww.b.trial=c(ww.b.trial,ww.b.append)
      ww.s.trial=c(ww.s.trial,ww.s.append)
      ww.p.trial=c(ww.p.trial,ww.p.append)
      
      rd.append=rd.posterior.sim(z, rd, tn, x, gamma, nu, phi, rho,
                                 match.b, match.s, match.p, 
                                 ww.b.append, ww.s.append, ww.p.append,
                                 seed.eval=((999+s)*seed.num*(iter+1)))
      rd.trial=c(rd.trial,rd.append)
    }
    
    length.ww.b=length(ww.b)
    length.ww.s=length(ww.s)
    length.ww.p=length(ww.p)
    match.b.trial=match.b
    match.s.trial=match.s
    match.p.trial=match.p
    for (s in 2:num.sim.e) {
      match.b.trial = c(match.b.trial, match.b+(s-1)*length.ww.b)
      match.s.trial = c(match.s.trial, match.s+(s-1)*length.ww.s)
      match.p.trial = c(match.p.trial, match.p+(s-1)*length.ww.p)
    }
    
    z.e.trial= recur.z(z.trial, rd.trial, tn.trial, x.trial, 
                       alpha, beta, gamma, nu, phi, rho,
                       match.b.trial, match.s.trial, match.p.trial, 
                       ww.b.trial, ww.s.trial, ww.p.trial)
    
    ## M-step
    update.logit=recur.logit(z.e.trial, rd.trial, tn.trial, x.trial, 
                             alpha, beta, gamma, nu, phi, rho,
                             match.b.trial, match.s.trial, match.p.trial, 
                             ww.b.trial, ww.s.trial, ww.p.trial, iter.eval=10, sig.pen)
    alpha=update.logit$alpha
    beta=update.logit$beta
    update.gamma=recur.gamma(z.e.trial, rd.trial, tn.trial, x.trial, 
                             alpha, beta, gamma, nu, phi, rho, 
                             match.b.trial, match.s.trial, match.p.trial, 
                             ww.b.trial, ww.s.trial, ww.p.trial, iter.eval=10, sig.pen)
    gamma=update.gamma$gamma
    nu=update.gamma$nu
    update.phi=recur.phi(z.e.trial, rd.trial, tn.trial, x.trial, 
                         alpha, beta, gamma, nu, phi,rho,
                         match.b.trial, match.s.trial, match.p.trial, 
                         ww.b.trial, ww.s.trial, ww.p.trial)
    phi=update.phi
    
    ww.sel.seq=num.burn+ceiling((num.sim-num.burn)*seq(num.sim.e)/num.sim.e)
    ww.b.multi.trial=as.vector(ww.matrix.b[1:ww.length.multi,ww.sel.seq])
    ww.s.multi.trial=as.vector(ww.matrix.s[1:ww.length.multi,ww.sel.seq])
    MS.bb=sum(ww.b.multi.trial^2)/length(ww.b.multi.trial)
    MS.ss=sum(ww.s.multi.trial^2)/length(ww.b.multi.trial)
    MS.bs=sum(ww.b.multi.trial*ww.s.multi.trial)/length(ww.b.multi.trial)
    rho.init=cor(ww.b.multi.trial,ww.s.multi.trial)
    update.rho=recur.rho(rho.init, MS.bb, MS.ss, MS.bs)
    rho=update.rho
    
    iter = iter + 1
    
    if(print) 
    {
      cat("===== Iteration", iter, "; ECM step =====", "\n")
      cat("alpha:", "\n")
      print(alpha)
      cat("beta:", "\n")
      print(beta)
      cat("gamma:", "\n")
      print(gamma)
      cat("nu:", "\n")
      print(nu)
      cat("phi:", c(phi), sep="\t", "\n")
      cat("rho:", c(rho), sep="\t", "\n")
      #cat("loglik:", loglik.em, sep="\t", "\n")
    }
    param.list[[iter]] = list(alpha=alpha, beta=beta, gamma=gamma, nu=nu, phi=phi, rho=rho)
  }
  return(param.list)
}


### Observed data log-likelihood (simulated using importance sampling after MH)
loglik.obs.sim.IS = function(z, rd, tn, x, alpha, beta, gamma, nu, phi, rho,
                             t.mat, match.b, match.s, match.p, ww.b, ww.s, ww.p, ww.length.multi,
                             num.sim=1000, seed.eval=1001, 
                             ww.b.mean.post, ww.b.sd.post, ww.s.mean.post, ww.s.sd.post,
                             ww.p.mean.post, ww.p.sd.post, ww.cor.post)
{
  set.seed(seed.eval)
  effect.fix.z=(x%*%alpha)[,1]
  effect.fix.rd=(x%*%gamma)[,1]
  length.ww.b=length(ww.b)
  length.ww.s=length(ww.s)
  length.ww.p=length(ww.p)
  ll.obs.ind.naive=numeric(num.sim)
  for (s in 1:num.sim)
  {
    ww.b.trial=rnorm(ww.length.multi)
    ww.s.trial=rnorm(ww.length.multi)
    ww.trial=rnorm(ww.length.multi)
    ww.b.std=sqrt(1-ww.cor.post)*ww.b.trial+sign(ww.cor.post)*sqrt(abs(ww.cor.post))*ww.trial
    ww.s.std=sqrt(1-ww.cor.post)*ww.s.trial+sign(ww.cor.post)*sqrt(abs(ww.cor.post))*ww.trial
    ww.b[1:ww.length.multi]=ww.b.mean.post[1:ww.length.multi]+ww.b.std*ww.b.sd.post[1:ww.length.multi]
    ww.s[1:ww.length.multi]=ww.s.mean.post[1:ww.length.multi]+ww.s.std*ww.s.sd.post[1:ww.length.multi]
    ww.b[-(1:ww.length.multi)]=ww.b.mean.post[-(1:ww.length.multi)]+rnorm(length.ww.b-ww.length.multi)*ww.b.sd.post[-(1:ww.length.multi)]
    ww.s[-(1:ww.length.multi)]=ww.s.mean.post[-(1:ww.length.multi)]+rnorm(length.ww.s-ww.length.multi)*ww.s.sd.post[-(1:ww.length.multi)]
    ww.p=ww.p.mean.post+rnorm(length.ww.p)*ww.p.sd.post
    loglik.obs.list=loglik(z, rd, tn, x, alpha, beta, gamma, nu, phi, rho,
                           match.b, match.s, match.p, ww.b, ww.s, ww.p, ww.length.multi,
                           effect.fix.z, effect.fix.rd)
    ll.ww.multi.post=log(as.numeric(dnorm2d(ww.b.std,ww.s.std,rho=ww.cor.post)))
    ll.ww.uni.b=dnorm(ww.b[-(1:ww.length.multi)],mean=ww.b.mean.post[-(1:ww.length.multi)],sd=ww.b.sd.post[-(1:ww.length.multi)],log=T)
    ll.ww.uni.s=dnorm(ww.s[-(1:ww.length.multi)],mean=ww.s.mean.post[-(1:ww.length.multi)],sd=ww.s.sd.post[-(1:ww.length.multi)],log=T)
    ll.ww.p=dnorm(ww.p,mean=ww.p.mean.post,sd=ww.p.sd.post,log=T)
    ll.obs.ind.naive[s]=loglik.obs.list$ll-sum(ll.ww.multi.post)-sum(ll.ww.uni.b)-sum(ll.ww.uni.s)-sum(ll.ww.p)
  }
  ll.obs.aggre.naive=logSumExp(ll.obs.ind.naive)-log(num.sim)
  list(ll.obs.ind.naive=ll.obs.ind.naive, ll.obs.aggre.naive=ll.obs.aggre.naive)
}

loglik.obs.sim.full.IS = function(z.full, rd.full, tn.full, x.full, 
                                  alpha, beta, gamma, nu, phi, rho,
                                  unique.b, unique.s, unique.p,
                                  unique.b.full, unique.s.full, unique.p.full,
                                  match.b.full, match.s.full, match.p.full, 
                                  ww.b.full, ww.s.full, ww.p.full, 
                                  ww.length.multi.full,
                                  num.sim=1000, seed.eval=1001, 
                                  ww.b.mean.post, ww.b.sd.post, 
                                  ww.s.mean.post, ww.s.sd.post,
                                  ww.p.mean.post, ww.p.sd.post)
{
  set.seed(seed.eval)
  effect.fix.z.full=(x.full%*%alpha)[,1]
  effect.fix.rd.full=(x.full%*%gamma)[,1]
  length.ww.b.full=length(ww.b.full)
  length.ww.s.full=length(ww.s.full)
  length.ww.p.full=length(ww.p.full)
  match.b.unique=match(unique.b.full,unique.b)
  match.s.unique=match(unique.s.full,unique.s)
  match.p.unique=match(unique.p.full,unique.p)
  ll.obs.ind.val=ll.obs.ind.adj=numeric(num.sim)
  for (s in 1:num.sim)
  {
    ww.b.mean.post.full=ww.b.mean.post[match.b.unique]
    ww.s.mean.post.full=ww.s.mean.post[match.s.unique]
    ww.p.mean.post.full=ww.p.mean.post[match.p.unique]
    ww.b.sd.post.full=ww.b.sd.post[match.b.unique]
    ww.s.sd.post.full=ww.s.sd.post[match.s.unique]
    ww.p.sd.post.full=ww.p.sd.post[match.p.unique]
    ww.b.mean.post.full[is.na(ww.b.mean.post.full)]=0
    ww.s.mean.post.full[is.na(ww.s.mean.post.full)]=0
    ww.p.mean.post.full[is.na(ww.p.mean.post.full)]=0
    ww.b.sd.post.full[is.na(ww.b.sd.post.full)]=1
    ww.s.sd.post.full[is.na(ww.s.sd.post.full)]=1
    ww.p.sd.post.full[is.na(ww.p.sd.post.full)]=1
    
    ww.b.full=ww.b.mean.post.full+rnorm(length.ww.b.full)*ww.b.sd.post.full
    ww.s.full=ww.s.mean.post.full+rnorm(length.ww.s.full)*ww.s.sd.post.full
    ww.p.full=ww.p.mean.post.full+rnorm(length.ww.p.full)*ww.p.sd.post.full
    
    loglik.obs.list=loglik(z=z.full, rd=rd.full, tn=tn.full, x=tn.full, 
                           alpha, beta, gamma, nu, phi, rho,
                           match.b=match.b.full, match.s=match.s.full, match.p=match.p.full, 
                           ww.b=ww.b.full, ww.s=ww.s.full, ww.p=ww.p.full, 
                           ww.length.multi=ww.length.multi.full,
                           effect.fix.z=effect.fix.z.full, effect.fix.rd=effect.fix.rd.full)
    ll.ww.b=dnorm(ww.b.full,mean=ww.b.mean.post.full,sd=ww.b.sd.post.full,log=T)
    ll.ww.s=dnorm(ww.s.full,mean=ww.s.mean.post.full,sd=ww.s.sd.post.full,log=T)
    ll.ww.p=dnorm(ww.p.full,mean=ww.p.mean.post.full,sd=ww.p.sd.post.full,log=T)
    ll.obs.ind.val[s]=loglik.obs.list$ll
    ll.obs.ind.adj[s]=-sum(ll.ww.b)-sum(ll.ww.s)-sum(ll.ww.p)
  }
  ll.obs.ind.naive=ll.obs.ind.val+ll.obs.ind.adj
  ll.obs.aggre.naive=logSumExp(ll.obs.ind.naive)-log(num.sim)
  list(ll.obs.ind.val=ll.obs.ind.val, ll.obs.ind.adj=ll.obs.ind.adj,
       ll.obs.ind.naive=ll.obs.ind.naive, ll.obs.aggre.naive=ll.obs.aggre.naive)
}

### Analysis of fit
reg.link.matrix = function(x, coef.fix, coef.random, 
                    match.b, match.s, match.p, ww.b.matrix, ww.s.matrix, ww.p.matrix, 
                    logit=F, effect.fix=NULL) 
{ # coef.fix=alpha(for p) or gamma(for mu); coef.random=beta(for p) or nu(for mu)
  # logit=F (for mu) or T(for p) 
  if(is.null(effect.fix)){effect.fix=(x%*%coef.fix)[,1]}
  effect.random.b = (coef.random[1]*ww.b.matrix)[match.b,]
  effect.random.s = (coef.random[2]*ww.s.matrix)[match.s,]
  effect.random.p = (coef.random[3]*ww.p.matrix)[match.p,]
  effect.random = effect.random.b + effect.random.s + effect.random.p
  effect.total = sweep(effect.random,1,effect.fix,FUN="+")
  if(logit==T){return(effect.total-log1pexp(effect.total))} else {
    return(effect.total)
  }
}

reg.link.matrix.OS = function(x, coef.fix, coef.random, 
                           w.b.matrix, w.s.matrix, w.p.matrix, 
                           logit=F, effect.fix=NULL) 
{ # coef.fix=alpha(for p) or gamma(for mu); coef.random=beta(for p) or nu(for mu)
  # logit=F (for mu) or T(for p) 
  if(is.null(effect.fix)){effect.fix=(x%*%coef.fix)[,1]}
  effect.random.b = coef.random[1]*w.b.matrix
  effect.random.s = coef.random[2]*w.s.matrix
  effect.random.p = coef.random[3]*w.p.matrix
  effect.random = effect.random.b + effect.random.s + effect.random.p
  effect.total = sweep(effect.random,1,effect.fix,FUN="+")
  if(logit==T){return(effect.total-log1pexp(effect.total))} else {
    return(effect.total)
  }
}

analysis.fit = function(z, z0, rd, tn, x, t.mat,
                        alpha, beta, gamma, nu, phi, rho,
                        match.b, match.s, match.p, 
                        ww.b.matrix, ww.s.matrix, ww.p.matrix,
                        z.valid, z0.valid, rd.valid, tn.valid, x.valid, t.mat.valid,
                        z.test, z0.test, rd.test, tn.test, x.test, t.mat.test, 
                        seed.eval=1001) ## unreported claim estimation
{
  match.b.valid=match(t.mat.valid[,1],unique.b)
  match.s.valid=match(t.mat.valid[,2],unique.s)
  match.p.valid=match(t.mat.valid[,3],unique.p)
  match.b.test=match(t.mat.test[,1],unique.b)
  match.s.test=match(t.mat.test[,2],unique.s)
  match.p.test=match(t.mat.test[,3],unique.p)
  
  set.seed(seed.eval)
  ## Training ##
  p.matrix=exp(reg.link.matrix(x, coef.fix=alpha, coef.random=beta, 
                               match.b, match.s, match.p, ww.b.matrix, ww.s.matrix, ww.p.matrix, 
                               logit=T))
  mu.matrix=exp(reg.link.matrix(x, coef.fix=gamma, coef.random=nu, 
                                match.b, match.s, match.p, ww.b.matrix, ww.s.matrix, ww.p.matrix,
                                logit=F))
  gamma.tn.bar=pgamma(tn,shape=1/phi,scale=mu.matrix*phi,lower.tail=FALSE,log.p=FALSE)
  gamma.tn=pgamma(tn,shape=1/phi,scale=mu.matrix*phi,lower.tail=TRUE,log.p=FALSE)
  
  z.e0=p.matrix*gamma.tn.bar/(p.matrix*gamma.tn.bar+(1-p.matrix))
  ur.e.matrix=array(0,dim=dim(p.matrix))
  ur.e.matrix[z==0,]=z.e0[z==0,]
  ur.e.ind=apply(ur.e.matrix,1,mean)
  ur.e.agg=sum(ur.e.ind)
  
  p.obs=pmax(apply(p.matrix*gamma.tn,1,mean),1e-100)
  p.ur=pmax(apply(ur.e.matrix,1,mean),1e-100)
  p.com=pmax(apply(p.matrix,1,mean),1e-100)
  
  ## Validation ##
  unique.bs.valid=unique(t.mat.valid[,c(1,2)])
  unique.bs.new.valid=unique.bs.valid[(!(unique.bs.valid%in%t.mat[,1]))&(!(unique.bs.valid%in%t.mat[,2]))]
  w.bs.new.valid=array(rmvnorm(length(unique.bs.new.valid)*ncol(p.matrix),mean=c(0,0),sigma=array(c(1,rho,rho,1),dim=c(2,2))),
                      dim=c(length(unique.bs.new.valid),ncol(p.matrix),2))
  match.bs.valid=match(t.mat.valid[,1],unique.s)
  match.sb.valid=match(t.mat.valid[,2],unique.b)
  w.b.matrix.valid=array(NA,dim=c(length(z.valid),ncol(p.matrix)))
  w.b.idx.1.valid=(!is.na(match.b.valid))
  w.b.idx.2.valid=(is.na(match.b.valid))&(!is.na(match.bs.valid))
  w.b.idx.3.valid=(is.na(match.b.valid))&(is.na(match.bs.valid))
  w.b.matrix.valid[w.b.idx.1.valid,]=ww.b.matrix[match.b.valid[w.b.idx.1.valid],]
  w.b.matrix.valid[w.b.idx.2.valid,]=sqrt(1-rho^2)*array(rnorm(sum(w.b.idx.2.valid)*ncol(p.matrix)),
                                                       dim=c(sum(w.b.idx.2.valid),ncol(p.matrix)))+
                                                       rho*ww.s.matrix[match.bs.valid[w.b.idx.2.valid],]
  w.b.matrix.valid[w.b.idx.3.valid,]=w.bs.new.valid[match(t.mat.valid[w.b.idx.3.valid,1],unique.bs.new.valid),,1]
  w.s.matrix.valid=array(NA,dim=c(length(z.valid),ncol(p.matrix)))
  w.s.idx.1.valid=(!is.na(match.s.valid))
  w.s.idx.2.valid=(is.na(match.s.valid))&(!is.na(match.sb.valid))
  w.s.idx.3.valid=(is.na(match.s.valid))&(is.na(match.sb.valid))
  w.s.matrix.valid[w.s.idx.1.valid,]=ww.s.matrix[match.s.valid[w.s.idx.1.valid],]
  w.s.matrix.valid[w.s.idx.2.valid,]=sqrt(1-rho^2)*array(rnorm(sum(w.s.idx.2.valid)*ncol(p.matrix)),
                                                       dim=c(sum(w.s.idx.2.valid),ncol(p.matrix)))+
                                                       rho*ww.b.matrix[match.sb.valid[w.s.idx.2.valid],]
  w.s.matrix.valid[w.s.idx.3.valid,]=w.bs.new.valid[match(t.mat.valid[w.s.idx.3.valid,2],unique.bs.new.valid),,2]
  w.p.matrix.valid=array(NA,dim=c(length(z.valid),ncol(p.matrix)))
  w.p.matrix.valid[!is.na(match.p.valid),]=ww.p.matrix[match.p.valid[!is.na(match.p.valid)],]
  w.p.matrix.valid[is.na(w.p.matrix.valid[,1]),]=rnorm(sum(is.na(w.p.matrix.valid)),mean=0,sd=1)
  
  p.matrix.valid=exp(reg.link.matrix.OS(x.valid, coef.fix=alpha, coef.random=beta, 
                                     w.b.matrix.valid, w.s.matrix.valid, w.p.matrix.valid, logit=T)) 
  mu.matrix.valid=exp(reg.link.matrix.OS(x.valid, coef.fix=gamma, coef.random=nu, 
                                      w.b.matrix.valid, w.s.matrix.valid, w.p.matrix.valid,logit=F)) 
  gamma.tn.bar.valid=pgamma(tn.valid,shape=1/phi,scale=mu.matrix.valid*phi,lower.tail=FALSE,log.p=FALSE)
  gamma.tn.valid=pgamma(tn.valid,shape=1/phi,scale=mu.matrix.valid*phi,lower.tail=TRUE,log.p=FALSE)
  
  z.e0.valid=p.matrix.valid*gamma.tn.bar.valid/(p.matrix.valid*gamma.tn.bar.valid+(1-p.matrix.valid))
  ur.e.matrix.valid=array(0,dim=dim(p.matrix.valid))
  ur.e.matrix.valid[z.valid==0,]=z.e0.valid[z.valid==0,]
  ur.e.ind.valid=apply(ur.e.matrix.valid,1,mean)
  ur.e.agg.valid=sum(ur.e.ind.valid)
  
  p.obs.valid=pmax(apply(p.matrix.valid*gamma.tn.valid,1,mean),1e-100)
  p.ur.valid=pmax(apply(ur.e.matrix.valid,1,mean),1e-100)
  p.com.valid=pmax(apply(p.matrix.valid,1,mean),1e-100)
  
  ## Testing ##
  unique.bs.test=unique(t.mat.test[,c(1,2)])
  unique.bs.new.test=unique.bs.test[(!(unique.bs.test%in%t.mat[,1]))&(!(unique.bs.test%in%t.mat[,2]))]
  w.bs.new.test=array(rmvnorm(length(unique.bs.new.test)*ncol(p.matrix),mean=c(0,0),sigma=array(c(1,rho,rho,1),dim=c(2,2))),
                      dim=c(length(unique.bs.new.test),ncol(p.matrix),2))
  match.bs.test=match(t.mat.test[,1],unique.s)
  match.sb.test=match(t.mat.test[,2],unique.b)
  w.b.matrix.test=array(NA,dim=c(length(z.test),ncol(p.matrix)))
  w.b.idx.1.test=(!is.na(match.b.test))
  w.b.idx.2.test=(is.na(match.b.test))&(!is.na(match.bs.test))
  w.b.idx.3.test=(is.na(match.b.test))&(is.na(match.bs.test))
  w.b.matrix.test[w.b.idx.1.test,]=ww.b.matrix[match.b.test[w.b.idx.1.test],]
  w.b.matrix.test[w.b.idx.2.test,]=sqrt(1-rho^2)*array(rnorm(sum(w.b.idx.2.test)*ncol(p.matrix)),
                                                       dim=c(sum(w.b.idx.2.test),ncol(p.matrix)))+
                                   rho*ww.s.matrix[match.bs.test[w.b.idx.2.test],]
  w.b.matrix.test[w.b.idx.3.test,]=w.bs.new.test[match(t.mat.test[w.b.idx.3.test,1],unique.bs.new.test),,1]
  w.s.matrix.test=array(NA,dim=c(length(z.test),ncol(p.matrix)))
  w.s.idx.1.test=(!is.na(match.s.test))
  w.s.idx.2.test=(is.na(match.s.test))&(!is.na(match.sb.test))
  w.s.idx.3.test=(is.na(match.s.test))&(is.na(match.sb.test))
  w.s.matrix.test[w.s.idx.1.test,]=ww.s.matrix[match.s.test[w.s.idx.1.test],]
  w.s.matrix.test[w.s.idx.2.test,]=sqrt(1-rho^2)*array(rnorm(sum(w.s.idx.2.test)*ncol(p.matrix)),
                                                       dim=c(sum(w.s.idx.2.test),ncol(p.matrix)))+
                                   rho*ww.b.matrix[match.sb.test[w.s.idx.2.test],]
  w.s.matrix.test[w.s.idx.3.test,]=w.bs.new.test[match(t.mat.test[w.s.idx.3.test,2],unique.bs.new.test),,2]
  w.p.matrix.test=array(NA,dim=c(length(z.test),ncol(p.matrix)))
  w.p.matrix.test[!is.na(match.p.test),]=ww.p.matrix[match.p.test[!is.na(match.p.test)],]
  w.p.matrix.test[is.na(w.p.matrix.test[,1]),]=rnorm(sum(is.na(w.p.matrix.test)),mean=0,sd=1)
  
  p.matrix.test=exp(reg.link.matrix.OS(x.test, coef.fix=alpha, coef.random=beta, 
                                           w.b.matrix.test, w.s.matrix.test, w.p.matrix.test, logit=T)) 
  
  p.com.test=pmax(apply(p.matrix.test,1,mean),1e-100)
  
  analysis.list=list(ur.e.ind=ur.e.ind, ur.e.agg=ur.e.agg,
                     p.obs=p.obs, p.ur=p.ur, p.com=p.com,
                     ur.e.ind.valid=ur.e.ind.valid, ur.e.agg.valid=ur.e.agg.valid,
                     p.obs.valid=p.obs.valid, p.ur.valid=p.ur.valid, p.com.valid=p.com.valid,
                     p.com.test=p.com.test)
  return(analysis.list)
}
