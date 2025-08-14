## Load package and data
library(doParallel)
library(matrixStats)
library(VGAM)
library(data.table)
library(mvtnorm)
library(fMultivar)
source("main_function.R")
load("x_syn.Rda")
load("tn_syn.Rda")
load("rd_syn.Rda")
load("z_syn.Rda")
load("tmat_syn.Rda")
load("fit_GLMM.Rda") # with DC network variables; run load("fit_GLMMb.Rda") without DC variables.
## Run x=x[,1:28] if evaluating without DC variables.
## Initial parameters
glm.z=glm.fit(x,z,family=binomial(link="logit"))
alpha=glm.z$coefficients
glm.rd=glm.fit(x[!is.na(rd),],rd[!is.na(rd)],family=Gamma(link="log"))
gamma=glm.rd$coefficients
phi=sigma(glm.rd)^2
beta=c(3,3,3)
nu=c(0.1,0.1,0.1)
rho=0
sig.pen=0.5*c(Inf+numeric(ncol(x)-1), 1, 1, 1) #total penalty=0.5~lambda=10^-5 in paper

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

length.ww.b=length(unique.b)
length.ww.s=length(unique.s)
length.ww.p=length(unique.p)
ww.b=numeric(length.ww.b)
ww.s=numeric(length.ww.s)
ww.p=numeric(length.ww.p)

### Retrieve ecm.list
iter.max=200
iter.eval=10
iter.min=(iter.max-iter.eval+1)
alpha.length=gamma.length=length(alpha)
alpha.matrix=gamma.matrix=array(0,dim=c(iter.max,alpha.length))
beta.matrix=nu.matrix=array(0,dim=c(iter.max,3))
phi.vec=rho.vec=as.numeric(iter.max)
for (i in 1:iter.max){
  alpha.matrix[i,]=ecm.list[[i]]$alpha
  beta.matrix[i,]=ecm.list[[i]]$beta
  gamma.matrix[i,]=ecm.list[[i]]$gamma
  nu.matrix[i,]=ecm.list[[i]]$nu
  phi.vec[i]=ecm.list[[i]]$phi
  rho.vec[i]=ecm.list[[i]]$rho
}

## retrieve fitted parameters
alpha.fit=apply(alpha.matrix[iter.min:iter.max,],2,mean)
beta.fit=apply(beta.matrix[iter.min:iter.max,],2,mean)
gamma.fit=apply(gamma.matrix[iter.min:iter.max,],2,mean)
nu.fit=apply(nu.matrix[iter.min:iter.max,],2,mean)
phi.fit=mean(phi.vec[iter.min:iter.max])
rho.fit=mean(rho.vec[iter.min:iter.max])

## Run ECM with parametric bootstrapping to find standard errors of parameters
ncore=60
n.ucty=300
cl = makePSOCKcluster(ncore)
registerDoParallel(cl)
ecm.list.ucty=foreach(j=1:n.ucty, .packages=c("MASS","matrixStats","VGAM","data.table",
                                             "mvtnorm","fMultivar")) %dopar% {
  tryCatch({
  set.seed(j)
  ww.multi=rmvnorm(ww.length.multi,mean=rep(0,2),sigma=array(c(1,rho.fit,rho.fit,1),dim=c(2,2)))
  ww.uni.b=rnorm(length(unique.uni.b), mean=0, sd=1)
  ww.uni.s=rnorm(length(unique.uni.s), mean=0, sd=1)
  ww.b.sim = c(ww.multi[,1],ww.uni.b)
  ww.s.sim = c(ww.multi[,2],ww.uni.s)
  ww.p.sim = rnorm(length.ww.p, mean=0, sd=1)                                             
  p.vec=exp(reg.link(x=x, coef.fix=alpha.fit, coef.random=beta.fit, 
                     match.b, match.s, match.p, 
                     ww.b=ww.b.sim, ww.s=ww.s.sim, ww.p=ww.p.sim, logit=T)) 
  z.sim=z0.sim=1*(runif(length(z))<=p.vec)
  mu.vec=exp(reg.link(x=x, coef.fix=gamma.fit, coef.random=nu.fit, 
                      match.b, match.s, match.p, 
                      ww.b=ww.b.sim, ww.s=ww.s.sim, ww.p=ww.p.sim, logit=F))
  rd.sim=rd0.sim=rgamma(runif(length(z)),shape=1/phi.fit,scale=mu.vec*phi.fit)
  rd.sim[rd0.sim>tn]=Inf
  rd.sim[z0.sim==0]=Inf
  z.sim[rd0.sim>tn]=0
  
  ecm(z.sim, rd.sim, tn, x, 
      alpha=alpha.fit, beta=beta.fit, gamma=gamma.fit, nu=nu.fit, phi=phi.fit, rho=rho.fit,
      t.mat, ww.b=ww.b.sim, ww.s=ww.s.sim, ww.p=ww.p.sim,
      MH.sd.ww=1, iter.max=40, num.sim=20, num.burn=10, num.sim.e=2, 
      sig.pen=NULL, seed.num=j, print=TRUE)                                             
  }, error=function(e){"Error!"})
  }
stopCluster(cl)

######## Retrieve uncertainty ##########
n.ucty=300
iter.max=40
iter.eval=10
iter.min=(iter.max-iter.eval+1)
alpha.length=gamma.length=length(alpha)
alpha.matrix=gamma.matrix=array(0,dim=c(iter.max,alpha.length))
beta.matrix=nu.matrix=array(0,dim=c(iter.max,3))
phi.vec=rho.vec=as.numeric(iter.max)
alpha.ucty.matrix=gamma.ucty.matrix=array(0,dim=c(n.ucty,alpha.length))
beta.ucty.matrix=nu.ucty.matrix=array(0,dim=c(n.ucty,3))
phi.ucty.vec=rho.ucty.vec=as.numeric(n.ucty)
for (j in 1:n.ucty){
  for (i in 1:iter.max){
    alpha.matrix[i,]=ecm.list.ucty[[j]][[i]]$alpha
    beta.matrix[i,]=ecm.list.ucty[[j]][[i]]$beta
    gamma.matrix[i,]=ecm.list.ucty[[j]][[i]]$gamma
    nu.matrix[i,]=ecm.list.ucty[[j]][[i]]$nu
    phi.vec[i]=ecm.list.ucty[[j]][[i]]$phi
    rho.vec[i]=ecm.list.ucty[[j]][[i]]$rho
  }
  alpha.ucty.matrix[j,]=apply(alpha.matrix[iter.min:iter.max,],2,mean)
  beta.ucty.matrix[j,]=apply(beta.matrix[iter.min:iter.max,],2,mean)
  gamma.ucty.matrix[j,]=apply(gamma.matrix[iter.min:iter.max,],2,mean)
  nu.ucty.matrix[j,]=apply(nu.matrix[iter.min:iter.max,],2,mean)
  phi.ucty.vec[j]=mean(phi.vec[iter.min:iter.max])
  rho.ucty.vec[j]=mean(rho.vec[iter.min:iter.max])
}
# Output the bootstrap standard error for each parameter
apply(alpha.ucty.matrix,2,sd)
apply(beta.ucty.matrix,2,sd)
apply(gamma.ucty.matrix,2,sd)
apply(nu.ucty.matrix,2,sd)
sd(phi.ucty.vec)
sd(rho.ucty.vec)
# Output the bootstrap z estimate for each parameter
alpha.fit/apply(alpha.ucty.matrix,2,sd)
beta.fit/apply(beta.ucty.matrix,2,sd)
gamma.fit/apply(gamma.ucty.matrix,2,sd)
nu.fit/apply(nu.ucty.matrix,2,sd)
phi.fit/sd(phi.ucty.vec)
rho.fit/sd(rho.ucty.vec)

