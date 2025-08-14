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
load("fit_GLMM.Rda")#ecm.list; 
#similarly do load("fit_GLMMb.Rda") for without DC
#but remember have x=x[,1:28]; x.valid=x.valid[,1:28]; x.test=x.test[,1:28] to exclude DC variables!

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

###################################
iter.max=200
iter.eval=10
iter.min=(iter.max-iter.eval+1)

### Retrieve ecm.list
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

###### In-sample observed data log-likelihood #####
num.sim.post=1000
num.sim.IS=5000
num.skip=2
num.burn=25
ncore=50
ww.MH.burn=ww.posterior.sim(z, rd, tn, x, alpha=alpha.fit, beta=beta.fit, gamma=gamma.fit, 
                            nu=nu.fit, phi=phi.fit, rho=rho.fit,
                            t.mat, match.b, match.s, match.p, 
                            ww.b=ww.b, ww.s=ww.s, ww.p=ww.p, 
                            unique.b, unique.s, unique.p, ww.length.multi,
                            MH.sd.ww=1, num.sim=num.burn, seed.eval=1001)


cl = makePSOCKcluster(ncore)
registerDoParallel(cl)
ww.MH.list <-
  foreach(j=1:ncore, .packages=c("MASS","matrixStats","VGAM","data.table",
                                 "mvtnorm","fMultivar")) %dopar% {
   ww.MH.list.trial=ww.posterior.sim(z, rd, tn, x, alpha=alpha.fit, 
                          beta=beta.fit, gamma=gamma.fit, 
                          nu=nu.fit, phi=phi.fit, rho=rho.fit,
                          t.mat, match.b, match.s, match.p, 
                          ww.b=ww.MH.burn$ww.MH.matrix.b[,num.burn], 
                          ww.s=ww.MH.burn$ww.MH.matrix.s[,num.burn], 
                          ww.p=ww.MH.burn$ww.MH.matrix.p[,num.burn], 
                          unique.b, unique.s, unique.p, ww.length.multi,
                          MH.sd.ww=1, num.sim=(num.skip*num.sim.post/ncore), seed.eval=(1001*j))
   ww.MH.matrix.b=ww.MH.list.trial$ww.MH.matrix.b[,num.skip*seq(num.sim.post/ncore)]
   ww.MH.matrix.s=ww.MH.list.trial$ww.MH.matrix.s[,num.skip*seq(num.sim.post/ncore)]
   ww.MH.matrix.p=ww.MH.list.trial$ww.MH.matrix.p[,num.skip*seq(num.sim.post/ncore)]
   list(ww.MH.matrix.b=ww.MH.matrix.b, ww.MH.matrix.s=ww.MH.matrix.s, ww.MH.matrix.p=ww.MH.matrix.p)
  }
stopCluster(cl)

ww.b.matrix=ww.MH.list[[1]]$ww.MH.matrix.b
ww.s.matrix=ww.MH.list[[1]]$ww.MH.matrix.s
ww.p.matrix=ww.MH.list[[1]]$ww.MH.matrix.p
for (j in 2:ncore){
  ww.b.matrix=cbind(ww.b.matrix,ww.MH.list[[j]]$ww.MH.matrix.b)
  ww.s.matrix=cbind(ww.s.matrix,ww.MH.list[[j]]$ww.MH.matrix.s)
  ww.p.matrix=cbind(ww.p.matrix,ww.MH.list[[j]]$ww.MH.matrix.p)
}

### save the posterior buyer, seller, and policy-level random effects
save(ww.b.matrix,file="ww_b_post.Rda") ## or save(ww.b.matrix,file="ww_b_postb.Rda") without DC
save(ww.s.matrix,file="ww_s_post.Rda") ## or save(ww.s.matrix,file="ww_s_postb.Rda") without DC
save(ww.p.matrix,file="ww_p_post.Rda") ## or save(ww.p.matrix,file="ww_p_postb.Rda") without DC

### in-sample observed likelihood evaluation
load("ww_b_post.Rda") ## or load("ww_b_postb.Rda") without DC
load("ww_s_post.Rda") ## or load("ww_s_postb.Rda") without DC
load("ww_p_post.Rda") ## or load("ww_p_postb.Rda") without DC
ww.b.mean.post=apply(ww.b.matrix,1,mean)
ww.s.mean.post=apply(ww.s.matrix,1,mean)
ww.p.mean.post=apply(ww.p.matrix,1,mean)
ww.b.sd.post=apply(ww.b.matrix,1,sd)
ww.s.sd.post=apply(ww.s.matrix,1,sd)
ww.p.sd.post=apply(ww.p.matrix,1,sd)
ww.MH.cor.array=array(0,dim=c(ww.length.multi,num.sim.post,2))
ww.MH.cor.array[,,1]=ww.b.matrix[1:ww.length.multi,]
ww.MH.cor.array[,,2]=ww.s.matrix[1:ww.length.multi,]
ww.cor.post=apply(ww.MH.cor.array,1,function(x)cor(x)[1,2])

cl = makePSOCKcluster(ncore)
registerDoParallel(cl)
loglik.list <-
  foreach(j=1:ncore, .packages=c("MASS","matrixStats","VGAM","data.table",
                                 "mvtnorm","fMultivar")) %dopar% {
                                   loglik.obs.sim.IS(z, rd, tn, x, alpha=alpha.fit, beta=beta.fit, gamma=gamma.fit,
                                                     nu=nu.fit, phi=phi.fit, rho=rho.fit,
                                                     t.mat, match.b, match.s, match.p,
                                                     ww.b=ww.b, ww.s=ww.s, ww.p=ww.p, ww.length.multi,
                                                     num.sim=(num.sim.IS/ncore), seed.eval=(1000*j),
                                                     ww.b.mean.post, ww.b.sd.post, ww.s.mean.post, ww.s.sd.post,
                                                     ww.p.mean.post, ww.p.sd.post, ww.cor.post)
                                 }
stopCluster(cl)

ll.obs.ind.naive=loglik.list[[1]]$ll.obs.ind.naive
for (j in 2:ncore) {ll.obs.ind.naive=c(ll.obs.ind.naive,loglik.list[[j]]$ll.obs.ind.naive)}

## output the observed log-likelihood of the proposed bivariate GLMM fitting results!
logSumExp(ll.obs.ind.naive)-log(length(ll.obs.ind.naive))
