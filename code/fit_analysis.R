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
load("rd0_syn.Rda")
load("z_syn.Rda")
load("z0_syn.Rda")
load("tmat_syn.Rda")

load("x_valid_syn.Rda")
load("tn_valid_syn.Rda")
load("rd_valid_syn.Rda")
load("z_valid_syn.Rda")
load("z0_valid_syn.Rda")
load("tmat_valid_syn.Rda")

load("x_test_syn.Rda")
load("tn_test_syn.Rda")
load("rd_test_syn.Rda")
load("z_test_syn.Rda")
load("z0_test_syn.Rda")
load("tmat_test_syn.Rda")

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

### Load fitting results ###
load("fit_GLMM.Rda") #ecm.list; 
#similarly do load("fit_GLMMb.Rda") for without DC
#but remember have x=x[,1:28]; x.valid=x.valid[,1:28]; x.test=x.test[,1:28] to exclude DC variables!

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

########## Performance analysis #########
# loading the posterior simulated samples of buyer/seller/policy random effects from each level
load("ww_b_post.Rda") ## or load("ww_b_postb.Rda") without DC
load("ww_s_post.Rda") ## or load("ww_s_postb.Rda") without DC
load("ww_p_post.Rda") ## or load("ww_p_postb.Rda") without DC

# This gives an analysis summary
analysis.list=analysis.fit(z, z0, rd, tn, x, t.mat,
                           alpha.fit, beta.fit, gamma.fit, nu.fit, phi.fit, rho.fit,
                           match.b, match.s, match.p, 
                           ww.b.matrix, ww.s.matrix, ww.p.matrix,
                           z.valid, z0.valid, rd.valid, tn.valid, x.valid, t.mat.valid,
                           z.test, z0.test, rd.test, tn.test, x.test, t.mat.test, 
                           seed.eval=1001)
# ADEV (Absolute deviance statistic)
# p.obs is the prob. of claim being observed before evaluation date (training set)
# p.ur is the prob of claim occur but not observed before evaluation date; reserve
# p.com=p.obs+p.ur is the prob of having claim anytime
# .valid represents validation set
# .test represents test set
sum(abs(analysis.list$p.obs-z))
sum(abs(analysis.list$p.ur-(z0-z)))
sum(abs(analysis.list$p.com-z0))
sum(abs(analysis.list$p.obs.valid-z.valid))
sum(abs(analysis.list$p.ur.valid-(z0.valid-z.valid)))
sum(abs(analysis.list$p.com.valid-z0.valid))
sum(abs(analysis.list$p.com.test-z0.test))
# Reserve (estimated vs actual)
c(sum(analysis.list$p.ur), sum(z0-z))
c(sum(analysis.list$p.ur.valid), sum(z0.valid-z.valid))

### Reporting delay distribution
set.seed(1001)
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
# Density plot
plot(density(rd.sim[is.finite(rd.sim)]),ylim=c(0,0.0028)) #model rd after truncation
lines(density(rd[!is.na(rd)]),lty=2) #actual observed rd (after truncation)
lines(density(rd0.sim[z0.sim==1]),col="red") #model rd before truncation
lines(density(rd0[!is.na(rd0)]),col="red",lty=2) #actual true rd
#lines(density(rd.valid[!is.na(rd.valid)]),col="blue")
# QQ-plot
# reporting delay truncated before evaluation date
qqplot(rd.sim[is.finite(rd.sim)],rd[!is.na(rd)])
abline(a=0,b=1,col="red")
# true reporting delay
qqplot(rd0.sim[z0.sim==1],rd0[!is.na(rd0)])
abline(a=0,b=1,col="red")

### Posterior classification
# Posterior mean of random effects for each level
ww.b.mean.post=apply(ww.b.matrix,1,mean)
ww.s.mean.post=apply(ww.s.matrix,1,mean)
ww.p.mean.post=apply(ww.p.matrix,1,mean)
ww.b.sd.post=apply(ww.b.matrix,1,sd)
ww.s.sd.post=apply(ww.s.matrix,1,sd)
ww.p.sd.post=apply(ww.p.matrix,1,sd)
# w.b.mean.post=ww.b.mean.post[t.mat[,1]]
#plot histogram
hist(ww.b.mean.post,breaks=100) #bimodal for buyer posterior effects
hist(ww.s.mean.post,breaks=100)
hist(ww.p.mean.post,breaks=100)

### Correlation to a trade connection
i=1
idx.b=t.mat[i,1]
idx.s=t.mat[i,2]
adj.idx=unique(c(which(t.mat[,1]==idx.b),which(t.mat[,2]==idx.b),
                 which(t.mat[,1]==idx.s),which(t.mat[,2]==idx.s)))
x.adj=x[adj.idx,]
tn.adj=tn[adj.idx]
t.mat.adj=t.mat[adj.idx,]

unique.b.adj=unique(t.mat.adj[,1])
unique.s.adj=unique(t.mat.adj[,2])
unique.multi.adj=unique.b.adj[unique.b.adj%in%unique.s.adj]
unique.uni.b.adj=unique.b.adj[!(unique.b.adj%in%unique.s.adj)]
unique.uni.s.adj=unique.s.adj[!(unique.s.adj%in%unique.b.adj)]
unique.b.adj=c(unique.multi.adj,unique.uni.b.adj)
unique.s.adj=c(unique.multi.adj,unique.uni.s.adj)
unique.p.adj=unique(t.mat.adj[,3])
match.b.adj=match(t.mat.adj[,1],unique.b.adj)
match.s.adj=match(t.mat.adj[,2],unique.s.adj)
match.p.adj=match(t.mat.adj[,3],unique.p.adj)
ww.length.multi.adj=length(unique.multi.adj)

n.sim.adj=100000
z0.sim.matrix=array(0,dim=c(n.sim.adj,nrow(x.adj)))
effect.fix.z.adj=(x.adj%*%alpha.fit)[,1]
effect.fix.rd.adj=(x.adj%*%gamma.fit)[,1]

for(i in 1:n.sim.adj){
  ww.multi.adj=rmvnorm(ww.length.multi.adj,mean=rep(0,2),sigma=array(c(1,rho.fit,rho.fit,1),dim=c(2,2)))
  ww.uni.b.adj=rnorm(length(unique.uni.b.adj), mean=0, sd=1)
  ww.uni.s.adj=rnorm(length(unique.uni.s.adj), mean=0, sd=1)
  ww.b.sim = c(ww.multi.adj[,1],ww.uni.b.adj)
  ww.s.sim = c(ww.multi.adj[,2],ww.uni.s.adj)
  ww.p.sim = rnorm(length(unique.p.adj), mean=0, sd=1)                                             
  p.vec=exp(reg.link(x=x.adj, coef.fix=alpha.fit, coef.random=beta.fit, 
                     match.b.adj, match.s.adj, match.p.adj, 
                     ww.b=ww.b.sim, ww.s=ww.s.sim, ww.p=ww.p.sim, 
                     logit=T, effect.fix = effect.fix.z.adj)) 
  z0.sim.matrix[i,]=1*(runif(nrow(x.adj))<=p.vec)
}
cor(z0.sim.matrix)[which(adj.idx==i),]
