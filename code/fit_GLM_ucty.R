## Load package and data
library(doParallel)
library(matrixStats)
library(VGAM)
library(data.table)
library(mvtnorm)
library(fMultivar)
source("main_function_GLM.R")
load("x_syn.Rda")
load("tn_syn.Rda")
load("rd_syn.Rda")
load("z_syn.Rda")
load("z0_syn.Rda")
load("tmat_syn.Rda")
load("fit_GLM.Rda") 

### Retrieve ecm.list.GLM
iter.max=100
iter.eval=10
iter.min=(iter.max-iter.eval+1)
alpha.length=gamma.length=ncol(x)
alpha.matrix=gamma.matrix=array(0,dim=c(iter.max,alpha.length))
phi.vec=as.numeric(iter.max)
for (i in 1:iter.max){
  alpha.matrix[i,]=ecm.list.GLM[[i]]$alpha
  gamma.matrix[i,]=ecm.list.GLM[[i]]$gamma
  phi.vec[i]=ecm.list.GLM[[i]]$phi
}
alpha=alpha.fit=apply(alpha.matrix[iter.min:iter.max,],2,mean)
gamma=gamma.fit=apply(gamma.matrix[iter.min:iter.max,],2,mean)
phi=phi.fit=mean(phi.vec[iter.min:iter.max])

## Run ECM
ncore=60
n.ucty=300
cl = makePSOCKcluster(ncore)
registerDoParallel(cl)
ecm.list.ucty=foreach(j=1:n.ucty, .packages=c("MASS","matrixStats","VGAM","data.table",
                                             "mvtnorm","fMultivar")) %dopar% {
     set.seed(j)                                             
     p.vec=exp(reg.link(x=x, coef.fix=alpha, logit=T)) 
     z.sim=z0.sim=1*(runif(length(z))<=p.vec)
     mu.vec=exp(reg.link(x=x, coef.fix=gamma, logit=F))
     rd.sim=rd0.sim=rgamma(runif(length(z)),shape=1/phi,scale=mu.vec*phi)
     rd.sim[rd0.sim>tn]=Inf
     rd.sim[z0.sim==0]=Inf
     z.sim[rd0.sim>tn]=0
     
     ecm(z.sim, rd.sim, tn, x, 
         alpha=alpha, gamma=gamma.fit, phi=phi.fit,
         iter.max=30, sig.pen=NULL, seed.num=(999*j), print=TRUE)                                             
   }
stopCluster(cl)

save(ecm.list.ucty,file="fit_GLM_ucty.Rda")

## evaluate uncertainty (bootstrap standard errors) [Need to load bootstrap results first!]
load("fit_GLM_ucty.Rda"); load("x_syn.Rda") 
#load("fit_GLMb_ucty.Rda"); load("x_syn.Rda"); x=x[,1:28]
n.ucty=300
iter.max=30
iter.eval=10
iter.min=(iter.max-iter.eval+1)
alpha.length=gamma.length=ncol(x)
alpha.matrix=gamma.matrix=array(0,dim=c(iter.max,alpha.length))
phi.vec=as.numeric(iter.max)
alpha.ucty.matrix=gamma.ucty.matrix=array(0,dim=c(n.ucty,alpha.length))
phi.ucty.vec=as.numeric(n.ucty)
for (j in 1:n.ucty){
  for (i in 1:iter.max){
    alpha.matrix[i,]=ecm.list.ucty[[j]][[i]]$alpha
    gamma.matrix[i,]=ecm.list.ucty[[j]][[i]]$gamma
    phi.vec[i]=ecm.list.ucty[[j]][[i]]$phi
  }
  alpha.ucty.matrix[j,]=apply(alpha.matrix[iter.min:iter.max,],2,mean)
  gamma.ucty.matrix[j,]=apply(gamma.matrix[iter.min:iter.max,],2,mean)
  phi.ucty.vec[j]=mean(phi.vec[iter.min:iter.max])
}
# Output the bootstrap standard error for each parameter
apply(alpha.ucty.matrix,2,sd)
apply(gamma.ucty.matrix,2,sd)
sd(phi.ucty.vec)
