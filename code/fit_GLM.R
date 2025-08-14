## Load package and data
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

################# With network variables ############
## Initial parameters
glm.z=glm.fit(x,z,family=binomial(link="logit"))
alpha=glm.z$coefficients
glm.rd=glm.fit(x[!is.na(rd),],rd[!is.na(rd)],family=Gamma(link="log"))
gamma=glm.rd$coefficients
phi=sigma(glm.rd)^2

ecm.list.GLM = ecm(z, rd, tn, x, alpha, gamma, phi, iter.max=100, 
                   sig.pen=NULL, seed.num=999, print=TRUE)

save(ecm.list.GLM,file="fit_GLM.Rda")

################# Without network variables##############
## Initial parameters
x=x[,1:28]
glm.z=glm.fit(x,z,family=binomial(link="logit"))
alpha=glm.z$coefficients
glm.rd=glm.fit(x[!is.na(rd),],rd[!is.na(rd)],family=Gamma(link="log"))
gamma=glm.rd$coefficients
phi=sigma(glm.rd)^2

ecm.list.GLM = ecm(z, rd, tn, x, alpha, gamma, phi, iter.max=100, 
                   sig.pen=NULL, seed.num=999, print=TRUE)

save(ecm.list.GLM,file="fit_GLMb.Rda")



#################### Analysis of fit ######################
##### Run one of the following two cases: ###########
# Case 1: GLM with DC
load("fit_GLM.Rda") 
load("x_syn.Rda") 
load("x_valid_syn.Rda")
load("x_test_syn.Rda")
# Case 2: GLM without DC
load("fit_GLMb.Rda")
load("x_syn.Rda") 
load("x_valid_syn.Rda")
load("x_test_syn.Rda")
x=x[,1:28]; x.valid=x.valid[,1:28]; x.test=x.test[,1:28]
####################################################

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

## log-likelihood
loglik(z, rd, tn, x, alpha, gamma, phi)$ll
## analysis
analysis.list=analysis.fit(z, z0, rd, tn, x, 
                           alpha, gamma, phi,
                           z.valid, z0.valid, rd.valid, tn.valid, x.valid,
                           z.test, z0.test, rd.test, tn.test, x.test,
                           seed.eval=1001)
# performance comparison: compute ADEV, reserve, etc.
c(sum(abs(analysis.list$p.obs-z)), sum(abs(analysis.list$p.ur-(z0-z))), sum(abs(analysis.list$p.com-z0)))
c(sum(abs(analysis.list$p.obs.valid-z.valid)), sum(abs(analysis.list$p.ur.valid-(z0.valid-z.valid))), 
  sum(abs(analysis.list$p.com.valid-z0.valid)))
c(sum(abs(analysis.list$p.com.test-z0.test)))
c(sum(analysis.list$p.obs),sum(z))
c(sum(analysis.list$p.ur),sum(z0-z))
c(sum(analysis.list$p.com),sum(z0))
c(sum(analysis.list$p.obs.valid),sum(z.valid))
c(sum(analysis.list$p.ur.valid),sum(z0.valid-z.valid))
c(sum(analysis.list$p.com.valid),sum(z0.valid))
c(sum(analysis.list$p.com.test),sum(z0.test))



