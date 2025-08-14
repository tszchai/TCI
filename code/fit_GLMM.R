## Load package and data
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
#### Initial parameters
## Note: Here we assume fitting with DC network variables ##
## If we fit without DC network variables, run the following additionally:
## x=x[,1:28]
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

## Run ECM
ecm.list=ecm(z, rd, tn, x, alpha, beta, gamma, nu, phi, rho,
             t.mat, ww.b, ww.s, ww.p,
             MH.sd.ww=1, iter.max=200, num.sim=20, num.burn=10, num.sim.e=2, 
             sig.pen, seed.num=999, print=T)

## Save fitting result (assume fitting with DC network variables)
save(ecm.list,file="fit_GLMM.Rda")
## If we assume fitting without DC network variables, save the following instead:
## save(ecm.list,file="fit_GLMMb.Rda")
