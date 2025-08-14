## load data ##
load("data_synthetic.Rda") ## for synthetic data
library(data.table)
library(lubridate)
library(Matrix)

##### DC calculation and append to the data #####
data$dc.o.s=data$dc.i.s=data$dc.oo.s=data$dc.ii.s=data$dc.oi.s=data$dc.io.s=NA
data$dc.o.b=data$dc.i.b=data$dc.oo.b=data$dc.ii.b=data$dc.oi.b=data$dc.io.b=NA
t.total=as.numeric(max(data$start_date)-min(data$start_date))+1
t.min=min(data$start_date)
Dt.dim=max(data$id.b,data$id.s)
for (s in 1:t.total)
{
  date.idx.eval=(data$start_date==(t.min+s-1))
  if(sum(date.idx.eval)>0){
    data.sel=data[(data$start_date<=(t.min+s-1))&(data$start_date>(t.min+s-365-1)),]
    data.eval=data[date.idx.eval,]
    Dt.mat=sparseMatrix(i=data.sel$id.s, j=data.sel$id.b, x = 1L, dims = c(Dt.dim,Dt.dim))
    id.s.eval=unique(data.eval$id.s)
    id.b.eval=unique(data.eval$id.b)
    dc.o.s.eval=rowSums(Dt.mat[id.s.eval,,drop=FALSE])
    dc.i.s.eval=colSums(Dt.mat[,id.s.eval,drop=FALSE])
    dc.o.b.eval=rowSums(Dt.mat[id.b.eval,,drop=FALSE])
    dc.i.b.eval=colSums(Dt.mat[,id.b.eval,drop=FALSE])
    Dt.oo.mat=Dt.mat%*%Dt.mat
    dc.oo.s.eval=rowSums(Dt.oo.mat[id.s.eval,,drop=FALSE])
    dc.ii.s.eval=colSums(Dt.oo.mat[,id.s.eval,drop=FALSE])
    dc.oo.b.eval=rowSums(Dt.oo.mat[id.b.eval,,drop=FALSE])
    dc.ii.b.eval=colSums(Dt.oo.mat[,id.b.eval,drop=FALSE])
    Dt.oi.mat=Dt.mat%*%t(Dt.mat)-Diagonal(x = as.numeric(diag(Dt.mat%*%t(Dt.mat))))
    dc.oi.s.eval=rowSums(Dt.oi.mat[id.s.eval,,drop=FALSE])
    dc.oi.b.eval=rowSums(Dt.oi.mat[id.b.eval,,drop=FALSE])
    Dt.io.mat=t(Dt.mat)%*%Dt.mat-Diagonal(x = as.numeric(diag(t(Dt.mat)%*%Dt.mat)))
    dc.io.s.eval=rowSums(Dt.io.mat[id.s.eval,,drop=FALSE])
    dc.io.b.eval=rowSums(Dt.io.mat[id.b.eval,,drop=FALSE])
    match.idx.s=match(data.eval$id.s,id.s.eval)
    match.idx.b=match(data.eval$id.b,id.b.eval)
    data$dc.o.s[date.idx.eval]=log(1+dc.o.s.eval[match.idx.s]-1)
    data$dc.i.s[date.idx.eval]=log(1+dc.i.s.eval[match.idx.s])
    data$dc.oo.s[date.idx.eval]=log(1+dc.oo.s.eval[match.idx.s])
    data$dc.ii.s[date.idx.eval]=log(1+dc.ii.s.eval[match.idx.s])
    data$dc.oi.s[date.idx.eval]=log(1+dc.oi.s.eval[match.idx.s])
    data$dc.io.s[date.idx.eval]=log(1+dc.io.s.eval[match.idx.s])
    data$dc.o.b[date.idx.eval]=log(1+dc.o.b.eval[match.idx.b])
    data$dc.i.b[date.idx.eval]=log(1+dc.i.b.eval[match.idx.b]-1)
    data$dc.oo.b[date.idx.eval]=log(1+dc.oo.b.eval[match.idx.b])
    data$dc.ii.b[date.idx.eval]=log(1+dc.ii.b.eval[match.idx.b])
    data$dc.oi.b[date.idx.eval]=log(1+dc.oi.b.eval[match.idx.b])
    data$dc.io.b[date.idx.eval]=log(1+dc.io.b.eval[match.idx.b])
  }
}

## save the full data with DC appended
data_full=data
save(data_full,file="data_synthetic_full.Rda")

## retrieve training, validation and testing datasets
date.eval=as_date("2019-12-31")
prob.train=0.8 ## percentage of training data points
set.seed(1)

data_full$rd=as.numeric(data_full$claim_date-data_full$start_date)
data_full$tn=as.numeric(date.eval-data_full$start_date+1)
data_full$tn[data_full$tn<0]=NA
data_full$z0=data_full$z
data_full$z[data_full$claim_date>date.eval]=0
data_full$rd0=data_full$rd
data_full$rd[data_full$claim_date>date.eval]=NA

data_eval=data_full[!is.na(data_full$tn),]
data_test=data_full[is.na(data_full$tn),]
train.sample=sample(nrow(data_eval))
data_train=data_eval[train.sample<=ceiling(prob.train*nrow(data_eval)),]
data_valid=data_eval[train.sample>ceiling(prob.train*nrow(data_eval)),]

#### training data feed to model ####
n=nrow(data_train)
rd=data_train$rd
rd0=data_train$rd0
tn=data_train$tn
z=data_train$z
z0=data_train$z0
t.mat=array(0,dim=c(n,3))
t.mat[,1]=data_train$id.b
t.mat[,2]=data_train$id.s
t.mat[,3]=data_train$id.p
x=as.matrix(data_train[,-c("start_date", "claim_date", "id.b", "id.s", "id.p", 
                           "z", "rd", "tn", "z0", "rd0")])
x=cbind(1,x)
colnames(x)=NULL
## Save data
save(x,file="x_syn.Rda")
save(tn,file="tn_syn.Rda")
save(rd,file="rd_syn.Rda")
save(rd0,file="rd0_syn.Rda")
save(z,file="z_syn.Rda")
save(z0,file="z0_syn.Rda")
save(t.mat,file="tmat_syn.Rda")

#### validation data####
n=nrow(data_valid)
rd.valid=data_valid$rd
tn.valid=data_valid$tn
z.valid=data_valid$z
z0.valid=data_valid$z0
t.mat.valid=array(0,dim=c(n,3))
t.mat.valid[,1]=data_valid$id.b
t.mat.valid[,2]=data_valid$id.s
t.mat.valid[,3]=data_valid$id.p
x.valid=as.matrix(data_valid[,-c("start_date", "claim_date", "id.b", "id.s", "id.p", 
                                 "z", "rd", "tn", "z0", "rd0")])
x.valid=cbind(1,x.valid)
colnames(x.valid)=NULL
## Save data
save(x.valid,file="x_valid_syn.Rda")
save(tn.valid,file="tn_valid_syn.Rda")
save(rd.valid,file="rd_valid_syn.Rda")
save(z.valid,file="z_valid_syn.Rda")
save(z0.valid,file="z0_valid_syn.Rda")
save(t.mat.valid,file="tmat_valid_syn.Rda")

#### testing data####
n=nrow(data_test)
rd.test=data_test$rd
tn.test=data_test$tn
z.test=data_test$z
z0.test=data_test$z0
t.mat.test=array(0,dim=c(n,3))
t.mat.test[,1]=data_test$id.b
t.mat.test[,2]=data_test$id.s
t.mat.test[,3]=data_test$id.p
x.test=as.matrix(data_test[,-c("start_date", "claim_date", "id.b", "id.s", "id.p", 
                               "z", "rd", "tn", "z0", "rd0")])
x.test=cbind(1,x.test)
colnames(x.test)=NULL
## Save data
save(x.test,file="x_test_syn.Rda")
save(tn.test,file="tn_test_syn.Rda")
save(rd.test,file="rd_test_syn.Rda")
save(z.test,file="z_test_syn.Rda")
save(z0.test,file="z0_test_syn.Rda")
save(t.mat.test,file="tmat_test_syn.Rda")