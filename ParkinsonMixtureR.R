##################################################
### 
### A Flexible Replication-Based Classification Approach for Parkinson's Disease Detection by Using Voice Recordings
###  
### Lizbeth Naranjo (1), Ruth Fuentes-Garcia (1), Carlos J. Perez (2).
###
### (1) Universidad Nacional Autonoma de Mexico, Mexico
### (2) Universidad de Extremadura, Spain
### 
###  
### 33rd National Forum of Statistics (FNE) and 13th Latin-American Congress of Statistical Societies (CLATSE)
### Springer Proceedings in Mathematics & Statistics
### 
##################################################

##################################################
### R packages
###
### Instructions: 
### Load the R libraries
##################################################
library(rjags)
library(cpca)
require(plyr)
library(pROC)

##################################################

##################################################
### READ DATA
###
### Instructions: 
### Change the address where the data and codes are located. 
### setwd("HERE")
##################################################

setwd("~/FileDirectory/")
getwd()

Park <- read.table(file="ReplicatedAcousticFeatures-ParkinsonDatabase.csv", dec=".", sep=",", header=TRUE)

summary(Park)
colnames(Park)
dim(Park)
attach(Park)

N = 80

### Standardized data
for(j in 5:48){
  Park[,j] = (Park[,j]-mean(Park[,j]))/sd(Park[,j])
}

Park1 = Park[Park$Recording==1,]
Park2 = Park[Park$Recording==2,]
Park3 = Park[Park$Recording==3,]

##################################################

##################################################
### ANALYSIS
##################################################

##################################################
### Common Principal Components (CPC)

### CPC Jitter
columjit = c("Recording","Jitter_rel","Jitter_abs","Jitter_RAP","Jitter_PPQ")
covjit <- daply(Park[,columjit], "Recording", function(x) cov(x[, -1]))
covjit <- aperm(covjit, c(2, 3, 1)) # put the 1st dimension to the end

modjit <- cpc(covjit)
cumsum(modjit$D[,1])/sum(modjit$D[,1],na.rm=T)
cumsum(modjit$D[,2])/sum(modjit$D[,2],na.rm=T)
cumsum(modjit$D[,3])/sum(modjit$D[,3],na.rm=T)

modjit <- cpc(covjit, k = 1)
jit1 = as.matrix(Park1[,columjit[-1]])%*%as.matrix(modjit$CPC)
jit2 = as.matrix(Park2[,columjit[-1]])%*%as.matrix(modjit$CPC)
jit3 = as.matrix(Park3[,columjit[-1]])%*%as.matrix(modjit$CPC)

### CPC Shimmer
columshi = c("Recording","Shim_loc","Shim_dB","Shim_APQ3","Shim_APQ5","Shi_APQ11")
covshi <- daply(Park[,columshi], "Recording", function(x) cov(x[, -1]))
covshi <- aperm(covshi, c(2, 3, 1)) # put the 1st dimension to the end

modshi <- cpc(covshi)
cumsum(modshi$D[,1])/sum(modshi$D[,1],na.rm=T)
cumsum(modshi$D[,2])/sum(modshi$D[,2],na.rm=T)
cumsum(modshi$D[,3])/sum(modshi$D[,3],na.rm=T)

modshi <- cpc(covshi, k = 1)
shi1 = as.matrix(Park1[,columshi[-1]])%*%as.matrix(modshi$CPC)
shi2 = as.matrix(Park2[,columshi[-1]])%*%as.matrix(modshi$CPC)
shi3 = as.matrix(Park3[,columshi[-1]])%*%as.matrix(modshi$CPC)

### CPC HNR
columhnr = c("Recording","HNR05","HNR15","HNR25","HNR35","HNR38")
covhnr <- daply(Park[,columhnr], "Recording", function(x) cov(x[, -1]))
covhnr <- aperm(covhnr, c(2, 3, 1)) # put the 1st dimension to the end

modhnr <- cpc(covhnr)
cumsum(modhnr$D[,1])/sum(modhnr$D[,1],na.rm=T)
cumsum(modhnr$D[,2])/sum(modhnr$D[,2],na.rm=T)
cumsum(modhnr$D[,3])/sum(modhnr$D[,3],na.rm=T)

modhnr <- cpc(covhnr, k = 1)
hnr1 = as.matrix(Park1[,columhnr[-1]])%*%as.matrix(modhnr$CPC)
hnr2 = as.matrix(Park2[,columhnr[-1]])%*%as.matrix(modhnr$CPC)
hnr3 = as.matrix(Park3[,columhnr[-1]])%*%as.matrix(modhnr$CPC)

### CPC MFCC & Delta
colummfcdel = c("Recording","MFCC0", "MFCC1","MFCC2","MFCC3","MFCC4", "MFCC5","MFCC6","MFCC7","MFCC8", "MFCC9","MFCC10","MFCC11","MFCC12",
                "Delta0", "Delta1","Delta2","Delta3","Delta4", "Delta5","Delta6","Delta7","Delta8", "Delta9","Delta10","Delta11","Delta12")
covmfcdel <- daply(Park[,colummfcdel], "Recording", function(x) cov(x[, -1]))
covmfcdel <- aperm(covmfcdel, c(2, 3, 1)) # put the 1st dimension to the end

modmfcdel <- cpc(covmfcdel)
cumsum(modmfcdel$D[,1])/sum(modmfcdel$D[,1],na.rm=T)
cumsum(modmfcdel$D[,2])/sum(modmfcdel$D[,2],na.rm=T)
cumsum(modmfcdel$D[,3])/sum(modmfcdel$D[,3],na.rm=T)
plot(modmfcdel$D[,1]/sum(modmfcdel$D[,1],na.rm=T))
plot(modmfcdel$D[,2]/sum(modmfcdel$D[,2],na.rm=T))
plot(modmfcdel$D[,3]/sum(modmfcdel$D[,3],na.rm=T))

modmfcdel <- cpc(covmfcdel, k = 3)
mfcdel1r1 = as.matrix(Park1[,colummfcdel[-1]])%*%as.matrix(modmfcdel$CPC[,1])
mfcdel1r2 = as.matrix(Park2[,colummfcdel[-1]])%*%as.matrix(modmfcdel$CPC[,1])
mfcdel1r3 = as.matrix(Park3[,colummfcdel[-1]])%*%as.matrix(modmfcdel$CPC[,1])
mfcdel2r1 = as.matrix(Park1[,colummfcdel[-1]])%*%as.matrix(modmfcdel$CPC[,2])
mfcdel2r2 = as.matrix(Park2[,colummfcdel[-1]])%*%as.matrix(modmfcdel$CPC[,2])
mfcdel2r3 = as.matrix(Park3[,colummfcdel[-1]])%*%as.matrix(modmfcdel$CPC[,2])
mfcdel3r1 = as.matrix(Park1[,colummfcdel[-1]])%*%as.matrix(modmfcdel$CPC[,3])
mfcdel3r2 = as.matrix(Park2[,colummfcdel[-1]])%*%as.matrix(modmfcdel$CPC[,3])
mfcdel3r3 = as.matrix(Park3[,colummfcdel[-1]])%*%as.matrix(modmfcdel$CPC[,3])

##################################################

##################################################
### Data 

gender = Park1$Gender
sano = status = Park1$Status

Rep1 = data.frame( status=Park1$Status, sex=Park1$Gender, 
                   jit=jit1, shi=shi1, hnr=hnr1, 
                   mfcdel1=mfcdel1r1, mfcdel2=mfcdel2r1, mfcdel3=mfcdel3r1, 
                   rpde=Park1$RPDE, dfa=Park1$DFA, ppe=Park1$PPE, gne=Park1$GNE )

Rep2 = data.frame( status=Park1$Status, sex=Park2$Gender, 
                   jit=jit2, shi=shi2, hnr=hnr2, 
                   mfcdel1=mfcdel1r2, mfcdel2=mfcdel2r2, mfcdel3=mfcdel3r2, 
                   rpde=Park2$RPDE, dfa=Park2$DFA, ppe=Park2$PPE, gne=Park2$GNE )

Rep3 = data.frame( status=Park1$Status, sex=Park3$Gender,
                   jit=jit3, shi=shi3, hnr=hnr3, 
                   mfcdel1=mfcdel1r3, mfcdel2=mfcdel2r3, mfcdel3=mfcdel3r3, 
                   rpde=Park3$RPDE, dfa=Park3$DFA, ppe=Park3$PPE, gne=Park3$GNE )

DATOS <- list( Rep1=Rep1, Rep2=Rep2, Rep3=Rep3)

attributes(DATOS)
attributes(DATOS$Rep1)
attributes(DATOS$Rep2)
attributes(DATOS$Rep3)

##################################################

##################################################
### Mixture model 
###
### Binary response 
### Mixture of Normal Distributions
### Replications
### Options: Normal, Mixture of 2 Normals, Mixture of 3 Normals
##################################################

### Data

N = 80
K = 10

X = array(NA,dim=c(N,3,K))
X[,,1] = cbind(Rep1$jit, Rep2$jit, Rep3$jit) 
X[,,2] = cbind(Rep1$shi, Rep2$shi, Rep3$shi) 
X[,,3] = cbind(Rep1$hnr, Rep2$hnr, Rep3$hnr) 
X[,,4] = cbind(Rep1$mfcdel1, Rep2$mfcdel1, Rep3$mfcdel1) 
X[,,5] = cbind(Rep1$mfcdel2, Rep2$mfcdel2, Rep3$mfcdel2) 
X[,,6] = cbind(Rep1$mfcdel3, Rep2$mfcdel3, Rep3$mfcdel3)  
X[,,7] = cbind(Rep1$rpde, Rep2$rpde, Rep3$rpde) 
X[,,8] = cbind(Rep1$dfa, Rep2$dfa, Rep3$dfa) 
X[,,9] = cbind(Rep1$ppe, Rep2$ppe, Rep3$ppe) 
X[,,10] = cbind(Rep1$gne, Rep2$gne, Rep3$gne) 

Z = Rep1$sex

Xmid = R2 = rep(NA,K)
for(k in 1:K){
  Xmid[k] = median(X[,,k])
  R2[k] = (max(X[,,k])-min(X[,,k]))^2	
}

pdBin.data <- list(
  Y = status , 
  X = X , 
  Z = gender ,
  N = 80 ,
  X.new = X , 
  Z.new = gender ,
  N.new = 80 ,
  Xmid = Xmid ,
  R2 = R2
)

### Parameters
pdBin.params <- c(
  "Beta" , "Gsex" , "Gintercept" 
#  , "sigmaX" , 
#  "muWclust" , "tauWclust" , "qClust" 
#  , "P" , "P.new" 
)

### Initial values
pdBin.inits <- function(){	list(
  "Beta" = rnorm(K,0,0.1) ,
  "Gsex" = rnorm(1,0,0.1) ,
  "Gintercept" = rnorm(1,0,0.1) ,
  "W" = matrix(rnorm(N*K,0,0.1),N,K) , 
  "W.new" = matrix(rnorm(N*K,0,0.1),N,K) , 
  "sigmaX" = rep(1,K) , 
  "tauWclust" = matrix(1,K,3) ,
  "muWclust00" = t(matrix(c(-0.1,0,0.1),3,K))
)	}

##################################################
### Binary regression
### Replications

### Mixture of Normal
pdBin.fit <- jags.model("ParkinsonMixtureJAGS.bug", pdBin.data, pdBin.inits, n.chains=3)

update(pdBin.fit,2000)

pdBin.sample <- coda.samples(pdBin.fit, pdBin.params, n.iter=2000, thin=10)

pdBin.post <- summary(pdBin.sample)
pdBin.post
plot(pdBin.sample)

par(mfrow=c(4,4))
traceplot(pdBin.sample)

##################################################

Pfit <- pdBin.post$statistics[13:92,1]
Yfit <- ifelse(Pfit>0.5,1,0)
table(status,Yfit)

Ppred <- pdBin.post$statistics[93:172,1]
Ypred <- ifelse(Ppred>0.5,1,0)
table(status,Ypred)

auc(status,Ppred)

##################################################

##################################################
