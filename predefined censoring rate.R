##########################################################################
###
### generate simulated survival data with predefiend censoring rate 

### Author: Zihang Zhong
### Date: Aug 2020

###  Specifications: 
### - the main aim here is to root-find the censoring parameter for censoring distribution and generate the right censoring data
### - As the simulated righted censored data requires two independent survival distributions,one for event time EVENT,
###    and the other for censoring time C. here, the baseline hazard function was set as from Weibull and the censoring 
###    time from uniform distribution. other distribution is also allowed but not discussed here. 
### - considering two kind right-censored survival data: 
###    A) data with no consideration of covariates
###    B) data incorporating with covariates
###       B.1) all covariates are binary variables
###       B.2) all covariates are normal variables -- default are standard normal
###       B.3) all covariates are mixed, the mix of normal, binary, uniform or possion 
###   
###
### Structure: 
### - 1) Support functions

###     1.1) data generation 
###     - CenDatNoCov
###     - CenDatBin
###     - CenDatNorm
###     - CenDatMixed

###     1.2) core code to calculate the censoring parameter 
###     - CensProbBin
###     - CensProbNorm
###     - CensProbMixed

### - 2) verify the average censoring rate of simulated data is coincident with the nominal censoring rate by simulations
###
###
### Parameters:

### - alpha: the shape parameter for hazard function(weibull distribuion)
### - lambda: the scale parameter for hazard function(weibull distribuion)
### - cens.p: the nominal censoring rate
### - beta: the predefined coefficients of the covariates
### - class: the distribution of the covariates
### - para: the distribution parameters of the covariates

### - theta: the censoring parameter
### - size: the sample size
### - nSims: the number of simulations

###################################################################################
###################################################################################
# library(simsurv)


#====================================================================
#A) data with  no covariates--CenDatNoCov
#====================================================================
CenDatNoCov <- function(alpha,lambda,cens.p,size){
  theta <- round(uniroot(function(x)lambda/(alpha*x)*pgamma((x/lambda)^alpha, 1/alpha) * gamma(1/alpha)-cens.p,c(0.00001,1000))$root,3)

  data<-data.frame(
    T = rweibull(size,alpha,lambda),
    C = runif(size,0,theta))
  cens.data <- data.frame(time = ifelse(data$T<=data$C,data$T,data$C),
                          status = ifelse(data$T<=data$C,1,0))
    return(list(theta = theta, cens.data = cens.data))
  }

CenDatNoCov(alpha = 2,lambda = 4,cens.p = 0.3,size =200)

# ---------------------------------------------------------------------------
# simulate 1000 times to verify the value of theta is ruboost.
# ---------------------------------------------------------------------------
SimCensP <- c()
set.seed(20200810)
for (i in 1:1000){
  data <- CensDataNoCov(alpha = 2,lambda = 4,cens.p = 0.3,size =200)$cens.data
  SimCensP[i] <- sum(data["status"])/200
}
mean(SimCensP)


#====================================================================
# B) incorporating with covariates
#====================================================================

# ====================================================================
# B.1) all covariates are binary variables: generate the censoring data -- CenDatBin
# ====================================================================
CenDatBin <- function(alpha,lambda,cens.p,beta,p,size,seed=20200812){
  alpha = alpha
  lambda =lambda
  cens.p = cens.p
  beta = beta
  p = p
  n = size
  # ---------------------------------------------------------------------------------
  #  the core code to calculate the censoring parameter : a more adaptive function -- censProbBin
  censProbBin <- function(theta){
   beta.0 <--alpha*log(lambda)
    
    LenCovar <- length(beta) 
    CombInd <-list(NULL)
    ncomb <- c()
    finalInd <- c()
    for (i in 1:LenCovar){
      ind<-combn(1:LenCovar,i)
      CombInd[[i]]<-ind
      ncomb[i]<-ncol(CombInd[[i]])
      for (j in 1:ncomb[i]) {
        location <- c(CombInd[[i]][, j],rep(0,LenCovar-i))
        finalInd<- cbind(finalInd,location)
      }
    }
    comb<-matrix(0,nrow = LenCovar, ncol = 2^LenCovar)
    for(i in 1: 2^LenCovar-1){
      comb[finalInd[,i],i]=1
    }
    lambda.i<-exp(-(beta.0+apply(beta*comb,2,sum))/alpha)
    cond.cens.Prob<-(lambda.i/(alpha*theta))*pgamma((theta/lambda.i)^alpha, 1/alpha) * gamma(1/alpha)
    pdf.lambda.i<-apply(matrix(p,nrow=LenCovar,ncol = 2^LenCovar)**comb*((1-matrix(p,nrow=LenCovar,ncol = 2^LenCovar))**(1-comb)),2,prod)
    
    return(t(cond.cens.Prob)%*%pdf.lambda.i)
  }
  # ---------------------------------------------------------------------------------
  
  theta <- round(uniroot(function(x)censProbBin(x)-cens.p,c(0.0000001,100))$root,3)
  
  set.seed(seed)
  cov<-mapply(rbinom,n,1,p)
  colnames(cov) <- paste("x",1:length(p),sep = "")
  EVENT<-rweibull(n,alpha,lambda*exp(-1/alpha*(cov%*%beta)))
  data <-as.data.frame(cbind(id = 1:n,cov = cov,EVENT = EVENT,C=runif(n,0,theta))) 
  
  data$time = ifelse(data$EVENT<=data$C,data$EVENT,data$C)
  data$status = ifelse(data$EVENT<=data$C,1,0)
  return(list(theta = theta, cens.data = data))
}

CenDatBin(alpha = 2,lambda = 4,cens.p = 0.3,beta = c(-0.1,0.2,-0.3,0.4),p = c(0.3,0.4,0.5,0.6),size = 200,seed=20200812)

# -----------------------------------------------------------------------------------
#  simulations to verify
# -----------------------------------------------------------------------------------
cens.p <- c()
set.seed(20200812)
n <- 200
for(i in 1:1000){
  data <-CenDatBin(alpha = 2,lambda = 4,cens.p = 0.3,beta = c(-0.1,0.2,-0.3,0.4),p = c(0.3,0.4,0.5,0.6),size = n,seed=2020+i)$cens.data
  cens.p[i]<-sum(data["status"])/n 
}
1-mean(cens.p)
# coxph(Surv(time,status))~x1+x2+x3+x4,data = data)


# ====================================================================
# B.2) all covariates are normal variables :generate the censoring data -- CenDatNorm
# default all are from standard normal distribution
# ====================================================================s
CenDatNorm <- function(alpha,lambda,cens.p,beta,size,seed=20200812){
  alpha = alpha
  lambda =lambda
  cens.p = cens.p
  beta = beta
  n = size
  # ---------------------------------------------------------------------------------
  #  2.1) the core code to calculate the censoring parameter : a more adaptive function -- censProbNorm
  censProbNorm<-function(theta){
    beta.0 <- -alpha*log(lambda)
    PdfLambdai<-function(u) dlnorm(u,-beta.0/alpha, beta%*%beta/alpha^2)
    CondCensProb<-function(u) (u/(alpha*theta))*pgamma((theta/u)^alpha, 1/alpha) * gamma(1/alpha)
    cens.Prob<-integrate(function(u){PdfLambdai(u)*CondCensProb(u)},-Inf,Inf)$value
    return(cens.Prob)
  }
  # ---------------------------------------------------------------------------------
  
  theta <- round(uniroot(function(x)censProbNorm(x)-cens.p,c(0.0001,1000))$root,3)
  
  set.seed(seed)
  cov<-mapply(rnorm,n,0,rep(1,length(beta)))
  colnames(cov) <- paste("x",1:length(beta),sep = "")
  EVENT<-rweibull(n,alpha,lambda*exp(-1/alpha*(cov%*%beta)))
  data <-as.data.frame(cbind(id = 1:n,cov = cov,EVENT = EVENT,C=runif(n,0,theta))) 
  
  data$time = ifelse(data$EVENT<=data$C,data$EVENT,data$C)
  data$status = ifelse(data$EVENT<=data$C,1,0)
  return(list(theta = theta, cens.data = data))
} 
CenDatNorm(alpha = 2,lambda = 4,cens.p = 0.3,beta = c(-0.1,0.2,-0.3,0.4),size = 200,seed=20200812)

# -----------------------------------------------------------------------------------
#  simulations to verify
# -----------------------------------------------------------------------------------
cens.p<-c()
set.seed(20200810)
n<-200
for(i in 1:1000){
  data <- CenDatNorm(alpha = 2,lambda = 4,cens.p = 0.3,beta = c(-0.1,0.2,-0.3,0.4),size = n, seed=202008+i)$cens.data
  cens.p[i]<-sum(data["status"])/n 
}
1-mean(cens.p)


# ====================================================================
# B.3) all covariates are mixed variables : generate the censoring data -- CenDatMixed
# - specify the distribution of all variables: "N" = normal,"B"  = binary, "P"= possion, "U"= uniform
# - specify the distribution parameters in list for "B", "P", "U";for "U", the start and end should be set in matrix with two colomns
# ====================================================================
CenDatMixed <- function(alpha,lambda,cens.p,beta,size,
                       class = c("N","B","P","U"),
                       para = list(B=c(0.5),P = c(5), U = matrix(c(0,1),ncol = 2,byrow = TRUE)),
                       seed=20200812){
  alpha = alpha
  lambda =lambda
  cens.p = cens.p
  beta = beta
  size = size
  seed = seed
  # ---------------------------------------------------------------------------------
  # the core code to calculate the censoring parameter : a more adaptive function -- CensProbMixed
  CensProbMixed <- function(class,para,theta){
    set.seed(20200810)

    n<-10000  
    nNorm <- length(class[which(class =="N")])
    nBin <- length(class[which(class =="B")])
    nPos <- length(class[which(class =="P")])
    nUni <- length(class[which(class =="U")])
    
    if (nNorm > 0){
      a<-as.matrix(mapply(rnorm,n,0,rep(1,nNorm)))
    }
    if (nBin > 0){
      b<-as.matrix(mapply(rbinom,n,1,para$B))
    }
    if (nPos > 0){
      c<-as.matrix(mapply(rpois,n,para$P))
    }
    
    if (nUni > 0){
      d<-as.matrix(mapply(runif,n,para$U[,1],para$U[,2]))
    }
    x <-cbind(a,b,c,d)
    
    beta.0<--alpha*log(lambda)
    lambda.i<-exp(-(beta.0+x%*%beta)/alpha)
    
    max.lambda.i<-max(lambda.i)
    min.lambda.i<-min(lambda.i)
    
    PdfLambdai<-function(u){
      dens<-density(lambda.i,bw = "nrd0",kernel="gaussian",na.rm=TRUE)
      y.loess<-loess(dens$y~dens$x,span = 0.1)
      pred.y<-predict(y.loess,newdata = u)
      return(pred.y)
    }
    
    CondCensProb<-function(u) (u/(alpha*theta))*pgamma((theta/u)^alpha, 1/alpha) * gamma(1/alpha)
    cens.Prob<-integrate(function(u){PdfLambdai(u)*CondCensProb(u)},min.lambda.i,max.lambda.i)$value
    return(cens.Prob)
  }
  # ---------------------------------------------------------------------------------  
  theta <- round(uniroot(function(y) CensProbMixed(class = class,para = para,y)-cens.p,c(0.001,100))$root,3)

  set.seed(seed)
  nNorm <- length(class[which(class =="N")])
  nBin <- length(class[which(class =="B")])
  nPos <- length(class[which(class =="P")])
  nUni <- length(class[which(class =="U")])
  
  if (nNorm > 0){
    a<-as.matrix(mapply(rnorm,size,0,rep(1,nNorm)))
  }
  if (nBin > 0){
    b<-as.matrix(mapply(rbinom,size,1,para$B))
  }
  if (nPos > 0){
    c<-as.matrix(mapply(rpois,size,para$P))
  }
  
  if (nUni > 0){
    d<-as.matrix(mapply(runif,size,para$U[,1],para$U[,2]))
  }
  cov <-cbind(a,b,c,d)
  colnames(cov) <- paste("x",1:length(beta),sep = "")
  EVENT<-rweibull(size,alpha,lambda*exp(-1/alpha*(cov%*%beta)))
  data <-as.data.frame(cbind(id = 1:size,cov = cov,EVENT = EVENT,C=runif(size,0,theta))) 
  
  data$time = ifelse(data$EVENT<=data$C,data$EVENT,data$C)
  data$status = ifelse(data$EVENT<=data$C,1,0)
  return(list(theta = theta, cens.data = data))
} 

CenDatMixed(alpha = 2,lambda = 4,cens.p = 0.3,beta = c(-0.1,0.2,-0.3,0.4),size=200,
                        class = c("N","B","P","U"),
                        para = list(B=c(0.5),P = c(5), U = matrix(c(0,1),ncol = 2,byrow = TRUE)),
                        seed=20200812)
  
# -----------------------------------------------------------------------------------
#  simulations to verify
# -----------------------------------------------------------------------------------
cens.p<-c()
for(i in 1:1000){
  data <- CenDatMixed(alpha = 2,lambda = 4,cens.p = 0.3,beta = c(-0.1,0.2,0.5,-0.6,-0.3,0.4),size=200,
              class = c("N","B","B","P","U","U"),
              para = list(B=c(0.5,0.4),P = c(5), U = matrix(c(0,1,0.2,0.6),ncol = 2,byrow = TRUE)),
              seed=202008+i)$cens.data
  cens.p[i]<-sum(data["status"])/200 
}
1-mean(cens.p)




