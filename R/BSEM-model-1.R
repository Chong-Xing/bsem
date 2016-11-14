##--------------------------------##
## Chong Xing                     ## 
## Computing Project -  BSEM Test ##
##--------------------------------##

library(mvtnorm)   #Load mvtnorm package
library(R2OpenBUGS) #Load R2WinBUGS package

N=500
Iters <- 5
##Sample size
# BZ=numeric(N)                 #Fixed covariate in structural equation
# XI=matrix(NA, nrow=N, ncol=2) #Explanatory latent variables
# Eta=numeric(N)                #Outcome latent variables
# Y=matrix(NA, nrow=N, ncol=10)  #Observed variables

#The covariance matrix of xi
phi=matrix(c(1, 0.3, 0.3, 1), nrow=2)

#Estimates and standard error estimates
Eu=matrix(NA, nrow=Iters, ncol=10);    SEu=matrix(NA, nrow=Iters, ncol=10)
Elam=matrix(NA, nrow=Iters, ncol=7);   SElam=matrix(NA, nrow=Iters, ncol=7)
Eb=numeric(Iters);                     SEb=numeric(Iters)
Egam=matrix(NA, nrow=Iters, ncol=5);   SEgam=matrix(NA, nrow=Iters, ncol=5)
Esgm=matrix(NA, nrow=Iters, ncol=10);  SEsgm=matrix(NA, nrow=Iters, ncol=10)
Esgd=numeric(Iters);                   SEsgd=numeric(Iters)
Ephx=matrix(NA, nrow=Iters, ncol=3);   SEphx=matrix(NA, nrow=Iters, ncol=3)

R=matrix(c(1.0, 0.3, 0.3, 1.0), nrow=2)

parameters=c("u", "lam", "b", "gam", "sgm", "sgd", "phx")

init1=list(u=rep(0,10), lam=rep(0,7), b=0, gam=rep(0,5), psi=rep(1,10),
           psd=1, phi=matrix(c(1, 0, 0, 1), nrow=2))

init2=list(u=rep(1,10), lam=rep(1,7), b=1, gam=rep(1,5), psi=rep(2,10),
           psd=2, phi=matrix(c(2, 0, 0, 2), nrow=2))

inits=list(init1, init2)

dataGen <- function(N, phi) {
    
    BZ <- rt(N, 5)
    XI <- rmvnorm(N, c(0,0), phi)

    eps=numeric(10)    
    eps[1:3]=rnorm(3, 0, sqrt(0.3))
    eps[4:7]=rnorm(4, 0, sqrt(0.5))
    eps[8:10]=rnorm(3, 0, sqrt(0.4))

    delta=rnorm(N, 0, sqrt(0.36))
    
    Eta=0.5*BZ+0.4*XI[,1]+0.4*XI[,2]+0.3*XI[,1]*XI[,2]
    +0.2*XI[,1]*XI[,1]+0.5*XI[,2]*XI[,2]+delta

    Y=matrix(NA, nrow=N, ncol=10)
    
    Y[, 1] = Eta + eps[1]
    Y[, 2] = 0.9 * Eta+eps[2]
    Y[, 3] = 0.7 * Eta+eps[3]
    
    Y[, 4] = XI[, 1] + eps[4]
    Y[, 5] = 0.9 * XI[, 1] + eps[5]
    Y[, 6] = 0.7 * XI[, 1] + eps[6]
    Y[, 7] = 0.5 * XI[, 1] + eps[7]

    Y[, 8] = XI[, 2] + eps[8]
    Y[, 9] = 0.9 * XI[, 2] + eps[9]
    Y[, 10] = 0.7 * XI[, 2] + eps[10]

    list(BZ = BZ, XI = XI, eps = eps, delta = delta, Eta = Eta, Y = Y)
}



run1dataset <- function(run, N, phi, R, inits, parameters, n.iter = 10, n.burnin = 5, n.thin = 1) {
    
    mydata <- dataGen(N, phi)
    
    data = list(N = N, zero = c(0, 0), z = mydata$BZ, R = R, y = mydata$Y)

    model <- bugs(data, inits, parameters,
                  model.file = "BayesianCFA.txt", n.chains = 2,
                  n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)
}
    

modellist <- lapply(1:Iters, run1dataset, N, phi, R, inits, parameters, n.iter = 10, n.burnin = 5, n.thin = 1)


for (t in 1:Iters) {
    
    model <- modellist[[t]]
    
    #Save Estimates
    Eu[t,]=model$mean$u;               SEu[t,]=model$sd$u
    Elam[t,]=model$mean$lam;           SElam[t,]=model$sd$lam
    Eb[t]=model$mean$b;                SEb[t]=model$sd$b
    Egam[t,]=model$mean$gam;           SEgam[t,]=model$sd$gam
    Esgm[t,]=model$mean$sgm;           SEsgm[t,]=model$sd$sgm
    Esgd[t]=model$mean$sgd;            SEsgd[t]=model$sd$sgd
    Ephx[t,1]=model$mean$phx[1,1];     SEphx[t,1]=model$sd$phx[1,1]
    Ephx[t,2]=model$mean$phx[1,2];     SEphx[t,2]=model$sd$phx[1,2]
    Ephx[t,3]=model$mean$phx[2,2];     SEphx[t,3]=model$sd$phx[2,2]
}

Eu
modellist
model          # posterior descriptives
class(model)   # "bug"
summary(model)
model$DIC
hist(Elam[1, ])
hist(Elam[2, ])
hist(Eb)
plot(model) 
modellist[[1]]$mean
