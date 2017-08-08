######################################################
###### Bootstrap of Functional Autoregression ########
######################################################

p=100
n=400
nlearning=200
ntest=50

x <- seq(0,1,length=p)     ## 100 by 1
X_1 <- t(cos(x)) # X_1(t)  ## 1 by 100

############### Type 1 Error Process #################
W <- t(cumsum(rnorm(p,mean=0,sd=1)))/sqrt(p)
EPS1 <- W - t(x) * W[100]

############### Type 2 Error Process #################
a1 <- rnorm(1,0,1)
a2 <- rnorm(1,0,1)
EPS2 <- a1*sqrt(2)*sin(2*pi*t(x)) + a2*cos(2*pi*t(x))

############### Type 1 Error Process #################
a=1
EPS3 <- EPS2 + a * EPS1

ERROR <- 0.1*EPS3

X = matrix(0,n,p)
X[1,] <- X_1 + ERROR

plot(x,X[1,],type="l")

##########################################################################################
######### Autoregression Model: X_i+1(t) = int_0^1 psi(t,s) X_i(u) du + e_i+1(t) #########
##########################################################################################

### Sloping plane(t) Kernel: psi(s,t) = C * t ####

C = 1.6
for (i in 2:nlearning){
	W <- t(cumsum(rnorm(p,mean=0,sd=1)))/sqrt(p)
	EPS1 <- W - t(x) * W[100]
	a1 <- rnorm(1,0,1)
	a2 <- rnorm(1,0,1)
	EPS2 <- a1*sqrt(2)*sin(2*pi*t(x)) + a2*cos(2*pi*t(x))
	a=1
	EPS3 <- EPS2 + a * EPS1

	X[i,] <- C*t(x)*cumsum(X[i-1,])/p + 0.1*EPS3
}

par(mfrow=c(2,3))
Index <- 101:105
for (i in Index) plot(x,X[i,],type="l")
matplot(x,t(X[Index,]),lwd=3,lty=1,type="l",xlab="",ylab="")


### Sloping plane(s) Kernel: psi(s,t) = C * s ####

### Indentity Kernel: psi(s,t) = C  ####

### Guassian Kernel ####




##################################
### uploading npfda R routines ###
##################################
##################################
file("http://www.math.univ-toulouse.fr/staph/npfda/npfda-routinesR.txt",open="r")
source("http://www.math.univ-toulouse.fr/staph/npfda/npfda-routinesR.txt")
#################################
#################################


##############################################################################
##############################################################################

CURVES1 <- X[1:200,]             # learning sample of curves
CURVES2 <- X[251:300,]           # testing sample of curves
RESP1 <- X[2:201,]               # learning sample of responses
RESP2 <- X[252:301,]             # testing sample of responses



#####################################################
############# Compute Empirical Basis ###############
#####################################################


Regmean <- apply(X,2,mean)
COVREG <- (t(X) - Regmean) %*% t(t(X) - Regmean)
COVREG <- COVREG/n

eigCOVREG <- eigen(COVREG)
Values <- eigCOVREG$values
VECTORS <- eigCOVREG$vectors



##############################################################################
##############################################################################
#  ASYMPTOTIC NORMALITY OF BOOTSTRAP (bootstrapped estimator minus estimator)
##############################################################################
##############################################################################
res <- ffunopare.knn.gcv(RESP1, CURVES1, CURVES2, 4, Knearest=6:14, semimetric="pca")
p <- ncol(CURVES1)
n1 <- nrow(CURVES1)
n2 <- nrow(CURVES2)
u <- sqrt(5)/10
pu <- 0.5+u
bwboot1 <- res$knearest.opt
bwboot2 <- bwboot1
rhat.bwboot1 <- ffunopare.knn(RESP1, CURVES1, CURVES2, bwboot1, 4, semimetric="pca")
##############################
# Computing centered residuals
##############################
RESID <- RESP1 - rhat.bwboot1$Estimated.values 
RESID <- t(t(RESID) - apply(RESID, 2, mean))
###################################
# Drawing "wild-bootstrapped errors
###################################
nb.boot <- 10
ERROR.BOOT <- matrix(0,n2*nb.boot,p)
for(boot in 1:nb.boot){
    cat(paste("iteration=",boot,"\n"))
    
    U <- runif(n1)
    U[U <= pu] <- 0.5 * (1 - sqrt(5))
    U[U > pu] <- 0.5 * (1 + sqrt(5))
    RESID.BOOT <- RESID * U
    RESID.BOOT <- t(t(RESID.BOOT) - apply(RESID.BOOT,2,mean))
    #######################
    # Drawing new responses
    #######################
    
    X.BOOT = matrix(0,nlearning+1,p)
    X.BOOT[1,]=X[1,]
    for (i in 2:(nlearning+1)){
    	     PSI.hat = ffunopare.knn(RESP1, CURVES1, X.BOOT[i-1,],bwboot1, 4, semimetric="pca")
    	     X.BOOT[i,] = PSI.hat$Predicted.values + RESID.BOOT[i-1,]
    	}
    
    # Computing predicted fat content
    #################################
    rhat.bwboot2 <- ffunopare.knn.single(X.BOOT[2:201,], X.BOOT[1:200,], CURVES2, bwboot2, 4, semimetric="pca")
    ERROR.BOOT[((boot-1)*n2+1):(boot*n2),] <- rhat.bwboot2$Predicted.values - rhat.bwboot1$Predicted.values
}


##############################################################################
##############################################################################
# ASYMPTOTIC NORMALITY OF ESTIMATOR MINUS TRUE VALUE  (MONTE-CARLO SCHEME)
##############################################################################
##############################################################################

nb.mc <- 200
ERROR.MC <-matrix(0,n2*nb.mc,p)
for(mc in 1:nb.mc){
###########################################################
# SIMULATING FAR MODEL
###########################################################
       #a.mc <- runif(1,min=0.5,max=1.5)
       X_1.MC <- t(cos(x))
       
       W <- t(cumsum(rnorm(p,mean=0,sd=1)))/sqrt(p)
       EPS1 <- W - t(x) * W[100]
       a1 <- rnorm(1,0,1)
       a2 <- rnorm(1,0,1)
       EPS2 <- a1*sqrt(2)*sin(2*pi*t(x)) + a2*cos(2*pi*t(x))
       a=1
       EPS3 <- EPS2 + a * EPS1
       ERROR <- EPS3
       
       X.MC = matrix(0,n,p)
       X.MC[1,] <- X_1.MC + ERROR

       C = 1.6
       for (i in 2:n){
       	   W <- t(cumsum(rnorm(p,mean=0,sd=1)))/sqrt(p)
       	   EPS1 <- W - t(x) * W[100]
       	   a1 <- rnorm(1,0,1)
       	   a2 <- rnorm(1,0,1)
       	   PS2 <- a1*sqrt(2)*sin(2*pi*t(x)) + a2*cos(2*pi*t(x))
       	   a=1
       	   EPS3.MC <- EPS2 + a * EPS1
       	   X.MC[i,] <- C*t(x)*cumsum(X.MC[i-1,])/p + EPS3.MC
       }
       	
       #Regmean.MC <- apply(X.MC,2,mean)
       #COVREG.MC <- (t(X) - Regmean.MC) %*% t(t(X) - Regmean.MC)
       #COVREG.MC <- COVREG.MC/n
       	
       #eigCOVREG.MC <- eigen(COVREG.MC)
       #Values.MC <- eigCOVREG.MC$values
       #VECTORS.MC <- eigCOVREG.MC$vectors
       
       
       X.TEST.MC <- X.MC[251:300,]
       X.TEST2.MC <- X.MC[252:301,]
       
       rhat.bwboot2.mc <- ffunopare.knn.single(X.MC[2:201,], X.MC[1:200,], X.TEST.MC, bwboot2, 4, semimetric="pca")
       ERROR.MC[((mc-1)*n2+1):(mc*n2),] <- rhat.bwboot2.mc$Predicted.values - X.TEST2.MC
       cat(paste("iteration=", mc, "\n"))


}



#####################################################################################
# DISPLAYING COMPONENTWISE COMPARISONS BETWEEN BOOTSTRAP AND MC DISTRIBUTION
#####################################################################################
Bern.mc <- sample(c(-1,1),nb.mc, replace = TRUE)
Weight.mc <- rep(Bern.mc, rep(n2,nb.mc))
ERROR.MCC <- ERROR.MC * Weight.mc
ERROR.MC.PROJ <- ERROR.MCC %*% VECTORS/p
Bern.boot <- sample(c(-1,1),nb.boot, replace = TRUE)
Weight.boot <- rep(Bern.boot, rep(n2,nb.boot))
ERROR.BOOTC <- ERROR.BOOT * Weight.boot
ERROR.BOOT.PROJ <- ERROR.BOOTC %*% VECTORS/p
nb.comp <- 4
par(mfrow=c(5,nb.comp), mar=c(2,2,2,1))
for(ind in 1:5){
	for(comp in 1:nb.comp){
		width=6
		Errors.boot <- ERROR.BOOT.PROJ[seq(ind,nb.boot*n2,n2),comp]
		#Errors.boot.std <- Errors.boot/sqrt(var(Errors.boot))
		Errors.boot.std <- Errors.boot
		bandwidth.boot <- bw.nrd(Errors.boot.std)
		Density.boot <- density(Errors.boot.std,bw=bandwidth.boot,from=-width,to=width)
		Errors.mc <- ERROR.MC.PROJ[seq(ind,nb.mc*n2,n2),comp]
		#Errors.mc.std <- Errors.mc/sqrt(var(Errors.mc))
		Errors.mc.std <- Errors.mc
		bandwidth.mc <- bw.nrd(Errors.mc.std)
		Density.mc <- density(Errors.mc.std,bw=bandwidth.mc,from=-width,to=width)
		xlimits=range(c(Density.boot$x,Density.mc$x))
		ylimits=range(c(Density.boot$y,Density.mc$y))
		plot(Density.mc$x, Density.mc$y, xlab="",ylab="",main=paste("curve", ind,"- comp",comp), type="l", xlim=xlimits, ylim=ylimits, lwd=1)
		abline(v=0)
		par(new=T)
		plot(Density.boot$x, Density.boot$y, xlab="",ylab="", type="l", xlim=xlimits, ylim=ylimits, lty=2, lwd=2)
		}
}




