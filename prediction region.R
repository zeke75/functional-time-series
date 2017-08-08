#################################################################
##### Prediction Region of Functional Autoregression (FAR1) #####
#################################################################
p=100
n=150
nlearning=100
ntest=50

x <- seq(0,1,length=p)     ## 100 by 1
#X_1 <- t(cos(x)) # X_1(t)  ## 1 by 100
X_1=rep(0,100)

############### Type 1 Error Process #################

W <- t(cumsum(rnorm(p,mean=0,sd=1)))/sqrt(p)
EPS1 <- W - t(x) * W[100]

#EPS1 <- 10*t((rnorm(p,mean=0,sd=1)))/sqrt(p)

ERROR <- EPS1
X = matrix(0,n,p)
X[1,] <- X_1 + ERROR
#plot(x,X[1,],type="l")

##########################################################################################
######### Autoregression Model: X_i+1(t) = int_0^1 psi(t,s) X_i(u) du + e_i+1(t) #########
##########################################################################################

###  Kernel: psi(s,t)  ####

C = 3
for (i in 2:n){
	W <- t(cumsum(rnorm(p,mean=0,sd=1)))/sqrt(p)
	EPS1 <- W - t(x) * W[100]
	
	#EPS1=10* t((rnorm(p,mean=0,sd=1)))/sqrt(p)
	
	X[i,] <- C*cumsum((x/p) * X[i-1,]) + EPS1
}
#par(mfrow=c(2,3))
#par(ask=T)
Index <- 101:105
#for (i in Index) plot(x,X[i,],type="l")
#### Plot of X_101,...X_105 #########
#matplot(x,t(X[Index,]),lwd=3,lty=1,type="l",xlab="",ylab="")
##################################
### uploading npfda R routines ###
##################################
##################################
#file("http://www.math.univ-toulouse.fr/staph/npfda/npfda-routinesR.txt",open="r")
#source("http://www.math.univ-toulouse.fr/staph/npfda/npfda-routinesR.txt")
#################################
#################################
CURVES1 <- X[1:(nlearning-1),]             # learning sample of curves
CURVES2 <- X[nlearning:(nlearning+49),]           # testing sample of curves
RESP1 <- X[2:nlearning,]               # learning sample of responses
RESP2 <- X[(nlearning+1):(nlearning+50),]             # testing sample of responses

n1 <- nrow(CURVES1)
n2 <- nrow(CURVES2)

Reg <- C*t(apply(t(t(X)*(x/p)),1,cumsum))
Reg2 <- Reg[(nlearning+1):(nlearning+50),]

Regmean <- apply(Reg,2,mean)
COVREG <- (t(Reg) - Regmean) %*% t(t(Reg) - Regmean)
COVREG <- COVREG/n

eigCOVREG <- eigen(COVREG)
Values <- eigCOVREG$values
VECTORS <- eigCOVREG$vectors
##############################################################################
##############################################################################
#  ASYMPTOTIC NORMALITY OF BOOTSTRAP (bootstrapped estimator minus estimator)
##############################################################################
##############################################################################


fit <- ffunopare.knn.gcv(RESP1, CURVES1, CURVES2, Knearest=2:20, q=2, kind.of.kernel="indicator", semimetric="pca")

h <- fit$knearest.opt
b <- 2*h

#b <- fit$knearest.opt
#h <- b

Psihat <- ffunopare.knn(RESP1, CURVES1, CURVES2, neighbour=b, q=2, kind.of.kernel="indicator", semimetric="pca")

Psihat2 <- ffunopare.knn(RESP1, CURVES1, CURVES2, neighbour=h, q=2, kind.of.kernel="indicator", semimetric="pca")

##################################
##### Compare Psi & Psihat #######
##################################

#readline(prompt="Press [enter] to continue:")

#dev.new()
#par(mfrow=c(2,2),ask=T)
for (i in 201:204){
	#plot(x,t(Psihat$Predicted.values[i-200,]),lty=2,type="l")
	#lines(x,t(C*t(x)*cumsum(X[i,])/p))
	#lines(x,C*cumsum((x/p)*X[i,]))
}
##############################
# Computing centered residuals
##############################
RESID <- RESP1 - Psihat$Estimated.values 
RESID <- t(t(RESID) - apply(RESID, 2, mean))

###################################
# Drawing bootstrap errors
###################################

nb.boot <- 1000
ERROR.boot <- matrix(0,n2*nb.boot,p)
error.norm1.boot <- rep(0,nb.boot)
error.norm2.boot <- rep(0,nb.boot)
error.norm3.boot <- rep(0,nb.boot)
error.norm4.boot <- rep(0,nb.boot)
error.norm5.boot <- rep(0,nb.boot)

#Y=matrix(0,B,p)
#norm=rep(0,B)


#plot(X[201,],type="l",ylim=c(-5,5),col="red")

for (boot in 1:nb.boot){
	#cat(paste("iteration=",boot,"\n"))
	RESID.boot <-RESID[sample(1:n1,replace=T),]
	#######################
    # Drawing new responses
    #######################
	X.boot = matrix(0,nlearning,p)
	X.boot[1,] <- X[1,]
	X.boot[2:nlearning,] <- Psihat$Estimated.values + RESID.boot
	
	Psihat.boot <- ffunopare.knn(X.boot[2:nlearning,], CURVES1, CURVES2, neighbour=h, q=2, kind.of.kernel="indicator", semimetric="pca")
	ERROR.boot[((boot-1)*n2+1):(boot*n2),] <- Psihat.boot$Predicted.values - Psihat$Predicted.values



observation.boot=Psihat$Predicted.values[1,]+RESID[sample(n1)[1],]
predictor.boot=Psihat.boot$Predicted.values[1,]

# l_1 norm
error.norm1.boot[boot] = (sum(abs(observation.boot-predictor.boot)))/p
# l_2 norm
error.norm2.boot[boot] = sqrt(sum((observation.boot-predictor.boot)^2)/p)
# l_inf norm
error.norm3.boot[boot] = max(abs(observation.boot-predictor.boot))
# 1st coordinate
error.norm4.boot[boot] = abs(observation.boot[1]-predictor.boot[1])
# 1st component
error.norm5.boot[boot] = abs((observation.boot-predictor.boot) %*% VECTORS[,1])

}

predictor=Psihat2$Predicted.values[1,]

print(quantile(error.norm1.boot,.95))
print(quantile(error.norm1.boot,.90))

print(quantile(error.norm2.boot,.95))
print(quantile(error.norm2.boot,.90))

B=1000
Y=matrix(0,B,p)
norm1=rep(0,B)
norm2=rep(0,B)
norm3=rep(0,B)
norm4=rep(0,B)
norm5=rep(0,B)


for  (k in 1:B){

W <- t(cumsum(rnorm(p,mean=0,sd=1)))/sqrt(p)
ERR <- (W - t(x) * W[100])


Y[k,]= C*cumsum((x/p) * X[nlearning,]) + ERR
norm1[k]=(sum(abs(Y[k,]-predictor)))/p
norm2[k]=sqrt(sum((Y[k,]-predictor)^2)/p)
norm3[k]=max(abs(Y[k,]-predictor))
norm4[k]=abs(Y[k,1]-predictor[1])
norm5[k]=abs((Y[k,]-predictor) %*% VECTORS[,1])


}
CVR1=sum(norm1<quantile(error.norm1.boot,.95))/B
CVR2=sum(norm1<quantile(error.norm1.boot,.90))/B
CVR3=sum(norm2<quantile(error.norm2.boot,.95))/B
CVR4=sum(norm2<quantile(error.norm2.boot,.90))/B
CVR5=sum(norm3<quantile(error.norm3.boot,.95))/B
CVR6=sum(norm3<quantile(error.norm3.boot,.90))/B
CVR7=sum(norm4<quantile(error.norm4.boot,.95))/B
CVR8=sum(norm4<quantile(error.norm4.boot,.90))/B
CVR9=sum(norm5<quantile(error.norm5.boot,.95))/B
CVR10=sum(norm5<quantile(error.norm5.boot,.90))/B

print(CVR1)
print(CVR2)
print(CVR3)
print(CVR4)
print(CVR5)
print(CVR6)
print(CVR7)
print(CVR8)
print(CVR9)
print(CVR10)
print(b)