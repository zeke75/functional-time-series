#######################################################
###############STUDY OF THE QUANTILES #################
#######################################################

m=100

q1.boot=array(0,c(5,nb.comp,m))
q2.boot=array(0,c(5,nb.comp,m))

q1.mc=array(0,c(5,nb.comp,m))
q2.mc=array(0,c(5,nb.comp,m))

for (k in 1:m){
	cat(paste("iteration=",k,"\n"))
	
	###################
	# BOOTSTRAP ERROR #
	###################
	nb.boot <- 200
	ERROR.boot <- matrix(0,n2*nb.boot,p)
	for (boot in 1:nb.boot){
		RESID.boot <-RESID[sample(1:n1,replace=T),]
		
		X.boot = matrix(0,nlearning,p)
		X.boot[1,] <- X[1,]
		X.boot[2:nlearning,] <- Psihat$Estimated.values + RESID.boot
		Psihat.boot <- ffunopare.knn(X.boot[2:nlearning,], CURVES1, CURVES2, neighbour=2, q=4, kind.of.kernel="quadratic", semimetric="pca")
		ERROR.boot[((boot-1)*n2+1):(boot*n2),] <- Psihat.boot$Predicted.values - Psihat$Predicted.values
	}
    
    	###################
	### TURE ERROR ####
	###################
    nb.mc <- 200
    ERROR.mc <-matrix(0,n2*nb.mc,p)
    for(mc in 1:nb.mc){
    	    X_1.mc <- t(cos(x))
    	    W <- t(cumsum(rnorm(p,mean=0,sd=1)))/sqrt(p)
    	    EPS1 <- W - t(x) * W[100]
    	    ERROR <- EPS1
    	    X.mc = matrix(0,n,p)
    	    X.mc[1,] <- X_1.mc + ERROR
    	
    	    C = 3
    	    for (i in 2:n){
    	    	    W <- t(cumsum(rnorm(p,mean=0,sd=1)))/sqrt(p)
    	    	    EPS1.mc <- W - t(x) * W[100]
    	    	    X.mc[i,] <- C*t(x)*cumsum(X.mc[i-1,])/p + EPS1.mc
    		}
    	    X.test.mc <- X.mc[251:300,]
    	    X.test2.mc <- X.mc[252:301,]
    	    Psihat.mc <- ffunopare.knn(X.mc[2:200,], X.mc[1:199,], X.test.mc, neighbour=2, q=4,kind.of.kernel="quadratic", semimetric="pca")
    	    ERROR.mc[((mc-1)*n2+1):(mc*n2),] <- Psihat.mc$Predicted.values - X.test2.mc
    }
    
    
    Bern.mc <- sample(c(-1,1),nb.mc, replace = TRUE)
    Weight.mc <- rep(Bern.mc, rep(n2,nb.mc))
    ERROR.mcc <- ERROR.mc * Weight.mc
    ERROR.mc.proj <- ERROR.mcc %*% VECTORS/p
    
    Bern.boot <- sample(c(-1,1),nb.boot, replace = TRUE)
    Weight.boot <- rep(Bern.boot, rep(n2,nb.boot))
    ERROR.bootc <- ERROR.boot * Weight.boot
    ERROR.boot.proj <- ERROR.bootc %*% VECTORS/p
    
    nb.comp <- 4
    for(ind in 1:5){
    	    for(comp in 1:nb.comp){
    	    	     Errors.boot <- ERROR.boot.proj[seq(ind,nb.boot*n2,n2),comp]
    	    	     Errors.mc <- ERROR.mc.proj[seq(ind,nb.mc*n2,n2),comp]
    	    	     
    		     ### 2.5% QUANTILE OF BOOTSTRAP ERROR
    		     q1.boot[ind,comp,k]=quantile(Errors.boot,prob=.025)
    		     ### 2.5% QUANTILE OF TRUE ERROR
    		     q1.mc[ind,comp,k]=quantile(Errors.mc,prob=.025)
    		     
    		     ### 97.5% QUANTILE OF BOOTSTRAP ERROR
    		     q2.boot[ind,comp,k]=quantile(Errors.boot,prob=.975)
    		     ### 97.5% QUANTILE OF TRUE ERROR
    		     q2.mc[ind,comp,k]=quantile(Errors.mc,prob=.975)
    		}
    	}
}

#### MEAN OF 2.5% QUANTILES #####
Q1.boot.mean = apply(q1.boot,3,mean)
Q1.mc.mean = apply(q1.mc,3,mean)

#### MEAN OF 97.5% QUANTILES #####
Q2.boot.mean = apply(q2.boot,3,mean)
Q2.mc.mean = apply(q2.mc,3,mean)



