############# functions

# The FindBifactor.orth() function is found in the syntax companion for
# Latent Variable Models (5th Edition) available from:
# http://sites.baylor.edu/lvm5

FindBifactor.orth <- function(A,reps=10,rotation="bifactor",solutions=1,round=8,maxit=1000,
seed=NA){
	require(GPArotation)
	m <- dim(A)[2]
	seed <- ifelse(!is.na(seed),round(seed),ceiling(runif(1, 0, 10^9)))
	set.seed(seed)
	ran.mat <- replicate(reps,Random.Start(m),simplify=FALSE)
	y <- lapply(ran.mat, function(z) GPForth(A,method=rotation,Tmat =z,maxit=maxit))
	y <- y[lapply(y, function(z) z$convergence) ==TRUE]
	criterion <- sapply(y, function(z) min(z$Table[,2]))
	results <- lapply(y, function(z) z$loadings)
	criterion <- round(criterion,round)
	results <- lapply(results, function(x) round(x,3))
	index.val <- 1:length(y)
	criterion.index <- data.frame(criterion,index.val)
	criterion.index <- criterion.index[order(criterion),]
	criterion.index.u <- criterion.index[match(unique(criterion.index$criterion), criterion.index$criterion),]
	index.keep <- criterion.index.u$index.val[1:solutions]
	output <- list(criterion=criterion[index.keep], loadings=results[index.keep])
	return(output)
}


HancockMueller<-function(x){
  loads<-inspect(x,"std")$lambda
  facs<-colnames(loads)
  out<-rep(NA,length(facs))
  names(out)<-facs
  
  for(i in facs){
    out[i]<-(1+(1/sum((loads[,i]^2)/(1-loads[,i]^2))))^-1
  }
  out
}

