# For HPC




# make data
makedata <- function(path, pattern="*.txt", TT, B, p, N, J, standardize=TRUE){
  temp <- list.files(path=path,  pattern=pattern) 
  eta <- eta.sd <- vector("list", N*J)
  
  for(n in 1:N){
    eta.temp <- as.matrix(read.table(paste(path, "\\", temp[n], sep=""), header=FALSE))
    
    if(standardize){
      
      eta.sd[[(n-1)*J+1]] <- sd1 <- apply(eta.temp[1:TT,], 2, sd)
      eta.sd[[(n-1)*J+2]] <- sd2 <- apply(eta.temp[(TT+1):(2*TT),], 2, sd)
      eta.sd[[(n-1)*J+3]] <- sd3 <- apply(eta.temp[(2*TT+1):(3*TT),], 2, sd)
      eta.sd[[(n-1)*J+4]] <- sd4 <- apply(eta.temp[(3*TT+1):(4*TT),], 2, sd)
      
      
      eta[[(n-1)*J+1]] <- t(t(eta.temp[1:TT,])/sd1)
      eta[[(n-1)*J+2]] <- t(t(eta.temp[(TT+1):(2*TT),])/sd2)
      eta[[(n-1)*J+3]] <- t(t(eta.temp[(2*TT+1):(3*TT),])/sd3)
      eta[[(n-1)*J+4]] <- t(t(eta.temp[(3*TT+1):(4*TT),])/sd4)
      
      
    } else {
      eta[[(n-1)*J+1]] <- eta.temp[1:TT,]
      eta[[(n-1)*J+2]] <- eta.temp[(TT+1):(2*TT),]
      eta[[(n-1)*J+3]] <- eta.temp[(2*TT+1):(3*TT),]
      eta[[(n-1)*J+4]] <- eta.temp[(3*TT+1):(4*TT),]
    }
    
    
    cat("n=", n, "\n")
  }
  H <- lapply(eta, function(x){H <- NULL; for(i in 1:(TT-p)){H <- cbind(H, c(t(x[(i+p-1):i,])))}; H})
  eta0 <- lapply(eta, function(x){t(x[(p+1):TT,])})
  HH <- lapply(H, tcrossprod)
  w_MLE <- lapply(1:(J*N), function(x){c(eta0[[x]] %*% t(H[[x]]) %*% chol2inv(chol(HH[[x]])))})
  return(list(eta=eta,
              H=H,
              eta0=eta0,
              HH=HH,
              w_MLE=w_MLE,
              eta.sd=eta.sd))
}
# save.image("HCP_all_p2.RData")



eta2H <- function(x){
  H <- NULL
  for(i in 1:(TT-p)){
    H <- cbind(H, c(t(x[(i+p-1):i,])))
  } 
  H
}


