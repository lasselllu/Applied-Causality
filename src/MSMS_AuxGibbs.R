#################################
# Fit SVAR by Gibbs Sampler     #
# Model 3                       #  
# multi subject multi session   #
#                               #
# using auxilliary variables    #
# to break variable dependence  #
#################################



library(MASS)
library(statmod)
library(mvnfast)




MSMSGibbs <- function(eta,                          # array of dim TT * B * J * N
                      H,
                      HH,
                      eta0,
                      w_MLE,
                      TT,
                      B,
                      p=2, 
                      N,
                      J,
                      seed=100, 
                      numIter=2000, 
                      thin=1, 
                      K_inv = diag(1/(B-1), B),      # inverse of scale matrix of Lambda prior
                      nu = 1,                        # df of Lambda prior
                      mu1 = 1,                       # mean of lambda1sq prior
                      nu1 = 0.001,                   # df of lambda1sq prior
                      mu2 = 1,                       # mean of lambda2 prior
                      nu2 = 0.01,                    # df of lambda2 prior
                      k1 = 200,                      # shape prior of omega_u
                      theta1 = 0.005,                # scale prior of omega_u
                      k2 = 200,                      # shape prior of omega_v
                      theta2 = 0.005,                # scale prior of omega_v
                      a = 1,
                      b = 1,
                      initial=NULL
                      # k3 = 200,
                      # theta3 = 0.005,
                      # k4 = 200,
                      # theta4 = 0.005
                      ){               
  
 

  # store values
  sim.u_j <- sim.v_n <- sim.w <- vector(mode="list", length=numIter/thin+1)
  sim.Lambda <- sim.lambda1sq <- sim.lambda2 <- sim.tausq <- vector(mode="list", length=numIter/thin+1)
  sim.U_u <- sim.U_v <- sim.alpha <- sim.beta <- vector(mode="list", length=numIter/thin+1)
  sim.v_n.star <- sim.u_j.star <- vector(mode="list", length=numIter/thin+1)
  sim.omega_u <- sim.omega_v <- vector(mode="list", length=numIter/thin+1)
  
  
  # initialize
  set.seed(seed)
  if(is.null(initial)){
    sim.w[[1]] <- w <- rnorm(B^2*p, 0, 1)
    sim.v_n[[1]] <- v_n <- lapply(vector("list", N), function(x){rnorm(B^2*p, 0, 1)})
    sim.u_j[[1]] <- u_j <- lapply(vector("list", J), function(x){rnorm(B^2*p, 0, 1)})
    sim.Lambda[[1]] <- Lambda <- Posdef(n=B, ev=runif(B, 0, 0.001))
    sim.lambda1sq[[1]] <- lambda1sq <- runif(B^2*p, 0.1, 1000)
    sim.lambda2[[1]] <- lambda2 <- runif(B^2*p, 0.1, 10)
    sim.tausq[[1]] <- tausq <- runif(B^2*p, 0.1, 10)
    sim.U_v[[1]] <- U_v <- runif(B^2*p, 0.1, 1)
    sim.U_u[[1]] <- U_u <- runif(B^2*p, 0.1, 1)
    sim.alpha[[1]] <- alpha <- runif(B^2*p, 0.1, 1)          
    sim.beta[[1]] <- beta <- runif(B^2*p, 0.1, 1)   
    sim.v_n.star[[1]] <- lapply(1:N, function(n){alpha*v_n[[n]]})
    sim.u_j.star[[1]] <- lapply(1:J, function(j){beta*u_j[[j]]})
    sim.omega_v[[1]] <- U_v/alpha^2
    sim.omega_u[[1]] <- U_u/beta^2
  } else {
    sim.w[[1]] <- w <- initial$w
    sim.v_n[[1]] <- v_n <- initial$v_n
    sim.u_j[[1]] <- u_j <- initial$u_j
    sim.Lambda[[1]] <- Lambda <- initial$Lambda
    sim.lambda1sq[[1]] <- lambda1sq <- initial$lambda1sq
    sim.lambda2[[1]] <- lambda2 <- initial$lambda2
    sim.tausq[[1]] <- tausq <- initial$tausq
    sim.U_v[[1]] <- U_v <- initial$U_v
    sim.U_u[[1]] <- U_u <- initial$U_u
    sim.alpha[[1]] <- alpha <- initial$alpha          
    sim.beta[[1]] <- beta <- initial$beta    
    sim.v_n.star[[1]] <- v_n.star <- initial$v_n.star
    sim.u_j.star[[1]] <- u_j.star <- initial$u_j.star
    sim.omega_v[[1]] <- omega_v <- initial$omega_v
    sim.omega_u[[1]] <- omega_u <- initial$omega_u
  }

  
  
  for(i in 1:numIter){
    
    # update parameters
    temp <- Gamma_recover(tausq=tausq, lambda1sq=lambda1sq, Lambda=Lambda, B=B, p=p)      
    G <- temp$G
    xisq <- temp$xisq
    prec <- prec_get(HH=HH, Lambda=Lambda)
    w <- w_update(prec=prec, w_MLE=w_MLE, v_n=v_n, u_j=u_j, alpha=alpha, beta=beta, tausq=tausq, lambda2, N=N, J=J)
    v_n <- v_n_update(prec=prec, w_MLE=w_MLE, w=w, u_j=u_j, alpha=alpha, beta=beta, U_v=U_v, N=N, J=J)
    u_j <- u_j_update(prec=prec, w_MLE=w_MLE, w=w, v_n=v_n, alpha=alpha, beta=beta, U_u=U_u, N=N, J=J)
    Lambda <- Lambda_update(eta0=eta0, w=w, v_n=v_n, u_j=u_j, alpha=alpha, beta=beta, H=H, G=G, K_inv=K_inv, nu=nu, TT=TT, B=B, p=p, N=N, J=J)
    lambda1sq <- lambda1sq_update(tausq=tausq, xisq=xisq, B=B, p=p, mu1=mu1, nu1=nu1)          
    lambda2 <- lambda2_update(w=w, B=B, p=p, mu2=mu2, nu2=nu2)                         
    tausq <- tausq_update(lambda1sq=lambda1sq, w=w, xisq=xisq, B=B, p=p) 
    U_v <- U_v_update(v_n=v_n, k1=k1, theta1=theta1, N=N, B=B, p=p)
    U_u <- U_u_update(u_j=u_j, k2=k2, theta2=theta2, J=J, B=B, p=p)
    alpha <- alpha_update(prec=prec, w_MLE=w_MLE, w=w, v_n=v_n, u_j=u_j, beta=beta, a=a, N=N, J=J, B=B, p=p)
    beta <- beta_update(prec=prec, w_MLE=w_MLE, w=w, v_n=v_n, u_j=u_j, alpha=alpha, b=b, N=N, J=J, B=B, p=p)
    v_n.star <- lapply(1:N, function(n){alpha*v_n[[n]]})
    u_j.star <- lapply(1:J, function(j){beta*u_j[[j]]})
    omega_v <- U_v/alpha^2
    omega_u <- U_u/beta^2
    
    
    # store thinned sims
    if(i %% thin == 0){  
      sim.w[[i/thin+1]] <- w
      sim.u_j[[i/thin+1]] <- u_j
      sim.v_n[[i/thin+1]] <- v_n
      sim.Lambda[[i/thin+1]] <- Lambda
      sim.lambda1sq[[i/thin+1]] <- lambda1sq
      sim.lambda2[[i/thin+1]] <- lambda2
      sim.tausq[[i/thin+1]] <- tausq
      sim.U_v[[i/thin+1]] <- U_v 
      sim.U_u[[i/thin+1]] <- U_u 
      sim.alpha[[i/thin+1]] <- alpha 
      sim.beta[[i/thin+1]] <- beta 
      sim.v_n.star[[i/thin+1]] <- v_n.star
      sim.u_j.star[[i/thin+1]] <- u_j.star
      sim.omega_v[[i/thin+1]] <- omega_v
      sim.omega_u[[i/thin+1]] <- omega_u
    }
    
    # report iter
    cat("i=", i, "\n", sep="")
  }
  
  
  return(list(sim.w=sim.w, 
              sim.u_j=sim.u_j,
              sim.v_n=sim.v_n,
              sim.Lambda=sim.Lambda, 
              sim.lambda1sq=sim.lambda1sq, 
              sim.lambda2=sim.lambda2, 
              sim.tausq=sim.tausq, 
              sim.U_v=sim.U_v,
              sim.U_u=sim.U_u,
              sim.alpha=sim.alpha,
              sim.beta=sim.beta,
              sim.v_n.star=sim.v_n.star,
              sim.u_j.star=sim.u_j.star,
              sim.omega_v=sim.omega_v,
              sim.omega_u=sim.omega_u))   
}








# generate n-by-n positive definite matrix
Posdef <- function (n, ev = runif(n, 0, 10)) {
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}







# (0.1) recover Gamma, xisq
Gamma_recover <- function(tausq, lambda1sq, Lambda, B, p){
  M <- kronecker(diag(1, B*p), chol2inv(chol(Lambda)))
  l <- B^2*p
  xisq <- c(rep(NA, B-1), M[B,B])                        # for xisq
  for(j in (B-1):1) xisq[j] <- M[j,j]- M[j,(j+1):B] %*% solve(M[(j+1):B,(j+1):B]) %*% M[(j+1):B,j]
  xisq <- rep(xisq, B*p)
  m.all.possible <- vector("list", B-1)                  # for m_j
  for(j in (B-1):1) m.all.possible[[j]] <- M[j,(j+1):B] %*% solve(M[(j+1):B,(j+1):B])
  G <- c(rep(NA, l-1), sqrt(tausq[l]*lambda1sq[l]))      # for gamma
  for(j in (l-1):1){
    k1 <- (l-j)%%B
    k2 <- floor((l-j)/B)
    if(k1>0){
      m <- m.all.possible[[B-k1]]
      G[j] <- sqrt(tausq[j]*lambda1sq[j]) - sum(m*G[(j+1):(j+k1)])
    } else {
      G[j] <- sqrt(tausq[j]*lambda1sq[j])
    }
  }
  return(list(G=G, xisq=xisq))
}



# (0.2) get prec: list of N*J
prec_get <- function(HH, Lambda){lapply(HH, function(x){kronecker(x, Lambda)})}




# (1) update w
w_update <- function(prec, w_MLE, v_n, u_j, alpha, beta, tausq, lambda2, N, J){
  D_vec <- 2*tausq/(2*lambda2*tausq+1)
  mu_nj <- lapply(1:(N*J), function(x){w_MLE[[x]] - alpha*v_n[[ceiling(x/J)]] - beta*u_j[[ifelse(x%%J, x%%J, J)]] })  
  Ssq_w <- chol2inv(chol(Reduce("+", prec) + diag(1/D_vec)))
  M_w <- Ssq_w %*% apply(mapply("%*%", prec, mu_nj), 1, sum)
  w <- as.vector(rmvn(n=1, mu=M_w, sigma=chol(Ssq_w), isChol=TRUE,  ncores = 4))
  return(w)
}




# (2) update v_n
v_n_update <- function(prec, w_MLE, w, u_j, alpha, beta, U_v, N, J){        
  v_n <- lapply(1:N, function(n){
    Ssq <- chol2inv(chol(t(t(Reduce("+", prec[((n-1)*J+1):(n*J)])*alpha)*alpha) + diag(U_v)))   
    mu <- lapply(1:J, function(j){w_MLE[[(n-1)*J+j]] - w - beta*u_j[[j]]})
    M <- Ssq %*% (alpha*apply(mapply("%*%", prec[((n-1)*J+1):(n*J)], mu), 1, sum))
    as.vector(rmvn(n=1, mu=M, sigma=chol(Ssq), isChol=TRUE, ncores=4))
  })
  return(v_n)
}




# (3) update u_j
u_j_update <- function(prec, w_MLE, w, v_n, alpha, beta, U_u, N, J){      
  u_j <- lapply(1:J, function(j){
    Ssq <- chol2inv(chol(t(t(Reduce("+", prec[(1:N-1)*J+j])*beta)*beta)+ diag(U_u)))       
    mu <- lapply(1:N, function(n){w_MLE[[(n-1)*J+j]] - w - alpha*v_n[[n]]})
    M <- Ssq %*% (beta*apply(mapply("%*%", prec[(1:N-1)*J+j], mu), 1, sum))
    as.vector(rmvn(n=1, mu=M, sigma=chol(Ssq), isChol=TRUE, ncores=4))
  })
  return(u_j)
}





# (4) update Lambda
Lambda_update <- function(eta0, w, v_n, u_j, alpha, beta, H, G, K_inv, nu, TT, B, p, N, J){
  S <- lapply(1:(J*N), function(x){
    W <- matrix(w + alpha*v_n[[ceiling(x/J)]] + beta*u_j[[ifelse(x%%J, x%%J, J)]], B)
    tcrossprod(eta0[[x]] - W%*%H[[x]])
    })
  S_Lambda <- chol2inv(chol(Reduce('+', S) + K_inv + 2*tcrossprod(matrix(G, nrow=B))))
  Sm <- min(eigen(S_Lambda)$values)
  if(Sm<=0) S_Lambda <- S_Lambda+diag(abs(Sm)+0.1, B)
  Lambda <- rWishart(1, df=(TT-p)*N*J+2*B*p+nu, Sigma=S_Lambda)[,,1]    
  while(sum(eigen(Lambda)$value<=0) | (!isSymmetric(Lambda))){                                    
    Lambda <- rWishart(1, df=(TT-p)*N*J+2*B*p+nu, Sigma=S_Lambda)[,,1]
  }
  return(Lambda)
}





# (5) update lambda1sq for j=1,...,B^2*p
lambda1sq_update <- function(tausq, xisq, B, p, mu1, nu1){
  if(length(mu1)==1) mu1 <- rep(mu1, B^2*p)
  if(length(nu1)==1) nu1 <- rep(nu1, B^2*p)
  nu_lambda1sq <- nu1+2
  mu_lambda1sq <- nu_lambda1sq*xisq*mu1/(2*tausq*mu1+nu1*xisq)
  lambda1sq <- abs(rgamma(n=rep(1, B^2*p), shape=nu_lambda1sq/2, scale=2*mu_lambda1sq/nu_lambda1sq))
  return(lambda1sq)
}





# (6) update lambda2 for j=1,...,B^2*p
lambda2_update <- function(w, B, p, mu2, nu2){
  if(length(mu2)==1) mu2 <- rep(mu2, B^2*p)
  if(length(nu2)==1) nu2 <- rep(nu2, B^2*p)
  nu_lambda2 <- nu2+2
  mu_lambda2 <- nu_lambda2*mu2/(w^2*mu2+nu2)
  lambda2 <- abs(rgamma(n=rep(1, B^2*p), shape=nu_lambda2/2, scale=2*mu_lambda2/nu_lambda2))
  return(lambda2)
}





# (7) update tausq for j=1,...,B^2*p
tausq_update <- function(lambda1sq, w, xisq, B, p){
  M <- sqrt(lambda1sq/(w^2)/xisq)
  S <- lambda1sq/xisq
  temp <- abs(rinvgauss(n=rep(1, B^2*p), mean=M, shape=S))
  tausq <- 1/2/temp
  return(tausq)
}




# (8) update U_v
U_v_update <- function(v_n, k1, theta1, N, B, p){
  d <- Reduce("+", lapply(v_n, function(x){x^2}))
  kv <- rep(N/2+k1, B^2*p)
  thetav <- 1/(d/2+1/theta1)     
  U_v <- abs(rgamma(n=rep(1, B^2*p), shape=kv, scale=thetav))  
  return(U_v)
}



# (9) update U_u
U_u_update <- function(u_j, k2, theta2, J, B, p){
  d <- Reduce("+", lapply(u_j, function(x){x^2}))
  ku <- rep(J/2+k2, B^2*p)
  thetau <- 1/(d/2+1/theta2)  
  U_u <- abs(rgamma(n=rep(1, B^2*p), shape=ku, scale=thetau))
  return(U_u)
}



# (10) update alpha
alpha_update <- function(prec, w_MLE, w, v_n, u_j, beta, a, N, J, B, p){
  Ssq_alpha <- chol2inv(chol(Reduce("+",lapply(1:(N*J), function(x){t(t(prec[[x]]*v_n[[ceiling(x/J)]])*v_n[[ceiling(x/J)]])})) + diag(a, B^2*p)))
  M_alpha <- Ssq_alpha %*%  Reduce("+", lapply(1:(N*J), function(x){(prec[[x]]%*%(w_MLE[[x]]-w-beta*u_j[[ifelse(x%%J, x%%J, J)]]))*v_n[[ceiling(x/J)]]}))
  alpha <- as.vector(rmvn(n=1, mu=M_alpha, sigma=chol(Ssq_alpha), isChol=TRUE, ncores=4))
  return(alpha)
}



# (11) update beta
beta_update <- function(prec, w_MLE, w, v_n, u_j, alpha, b, N, J, B, p){
  Ssq_beta <- chol2inv(chol(Reduce("+",lapply(1:(N*J), function(x){t(t(prec[[x]]*u_j[[ifelse(x%%J, x%%J, J)]])*u_j[[ifelse(x%%J, x%%J, J)]])})) + diag(b, B^2*p)))
  M_beta <- Ssq_beta %*%  Reduce("+", lapply(1:(N*J), function(x){(prec[[x]]%*%(w_MLE[[x]]-w-alpha*v_n[[ceiling(x/J)]]))*u_j[[ifelse(x%%J, x%%J, J)]]}))
  beta <- as.vector(rmvn(n=1, mu=M_beta, sigma=chol(Ssq_beta), isChol=TRUE, ncores=4))
  return(beta)
}





#   
# # (9) update omega_v
# omega_v_update <- function(v_n, k2, theta2, N, B, p){
#   d <- Reduce("+", lapply(v_n, function(x){x^2}))
#   kv <- rep(N/2+k2, B^2*p)
#   thetav <- 1/(d/2+1/theta1)
#   omega_v <- rgamma(n=rep(1, B^2*p), shape=kv, scale=thetav)
#   return(omega_v)
# }
# 
# 
# 
# # (8) update omega_u
# omega_u_update <- function(u_j, k1, theta1, J, B, p){
#   d <- Reduce("+", lapply(u_j, function(x){x^2}))
#   ku <- rep(J/2+k1, B^2*p)
#   thetau <- 1/(d/2+1/theta1)
#   omega_u <- rgamma(n=rep(1, B^2*p), shape=ku, scale=thetau)
#   return(omega_u)
# }





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








# save last iteration for more iterations
save_init <- function(chain1, chain2, chain3, chain4){
  
  n <- length(chain1)
  init1 <- list(w=chain1$sim.w[[n]],
                v_n=chain1$sim.v_n[[n]],
                u_j=chain1$sim.u_j[[n]],
                Lambda=chain1$sim.Lambda[[n]],
                lambda1sq=chain1$sim.lambda1sq[[n]],
                lambda2=chain1$sim.lambda2[[n]],
                tausq=chain1$sim.tausq[[n]],
                U_v=chain1$sim.U_v[[n]],
                U_u=chain1$sim.U_u[[n]],
                alpha=chain1$sim.alpha[[n]],
                beta=chain1$sim.beta[[n]],
                v_n.star=chain1$sim.v_n.star[[n]],
                u_j.star=chain1$sim.u_j.star[[n]],
                omega_v=chain1$sim.omega_v[[n]],
                omega_u=chain1$sim.omega_u[[n]])
  
  init2 <- list(w=chain2$sim.w[[n]],
                v_n=chain2$sim.v_n[[n]],
                u_j=chain2$sim.u_j[[n]],
                Lambda=chain2$sim.Lambda[[n]],
                lambda1sq=chain2$sim.lambda1sq[[n]],
                lambda2=chain2$sim.lambda2[[n]],
                tausq=chain2$sim.tausq[[n]],
                U_v=chain2$sim.U_v[[n]],
                U_u=chain2$sim.U_u[[n]],
                alpha=chain2$sim.alpha[[n]],
                beta=chain2$sim.beta[[n]],
                v_n.star=chain2$sim.v_n.star[[n]],
                u_j.star=chain2$sim.u_j.star[[n]],
                omega_v=chain2$sim.omega_v[[n]],
                omega_u=chain2$sim.omega_u[[n]])
  
  
  init3 <- list(w=chain3$sim.w[[n]],
                v_n=chain3$sim.v_n[[n]],
                u_j=chain3$sim.u_j[[n]],
                Lambda=chain3$sim.Lambda[[n]],
                lambda1sq=chain3$sim.lambda1sq[[n]],
                lambda2=chain3$sim.lambda2[[n]],
                tausq=chain3$sim.tausq[[n]],
                U_v=chain3$sim.U_v[[n]],
                U_u=chain3$sim.U_u[[n]],
                alpha=chain3$sim.alpha[[n]],
                beta=chain3$sim.beta[[n]],
                v_n.star=chain3$sim.v_n.star[[n]],
                u_j.star=chain3$sim.u_j.star[[n]],
                omega_v=chain3$sim.omega_v[[n]],
                omega_u=chain3$sim.omega_u[[n]])
  
  
  init4 <- list(w=chain4$sim.w[[n]],
                v_n=chain4$sim.v_n[[n]],
                u_j=chain4$sim.u_j[[n]],
                Lambda=chain4$sim.Lambda[[n]],
                lambda1sq=chain4$sim.lambda1sq[[n]],
                lambda2=chain4$sim.lambda2[[n]],
                tausq=chain4$sim.tausq[[n]],
                U_v=chain4$sim.U_v[[n]],
                U_u=chain4$sim.U_u[[n]],
                alpha=chain4$sim.alpha[[n]],
                beta=chain4$sim.beta[[n]],
                v_n.star=chain4$sim.v_n.star[[n]],
                u_j.star=chain4$sim.u_j.star[[n]],
                omega_v=chain4$sim.omega_v[[n]],
                omega_u=chain4$sim.omega_u[[n]])
  
  return(list(init1=init1, init2=init2, init3=init3, init4=init4))
  
}



# anaylize results
w_inference <- function(chain1, chain2, chain3, chain4, B, p, n.sim, filename, thin){
  
  l <- floor(n.sim/2/thin)
  sel <- seq(n.sim/2+1, n.sim, thin)
  sw <- array(NA, dim=c(l, 4, B^2*p))
  sw[,1,] <- do.call("rbind", chain1$sim.w)[sel,]
  sw[,2,] <- do.call("rbind", chain2$sim.w)[sel,]
  sw[,3,] <- do.call("rbind", chain3$sim.w)[sel,]
  sw[,4,] <- do.call("rbind", chain4$sim.w)[sel,]
  mon.w <- monitor(sw)
  
  pdf(filename)
  for(i in 1:(B^2*p)){
    plot(sw[,1,i], type="l", ylim=range(c(sw[,1,i], sw[,2,i], sw[,3,i], sw[,4,i])), main=paste("w", i, sep=""), ylab=paste("w", i, sep=""))
    lines(sw[,2,i], col=2)
    lines(sw[,3,i],  col=3)
    lines(sw[,4,i], col=4)
  }
  dev.off()
  return(list(sw=sw, mon.w=mon.w))
}











