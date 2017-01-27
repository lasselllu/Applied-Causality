#################################
# Fit SVAR by VB                #
# Model 3                       #  
# multi subject multi session   #
#################################






library(statmod)



MSMS_VB <- function(H, HH, eta0, w_MLE, 
                    N, J, B, TT, p, 
                    initial=NULL, 
                    nu=1, mu1=1, 
                    nu1=0.001, mu2=1, 
                    nu2=0.01, 
                    k1 = 200, theta1 = 0.005, 
                    k2 = 200, theta2 = 0.005,  
                    K_inv=diag(1/(B-1), B),
                    a = 1,
                    b = 1,
                    seed=2000, 
                    numIter=100,
                    thin=1){
  
  q <- B^2*p
  k_lambda1sq <- nu1/2+1
  k_lambda2 <- nu2/2+1/2
  nu_Lambda <- (TT-p)*N*J + nu + 2*B*p
  
  
  
  # store values
  sim.M_w <- sim.Ssq_w <- vector(mode="list", length=numIter/thin+1)
  sim.M_v <- sim.Ssq_v <- vector(mode="list", length=numIter/thin+1)
  sim.M_u <- sim.Ssq_u <- vector(mode="list", length=numIter/thin+1)
  sim.M_U_v <- sim.M_U_u <- sim.M_alpha <- sim.M_beta <- vector(mode="list", length=numIter/thin+1)
  sim.Ssq_alpha <- sim.Ssq_beta <- vector(mode="list", length=numIter/thin+1)
  sim.S_Lambda <- sim.mu_lambda1sq <- sim.mu_lambda2 <- vector(mode="list", length=numIter/thin+1)
  sim.mu_2tausqinv <- sim.lambda_2tausqinv <- vector(mode="list", length=numIter/thin+1)
  sim.M_v_star <- sim.M_u_star <- sim.omega_u <- sim.omega_v  <- vector(mode="list", length=numIter/thin+1)
  
  
  set.seed(seed)
  if(!is.null(initial)){
    sim.Ssq_w[[1]] <- Ssq_w <- initial$Ssq_w
    sim.M_w[[1]] <- M_w <- initial$M_w
    sim.Ssq_v[[1]] <- Ssq_v <- initial$Ssq_v      # list of N
    sim.M_v[[1]] <- M_v <- initial$M_v          # list of N
    sim.Ssq_u[[1]] <- Ssq_u <- initial$Ssq_u      # list of J
    sim.M_u[[1]] <- M_u <- initial$M_u          # list of J
    sim.M_U_v[[1]] <- M_U_v <- initial$M_U_v
    sim.M_U_u[[1]] <- M_U_u <- initial$M_U_u
    sim.M_alpha[[1]] <- M_alpha <- initial$M_alpha
    sim.M_beta[[1]] <- M_beta <- initial$M_beta
    sim.Ssq_alpha[[1]] <- Ssq_alpha <- initial$Ssq_alpha
    sim.Ssq_beta[[1]] <- Ssq_beta <- initial$Ssq_beta
    sim.S_Lambda[[1]] <- S_Lambda <- initial$S_Lambda
    sim.mu_lambda1sq[[1]] <- mu_lambda1sq <- initial$mu_lambda1sq
    sim.mu_lambda2[[1]] <- mu_lambda2 <- initial$mu_lambda2
    sim.mu_2tausqinv[[1]] <- mu_2tausqinv <- initial$mu_2tausqinv
    sim.lambda_2tausqinv[[1]] <- lambda_2tausqinv <- initial$lambda_2tausqinv
    
  } else {
    sim.Ssq_w[[1]] <- Ssq_w <- crossprod(matrix(rnorm(q^2), q, q)+diag(1, q))
    sim.M_w[[1]] <- M_w <- rnorm(q)
    sim.Ssq_v[[1]] <- Ssq_v <- lapply(vector("list", N), function(x){crossprod(matrix(rnorm(q^2, 0, 1), q, q))+diag(1, q)})
    sim.M_v[[1]] <- M_v <- lapply(vector("list", N), function(x){rnorm(q, 0, 1)})
    sim.Ssq_u[[1]] <- Ssq_u <- lapply(vector("list", J), function(x){crossprod(matrix(rnorm(q^2, 0, 1), q, q))+diag(1, q)})
    sim.M_u[[1]] <- M_u <- lapply(vector("list", J), function(x){rnorm(q, 0, 1)})
    sim.M_U_v[[1]] <- M_U_v <-  runif(q, 0.1, 1)
    sim.M_U_u[[1]] <- M_U_u <-  runif(q, 0.1, 1)
    sim.M_alpha[[1]] <- M_alpha <- runif(q, 0.1, 1)
    sim.M_beta[[1]] <- M_beta <- runif(q, 0.1, 1)
    sim.Ssq_alpha[[1]] <- Ssq_alpha <- crossprod(matrix(rnorm(q^2), q, q)+diag(1, q))
    sim.Ssq_beta[[1]] <- Ssq_beta <- crossprod(matrix(rnorm(q^2), q, q)+diag(1, q))
    sim.S_Lambda[[1]] <- S_Lambda <- crossprod(matrix(rnorm(B^2), B, B)+diag(1, B))
    sim.mu_lambda1sq[[1]] <- mu_lambda1sq <- runif(q, 0.1, 1)
    sim.mu_lambda2[[1]] <- mu_lambda2 <- runif(q, 0.1, 1)
    sim.mu_2tausqinv[[1]] <- mu_2tausqinv <- runif(q, 0.1, 1)
    sim.lambda_2tausqinv[[1]] <- lambda_2tausqinv <- runif(q, 0.1, 1)
  }
  
  sim.M_v_star[[1]] <- lapply(1:N, function(n){M_alpha*M_v[[n]]})
  sim.M_u_star[[1]] <- lapply(1:J, function(j){M_beta*M_u[[j]]})
  sim.omega_u[[1]] <- M_U_v / (M_alpha^2)
  sim.omega_v[[1]] <- M_U_u / (M_beta^2)  
  
  
  
  # update parameters
  for(i in 1:numIter){
    
    #temp <- Gamma_recover(tausq=0.5/mu_2tausqinv, lambda1sq=mu_lambda1sq, Lambda=nu*S_Lambda, B, p)
    #G <- temp$G
    #xisq <- temp$xisq
    temp <- EQQ_Exisq_recover(mu_2tausqinv, lambda_2tausqinv, mu_lambda1sq, k_lambda1sq, S_Lambda, nu, N, J, B, TT, p, q)
    EQQ <- temp$EQQ
    Exisq <- temp$Exisq
    
    prec <- lapply(HH, function(x){kronecker(x, nu_Lambda*S_Lambda)})
    
    Ssq_w <- chol2inv(chol(Reduce("+", prec) + diag(mu_lambda2+mu_2tausqinv)))
    mu_nj <- lapply(1:(N*J), function(x){w_MLE[[x]] - M_alpha * M_v[[ceiling(x/J)]] - M_beta * M_u[[ifelse(x%%J, x%%J, J)]] })  
    M_w <- c(Ssq_w %*% apply(mapply("%*%", prec, mu_nj), 1, sum))

    Ssq_v <- lapply(1:N, function(n){chol2inv(chol(t(t(Reduce("+", prec[((n-1)*J+1):(n*J)])*M_alpha)*M_alpha) + diag(M_U_v)))}) 
    M_v <- lapply(1:N, function(n){c(Ssq_v[[n]] %*% (M_alpha*apply(mapply("%*%", prec[((n-1)*J+1):(n*J)], lapply(1:J, function(j){w_MLE[[(n-1)*J+j]] - M_w - M_beta*M_u[[j]]})), 1, sum)))})

    Ssq_u <- lapply(1:J, function(j){chol2inv(chol(t(t(Reduce("+", prec[(1:N-1)*J+j])*M_beta)*M_beta) + diag(M_U_u)))})
    M_u <- lapply(1:J, function(j){c(Ssq_u[[j]] %*% (M_beta*apply(mapply("%*%", prec[(1:N-1)*J+j], lapply(1:N, function(n){w_MLE[[(n-1)*J+j]] - M_w - M_alpha*M_v[[n]]})), 1, sum)))})
    
    Ssq_alpha <- chol2inv(chol(Reduce("+",lapply(1:(N*J), function(x){t(t(prec[[x]]*M_v[[ceiling(x/J)]])*M_v[[ceiling(x/J)]])})) + diag(a, q)))
    M_alpha <- c(Ssq_alpha %*%  Reduce("+", lapply(1:(N*J), function(x){(prec[[x]]%*%(w_MLE[[x]]-M_w-M_beta*M_u[[ifelse(x%%J, x%%J, J)]]))*M_v[[ceiling(x/J)]]})))
    
    Ssq_beta <- chol2inv(chol(Reduce("+",lapply(1:(N*J), function(x){t(t(prec[[x]]*M_u[[ifelse(x%%J, x%%J, J)]])*M_u[[ifelse(x%%J, x%%J, J)]])})) + diag(b, q)))
    M_beta <- c(Ssq_beta %*%  Reduce("+", lapply(1:(N*J), function(x){(prec[[x]]%*%(w_MLE[[x]]-M_w-M_alpha*M_v[[ceiling(x/J)]]))*M_u[[ifelse(x%%J, x%%J, J)]]})))
    
    M_U_v <- (N/2+k1)/(0.5*Reduce("+", lapply(1:N, function(n){diag(Ssq_v[[n]])+M_v[[n]]^2}))+1/theta1)
    M_U_u <- (J/2+k2)/(0.5*Reduce("+", lapply(1:J, function(j){diag(Ssq_u[[j]])+M_u[[j]]^2}))+1/theta2)
    
    ES_nj <- ES_nj_recover(N, J, B, TT, p, M_w, M_v, M_u, M_alpha, M_beta, eta0, H, Ssq_w, Ssq_v, Ssq_u, Ssq_alpha, Ssq_beta)
    #S_Lambda <- chol2inv(chol(Reduce("+", ES_nj) + K_inv + 2*tcrossprod(matrix(G, nrow=B))))
    #mu_lambda1sq <- (nu1+2)*mu1/(nu1+mu1/xisq/mu_2tausqinv)
    #mu_2tausqinv <- sqrt(mu_lambda1sq/(diag(Ssq_w)+M_w^2)/xisq)
    S_Lambda <- chol2inv(chol(Reduce("+", ES_nj) + K_inv + 2*EQQ))
    mu_lambda1sq <- (nu1+2)*mu1/(nu1+mu1/Exisq/mu_2tausqinv)
    mu_lambda2 <- c((nu2+1)*mu2/(nu2+mu2*(diag(Ssq_w)+M_w^2)))
    mu_2tausqinv <- sqrt(mu_lambda1sq/(diag(Ssq_w)+M_w^2)/Exisq)
    lambda_2tausqinv <- mu_lambda1sq/Exisq

    

    
    
    
    # store  sims
    if(i %% thin == 0){  
      sim.Ssq_w[[i/thin+1]] <- Ssq_w 
      sim.M_w[[i/thin+1]] <- M_w 
      sim.Ssq_v[[i/thin+1]] <- Ssq_v 
      sim.M_v[[i/thin+1]] <- M_v 
      sim.Ssq_u[[i/thin+1]] <- Ssq_u 
      sim.M_u[[i/thin+1]] <- M_u 
      sim.M_U_v[[i/thin+1]] <- M_U_v 
      sim.M_U_u[[i/thin+1]] <- M_U_u 
      sim.M_alpha[[i/thin+1]] <- M_alpha 
      sim.M_beta[[i/thin+1]] <- M_beta 
      sim.Ssq_alpha[[i/thin+1]] <- Ssq_alpha 
      sim.Ssq_beta[[i/thin+1]] <- Ssq_beta 
      sim.S_Lambda[[i/thin+1]] <- S_Lambda 
      sim.mu_lambda1sq[[i/thin+1]] <- mu_lambda1sq 
      sim.mu_lambda2[[i/thin+1]] <- mu_lambda2 
      sim.mu_2tausqinv[[i/thin+1]] <- mu_2tausqinv 
      sim.lambda_2tausqinv[[i/thin+1]] <- lambda_2tausqinv 
      sim.M_v_star[[i/thin+1]] <- lapply(1:N, function(n){M_alpha*M_v[[n]]})
      sim.M_u_star[[i/thin+1]] <- lapply(1:J, function(j){M_beta*M_u[[j]]})
      sim.omega_u[[i/thin+1]] <- M_U_v / (M_alpha^2)
      sim.omega_v[[i/thin+1]] <- M_U_u / (M_beta^2)  
    }
    
    # report iter
    cat("i=", i, "\n", sep="")
   
    
  }
  
  
  return(list(
    sim.Ssq_w = sim.Ssq_w,
    sim.M_w = sim.M_w,
    sim.Ssq_v =sim.Ssq_v,
    sim.M_v = sim.M_v,
    sim.Ssq_u = sim.Ssq_u,
    sim.M_u = sim.M_u,
    sim.M_U_v=sim.M_U_v,
    sim.M_U_u=sim.M_U_u,
    sim.M_alpha=sim.M_alpha,
    sim.M_beta=sim.M_beta,
    sim.Ssq_alpha=sim.Ssq_alpha,
    sim.Ssq_beta=sim.Ssq_beta,
    sim.S_Lambda = sim.S_Lambda,
    sim.mu_lambda1sq = sim.mu_lambda1sq,
    sim.mu_lambda2 = sim.mu_lambda2,
    sim.mu_2tausqinv = sim.mu_2tausqinv,
    sim.lambda_2tausqinv = sim.lambda_2tausqinv,
    sim.M_v_star=sim.M_v_star,
    sim.M_u_star=sim.M_u_star,
    sim.omega_u=sim.omega_u,
    sim.omega_v=sim.omega_v,
    
    res = list(                    # final value
      Ssq_w = Ssq_w,
      M_w  = M_w, 
      Ssq_v = Ssq_v, 
      M_v = M_v, 
      Ssq_u = Ssq_u, 
      M_u = M_u, 
      M_U_v = M_U_v, 
      M_U_u = M_U_u, 
      M_alpha = M_alpha, 
      M_beta = M_beta, 
      Ssq_alpha = Ssq_alpha, 
      Ssq_beta = Ssq_beta, 
      S_Lambda = S_Lambda, 
      mu_lambda1sq  = mu_lambda1sq, 
      mu_lambda2 = mu_lambda2, 
      mu_2tausqinv = mu_2tausqinv, 
      lambda_2tausqinv = lambda_2tausqinv, 
      M_v_star = lapply(1:N, function(n){M_alpha*M_v[[n]]}),
      M_u_star = lapply(1:J, function(j){M_beta*M_u[[j]]}),
      omega_u = M_U_v / (M_alpha^2),
      omega_v = M_U_u / (M_beta^2)  
    )
    
    ))
  
}











ES_nj_recover <- function(N, J, B, TT, p, M_w, M_v, M_u, M_alpha, M_beta, eta0, H, Ssq_w, Ssq_v, Ssq_u, Ssq_alpha, Ssq_beta){
  lapply(1:(J*N), function(x){
    n <- ceiling(x/J)
    j <- ifelse(x%%J, x%%J, J)
    EW <- matrix(M_w + M_alpha* M_v[[n]] + M_beta*M_u[[j]], B)
    VW <- Ssq_w + (Ssq_alpha + tcrossprod(M_alpha)) * (Ssq_v[[n]] + tcrossprod(M_v[[n]])) - tcrossprod(M_alpha*M_v[[n]])
                + (Ssq_beta + tcrossprod(M_beta)) * (Ssq_u[[j]] + tcrossprod(M_u[[j]])) - tcrossprod(M_beta*M_u[[j]]) 
    # if(min(eigen(VW)$value) <= 0) VW <- VW + diag(abs(min(eigen(VW)$value))+0.001, B^2*p)
    S <- tcrossprod(eta0[[x]] - EW%*%H[[x]])  
    + Reduce("+", mapply(function(k, l){sum(H[[x]][k,]*H[[x]][l,])*VW[((k-1)*B+1):(k*B), ((l-1)*B+1):(l*B)]}, rep(1:(B*p), B*p), sort(rep(1:(B*p), B*p)), SIMPLIFY = FALSE))
    S
  })
}





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



EQQ_Exisq_recover <- function(mu_2tausqinv, lambda_2tausqinv, mu_lambda1sq, k_lambda1sq, S_Lambda, nu, N, J, B, TT, p, q){
  temp_tausq <- matrix(1/2/abs(rinvgauss(n=q*100, mean=mu_2tausqinv, shape=lambda_2tausqinv)), nrow=q)       # q-100
  temp_Lambda <- rWishart(100, df=(TT-p)*N*J+2*B*p+nu, Sigma=S_Lambda)                                       # B-B-100 array
  temp_lambda1sq <- matrix(abs(rgamma(n=q*100, shape=k_lambda1sq, scale=mu_lambda1sq/k_lambda1sq)), nrow=q)  # q-100
  G_mat <- xisq_mat <- matrix(NA, q, 100)
  for(i in 1:100){
    temp <- Gamma_recover(temp_tausq[,i], temp_lambda1sq[,i], temp_Lambda[,,i], B, p)
    G_mat[,i] <- temp$G
    xisq_mat[,i] <- temp$xisq
  }
  Ssq_G <- cov(t(G_mat))
  EQQ <- tcrossprod(matrix(rowMeans(G_mat), B)) + Reduce("+", lapply(1:(B*p), function(i){Ssq_G[((i-1)*B+1):(i*B), ((i-1)*B+1):(i*B)]}))
  return(list(EQQ=EQQ, Exisq=rowMeans(xisq_mat)))
}






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



