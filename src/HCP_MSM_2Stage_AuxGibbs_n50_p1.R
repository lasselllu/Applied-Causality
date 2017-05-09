################################
# Application 2. HCP           #
# Model 3. MSMS                #
# 2-stage + AuxGibbs           #
################################

library(fields)
library(gplots)
library(rstan)

source("C:\\Columbia\\projects\\Michael\\fmri\\code\\MSMS_AuxGibbs.R")
source("C:\\Columbia\\projects\\Michael\\fmri\\code\\MSMS_AnalyzeResults.R")
source("C:\\Columbia\\projects\\Michael\\fmri\\code\\makedata.R")


setwd("C:\\Columbia\\projects\\Michael\\fmri\\application\\HCP\\HCP_n50_p1_2Stage_AuxGibbs")

TT=1200
B=15
p=1
N=50
J=4

dat <- makedata(path="C:\\Columbia\\projects\\Michael\\fmri\\application\\HCP\\3T_HCP820_MSMAll_d15_ts2", TT=TT, B=B, p=p, N=N, J=J, standardize=FALSE)


w_hat_2S <- Reduce("+", dat$w_MLE)/N/J
v_hat_2S <- vector("list", N)
for(n in 1:N) v_hat_2S[[n]] <- Reduce("+", dat$w_MLE[((n-1)*J+1):(n*J)])/J - w_hat_2S
v_hat_2S <- unlist(v_hat_2S)
u_hat_2S <- vector("list", J)
for(j in 1:J) u_hat_2S[[j]] <- Reduce("+", dat$w_MLE[(1:N-1)*J+j])/N - w_hat_2S
u_hat_2S <- unlist(u_hat_2S)


# save estimates
rm(list=ls()[!ls() %in% c("w_hat_2S", "u_hat_2S", "v_hat_2S")])
save.image("HCP_MSMS_n50_p1_2Stage_estimates.RData")



######################################  
# 1.  MSMS Gibbs                     #
######################################

# (1) first find initials, using Naive 2-stage
load("HCP_MSMS_n50_p1_2stage_estimates.RData")
initial = list(w=w_hat_2S, v_n=v_hat_2S, u_j=u_hat_2S, 
               Lambda=Posdef(n=B, ev=runif(B, 0, 0.001)),
               lambda1sq = runif(B^2*p, 0.1, 1000),
               lambda2=runif(B^2*p, 0.1, 10),
               tausq =runif(B^2*p, 0.1, 10),
               U_v = runif(B^2*p, 0.1, 1),
               U_u = runif(B^2*p, 0.1, 1),
               alpha = runif(B^2*p, 0.1, 1),        
               beta = runif(B^2*p, 0.1, 1))
initial$v_n.star <- lapply(1:N, function(n){initial$alpha*initial$v_n[[n]]})
initial$u_j.star <- lapply(1:J, function(j){initial$beta*initial$u_j[[j]]})
initial$omega_v <- initial$U_v/(initial$alpha^2)
initial$omega_u <- initial$U_u/(initial$beta^2)



# (2) run 3 chains
n.iter <- 50000
n.thin <- 25
n.sim <- n.iter/n.thin

sims1 <- MSMSGibbs(eta=dat$eta, H=dat$H, HH=dat$HH, eta0=dat$eta0, w_MLE=dat$w_MLE, TT=TT, B=B, p=p, N=N, J=J, seed=1800, numIter=n.iter, thin=n.thin, initial=initial)
save.image("HCP_MSMS_n50_p1_2stage_Gibbs_sim1.RData")

system.time(sims2 <- MSMSGibbs(eta=dat$eta, H=dat$H, HH=dat$HH, eta0=dat$eta0, w_MLE=dat$w_MLE, TT=TT, B=B, p=p, N=N, J=J, seed=2000, numIter=n.iter, thin=n.thin, initial=initial))
save.image("HCP_MSMS_n50_p1_2stage_Gibbs_sim2.RData")

sims3 <- MSMSGibbs(eta=dat$eta, H=dat$H, HH=dat$HH, eta0=dat$eta0, w_MLE=dat$w_MLE, TT=TT, B=B, p=p, N=N, J=J, seed=4000, numIter=n.iter, thin=n.thin, initial=initial)
save.image("HCP_MSMS_n50_p1_2stage_Gibbs_sim3.RData")





######################################  
# 2.  posterior inference            #
######################################

rm(list=ls())

load("HCP_MSMS_n50_p1_2Stage_Gibbs_sim1.RData")
load("HCP_MSMS_n50_p1_2Stage_Gibbs_sim2.RData")
load("HCP_MSMS_n50_p1_2Stage_Gibbs_sim3.RData")






# (1) w
sw <- w_inference(chains=list(sims1, sims2, sims3), B=B, p=p, n.sim=2001, warmup=1001, thin=1, 
                  filename="HCP_MSMS_n50_p1_2Stage_Gibbs_w_trace_hist.pdf", plot=FALSE, d=3)
sw$mon.w[order(sw$mon.w[,10], decreasing=TRUE),][1:10,]
plot(sw$w.median, sw$w.mode, xlab="median", ylab="mode")   # median and mode are quite the same
w_hat_2S_Gibbs <- sw$w.mode
w.significant <- which(sw$mon.w[,4]*sw$mon.w[,8]>0)   # significantly non-zero
length(w.significant)



# (2) v
sv <- v_inference(chains=list(sims1, sims2, sims3), B=B, p=p, N=N, n.sim=2001, warmup=1001, thin=1, 
                  filename="HCP_MSMS_n50_p1_2Stage_Gibbs_v_trace_hist.pdf", plot=TRUE, d=3)
sv$mon.v[order(sv$mon.v[,10], decreasing=TRUE),][1:10,]
plot(sv$v.median, sv$v.mode, xlab="median", ylab="mode")
v_hat_2S_Gibbs <- sv$v.mode
v.significant <- which(sv$mon.v[,4]*sv$mon.v[,8]>0)   # significantly non-zero
length(v.significant)




# (3) u
su <- u_inference(chains=list(sims1, sims2, sims3), B=B, p=p, J=J, n.sim=2001, warmup=1001, thin=1, 
                  filename="HCP_MSMS_n50_p1_2Stage_Gibbs_u_trace_hist.pdf", plot=TRUE, d=3)
su$mon.u[order(su$mon.u[,10], decreasing=TRUE),][1:10,]
plot(su$u.median, su$u.mode, xlab="median", ylab="mode")
u_hat_2S_Gibbs <- su$u.mode
u.significant <- which(su$mon.u[,4]*su$mon.u[,8]>0)   # significantly non-zero
length(u.significant)






# (4) L
sL <- L_inference(chains=list(sims1, sims2, sims3), B=B, p=p, J=J, n.sim=2001, warmup=1001, thin=1, 
                  filename="HCP_MSMS_n50_p1_2Stage_Gibbs_Lambda_trace_hist.pdf", plot=T, d=6)
sL$mon.L[order(sL$mon.L[,10], decreasing=TRUE),][1:10,]
plot(sL$L.median, sL$L.mode, xlab="median", ylab="mode")
L_hat_2S_Gibbs <- sL$L.mode 




# (5) Omega_v
somega_v <- omega_v_inference(chains=list(sims1, sims2, sims3), B=B, p=p, n.sim=2001, warmup=1001,  thin=1, 
                               filename="HCP_MSMS_n50_p1_2Stage_Gibbs_omega_v_trace_hist.pdf", plot=TRUE)
somega_v$mon.omega_v[order(somega_v$mon.omega_v[,10], decreasing=TRUE),][1:10,]    # all mixed!
omega_v_hat_2S_Gibbs <- somega_v$omega_v.median




# (6) Omega_u
somega_u <- omega_u_inference(chains=list(sims1, sims2, sims3), B=B, p=p, n.sim=2001, warmup=1001,  thin=1, 
                               filename="HCP_MSMS_n50_p1_2Stage_Gibbs_omega_u_trace_hist.pdf", plot=TRUE)
somega_u$mon.omega_u[order(somega_u$mon.omega_u[,10], decreasing=TRUE),][1:10,]    # all mixed!
omega_u_hat_2S_Gibbs <- somega_u$omega_u.median




# (7) lambda1sq
slambda1sq <- lambda1sq_inference(chains=list(sims1, sims2, sims3), B=B, p=p, n.sim=2001, warmup=1001,  thin=1, 
                              filename="HCP_MSMS_n50_p1_2Stage_Gibbs_lambda1sq_trace_hist.pdf", plot=TRUE)
slambda1sq$mon.lambda1sq[order(slambda1sq$mon.lambda1sq[,10], decreasing=TRUE),][1:10,]    # all mixed!
lambda1sq_hat_2S_Gibbs <- slambda1sq$lambda1sq.median






# (8) lambda2
slambda2 <- lambda2_inference(chains=list(sims1, sims2, sims3), B=B, p=p, n.sim=2001, warmup=1001,  thin=1, 
                                  filename="HCP_MSMS_n50_p1_2Stage_Gibbs_lambda2_trace_hist.pdf", plot=TRUE)
slambda2$mon.lambda2[order(slambda2$mon.lambda2[,10], decreasing=TRUE),][1:10,]     # all mixed!
lambda2_hat_2S_Gibbs <- slambda2$lambda2.median





# (9) tausq
stausq <- tausq_inference(chains=list(sims1, sims2, sims3), B=B, p=p, n.sim=2001, warmup=1001,  thin=1, 
                                  filename="HCP_MSMS_n50_p1_2Stage_Gibbs_tausq_trace_hist.pdf", plot=TRUE)
stausq$mon.tausq[order(stausq$mon.tausq[,10], decreasing=TRUE),][1:10,]     # all mixed!
tausq_hat_2S_Gibbs <- stausq$tausq.median






# save estimates
save.image("HCP_MSMS_n50_p1_2Stage_Gibbs_estimates.RData")






################################
# 3. visualization             #
################################


load("HCP_MSMS_n50_p1_2Stage_Gibbs_estimates.RData")
# load("HCP_MSMS_n50_p1_2Stage_estimates.RData")
TT=1200
B=15
p=1
N=50
J=4

ord <- c(1, 3, 4, 8, 13, 9, 10, 11, 7, 5, 6, 14, 15, 2, 12)
reg.num <- c(4, 1, 3, 1, 4, 2)
reg <- c(1, 2, 3, 4, 6, 7)





# (1) w
pdf("HCP_MSMS_n50_p1_2Stage_Gibbs_w_image.pdf", height=5.3, width=15)
par(oma=c(1, 3, 4, 2), mar=c(2, 4, 6, 7), las=1,  mfrow=c(1,3))


d <- 1 # max(abs(w_hat_2S_Gibbs))
image(t(matrix(w_hat_2S_Gibbs, B, B)[ord[B:1],ord]), col="transparent", main="", axes=F)
#axis(side=3, at=seq(0, 1, l=15)[seq(1, 15, 2)], labels=c(1:15)[seq(1, 15, 2)], line=0)
#axis(side=4, at=seq(0, 1, l=15)[seq(1, 15, 2)], labels=c(15:1)[seq(1, 15, 2)], line=0)
axis(side=3, at=(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""), tick=FALSE)
axis(side=2, at=1-(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""),  tick=FALSE)
title(main="population (w)", line=4)
image.plot(t(matrix(w_hat_2S_Gibbs, B, B)[ord[B:1],ord]), add=T,
             col=colorpanel(20, "blue", "white", "red"), 
             breaks=sort(c(seq(-d, d, length.out=20), 0)), xlab="", ylab="", yaxt="n")
polygon(x=c(-0.5/14, 0.5/14+1, 0.5/14+1, -0.5/14), y=c(-0.5/14, -0.5/14, 0.5/14+1, 0.5/14+1))

abline(v=(cumsum(reg.num)+0.5-1)/(15-1))
abline(h=(cumsum(reg.num[6:1])+0.5-1)/(15-1))




# plot off-diagonal elements to compare with v_n and u_j
w_hat_2S_Gibbs_offdiag <- w_hat_2S_Gibbs
w_hat_2S_Gibbs_offdiag[seq(1, 225, 16)] <- 0   # make diagonal elements 0
d <- 0.22
#d <- max(abs(w_hat_2S_Gibbs_offdiag), c(abs(u_hat_2S_Gibbs), abs(u_hat_2S), abs(v_hat_2S_Gibbs), abs(v_hat_2S)))   # off diagonal elements
image(matrix(w_hat_2S_Gibbs_offdiag, B, B)[ord[B:1],ord], col="transparent", main="", axes=F)
axis(side=3, at=(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""), tick=FALSE)
axis(side=2, at=1-(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""),  tick=FALSE)
title(main="off-diagonal", line=4)
image.plot(t(matrix(w_hat_2S_Gibbs_offdiag, B, B)[ord[B:1],ord]), add=T,
           col=colorpanel(20, "blue", "white", "red"), 
           breaks=sort(c(seq(-d, d, length.out=20), 0)))
polygon(x=c(-0.5/14, 0.5/14+1, 0.5/14+1, -0.5/14), y=c(-0.5/14, -0.5/14, 0.5/14+1, 0.5/14+1))
abline(v=(cumsum(reg.num)+0.5-1)/(15-1))
abline(h=(cumsum(reg.num[6:1])+0.5-1)/(15-1))





# plot off-diagonal elements significantly non-zero
w_hat_2S_Gibbs_offdiag_significant <- w_hat_2S_Gibbs_offdiag
w_hat_2S_Gibbs_offdiag_significant[-w.significant] <- 0 
d <- 0.22
#d <- max(abs(w_hat_2S_Gibbs_offdiag), c(abs(u_hat_2S_Gibbs), abs(u_hat_2S), abs(v_hat_2S_Gibbs), abs(v_hat_2S)))
image(t(matrix(w_hat_2S_Gibbs_offdiag_significant, B, B)[ord[B:1],ord]), col="transparent", main="", axes=F)
axis(side=3, at=(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""), tick=FALSE)
axis(side=2, at=1-(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""),  tick=FALSE)
title(main="off-diagonal significantly non-zero", line=4)
image.plot(t(matrix(w_hat_2S_Gibbs_offdiag_significant, B, B)[ord[B:1],ord]), add=T,
           col=colorpanel(20, "blue", "white", "red"), 
           breaks=sort(c(seq(-d, d, length.out=20), 0)))
polygon(x=c(-0.5/14, 0.5/14+1, 0.5/14+1, -0.5/14), y=c(-0.5/14, -0.5/14, 0.5/14+1, 0.5/14+1))
abline(v=(cumsum(reg.num)+0.5-1)/(15-1))
abline(h=(cumsum(reg.num[6:1])+0.5-1)/(15-1))

mtext("N=50, p=1", side=3, outer=T, line=1, font=2, cex=1.3)

dev.off()




# (2) u
pdf("HCP_MSMS_n50_p1_2Stage_Gibbs_u_image.pdf", height=11, width=11.5)
# d <- max(c(abs(u_hat_2S_Gibbs), abs(u_hat_2S), abs(v_hat_2S_Gibbs), abs(v_hat_2S)))
d <- 0.22
par(mfcol=c(2, 2), oma=c(1, 3, 4, 2), mar=c(2, 4, 6, 7), las=1)
for(j in 1:J){
 image(t(matrix(u_hat_2S_Gibbs[(B^2*(j-1)+1):(B^2*j)], B, B)[ord[B:1],ord]), col="transparent", main="", axes=F)
 axis(side=3, at=(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""), tick=FALSE)
 axis(side=2, at=1-(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""),  tick=FALSE)
 title(main=paste("session", j, " (u)", sep=""), line=4)
 image.plot(t(matrix(u_hat_2S_Gibbs[(B^2*(j-1)+1):(B^2*j)], B, B)[ord[B:1],ord]), add=T, 
               col=colorpanel(20, "blue", "white", "red"), 
               breaks=sort(c(seq(-d, d, length.out=20), 0)))
 polygon(x=c(-0.5/14, 0.5/14+1, 0.5/14+1, -0.5/14), y=c(-0.5/14, -0.5/14, 0.5/14+1, 0.5/14+1))
 abline(v=(cumsum(reg.num)+0.5-1)/(15-1))
 abline(h=(cumsum(reg.num[6:1])+0.5-1)/(15-1))
 
 if(j %% 4 == 1) mtext("N=50, p=1", side=3, outer=T, line=1, font=2, cex=1.3)
}

dev.off()
u.significant





# (3) v
pdf("HCP_MSMS_n50_p1_2Stage_Gibbs_v_image.pdf", height=16.2, width=11.2)
# d <- max(c(abs(u_hat_2S_Gibbs), abs(u_hat_2S), abs(v_hat_2S_Gibbs), abs(v_hat_2S)))
d <- 0.22
par(mfcol=c(3, 2), oma=c(1, 1, 4, 2), mar=c(2, 4, 6, 7), las=1)
for(n in 1:N){
  image(t(matrix(v_hat_2S_Gibbs[(B^2*(n-1)+1):(B^2*n)], B, B)[ord[B:1],ord]), col="transparent", main="", axes=F)
  axis(side=3, at=(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""), tick=FALSE)
  axis(side=2, at=1-(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""),  tick=FALSE)
  title(main=paste("subject", n, " (v)", sep=""), line=4)
  image.plot(t(matrix(v_hat_2S_Gibbs[(B^2*(n-1)+1):(B^2*n)], B, B)[ord[B:1],ord]), add=T,
               col=colorpanel(20, "blue", "white", "red"), 
               breaks=sort(c(seq(-d, d, length.out=20), 0)))
  polygon(x=c(-0.5/14, 0.5/14+1, 0.5/14+1, -0.5/14), y=c(-0.5/14, -0.5/14, 0.5/14+1, 0.5/14+1))
  abline(v=(cumsum(reg.num)+0.5-1)/(15-1))
  abline(h=(cumsum(reg.num[6:1])+0.5-1)/(15-1))
  if(n %% 6 == 1)  mtext("N=50, p=1", side=3, outer=T, line=1, font=2, cex=1.3)
}



v_hat_2S_Gibbs_significant <- v_hat_2S_Gibbs
v_hat_2S_Gibbs_significant[-v.significant] <- 0
d <- 0.22
par(mfcol=c(3, 2), oma=c(1, 1, 4, 2), mar=c(2, 4, 6, 7), las=1)
for(n in 1:N){
  
  image(t(matrix(v_hat_2S_Gibbs_significant[(B^2*(n-1)+1):(B^2*n)], B, B)[ord[B:1],ord]), col="transparent", main="", axes=F)
  axis(side=3, at=(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""), tick=FALSE)
  axis(side=2, at=1-(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""),  tick=FALSE)
  title(main=paste("non-zero subject", n, " (v)", sep=""), line=4)
  
  image.plot(t(matrix(v_hat_2S_Gibbs_significant[(B^2*(n-1)+1):(B^2*n)], B, B)[ord[B:1],ord]), add=T,
             col=colorpanel(20, "blue", "white", "red"), 
             breaks=sort(c(seq(-d, d, length.out=20), 0)))
  polygon(x=c(-0.5/14, 0.5/14+1, 0.5/14+1, -0.5/14), y=c(-0.5/14, -0.5/14, 0.5/14+1, 0.5/14+1))
  abline(v=(cumsum(reg.num)+0.5-1)/(15-1))
  abline(h=(cumsum(reg.num[6:1])+0.5-1)/(15-1))
  if(n %% 6 == 1)  mtext("N=50, p=1", side=3, outer=T, line=1, font=2, cex=1.3)
  
}
dev.off()



# (4) Lambda
pdf("HCP_MSMS_n50_p1_2Stage_Gibbs_Lambda_image.pdf", height=5.2, width=5.7)
par(mfcol=c(1, 1), oma=c(1, 1, 2, 2), mar=c(2, 4, 6, 7), las=1)
d <- max(abs(L_hat_2S_Gibbs))
image(t(matrix(L_hat_2S_Gibbs, B, B)[ord[B:1],ord]), col="transparent", main="", axes=F)
axis(side=3, at=(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""), tick=FALSE)
axis(side=2, at=1-(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""),  tick=FALSE)
title(main="Lambda", line=4)
image.plot(t(matrix(L_hat_2S_Gibbs, B, B)[ord[B:1],ord]), add=T,
           col=colorpanel(20, "blue", "white", "red"), 
           breaks=sort(c(seq(-d, d, length.out=20), 0)))
polygon(x=c(-0.5/14, 0.5/14+1, 0.5/14+1, -0.5/14), y=c(-0.5/14, -0.5/14, 0.5/14+1, 0.5/14+1))
abline(v=(cumsum(reg.num)+0.5-1)/(15-1))
abline(h=(cumsum(reg.num[6:1])+0.5-1)/(15-1))
mtext("N=50, p=1", side=3, outer=T, line=0, font=2, cex=1.3)
dev.off()





# (5) omega_v
pdf("HCP_MSMS_n50_p1_2Stage_Gibbs_omega_v_image.pdf", height=5.2, width=5.7)
par(mfcol=c(1, p), oma=c(1, 1, 2, 2), mar=c(2, 4, 6, 7), las=1)
for(pp in 1:p){
  image(t(matrix(sqrt(1/omega_v_hat_2S_Gibbs)[((pp-1)*B^2+1) :(pp*B^2)], B, B)[ord[B:1],ord]), col="transparent", main="", axes=F)
  axis(side=3, at=(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""), tick=FALSE)
  axis(side=2, at=1-(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""),  tick=FALSE)
  title(main="SD across individual", line=4)
  image.plot(t(matrix(sqrt(1/omega_v_hat_2S_Gibbs)[((pp-1)*B^2+1) :(pp*B^2)], B, B)[ord[B:1],ord]), add=T,
             col=colorpanel(20, "white", "red"), 
             breaks=seq(0, 0.3, length.out=21))
  polygon(x=c(-0.5/14, 0.5/14+1, 0.5/14+1, -0.5/14), y=c(-0.5/14, -0.5/14, 0.5/14+1, 0.5/14+1))
  abline(v=(cumsum(reg.num)+0.5-1)/(15-1))
  abline(h=(cumsum(reg.num[6:1])+0.5-1)/(15-1))
}
mtext("N=50, p=1", side=3, outer=T, line=0, font=2, cex=1.3)
dev.off()




# (6) omega_u
pdf("HCP_MSMS_n50_p1_2Stage_Gibbs_omega_u_image.pdf", height=5.2, width=5.7)
par(mfcol=c(1, p), oma=c(1, 1, 2, 2), mar=c(2, 4, 6, 7), las=1)
for(pp in 1:p){
  image(t(matrix(sqrt(1/omega_u_hat_2S_Gibbs)[((pp-1)*B^2+1) :(pp*B^2)], B, B)[ord[B:1],ord]), col="transparent", main="", axes=F)
  axis(side=3, at=(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""), tick=FALSE)
  axis(side=2, at=1-(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""),  tick=FALSE)
  title(main="SD across sessions", line=4)
  image.plot(t(matrix(sqrt(1/omega_u_hat_2S_Gibbs)[((pp-1)*B^2+1) :(pp*B^2)], B, B)[ord[B:1], ord]), add=T,
             col=colorpanel(20, "white", "red"), 
             breaks=seq(0, 0.3, length.out=21))
  polygon(x=c(-0.5/14, 0.5/14+1, 0.5/14+1, -0.5/14), y=c(-0.5/14, -0.5/14, 0.5/14+1, 0.5/14+1))
  abline(v=(cumsum(reg.num)+0.5-1)/(15-1))
  abline(h=(cumsum(reg.num[6:1])+0.5-1)/(15-1))
}
mtext("N=50, p=1", side=3, outer=T, line=0, font=2, cex=1.3)
dev.off()





# (7) lambda1sq
pdf("HCP_MSMS_n50_p1_2Stage_Gibbs_lambda1sq_image.pdf", height=5.2, width=5.7)
par(mfcol=c(1, p), oma=c(1, 1, 2, 2), mar=c(2, 4, 6, 7), las=1)
for(pp in 1:p){
  image(t(matrix(lambda1sq_hat_2S_Gibbs[((pp-1)*B^2+1) :(pp*B^2)], B, B)[ord[B:1],ord]), col="transparent", main="", axes=F)
  axis(side=3, at=(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""), tick=FALSE)
  axis(side=2, at=1-(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""),  tick=FALSE)
  title(main="lambda1sq", line=4)
  
  image.plot(t(matrix(lambda1sq_hat_2S_Gibbs[((pp-1)*B^2+1) :(pp*B^2)], B, B)[ord[B:1],ord]), add=T,
             col=colorpanel(20, "white", "red"), 
             breaks=seq(min(lambda1sq_hat_2S_Gibbs),  max(lambda1sq_hat_2S_Gibbs), length.out=21))
  polygon(x=c(-0.5/14, 0.5/14+1, 0.5/14+1, -0.5/14), y=c(-0.5/14, -0.5/14, 0.5/14+1, 0.5/14+1))
  abline(v=(cumsum(reg.num)+0.5-1)/(15-1))
  abline(h=(cumsum(reg.num[6:1])+0.5-1)/(15-1))
}
mtext("N=50, p=1", side=3, outer=T, line=0, font=2, cex=1.3)
dev.off()





# (8) lambda2
pdf("HCP_MSMS_n50_p1_2Stage_Gibbs_lambda2_image.pdf", height=5.2, width=5.7)
par(mfcol=c(1, p), oma=c(1, 1, 2, 2), mar=c(2, 4, 6, 7), las=1)
for(pp in 1:p){
  image(t(matrix(lambda2_hat_2S_Gibbs[((pp-1)*B^2+1) :(pp*B^2)], B, B)[ord[B:1],ord]), col="transparent", main="", axes=F)
  axis(side=3, at=(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""), tick=FALSE)
  axis(side=2, at=1-(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""),  tick=FALSE)
  title(main="lambda2", line=4)
  
  image.plot(t(matrix(lambda2_hat_2S_Gibbs[((pp-1)*B^2+1) :(pp*B^2)], B, B)[ord[B:1],ord]), add=T,
             col=colorpanel(20, "white", "red"), 
             breaks=seq(min(lambda2_hat_2S_Gibbs),  max(lambda2_hat_2S_Gibbs), length.out=21))
  polygon(x=c(-0.5/14, 0.5/14+1, 0.5/14+1, -0.5/14), y=c(-0.5/14, -0.5/14, 0.5/14+1, 0.5/14+1))
  abline(v=(cumsum(reg.num)+0.5-1)/(15-1))
  abline(h=(cumsum(reg.num[6:1])+0.5-1)/(15-1))
}
mtext("N=50, p=1", side=3, outer=T, line=0, font=2, cex=1.3)
dev.off()








# (9) tausq
pdf("HCP_MSMS_n50_p1_2Stage_Gibbs_tausq_image.pdf",  height=5.2, width=5.7)
par(mfcol=c(1, p), oma=c(1, 1, 2, 2), mar=c(2, 4, 6, 7), las=1)
for(pp in 1:p){
  image(t(matrix(tausq_hat_2S_Gibbs[((pp-1)*B^2+1) :(pp*B^2)], B, B)[ord[B:1],ord]), col="transparent", main="", axes=F)
  axis(side=3, at=(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""), tick=FALSE)
  axis(side=2, at=1-(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""),  tick=FALSE)
  title(main="tausq", line=4)
  image.plot(t(matrix(tausq_hat_2S_Gibbs[((pp-1)*B^2+1) :(pp*B^2)], B, B)[ord[B:1],ord]), add=T,
             col=colorpanel(20, "white", "red"), 
             breaks=seq(min(tausq_hat_2S_Gibbs),  max(tausq_hat_2S_Gibbs), length.out=21), main="tausq", xlab="", ylab="")
  polygon(x=c(-0.5/14, 0.5/14+1, 0.5/14+1, -0.5/14), y=c(-0.5/14, -0.5/14, 0.5/14+1, 0.5/14+1))
  abline(v=(cumsum(reg.num)+0.5-1)/(15-1))
  abline(h=(cumsum(reg.num[6:1])+0.5-1)/(15-1))
}
mtext("N=50, p=1", side=3, outer=T, line=0, font=2, cex=1.3)
dev.off()






# (10) omega_u / (omega_u + omega_v)
bet.var <- 1/omega_u_hat_2S_Gibbs    # between-session variance
with.var <- 1/omega_v_hat_2S_Gibbs    # within-session variance
total.var <- bet.var+with.var         # total variance
pdf("HCP_MSMS_n50_p1_2Stage_Gibbs_test-retest_image.pdf",  height=5.2, width=5.7)
par(mfcol=c(1, p), oma=c(1, 1, 2, 2), mar=c(2, 4, 6, 7), las=1)
for(pp in 1:p){

  image(t(matrix(c(bet.var/total.var)[((pp-1)*B^2+1) :(pp*B^2)], B, B)[ord[B:1],ord]), col="transparent", main="", axes=F)
  axis(side=3, at=(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""), tick=FALSE)
  axis(side=2, at=1-(c(2.5, 5, 7, 9, 11.5, 14.5)-1)/(B-1), labels=paste("R", reg, sep=""),  tick=FALSE)
  title(main="between-sessions/total", line=4)
  
  image.plot(t(matrix(c(bet.var/total.var)[((pp-1)*B^2+1) :(pp*B^2)], B, B)[ord[B:1],ord]), add=T,
             col=colorpanel(20, "blue","white", "red"), 
             breaks=seq(0, 1, length.out=21))
  polygon(x=c(-0.5/14, 0.5/14+1, 0.5/14+1, -0.5/14), y=c(-0.5/14, -0.5/14, 0.5/14+1, 0.5/14+1))
  abline(v=(cumsum(reg.num)+0.5-1)/(15-1))
  abline(h=(cumsum(reg.num[6:1])+0.5-1)/(15-1))
}
mtext("N=50, p=1", side=3, outer=T, line=0, font=2, cex=1.3)
dev.off()





############### ------------------- compare with Naive 2Stage -------------------------- ##################

pdf("HCP_MSMS_n50_p1_2StageGibbs_vs_2Stage.pdf", width=13, height=5)
par(mfrow=c(1, 3))
plot(x=w_hat_2S, y=w_hat_2S_Gibbs)
abline(0, 1, col=2)
plot(x=v_hat_2S, y=v_hat_2S_Gibbs)
abline(0, 1, col=2)
plot(x=u_hat_2S, y=u_hat_2S_Gibbs)
abline(0, 1, col=2)
dev.off()

# (11) w_hat_2S
pdf("HCP_MSMS_n50_p1_2Stage_w_image.pdf", height=5, width=5.3)
d <- max(abs(w_hat_2S))
par(oma=c(1, 1, 1, 2))
image.plot(x=1:15, y=1:15, z=matrix(w_hat_2S, B, B)[B:1, ], 
             col=colorpanel(20, "blue", "white", "red"), 
             breaks=sort(c(seq(-d, d, length.out=20), 0)), main="population_2Stage (w)", xlab="", ylab="")
dev.off()





# (12) u_hat_2S
pdf("HCP_MSMS_n50_p1_2Stage_u_image.pdf", height=10, width=10)
d <- max(c(abs(u_hat_2S_Gibbs), abs(u_hat_2S), abs(v_hat_2S_Gibbs), abs(v_hat_2S)))
par(mfcol=c(2,2), oma=c(1, 1, 1, 2))
for(j in 1:J){
  image.plot(x=1:15, y=1:15, z=matrix(u_hat_2S[(B^2*(j-1)+1):(B^2*j)], B, B)[B:1, ], 
               col=colorpanel(20, "blue", "white", "red"), 
               breaks=sort(c(seq(-d, d, length.out=20), 0)), main=paste("session", j,  "_2Stage (u)", sep=""), xlab="", ylab="")
}
dev.off()



# (13) v
pdf("HCP_MSMS_n50_p1_2Stage_v_image.pdf", height=15, width=10)
d <- max(c(abs(u_hat_2S_Gibbs), abs(u_hat_2S), abs(v_hat_2S_Gibbs), abs(v_hat_2S)))
par(mfcol=c(3, 2), oma=c(1, 1, 1, 2))
for(n in 1:N){
  for(pp in 1:p){
    image.plot(x=1:15, y=1:15, z=matrix(v_hat_2S[(B^2*p*(n-1)+B^2*(pp-1)+1):(B^2*p*(n-1)+ B^2*pp)], B, B)[B:1, ], 
               col=colorpanel(20, "blue", "white", "red"), 
               breaks=sort(c(seq(-d, d, length.out=20), 0)), main=paste("subject", n,  "_2Stage (v)", sep=""), xlab="", ylab="")
  }
  
}
dev.off()

