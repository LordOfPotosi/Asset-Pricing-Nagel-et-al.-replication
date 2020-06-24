rm(list = ls())

library('dplyr')
library('lmtest')
library('MASS')
library('foreach')
library('parallel')
library('gmm')
library('xtable')

setwd("C:/Users/mayur/Dropbox/AP Replication Project/Replication Project Code files")

fac.dat <- read.csv('FF3q.csv', header = TRUE)
port.dat <- read.csv('25q_portfolios.csv', header = TRUE)


T <- nrow(port.dat)
N <- ncol(port.dat) - 1
K <- ncol(fac.dat) -2 ## Exlude RF 
n.sim <- 1000

fin.dat <- merge(port.dat, fac.dat, by = 'Date')

retex.dat <- cbind((fin.dat[,2:ncol(port.dat)] - fin.dat$RF), fin.dat[,(ncol(port.dat)+1):ncol(fin.dat)]) ### Get excess returns 
fac <- as.matrix(fac.dat[,2:(K+1)]) 

### GRS-F statistic Testing CAPM ###                            
P1 <- as.matrix(fac[,1:K])   ### I am using the FF 3 factors Mkt-RF factor 

### Calcualting the maximal Sharpe ratio: Write the general code for K factors. Here K = 1 ###
P1.mean <- matrix(colMeans(P1),K,1) 
P1.mean.mat <- diag(as.numeric(P1.mean), K, K)
P1.demeaned <- P1 - matrix(1,T,K) %*% P1.mean.mat

P1.var <- (t(P1.demeaned) %*% P1.demeaned)/T
P1.SR <- t(P1.mean) %*% solve(P1.var) %*% P1.mean

#### First stage regressions to get C : R = a + P C + e ####

P <- cbind(matrix(1,T,1),P1)   

R <- as.matrix(retex.dat[,1:N])

C1 <- matrix(solve(t(P) %*% P) %*% t(P) %*% R,(K+1),N) 
C <- matrix(C1[2:(K+1),],K,N)


e <- R - P%*% C1
Sigma.e <- (t(e) %*% e)/T 



#### Table 1 Results: 2nd stage regressions  ####

### OLS Results ####
mu.R <- as.matrix(colMeans(R))

C.t <- cbind(matrix(1,N,1),t(C))

lambda.capm <- solve(t(C.t) %*% C.t) %*% t(C.t) %*% mu.R

M.c <- diag(1,N,N) - C.t %*% solve(t(C.t) %*% C.t) %*% t(C.t)
capm.Rsq <- 1 - t(mu.R) %*% M.c %*% mu.R/(t(mu.R - mean(mu.R)) %*% (mu.R - mean(mu.R)))*N/(N-K)

Sigma.f <- cbind(matrix(0,K+1,1),rbind(matrix(0,1,K),P1.var)) ### To incorporate the constant 

shanken <- (1 + t(matrix(lambda.capm)) %*% ginv(Sigma.f) %*% matrix(lambda.capm))

Sigma.lambda <- ((ginv( t(C.t) %*% C.t) %*% t(C.t)%*% Sigma.e%*% C.t %*% ginv(t(C.t)%*% C.t))*shanken[1] + Sigma.f )/T

sd.lambda <- sqrt(diag(Sigma.lambda))

lambda.tstat <- lambda.capm/sd.lambda


### Tsq and q stat calculations ###

Sigma.a <- (M.c %*% Sigma.e %*% M.c)/T * shanken[1]
a.hat <- mu.R - C.t %*% lambda.capm 
capm.tsq <- t(a.hat) %*% ginv(Sigma.a) %*% a.hat 

capm.q <- t(a.hat) %*% ginv(M.c %*% Sigma.e %*% M.c) %*% a.hat 


capm.tsq.pval <- 1 - pchisq(capm.tsq, df = (N-K-1), ncp = capm.q )



q.min <- 0
q.max <- 2
incr <- 0.05
q.grid <- as.matrix(seq(q.min, q.max, incr))

tsq.quantiles <- matrix(NA,length(q.grid),3)

for (q.ind in 1:length(q.grid)){
  q.ncp <- q.grid[q.ind]*T/shanken[1]
  tsq.quantiles[q.ind,] <- qchisq(c(0.05,0.5,0.95),df = (N-K-1), ncp =q.ncp)
}

plot(q.grid, tsq.quantiles[,1], ylab = expression(T^2 ~ 'statistic'), xlab = expression(q), main = expression('FF 25 (Conf.Int. for'~T^2~ 'and q)'), type = 'l', xlim = c(0,max(q.grid)), ylim = c(min(tsq.quantiles),max(tsq.quantiles)))
lines(q.grid, tsq.quantiles[,2], col = 'red')
lines(q.grid, tsq.quantiles[,3], col = 'blue')
lines(q.grid, matrix(capm.tsq,length(q.grid),1),lwd = 2, col = 'green', lty = 2)
legend('topleft', c('5%-lie','Median','95%-ile','Tsq CAPM'), lty = c(1,1,1,2), col = c('black','red','blue','green'),cex = 0.7, bg = 'transparent')
locator()       ### You will have to select the points of intersection of the confidence interval lines with capm.tsq, 
### Then hit finish to retrieve the points. I have reported them below, 

q.conf.int <- matrix(c(0.06,0.31),2,1)
q.out <- matrix(c(round(capm.q,2), paste('(', as.character(q.conf.int[1]), ',', as.character(q.conf.int[2]) ,')')),2,1)
capm.tsq.out <- c(round(capm.tsq,2), round(capm.tsq.pval,2))

### GLS Rsq ###
# 
# lambda.gls.capm <- solve(t(C.t) %*% ginv(Sigma.e) %*% C.t) %*% t(C.t) %*% ginv(Sigma.e) %*% mu.R
# e.gls <- mu.R - C.t %*% lambda.gls.capm
# capm.gls.Rsq <- 1 - var(e.gls)/var(mu.R)*N/(N-K)

R.eig <- eigen(var(R))
R.Q <- R.eig$vectors
R.lambda <- R.eig$values

R.chol <- R.Q %*% diag(sqrt(R.lambda),N,N) %*% solve(R.Q)

C.gls.t <- solve(t(R.chol)) %*% C.t
mu.gls.R <- solve(t(R.chol)) %*% mu.R

M.gls.c <- diag(1,N,N) - C.gls.t %*% solve(t(C.gls.t) %*% C.gls.t) %*% t(C.gls.t)
capm.gls.Rsq <- 1 - (t(mu.gls.R) %*% M.gls.c %*% mu.gls.R/(t(mu.gls.R - mean(mu.gls.R)) %*% (mu.gls.R - mean(mu.gls.R))))*N/(N-K-1)


capm.out <- cbind(matrix(round(lambda.capm,2),(K+1),1),matrix(round(lambda.tstat,2),(K+1),1),matrix(c(round(capm.Rsq,2),'-','-','-'),(K+1),1),matrix(c(round(capm.gls.Rsq,2),'-','-','-'),(K+1),1),c(capm.tsq.out,'-','-'), c(q.out,'-','-'))

colnames(capm.out) <- c('Coefficients', 't-stat','Rsq OLS', 'Rsq GLS','Tsq Stat (p-val)', 'q (confidence interval)')
rownames(capm.out) <- c('Const.', 'Rm','SMB','HML')

xtable(capm.out, type = "latex", file = "tab1_ff3_a.tex")

#######################################
#### FF25 + 30 Industry Portfolios #### 
#######################################
ind.dat <- read.csv('30ind_q.csv', header = TRUE)
N.ind <- ncol(port.dat) + ncol(ind.dat) - 2
test.dat <- merge(port.dat, ind.dat,by = 'Date')
fin.ind.dat <- merge(test.dat, fac.dat, by = 'Date')
retex.ind.dat <- cbind((fin.ind.dat[,2:ncol(port.dat)] - fin.ind.dat$RF), (fin.ind.dat[,(ncol(port.dat)+1):(ncol(port.dat) + ncol(ind.dat))] - fin.ind.dat$RF), fin.ind.dat[,(ncol(port.dat) +ncol(ind.dat)+1):ncol(fin.ind.dat)]) ### Get excess returns 

R.ind <- as.matrix(retex.ind.dat[,1:N.ind])
mu.ind.R <- as.matrix(colMeans(R.ind))

C1.ind <-  matrix(solve(t(P) %*% P) %*% t(P) %*% R.ind,(K+1),N.ind) 
C.ind <- matrix(C1.ind[2:(K+1),],K,N.ind)
C.ind.t <-  cbind(matrix(1,N.ind,1),t(C.ind))

e.ind <- R.ind - P %*% C1.ind
Sigma.ind.e <- (t(e.ind) %*% e.ind)/T 


lambda.ind.capm <- solve(t(C.ind.t) %*% C.ind.t) %*% t(C.ind.t) %*% mu.ind.R
shanken.ind <- (1 + t(matrix(lambda.ind.capm)) %*% ginv(Sigma.f) %*% matrix(lambda.ind.capm))

Sigma.ind.lambda <- ((ginv( t(C.ind.t) %*% C.ind.t) %*% t(C.ind.t)%*% Sigma.ind.e%*% C.ind.t %*% ginv(t(C.ind.t)%*% C.ind.t))*shanken.ind[1] + Sigma.f )/T

sd.ind.lambda <- sqrt(diag(Sigma.ind.lambda))

lambda.ind.tstat <- lambda.ind.capm/sd.ind.lambda

e.ind.cs <- mu.ind.R - C.ind.t %*% lambda.ind.capm
capm.ind.Rsq <- 1 - var(e.ind.cs)/var(mu.ind.R)*(N)/(N-K)


### GLS Rsq Industry portfolios ###

R.ind.eig <- eigen(var(R.ind))
R.ind.Q <- R.ind.eig$vectors
R.ind.lambda <- R.ind.eig$values

R.chol <- R.ind.Q %*% diag(sqrt(R.ind.lambda),N.ind,N.ind) %*% solve(R.ind.Q)

C.ind.gls.t <- solve(t(R.chol)) %*% C.ind.t
mu.ind.gls.R <- solve(t(R.chol)) %*% mu.ind.R

M.ind.gls.c <- diag(1,N.ind,N.ind) - C.ind.gls.t %*% solve(t(C.ind.gls.t) %*% C.ind.gls.t) %*% t(C.ind.gls.t)
capm.ind.gls.Rsq <- 1 - (t(mu.ind.gls.R) %*% M.ind.gls.c %*% mu.ind.gls.R/(t(mu.ind.gls.R - mean(mu.ind.gls.R)) %*% (mu.ind.gls.R - mean(mu.ind.gls.R))))*N.ind/(N.ind-K-1)



### FF25 + Ind30: Tsq and q stat calculations ###
M.ind.c <- diag(1,N.ind,N.ind) - C.ind.t %*% solve(t(C.ind.t) %*% C.ind.t) %*% t(C.ind.t)

Sigma.ind.a <- (M.ind.c %*% Sigma.ind.e %*% M.ind.c)/T * shanken[1]
a.ind.hat <- mu.ind.R - C.ind.t %*% lambda.ind.capm 
capm.ind.tsq <- t(a.ind.hat) %*% ginv(Sigma.ind.a) %*% a.ind.hat 

capm.ind.q <- t(a.ind.hat) %*% ginv(M.ind.c %*% Sigma.ind.e %*% M.ind.c) %*% a.ind.hat 

capm.ind.tsq.pval <- 1 - pchisq(capm.ind.tsq, df = (N.ind-K-1), ncp = capm.q )



q.min <- 0
q.max <- 2
incr <- 0.05
q.grid <- as.matrix(seq(q.min, q.max, incr))

tsq.ind.quantiles <- matrix(NA,length(q.grid),3)

for (q.ind in 1:length(q.grid)){
  q.ncp <- q.grid[q.ind]*T/shanken.ind[1]
  tsq.ind.quantiles[q.ind,] <- qchisq(c(0.05,0.5,0.95),df = (N.ind-K-1), ncp =q.ncp)
}

plot(q.grid, tsq.ind.quantiles[,1], ylab = expression(T^2 ~ 'statistic'), xlab = expression(q), main = expression('FF 25 + 30 Ind. (Conf.Int. for'~T^2~ 'and q)'), type = 'l', xlim = c(0,max(q.grid)), ylim = c(min(tsq.quantiles),max(tsq.quantiles)))
lines(q.grid, tsq.ind.quantiles[,2], col = 'red')
lines(q.grid, tsq.ind.quantiles[,3], col = 'blue')
lines(q.grid, matrix(capm.ind.tsq,length(q.grid),1),lwd = 2, col = 'green', lty = 2)
legend('topleft', c('5%-lie','Median','95%-ile','Tsq CAPM'), lty = c(1,1,1,2), col = c('black','red','blue','green'),cex = 0.7, bg = 'transparent')
locator()       ### You will have to select the points of intersection of the confidence interval lines with capm.tsq, 
### Then hit finish to retrieve the points. I have reported them below, 

q.ind.conf.int <- matrix(c(0.49,0.98),2,1)
q.ind.out <- matrix(c(round(capm.ind.q,2), paste('(', as.character(q.ind.conf.int[1]), ',', as.character(q.ind.conf.int[2]) ,')')),2,1)
capm.ind.tsq.out <- c(round(capm.ind.tsq,2), round(capm.ind.tsq.pval,2))




capm.ind.out <- cbind(matrix(round(lambda.ind.capm,2),(K+1),1),matrix(round(lambda.ind.tstat,2),(K+1),1),matrix(c(round(capm.ind.Rsq,2),'-','-','-'),(K+1),1),matrix(c(round(capm.ind.gls.Rsq,2),'-','-','-'),(K+1),1),matrix(c(capm.ind.tsq.out,'-','-'),(K+1),1), matrix(c(q.ind.out,'-','-'),(K+1),1))
colnames(capm.ind.out) <- c('Coefficients ', 't-stat','R^2 OLS', '#Rsq GLS','Tsq Stat (p-val)', 'q (confidence interval)')
rownames(capm.ind.out) <- c('Const.', '$R_m$', 'SMB','HML')

xtable(capm.ind.out, type = "latex", file = "tab1_ff3_b.tex")



########### Code finish ########################
################################################
################################################
### Alternative way to calculate the GLS Rsq ### 

# lambda.ind.gls.capm <- solve(t(C.ind.t) %*% ginv(Sigma.ind.e) %*% C.ind.t) %*% t(C.ind.t) %*% ginv(Sigma.ind.e) %*% mu.ind.R
# e.ind.gls <- mu.ind.R - C.ind.t %*% lambda.gls.capm
# capm.ind.gls.Rsq1 <- 1 - var(e.ind.gls)/var(mu.ind.R)
