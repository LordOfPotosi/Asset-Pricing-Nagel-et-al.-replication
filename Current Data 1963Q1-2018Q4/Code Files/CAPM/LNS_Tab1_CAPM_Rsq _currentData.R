rm(list = ls())

library('dplyr')
library('lmtest')
library('MASS')
library('foreach')
library('parallel')
library('gmm')
library('xtable')

#### FF25 Results ####

setwd("C:/Users/mayur/Dropbox/AP Replication Project/Replication Project Code files")
#setwd("D:/Dropbox/AP Replication Project/Replication Project Code files")

fac.dat <- read.csv('FF_current.csv', header = TRUE)
port.dat <- read.csv('25_current.csv', header = TRUE)


T <- nrow(port.dat)
N <- ncol(port.dat) - 1
#K <- ncol(fac.dat) -2 ## Exlude RF 
n.sim <- 1000

fin.dat <- merge(port.dat, fac.dat, by = 'Date')

retex.dat <- cbind((fin.dat[,2:ncol(port.dat)] - fin.dat$RF), fin.dat[,(ncol(port.dat)+1):ncol(fin.dat)]) ### Get excess returns 
K <- 1                            ## Mkt-Rf
fac <- as.matrix(fac.dat[,1:(K+1)]) 

P1 <- as.matrix(fac[,2:(K+1)])   ### I am using the FF 3 factors Mkt-RF factor 

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
C <- matrix(C1[2,],1,N)


e <- R - P%*% C1
Sigma.e <- (t(e) %*% e)/T 



#### Table 1 Results: 2nd stage regressions  ####

### OLS Results ####
mu.R <- as.matrix(colMeans(R))

C.t <- cbind(matrix(1,N,1),t(C))

lambda.capm <- solve(t(C.t) %*% C.t) %*% t(C.t) %*% mu.R

M.c <- diag(1,N,N) - C.t %*% solve(t(C.t) %*% C.t) %*% t(C.t)
capm.Rsq <- 1 - t(mu.R) %*% M.c %*% mu.R/(t(mu.R - mean(mu.R)) %*% (mu.R - mean(mu.R)))*N/(N-K)



### GLS Rsq ###
### GLS Rsq calculations + Confidence interval ###  
R.eig <- eigen(var(R))
R.Q <- R.eig$vectors
R.lambda <- R.eig$values

R.chol <- R.Q %*% diag(sqrt(R.lambda),N,N) %*% solve(R.Q)

C.gls.t <- solve(t(R.chol)) %*% C.t
mu.gls.R <- solve(t(R.chol)) %*% mu.R

lambda.gls <- solve(t(C.gls.t) %*% C.gls.t) %*% t(C.gls.t) %*% mu.gls.R

e.gls <- mu.R - C.t %*% lambda.gls
capm.gls.Rsq <- 1 - var(e.gls)/var(mu.R)

###################################################
##### Rsq Confidence interval computations: FF25###
###################################################

## Form a grid for true Rsq ## 
n.sim <- 10

rsq.min <- 0
rsq.max <- 1
incr <- 0.1
rsq.grid <- as.matrix(seq(rsq.min, rsq.max, incr))

lambda.factors <- lambda.capm[2:(K+1),1]
mu.fit.R <- C.t %*% lambda.capm 

mu.fit.gls.R <- C.gls.t %*% lambda.gls


Rsq.quantiles <- matrix(NA,(length(rsq.grid)),3) 
colnames(Rsq.quantiles) <- c('5%','Median','95%')

Rsq.gls.quantiles <- matrix(NA,(length(rsq.grid)),3) 
colnames(Rsq.gls.quantiles) <- c('5%','Median','95%')




### demean factors ###
P1.demeaned <- P1 - matrix(1,T,N) %*% diag(colMeans(P1),N,K)


for(r2.ind in 1:nrow(Rsq.quantiles)) {
  ### Get the primitives of the first stage: alpha, h, sigma.alpha ###
  r2 <- rsq.grid[r2.ind] 
  h <- sqrt(r2*var(mu.R)/var(mu.fit.R)) 
  sigma.alpha <- sqrt((1-r2)*var(mu.R))
  
  ### GLS alpha and h  ### 
  h.gls <- sqrt(r2*var(mu.gls.R)/var(mu.fit.gls.R)) 
  sigma.gls.alpha <- sqrt((1-r2)*var(mu.fit.gls.R))
  
  if(r2 == 1) { alpha.sim <- matrix(rnorm(N,0,sigma.alpha),N,1)} else {
    alpha.sim <- matrix(rnorm(N,0,sigma.alpha),N,1)
    alpha.sim <- alpha.sim -  mu.fit.R %*% as.matrix(cov(mu.fit.R, alpha.sim)/var(mu.fit.R))
    alpha.sim <- alpha.sim %*% as.matrix(sigma.alpha/sqrt(var(alpha.sim)))
    alpha.sim <- alpha.sim - mean(alpha.sim)
  }
  
  if (r2 ==1){alpha.gls.sim <- matrix(rnorm(N,0,sigma.gls.alpha),N,1)} else {
    alpha.gls.sim <- matrix(rnorm(N,0,sigma.gls.alpha),N,1)
    alpha.gls.sim <- alpha.gls.sim -  mu.fit.gls.R %*% as.matrix(cov(mu.fit.gls.R, alpha.gls.sim)/var(mu.fit.gls.R))
    alpha.gls.sim <- alpha.gls.sim %*% as.matrix(sigma.gls.alpha/sqrt(var(alpha.gls.sim)))
    alpha.gls.sim <- alpha.gls.sim - mean(alpha.gls.sim)
  }
  
  Rsq.sim <- matrix(NA,n.sim,1)
  Rsq.gls.sim <- matrix(NA,n.sim,1)
  
  for(i in 1:n.sim){
    
    #### Sample the factors ### 
    x <- c(1:T)
    sample.index <- matrix(sample(x, T*K,replace = TRUE),T,K)
    for ( j in 1 : K){
      if (j==1){P.temp <- P1.demeaned[sample.index[,j],j] } 
      else {P.temp <- cbind(P.temp, P1.demeaned[sample.index[,j],j])}
    }
    P.sample <- matrix(P.temp, T, K)
    P.gls.sample <- cbind(matrix(1,T,1),P.sample)
    
    ###  means for OLS ### 
    R.sim.mean <- diag(as.numeric(mu.fit.R %*% h + alpha.sim),N,N) %*% matrix(1, N, T) + t(P.sample %*% C)
    R.sim <- matrix(NA,T,N) ## length(sample.index)
    
    ### now simulate the returns for OLS #### 
    for (t in 1:T){  ## length(sample.index)
      R.sim[t,] <- mvrnorm(1,R.sim.mean[,t],Sigma.e) 
    }
    
    ### means for GLS ###
    R.sim.gls.mean <- R.chol %*% (diag(as.numeric(mu.fit.gls.R %*% h + alpha.gls.sim),N,N) %*% matrix(1, N, T) + t(P.gls.sample %*% t(C.gls.t)))
    R.sim.gls <- matrix(NA,T,N)
    
    ### now simulate the returns for GLS #### 
    for (t in 1:T){  ## length(sample.index)
      R.sim.gls[t,] <- mvrnorm(1,R.sim.gls.mean[,t],Sigma.e) 
    }
    
    #### Compute Rsq #### 
    P.sim <- cbind(matrix(1,T,1),P.sample)
    C1.sim <- matrix(solve(t(P.sim) %*% P.sim) %*% t(P.sim) %*% R.sim,(K+1),N) ### Time series betas 
    C.sim <- matrix(C1.sim[2:(K+1),],K,N)                          ### Cross section regressor matrix          
    
    
    #### Table 1 Results: 2nd stage regressions  ####
    
    ### OLS Results ####
    mu.R.sim <- as.matrix(colMeans(R.sim))
    C.t.sim <- cbind(matrix(1,N,1),t(C.sim))
    
    lambda.capm.sim <- ginv(t(C.t.sim) %*% C.t.sim) %*% t(C.t.sim) %*% mu.R.sim
    
    M.c.sim <- diag(1,N,N) - C.t.sim %*% ginv(t(C.t.sim) %*% C.t.sim) %*% t(C.t.sim)
    Rsq.sim[i] <- 1 - (t(mu.R.sim) %*% M.c.sim %*% mu.R.sim)/(t(mu.R.sim - mean(mu.R.sim)) %*% (mu.R.sim - mean(mu.R.sim)))*N/(N-K-1)
    
    ### GLS Rsq results ### 
    mu.R.sim <- as.matrix(colMeans(R.sim.gls))
    C1.ind.sim <- matrix(solve(t(P.sim) %*% P.sim) %*% t(P.sim) %*% R.sim.gls ,(K+1),N) ### Time series betas 
    C.ind.sim <- matrix(C1.ind.sim[2:(K+1),],K,N)                                       ### Cross section regressor matrix         
    C.sim.gls.t <- cbind(matrix(1,N,1), t(C.ind.sim))
    
    R.ind.eig <- eigen(var(R.sim.gls))
    R.ind.Q <- R.ind.eig$vectors
    R.ind.lambda <- R.ind.eig$values
    
    R.gls.chol <- R.ind.Q %*% diag(sqrt(R.ind.lambda),N,N) %*% solve(R.ind.Q)
    
    C.ind.gls.t <- solve(t(R.gls.chol)) %*% C.sim.gls.t
    mu.ind.gls.R <- solve(t(R.gls.chol)) %*% mu.R.sim
    
    M.ind.gls.c <- diag(1,N,N) - C.ind.gls.t %*% solve(t(C.ind.gls.t) %*% C.ind.gls.t) %*% t(C.ind.gls.t)
    Rsq.gls.sim[i,] <- 1 - (t(mu.ind.gls.R) %*% M.ind.gls.c %*% mu.ind.gls.R/(t(mu.ind.gls.R - mean(mu.ind.gls.R)) %*% (mu.ind.gls.R - mean(mu.ind.gls.R))))*N/(N-K-1)
    
  }
  Rsq.quantiles[r2.ind,] <- matrix(quantile(Rsq.sim, probs = c(0.05,0.5,0.95),na.rm = TRUE),nrow = 1, ncol = 3)
  Rsq.gls.quantiles[r2.ind,] <- matrix(quantile(Rsq.gls.sim, probs = c(0.05,0.5,0.95),na.rm = TRUE),nrow = 1, ncol = 3)
  
}

### Plot the results ### 
rsq.grid.mod <- rsq.grid[1:r2.ind]
plot(rsq.grid.mod, Rsq.quantiles[,1], ylab = expression('Simulated' ~ R^2), xlab = expression('True' ~ R^2), main = expression(R^2~ 'OLS Confidence Intervals: FF25, FF3 factors'), type = 'l', xlim = c(0,max(rsq.grid)), ylim = c(min(Rsq.quantiles),max(Rsq.quantiles)))
lines(rsq.grid.mod, Rsq.quantiles[,2], col = 'red')
lines(rsq.grid.mod, Rsq.quantiles[,3], col = 'blue')
abline(h = capm.Rsq,lwd = 2, col = 'green', lty = 2)
legend('topleft', c('5%-lie','Median','95%-ile','Rsq'), lty = c(1,1,1,2), col = c('black','red','blue','green'),cex = 0.9, bg = 'transparent')
#locator() 

rsq.conf.int <- matrix(c(0.0,0.5),2,1)
rsq.out <- matrix(c(round(capm.Rsq,2), paste('(', as.character(rsq.conf.int[1]), ',', as.character(rsq.conf.int[2]) ,')')),2,1)


### Plot the GLS results ### 

plot(rsq.grid, Rsq.gls.quantiles[,1], ylab = expression('Simulated' ~ R^2), xlab = expression('True' ~ R^2), main = expression(R^2~ 'GLS Confidence Intervals: FF25, FF3 factors'), type = 'l', xlim = c(0,max(rsq.grid)), ylim = c(min(Rsq.gls.quantiles),max(Rsq.gls.quantiles)))
lines(rsq.grid, Rsq.gls.quantiles[,2], col = 'red')
lines(rsq.grid, Rsq.gls.quantiles[,3], col = 'blue')
abline(h = capm.gls.Rsq,lwd = 2, col = 'green', lty = 2)
legend('topleft', c('5%-lie','Median','95%-ile','Rsq'), lty = c(1,1,1,2), col = c('black','red','blue','green'),cex = 0.7, bg = 'transparent')
#locator() 

rsq.conf.gls.int <- matrix(c(0.00,0.01),2,1)
rsq.gls.out <- matrix(c(round(capm.gls.Rsq,2), paste('(', as.character(rsq.conf.gls.int[1]), ',', as.character(rsq.conf.gls.int[2]) ,')')),2,1)



#######################################
#### FF25 + 30 Industry Portfolios #### 
#######################################


ind.dat <- read.csv('30ind_current.csv', header = TRUE)
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

e.ind.cs <- mu.ind.R - C.ind.t %*% lambda.ind.capm
capm.ind.Rsq <- 1 - var(e.ind.cs)/var(mu.ind.R)*(N)/(N-K)


### GLS Rsq Industry portfolios ###
R.ind.eig <- eigen(var(R.ind))
R.ind.Q <- R.ind.eig$vectors
R.ind.lambda <- R.ind.eig$values

R.chol <- R.ind.Q %*% diag(sqrt(R.ind.lambda),N.ind,N.ind) %*% solve(R.ind.Q)

C.ind.gls.t <- solve(t(R.chol)) %*% C.ind.t

mu.ind.gls.R <- solve(t(R.chol)) %*% mu.ind.R

lambda.ind.gls <- solve(t(C.ind.gls.t) %*% C.ind.gls.t) %*% t(C.ind.gls.t) %*% mu.ind.gls.R

M.ind.gls.c <- diag(1,N.ind,N.ind) - C.ind.gls.t %*% solve(t(C.ind.gls.t) %*% C.ind.gls.t) %*% t(C.ind.gls.t)
capm.ind.gls.Rsq <- 1 - (t(mu.ind.gls.R) %*% M.ind.gls.c %*% mu.ind.gls.R/(t(mu.ind.gls.R - mean(mu.ind.gls.R)) %*% (mu.ind.gls.R - mean(mu.ind.gls.R))))*N.ind/(N.ind-K-1)


########################################################
##### FF25 + 30 Industry Rsq  Confidence intervals######
########################################################


lambda.ind.factors <- lambda.ind.capm[2:(K+1),1]
mu.fit.ind.R <- t(C.ind) %*% lambda.ind.factors

mu.fit.ind.gls.R <- C.ind.gls.t %*% lambda.ind.gls

Rsq.ind.quantiles <- matrix(NA,(length(rsq.grid)),3) 
colnames(Rsq.ind.quantiles) <- c('5%','Median','95%')

Rsq.ind.gls.quantiles <- matrix(NA,(length(rsq.grid)),3) 
colnames(Rsq.ind.gls.quantiles) <- c('5%','Median','95%')


for(r2.ind in 1:nrow(Rsq.ind.quantiles)) {
  ### Get the primitives of the first step: alpha, h, sigma.alpha ###
  r2 <- rsq.grid[r2.ind] 
  h <- sqrt(r2*var(mu.ind.R)/var(mu.fit.ind.R)) 
  sigma.alpha <- sqrt((1-r2)*var(mu.ind.R))
  
  h.gls <- sqrt(r2*var(mu.ind.gls.R)/var(mu.fit.ind.gls.R)) 
  sigma.gls.alpha <- sqrt((1-r2)*var(mu.ind.gls.R))
  
  #### alpha, h for OLS #### 
  if (r2 ==1){alpha.sim <- matrix(rnorm(N.ind,0,sigma.alpha),N.ind,1)} else {
    alpha.sim <- matrix(rnorm(N.ind,0,sigma.alpha),N.ind,1)
    alpha.sim <- alpha.sim -  mu.fit.ind.R %*% as.matrix(cov(mu.fit.ind.R, alpha.sim)/var(mu.fit.ind.R))
    alpha.sim <- alpha.sim %*% as.matrix(sigma.alpha/sqrt(var(alpha.sim)))
    alpha.sim <- alpha.sim - mean(alpha.sim)
  }
  
  #### alpha, h for GLS #### 
  if (r2 ==1){alpha.gls.sim <- matrix(rnorm(N.ind,0,sigma.gls.alpha),N.ind,1)} else {
    alpha.gls.sim <- matrix(rnorm(N.ind,0,sigma.gls.alpha),N.ind,1)
    alpha.gls.sim <- alpha.gls.sim -  mu.fit.ind.gls.R %*% as.matrix(cov(mu.fit.ind.gls.R, alpha.gls.sim)/var(mu.fit.ind.gls.R))
    alpha.gls.sim <- alpha.gls.sim %*% as.matrix(sigma.gls.alpha/sqrt(var(alpha.gls.sim)))
    alpha.gls.sim <- alpha.gls.sim - mean(alpha.gls.sim)
  }
  
  Rsq.sim <- matrix(NA,n.sim,1)
  Rsq.gls.sim <- matrix(NA,n.sim,1)
  
  for(i in 1:n.sim){
    
    ### Sample the Factors ### 
    x <- c(1:T)
    sample.index <- matrix(sample(x, T*K,replace = TRUE),T,K)
    for ( j in 1 : K){
      if (j==1){P.temp <- P1.demeaned[sample.index[,j],j] } 
      else {P.temp <- cbind(P.temp, P1.demeaned[sample.index[,j],j])}
    }
    P.sample <- matrix(P.temp, T, K)
    P.gls.sample <- cbind(matrix(1,T,1),P.sample)
    
    ### means for OLS 
    R.sim.mean <- diag(as.numeric(mu.fit.ind.R %*% h + alpha.sim),N.ind,N.ind) %*% matrix(1, N.ind, T) + t(P.sample %*% C.ind)
    R.sim <- matrix(NA,T,N.ind)
    
    ### now simulate the returns for OLS #### 
    for (t in 1:T){
      R.sim[t,] <- mvrnorm(1,R.sim.mean[,t],Sigma.ind.e) 
    }
    
    ### means for GLS 
    R.sim.gls.mean <- R.chol %*% (diag(as.numeric(mu.fit.ind.gls.R %*% h + alpha.gls.sim),N.ind,N.ind) %*% matrix(1, N.ind, T) + t(P.gls.sample %*% t(C.ind.gls.t)))
    R.sim.gls <- matrix(NA,T,N.ind)
    
    ### now simulate the returns for OLS #### 
    for (t in 1:T){
      R.sim.gls[t,] <- mvrnorm(1,R.sim.gls.mean[,t],Sigma.ind.e) 
    }
    
    
    #### Compute Rsq OLS #### 
    P.sim <- cbind(matrix(1,T,1),P.sample)
    C1.ind.sim <- matrix(solve(t(P.sim) %*% P.sim) %*% t(P.sim) %*% R.sim ,(K+1),N.ind) ### Time series betas 
    C.ind.sim <- matrix(C1.ind.sim[2:(K+1),],K,N.ind)                          ### Cross section regressor matrix          
    
    
    #### Table 1 Results: 2nd stage regressions  ####
    
    ### OLS Rsq Results ####
    mu.R.sim <- as.matrix(colMeans(R.sim))
    C.t.sim <- cbind(matrix(1,N.ind,1),t(C.ind.sim))
    
    lambda.capm.sim <- solve(t(C.t.sim) %*% C.t.sim) %*% t(C.t.sim) %*% mu.R.sim
    
    M.c.sim <- diag(1,N.ind,N.ind) - C.t.sim %*% solve(t(C.t.sim) %*% C.t.sim) %*% t(C.t.sim)
    
    Rsq.sim[i] <- 1 - (t(mu.R.sim) %*% M.c.sim %*% mu.R.sim)/(t(mu.R.sim - mean(mu.R.sim)) %*% (mu.R.sim - mean(mu.R.sim)))*N.ind/(N.ind-K)
    
    ### GLS Rsq results ### 
    mu.R.sim <- as.matrix(colMeans(R.sim.gls))
    C1.ind.sim <- matrix(solve(t(P.sim) %*% P.sim) %*% t(P.sim) %*% R.sim.gls ,(K+1),N.ind) ### Time series betas 
    C.ind.sim <- matrix(C1.ind.sim[2:(K+1),],K,N.ind)                                       ### Cross section regressor matrix         
    C.sim.gls.t <- cbind(matrix(1,N.ind,1), t(C.ind.sim))
    
    R.ind.eig <- eigen(var(R.sim.gls))
    R.ind.Q <- R.ind.eig$vectors
    R.ind.lambda <- R.ind.eig$values
    
    R.gls.chol <- R.ind.Q %*% diag(sqrt(R.ind.lambda),N.ind,N.ind) %*% solve(R.ind.Q)
    
    C.ind.gls.t <- solve(t(R.gls.chol)) %*% C.sim.gls.t
    mu.ind.gls.R <- solve(t(R.gls.chol)) %*% mu.R.sim
    
    M.ind.gls.c <- diag(1,N.ind,N.ind) - C.ind.gls.t %*% solve(t(C.ind.gls.t) %*% C.ind.gls.t) %*% t(C.ind.gls.t)
    Rsq.gls.sim[i,] <- 1 - (t(mu.ind.gls.R) %*% M.ind.gls.c %*% mu.ind.gls.R/(t(mu.ind.gls.R - mean(mu.ind.gls.R)) %*% (mu.ind.gls.R - mean(mu.ind.gls.R))))*N.ind/(N.ind-K-1)
    
  }
  
  Rsq.ind.quantiles[r2.ind,] <- matrix(quantile(Rsq.sim, probs = c(0.05,0.5,0.95),na.rm = TRUE),nrow = 1, ncol = 3)
  Rsq.ind.gls.quantiles[r2.ind,] <- matrix(quantile(Rsq.gls.sim, probs = c(0.05,0.5,0.95),na.rm = TRUE),nrow = 1, ncol = 3)
  
}

### Plot the results: First OLS ### 

plot(rsq.grid, Rsq.ind.quantiles[,1], ylab = expression('Simulated' ~ R^2), xlab = expression('True' ~ R^2), main = expression(R^2~ 'OLS Confidence Intervals: FF25 + Ind30, FF3 factors'), type = 'l', xlim = c(0,max(rsq.grid)), ylim = c(min(Rsq.ind.quantiles),max(Rsq.ind.quantiles)))
lines(rsq.grid, Rsq.ind.quantiles[,2], col = 'red')
lines(rsq.grid, Rsq.ind.quantiles[,3], col = 'blue')
abline(h = capm.ind.Rsq,lwd = 2, col = 'green', lty = 2)
legend('topleft', c('5%-lie','Median','95%-ile','Rsq'), lty = c(1,1,1,2), col = c('black','red','blue','green'),cex = 0.7, bg = 'transparent')
#locator() 

rsq.conf.ind.int <- matrix(c(0.00,0.00),2,1)
rsq.ind.out <- matrix(c(round(capm.ind.Rsq,2), paste('(', as.character(rsq.conf.ind.int[1]), ',', as.character(rsq.conf.ind.int[2]) ,')')),2,1)

### Plot the GLS Rsq confidence interval results ### 
plot(rsq.grid, Rsq.ind.gls.quantiles[,1], ylab = expression('Simulated' ~ R^2), xlab = expression('True' ~ R^2), main = expression(R^2~ 'GLS Confidence Intervals: FF25 + Ind30, FF3 factors'), type = 'l', xlim = c(0,max(rsq.grid)), ylim = c(min(Rsq.ind.gls.quantiles),max(Rsq.ind.gls.quantiles)))
lines(rsq.grid, Rsq.ind.gls.quantiles[,2], col = 'red')
lines(rsq.grid, Rsq.ind.gls.quantiles[,3], col = 'blue')
abline(h = capm.ind.gls.Rsq,lwd = 2, col = 'green', lty = 2)
legend('topleft', c('5%-lie','Median','95%-ile','Rsq'), lty = c(1,1,1,2), col = c('black','red','blue','green'),cex = 0.7, bg = 'transparent')
#locator() 

rsq.conf.ind.gls.int <- matrix(c(0.00,0.04),2,1)
rsq.ind.gls.out <- matrix(c(round(capm.ind.gls.Rsq,2), paste('(', as.character(rsq.conf.ind.gls.int[1]), ',', as.character(rsq.conf.ind.gls.int[2]) ,')')),2,1)


