
rm(list = ls())

library('dplyr')
library('lmtest')
library('MASS')
library('foreach')
library('parallel')
library('gmm')

setwd("C:/Users/mayur/Dropbox/AP Replication Project/Replication Project Code files")

fac.dat <- read.csv('FF_current.csv', header = TRUE)
port.dat <- read.csv('25_current.csv', header = TRUE)
ind.dat <- read.csv('30ind_current.csv', header = TRUE)

T <- nrow(port.dat)
N <- ncol(port.dat) - 1
N.ind <- ncol(port.dat) + ncol(ind.dat) - 2
K <- ncol(fac.dat) -2 ## Exlude RF 
n.sim <- 1000

fin.dat <- merge(port.dat, fac.dat,by = 'Date')
test.dat <- merge(port.dat, ind.dat,by = 'Date')
fin.ind.dat <- merge(test.dat, fac.dat, by = 'Date')

retex.dat <- cbind((fin.dat[,2:ncol(port.dat)] - fin.dat$RF), fin.dat[,(ncol(port.dat)+1):ncol(fin.dat)]) ### Get excess returns 
retex.ind.dat <- cbind((fin.ind.dat[,2:ncol(port.dat)] - fin.ind.dat$RF), (fin.ind.dat[,(ncol(port.dat)+1):(ncol(port.dat) + ncol(ind.dat))] - fin.ind.dat$RF), fin.ind.dat[,(ncol(port.dat) +ncol(ind.dat)+1):ncol(fin.ind.dat)]) ### Get excess returns 


fac <- as.matrix(fac.dat[,2:(K+1)]) 

w.sigma.sim <- diag(1, K , K)
w.mu.sim <- c(0,0,0) 

Rsq.a <- matrix(NA,n.sim,K) 
Rsq.ind.a <- matrix(NA,n.sim,K) 

R <- as.matrix(retex.dat[,1:N])
R.ind <- as.matrix(retex.ind.dat[,1:N.ind])

mu.R <- as.matrix(colMeans(R))
mu.R.ind <- as.matrix(colMeans(R.ind))

for(j in 1:n.sim){
  
  for (i in 1:K){  
    w <- matrix(rnorm(K*K,0,1),K,K)   #mvrnorm(K,w.mu.sim, w.sigma.sim)
    v <- rnorm(T,0,1)
    W <- as.matrix(w[,(1:i)])
    P1 <- fac %*% W + v
    #### First stage regressions to get C : R = a + P C + e ####
    P <- cbind(matrix(1,T,1),P1)
    
    C1 <- matrix(solve(t(P) %*% P) %*% t(P) %*% R,(i+1),N) 
    C1.ind <-  matrix(solve(t(P) %*% P) %*% t(P) %*% R.ind,(i+1),N.ind) 
    
    C <- matrix(C1[2:(i+1),],i,N)
    C.ind <- matrix(C1.ind[2:(i+1),],i,N.ind)
    
    C.t <- cbind(matrix(1,N,1),t(C))
    C.ind.t <-  cbind(matrix(1,N.ind,1),t(C.ind))
    
    M.c <- diag(1,N,N) - C.t %*% solve(t(C.t) %*% C.t) %*% t(C.t)
    M.ind.c <- diag(1,N.ind,N.ind) - C.ind.t %*% solve(t(C.ind.t) %*% C.ind.t) %*% t(C.ind.t)
    
    Rsq.a[j,i] <- 1 - t(mu.R) %*% M.c %*% mu.R/(t(mu.R - mean(mu.R)) %*% (mu.R - mean(mu.R)))
    Rsq.ind.a[j,i] <- 1 - t(mu.R.ind) %*% M.ind.c %*% mu.R.ind/(t(mu.R.ind - mean(mu.R.ind)) %*% (mu.R.ind - mean(mu.R.ind)))
  }
  
}

Rsq.a.out <- matrix(NA, nrow = 3, ncol = K)
Rsq.ind.a.out <- matrix(NA, nrow = 3, ncol = K)

for (i in 1:K){
  Rsq.a.out[,i] <-  matrix(quantile(Rsq.a[,i],probs = c(0.05,0.5,0.95),na.rm = TRUE),nrow = 3, ncol = 1)
  Rsq.ind.a.out[,i] <-  matrix(quantile(Rsq.ind.a[,i],probs = c(0.05,0.5,0.95),na.rm = TRUE),nrow = 3, ncol = 1)
  
}


colnames(Rsq.a.out) <- c('1','2','3')
rownames(Rsq.a.out) <- c('5%','50%','95%')

colnames(Rsq.ind.a.out) <- c('1','2','3')
rownames(Rsq.ind.a.out) <- c('5%','50%','95%')

plot(c(1,2,3), Rsq.a.out[2,], type = 'l', col = 'red', main = 'Fig 3 Panel A', xlab = '# of Factors', ylab = expression(R^2), ylim = c(0,1), xlim = c(1,5),lty = 2)
lines(c(1,2,3),  Rsq.a.out[3,], col = 'blue', lty = 4)
lines(c(1,2,3), Rsq.ind.a.out[2,], col = 'red', lty = 6, lwd = 2)
lines(c(1,2,3), Rsq.ind.a.out[3,], col = 'blue', lty = 8, lwd = 2)
legend('topright', c( 'FF 25 Median', 'FF 25 95%-ile','FF25 + Ind30 Median','FF 25 + Ind30 95%-ile'), lty = c(2,4,6,8), col = c('red','blue','red','blue'), lwd = c(1,1,2,2), cex = 0.5)



#### Figure 1: Panel B #### 

K <- 5

Rsq2.b <- matrix(NA, n.sim, K)
Rsq2.ind.b <- matrix(NA, n.sim, K)

for (k in 1:n.sim) {
  for (i in 1:K) {
    
    I.N <- diag(1,N,N)
    I.N.ind <- diag(1,N.ind, N.ind) 
    
    one.mat <- matrix(1,N,1)
    one.ind.mat <- matrix(1,N.ind,1)
    
    w.p <- matrix(rnorm(N*i,0,1), N, i)
    w.ind.p <- matrix(rnorm(N.ind*i,0,1),N.ind,i)
    
    w.p.demeaned <- t( t(w.p) %*% (I.N - (one.mat %*% t(one.mat))/N))
    w.ind.p.demeaned <- t( t(w.ind.p) %*% (I.N.ind - (one.ind.mat %*% t(one.ind.mat))/N.ind))
      
    for (j in 1:i) {  
      w.p.demeaned[w.p.demeaned[,j] >= 0,j] <- w.p.demeaned[w.p.demeaned[,j] >= 0,j]/sum(w.p.demeaned[w.p.demeaned[,j] >= 0,j])
      w.p.demeaned[w.p.demeaned[,j] < 0,j] <- -w.p.demeaned[w.p.demeaned[,j] < 0,j]/sum(w.p.demeaned[w.p.demeaned[,j] < 0,j])
      w.ind.p.demeaned[w.ind.p.demeaned[,j] >= 0,j] <- w.ind.p.demeaned[w.ind.p.demeaned[,j] >= 0,j]/sum(w.ind.p.demeaned[w.ind.p.demeaned[,j] >= 0,j])
      w.ind.p.demeaned[w.ind.p.demeaned[,j] < 0,j] <- -w.ind.p.demeaned[w.ind.p.demeaned[,j] < 0,j]/sum(w.ind.p.demeaned[w.ind.p.demeaned[,j] < 0,j])
    }
    
    P1 <- R %*% w.p.demeaned  
    P1.ind <- R.ind %*% w.ind.p.demeaned 
    
    P <- cbind(matrix(1,T,1),P1)
    P.ind <- cbind(matrix(1,T,1), P1.ind)
    
    C1 <- matrix(solve(t(P) %*% P) %*% t(P) %*% R,(i+1),N) 
    C1.ind <-  matrix(solve(t(P) %*% P) %*% t(P) %*% R.ind,(i+1),N.ind)
    
    C <- matrix(C1[2:(i+1),],i,N)
    C.ind <- matrix(C1.ind[2:(i+1),],i,N.ind)
    
    
    C.t <- cbind(matrix(1,N,1),t(C))
    C.ind.t <-  cbind(matrix(1,N.ind,1),t(C.ind))
    
    M.c <- diag(1,N,N) - C.t %*% solve(t(C.t) %*% C.t) %*% t(C.t)
    M.ind.c <- diag(1,N.ind,N.ind) - C.ind.t %*% solve(t(C.ind.t) %*% C.ind.t) %*% t(C.ind.t)
    
    Rsq2.b[k,i] <- 1 - t(mu.R) %*% M.c %*% mu.R/(t(mu.R - mean(mu.R)) %*% (mu.R - mean(mu.R)))
    Rsq2.ind.b[k,i] <- 1 - t(mu.R.ind) %*% M.ind.c %*% mu.R.ind/(t(mu.R.ind - mean(mu.R.ind)) %*% (mu.R.ind - mean(mu.R.ind)))
    
  }
}

Rsq2.b.out <- matrix(NA, nrow = 3, ncol = K)
Rsq2.ind.b.out <- matrix(NA, nrow = 3, ncol = K)

for (i in 1:K){
  Rsq2.b.out[,i] <-  matrix(quantile(Rsq2.b[,i],probs = c(0.05,0.5,0.95),na.rm = TRUE),nrow = 3, ncol = 1)
  Rsq2.ind.b.out[,i] <-  matrix(quantile(Rsq2.ind.b[,i],probs = c(0.05,0.5,0.95),na.rm = TRUE),nrow = 3, ncol = 1)
  }


colnames(Rsq2.b.out) <- c('1','2','3','4','5')
rownames(Rsq2.b.out) <- c('5%','50%','95%')

colnames(Rsq2.ind.b.out) <- c('1','2','3','4','5')
rownames(Rsq2.ind.b.out) <- c('5%','50%','95%')


plot(c(1,2,3,4,5), Rsq2.b.out[2,], type = 'l', col = 'red', main = 'Fig 3 Panel B', xlab = '# of Factors', ylab = expression(R^2), ylim = c(0,1), xlim = c(1,5),lty = 2)
lines(c(1,2,3,4,5),  Rsq2.b.out[3,], col = 'blue', lty = 4)
lines(c(1,2,3,4,5), Rsq2.ind.b.out[2,], col = 'red', lty = 6, lwd = 2)
lines(c(1,2,3,4,5), Rsq2.ind.b.out[3,], col = 'blue', lty = 8, lwd = 2)
legend('topright', c( 'FF 25 Median', 'FF 25 95%-ile','FF25 + Ind30 Median','FF 25 + Ind30 95%-ile'), lty = c(2,4,6,8), col = c('red','blue','red','blue'), lwd = c(1,1,2,2), cex = 0.5, bg = 'transparent')


#### Figure 1 Panel C: Mean zero factors retained #### 
K <- 5

Rsq2.c <- matrix(NA, n.sim, K)
Rsq2.ind.c <- matrix(NA, n.sim, K)

for (k in 1:n.sim) {
  for (i in 1:K) {
    while(TRUE){
      I.N <- diag(1,N,N)
      one.mat <- matrix(1,N,1)
      w.p <- matrix(rnorm(N*i,0,1), N, i)
      
      w.p.demeaned <- t( t(w.p) %*% (I.N - (one.mat %*% t(one.mat))/N))
      
      for (j in 1:i) {  
        w.p.demeaned[w.p.demeaned[,j] >= 0,j] <- w.p.demeaned[w.p.demeaned[,j] >= 0,j]/sum(w.p.demeaned[w.p.demeaned[,j] >= 0,j])
        w.p.demeaned[w.p.demeaned[,j] < 0,j] <- -w.p.demeaned[w.p.demeaned[,j] < 0,j]/sum(w.p.demeaned[w.p.demeaned[,j] < 0,j])
      }
      
      P1 <- R %*% w.p.demeaned
      if(sum((colMeans(P1))^2) <= 10^(-2)){ break }
    }
    
    while(TRUE){
    
      I.N.ind <- diag(1,N.ind, N.ind) 
      one.ind.mat <- matrix(1,N.ind,1)
      w.ind.p <- matrix(rnorm(N.ind*i,0,1),N.ind,i)
      w.ind.p.demeaned <- t( t(w.ind.p) %*% (I.N.ind - (one.ind.mat %*% t(one.ind.mat))/N.ind))
      
      
      for (j in 1:i) {  
        w.ind.p.demeaned[w.ind.p.demeaned[,j] >= 0,j] <- w.ind.p.demeaned[w.ind.p.demeaned[,j] >= 0,j]/sum(w.ind.p.demeaned[w.ind.p.demeaned[,j] >= 0,j])
        w.ind.p.demeaned[w.ind.p.demeaned[,j] < 0,j] <- -w.ind.p.demeaned[w.ind.p.demeaned[,j] < 0,j]/sum(w.ind.p.demeaned[w.ind.p.demeaned[,j] < 0,j])
      }
      
      P1.ind <- R.ind %*% w.ind.p.demeaned
      if(sum((colMeans(P1.ind))^2) <= 10^(-2)){ break }
    }
    
    
    P <- cbind(matrix(1,T,1),P1)
    P.ind <- cbind(matrix(1,T,1), P1.ind)
    
    C1 <- matrix(solve(t(P) %*% P) %*% t(P) %*% R,(i+1),N) 
    C1.ind <-  matrix(solve(t(P) %*% P) %*% t(P) %*% R.ind,(i+1),N.ind)
    
    C <- matrix(C1[2:(i+1),],i,N)
    C.ind <- matrix(C1.ind[2:(i+1),],i,N.ind)
    
    C.t <- cbind(matrix(1,N,1),t(C))
    C.ind.t <-  cbind(matrix(1,N.ind,1),t(C.ind))
    
    M.c <- diag(1,N,N) - C.t %*% solve(t(C.t) %*% C.t) %*% t(C.t)
    M.ind.c <- diag(1,N.ind,N.ind) - C.ind.t %*% solve(t(C.ind.t) %*% C.ind.t) %*% t(C.ind.t)
    
    Rsq2.c[k,i] <- 1 - t(mu.R) %*% M.c %*% mu.R/(t(mu.R - mean(mu.R)) %*% (mu.R - mean(mu.R)))
    Rsq2.ind.c[k,i] <- 1 - t(mu.R.ind) %*% M.ind.c %*% mu.R.ind/(t(mu.R.ind - mean(mu.R.ind)) %*% (mu.R.ind - mean(mu.R.ind)))
    
  }
}

Rsq2.c.out <- matrix(NA, nrow = 3, ncol = K)
Rsq2.ind.c.out <- matrix(NA, nrow = 3, ncol = K)


for (i in 1:K){
  Rsq2.c.out[,i] <-  matrix(quantile(Rsq2.c[,i],probs = c(0.05,0.5,0.95),na.rm = TRUE),nrow = 3, ncol = 1)
  Rsq2.ind.c.out[,i] <-  matrix(quantile(Rsq2.ind.c[,i],probs = c(0.05,0.5,0.95),na.rm = TRUE),nrow = 3, ncol = 1)
  
  }


colnames(Rsq2.c.out) <- c('1','2','3','4','5')
rownames(Rsq2.c.out) <- c('5%','50%','95%')

colnames(Rsq2.ind.c.out) <- c('1','2','3','4','5')
rownames(Rsq2.ind.c.out) <- c('5%','50%','95%')


plot(c(1,2,3,4,5), Rsq2.c.out[2,], type = 'l', col = 'red', main = 'Fig 3 Panel C', xlab = '# of Factors', ylab = expression(R^2), ylim = c(0,1), xlim = c(1,5),lty = 2)
lines(c(1,2,3,4,5),  Rsq2.c.out[3,], col = 'blue', lty = 4)
lines(c(1,2,3,4,5), Rsq2.ind.c.out[2,], col = 'red', lty = 6, lwd = 2)
lines(c(1,2,3,4,5), Rsq2.ind.c.out[3,], col = 'blue', lty = 8, lwd = 2)
legend('topright', c( 'FF 25 Median', 'FF 25 95%-ile','FF25 + Ind30 Median','FF 25 + Ind30 95%-ile'), lty = c(2,4,6,8), col = c('red','blue','red','blue'), lwd = c(1,1,2,2), cex = 0.5, bg = 'transparent')


