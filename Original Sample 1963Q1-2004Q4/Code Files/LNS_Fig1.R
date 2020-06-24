
rm(list = ls())

library('dplyr')
library('lmtest')
library('MASS')
library('foreach')
library('parallel')
library('gmm')

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

w.sigma.sim <- diag(1, K , K)
w.mu.sim <- c(0,0,0) 

Rsq.a <- matrix(NA,n.sim,K) 
 

      for(j in 1:n.sim){
          
          for (i in 1:K){  
              w <- matrix(rnorm(K*K,0,1),K,K)   #mvrnorm(K,w.mu.sim, w.sigma.sim)
              v <- rnorm(T,0,1)
              W <- as.matrix(w[,(1:i)])
              P1 <- fac %*% W + v
              #### First stage regressions to get C : R = a + P C + e ####
              P <- cbind(matrix(1,T,1),P1)
              R <- as.matrix(retex.dat[,1:N])
              
              C1 <- matrix(solve(t(P) %*% P) %*% t(P) %*% R,(i+1),N) 
              
              C <- matrix(C1[2:(i+1),],i,N)
              mu.R <- as.matrix(colMeans(R))
              
              C.t <- cbind(matrix(1,N,1),t(C))
              
              M.c <- diag(1,N,N) - C.t %*% solve(t(C.t) %*% C.t) %*% t(C.t)
              
              Rsq.a[j,i] <- 1 - t(mu.R) %*% M.c %*% mu.R/(t(mu.R - mean(mu.R)) %*% (mu.R - mean(mu.R)))
              
            }
          
        }
  
Rsq.a.out <- matrix(NA, nrow = 3, ncol = K)
   
    for (i in 1:K){
      Rsq.a.out[,i] <-  matrix(quantile(Rsq.a[,i],probs = c(0.05,0.5,0.95),na.rm = TRUE),nrow = 3, ncol = 1)
    }


colnames(Rsq.a.out) <- c('1','2','3')
rownames(Rsq.a.out) <- c('5%','50%','95%')

plot(c(1,2,3), Rsq.a.out[1,], type = 'l', col = 'red', main = 'Fig 1 Panel A', xlab = '# of Factors', ylab = expression(R^2), ylim = c(0,1), xlim = c(1,5))
lines(c(1,2,3),  Rsq.a.out[2,], col = 'blue', lty = 2)
lines(c(1,2,3), Rsq.a.out[3,], col = 'green', lty = 6)

legend('topright', c('5%-ile', 'Median', '95%-ile'), lty = c(1,2,6), col = c('red','blue','green'), cex = 0.5)



#### Figure 1: Panel B #### 

K <- 5

Rsq2.b <- matrix(NA, n.sim, K)

  for (k in 1:n.sim) {
     for (i in 1:K) {
        
           I.N <- diag(1,N,N)
           one.mat <- matrix(1,N,1)
           w.p <- matrix(rnorm(N*i,0,1), N, i)
           
           w.p.demeaned <- t( t(w.p) %*% (I.N - (one.mat %*% t(one.mat))/N))
          
          for (j in 1:i) {  
             w.p.demeaned[w.p.demeaned[,j] >= 0,j] <- w.p.demeaned[w.p.demeaned[,j] >= 0,j]/sum(w.p.demeaned[w.p.demeaned[,j] >= 0,j])
             w.p.demeaned[w.p.demeaned[,j] < 0,j] <- -w.p.demeaned[w.p.demeaned[,j] < 0,j]/sum(w.p.demeaned[w.p.demeaned[,j] < 0,j])
            }
        
          P1 <- R %*% w.p.demeaned  
          P <- cbind(matrix(1,T,1),P1) 
          
          R <- as.matrix(retex.dat[,1:N])
          
          C1 <- matrix(solve(t(P) %*% P) %*% t(P) %*% R,(i+1),N) 
          
          C <- matrix(C1[2:(i+1),],i,N)
          mu.R <- as.matrix(colMeans(R))
          
          C.t <- cbind(matrix(1,N,1),t(C))
          
          M.c <- diag(1,N,N) - C.t %*% solve(t(C.t) %*% C.t) %*% t(C.t)
          
          Rsq2.b[k,i] <- 1 - t(mu.R) %*% M.c %*% mu.R/(t(mu.R - mean(mu.R)) %*% (mu.R - mean(mu.R)))
          
        }
    }

Rsq2.b.out <- matrix(NA, nrow = 3, ncol = K)

for (i in 1:K){
  Rsq2.b.out[,i] <-  matrix(quantile(Rsq2.b[,i],probs = c(0.05,0.5,0.95),na.rm = TRUE),nrow = 3, ncol = 1)
}


colnames(Rsq2.b.out) <- c('1','2','3','4','5')
rownames(Rsq2.b.out) <- c('5%','50%','95%')

plot(c(1,2,3,4,5), Rsq2.b.out[1,], type = 'l', col = 'red', main = 'Fig 1 Panel B', xlab = '# of Factors', ylab = expression(R^2), ylim = c(0,1), xlim = c(1,5))
lines(c(1,2,3,4,5),  Rsq2.b.out[2,], col = 'blue', lty = 2)
lines(c(1,2,3,4,5), Rsq2.b.out[3,], col = 'green', lty = 6)

legend('bottomright', c('5%-ile', 'Median', '95%-ile'), lty = c(1,2,6), col = c('red','blue','green'),cex = 0.5)




#### Figure 1 Panel C: Mean zero factors retained #### 


K <- 5

Rsq2.c <- matrix(NA, n.sim, K)

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
    
        P <- cbind(matrix(1,T,1),P1) 
        
        R <- as.matrix(retex.dat[,1:N])
        
        C1 <- matrix(solve(t(P) %*% P) %*% t(P) %*% R,(i+1),N) 
        
        C <- matrix(C1[2:(i+1),],i,N)
        mu.R <- as.matrix(colMeans(R))
        
        C.t <- cbind(matrix(1,N,1),t(C))
        
        M.c <- diag(1,N,N) - C.t %*% solve(t(C.t) %*% C.t) %*% t(C.t)
        
        Rsq2.c[k,i] <- 1 - t(mu.R) %*% M.c %*% mu.R/(t(mu.R - mean(mu.R)) %*% (mu.R - mean(mu.R)))
        
  }
}

Rsq2.c.out <- matrix(NA, nrow = 3, ncol = K)

for (i in 1:K){
  Rsq2.c.out[,i] <-  matrix(quantile(Rsq2.c[,i],probs = c(0.05,0.5,0.95),na.rm = TRUE),nrow = 3, ncol = 1)
}


colnames(Rsq2.c.out) <- c('1','2','3','4','5')
rownames(Rsq2.c.out) <- c('5%','50%','95%')

plot(c(1,2,3,4,5), Rsq2.c.out[1,], type = 'l', col = 'red', main = 'Fig 1 Panel C', xlab = '# of Factors', ylab = expression(R^2), ylim = c(0,1), xlim = c(1,5))
lines(c(1,2,3,4,5),  Rsq2.c.out[2,], col = 'blue', lty = 2)
lines(c(1,2,3,4,5), Rsq2.c.out[3,], col = 'green', lty = 6)

legend('bottomright', c('5%-ile', 'Median', '95%-ile'), lty = c(1,2,6), col = c('red','blue','green'), cex = 0.5)





