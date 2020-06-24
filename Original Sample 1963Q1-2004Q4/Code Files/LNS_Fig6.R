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

### GRS-F statistic Testing CAPM ### 
      K <- 1                            ## Mkt-Rf
      P1 <- as.matrix(fac[,'Mkt.RF'])   ### I am using the FF 3 factors Mkt-RF factor 
      
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
      
      e <- R - P%*% C1
      Sigma.e <- (t(e) %*% e)/T 
      
      
      a.hat <- t(matrix(C1[1,],1,N))
    
      theta.z <- t(a.hat) %*% solve(Sigma.e) %*% a.hat 
      c <- (1+P1.SR)/T
      
      F.grs <- (T-N-K)/(T*N)*c^(-1)*theta.z
      
      F.grs.lns <- c^(-1)*theta.z*(T-N-K)/(N*(T-K-1))
      
      C <- matrix(C1[2,],1,N)
      
      
    ### Simulations for F GRS confidence interval ### 
      
      theta.z.min <- 0
      theta.z.max <- 1
      incr <- 0.05
      theta.z.grid <- as.matrix(seq(theta.z.min, theta.z.max, incr))
      
      F.grs.quantiles <- matrix(NA,length(theta.z.grid),3)
      
      for (theta.ind in 1:length(theta.z.grid)){
              F.theta.z <- c^(-1)*theta.z.grid[theta.ind]*(T-N-K)/(N*(T-K-1))
              F.ncp <- c^(-1)*theta.z.grid[theta.ind]
              
              F.grs.quantiles[theta.ind,] <- qf(c(0.05,0.5,0.95),df1 = N, df2 = T - N - K, ncp = F.ncp )
        }
      
    ### Plot the results ###  
      plot(theta.z.grid, F.grs.quantiles[,1], ylab = 'GRS F-statistic', xlab = expression(theta[z]^{2}), main = 'Confidence intervals for GRS F', type = 'l', xlim = c(0,1), ylim = c(0,max(F.grs.quantiles)))
      lines(theta.z.grid, F.grs.quantiles[,2], col = 'red')
      lines(theta.z.grid, F.grs.quantiles[,3], col = 'blue')
      lines(theta.z.grid, matrix(F.grs.lns,length(theta.z.grid),1),lty = 2)
      lines(theta.z.grid, matrix(3.49,length(theta.z.grid),1),lty = 6,col = 'gray') ### As reported in the paper, I get a different value since I am using a different market portfolio. 
      legend('topleft', c('5%-lie','Median','95%-ile','GRS F statistic (observed)','GRS F statistic (LNS)'), lty = c(1,1,1,2,6), col = c('black','red','blue','black','gray'),cex = 0.9, bg = 'transparent')
      
 
      