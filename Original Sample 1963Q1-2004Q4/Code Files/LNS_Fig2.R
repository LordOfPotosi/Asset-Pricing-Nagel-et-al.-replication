
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
n.sim <- 4000
n.fac <- 10
f <- c(1,3,5)

rsq.min <- 0
rsq.max <- 1
incr <- 0.1
rsq.grid <- as.matrix(seq(rsq.min, rsq.max, incr))

fin.dat <- merge(port.dat, fac.dat, by = 'Date')

retex.dat <- cbind((fin.dat[,2:ncol(port.dat)] - fin.dat$RF), fin.dat[,(ncol(port.dat)+1):ncol(fin.dat)]) ### Get excess returns 
fac <- as.matrix(fac.dat[,2:(K+1)]) 


Rsq.fig2 <- array(NA,c(length(rsq.grid),(n.sim*n.fac),K)) 
Rsq.fig2.noise <- array(NA,c(length(rsq.grid),(n.sim*n.fac),K)) 

### Create the C vector using the procedure described 
   for (r2.ind in 1:length(rsq.grid)){ 
     r2 <- rsq.grid[r2.ind]
     
      for (k in 1 : n.sim){
        for (i in f){
            
          
            c <- sqrt(r2/i) 
            
            G <- matrix(NA,N,i)
            R <- as.matrix(retex.dat[,1:N])
            mu.R <- as.matrix(colMeans(R))
            
            mu.S.num <- mu.R - mean(mu.R)
            mu.S.den <- sqrt(diag(var(mu.R)))
            mu.S <- mu.S.num/mu.S.den
              
              for ( j in 1 : i){
                        ej <- matrix(rnorm(N,0,1),N,1)
                        G[,j] <- c*mu.S + sqrt(1-c^2)*ej  
                   }
            
           
            
            W = solve(var(R)) %*% G   #### From G get the implied weights on the portfolio 
                 
               for (count in 1:n.fac){                                ### Sample the factors and the data n.fac times 
                      x <- 1:T 
                      sample.index <- sample(x, T, replace = TRUE)
                      R.sample <- as.matrix(retex.dat[sample.index,1:N])
                      mu.R.sample <- as.matrix(colMeans(R.sample))
                      
                      P1 <- R.sample %*% W 
                      ### Create noise and then transform it to have mean 0 and variance 3 times the variance of the factors ###
                      mu.noise <- matrix(0, i, 1)
                      sigma.noise <- 3*diag(diag(var(P1)),i,i)
                      
                      e.noise <- mvrnorm(T, mu.noise, sigma.noise)  ### This noise will not have a 0 sample mean and the sample variance as desired, hence will need to be tranformed
                      e.noise.demeaned <- t(t(e.noise) %*% (diag(T) - (matrix(1,T,1) %*% t(matrix(1,T,1)))/T)) %*% diag(sqrt(diag(sigma.noise)/diag(var(e.noise))),i,i)
                      
                      P1.noise <- P1 + e.noise.demeaned
                        
                      P <- cbind(matrix(1,T,1),P1) 
                      P.noise <- cbind(matrix(1,T,1),P1.noise) 
                      
                      C1 <- matrix(solve(t(P) %*% P) %*% t(P) %*% R.sample,(i+1),N) 
                      C1.noise <- matrix(solve(t(P.noise) %*% P.noise) %*% t(P.noise) %*% R.sample,(i+1),N) 
                      
                      C <- matrix(C1[2:(i+1),],i,N)
                      C.noise <- matrix(C1.noise[2:(i+1),],i,N)
                      
                      C.t <- cbind(matrix(1,N,1),t(C))
                      C.t.noise <- cbind(matrix(1,N,1),t(C.noise))
                      
                      M.c <- diag(1,N,N) - C.t %*% solve(t(C.t) %*% C.t) %*% t(C.t)
                      M.c.noise <- diag(1,N,N) - C.t.noise %*% solve(t(C.t.noise) %*% C.t.noise) %*% t(C.t.noise)
                      
                      Rsq.fig2[r2.ind, (n.fac*(k-1) + count),((i+1)/2)] <- 1 - (t(mu.R.sample) %*% M.c %*% mu.R.sample/(t(mu.R.sample - mean(mu.R.sample)) %*% (mu.R.sample - mean(mu.R.sample))))*(N/(N-(i+1)))
                      Rsq.fig2.noise[r2.ind, (n.fac*(k-1) + count),((i+1)/2)] <- 1 - (t(mu.R.sample) %*% M.c.noise %*% mu.R.sample/(t(mu.R.sample - mean(mu.R.sample)) %*% (mu.R.sample - mean(mu.R.sample))))*(N/(N-(i+1)))
                      
                  }
                      
                      
           }   
        }

   }

Rsq.fig2.out <- array(NA, c(length(rsq.grid),3, ncol = K))
Rsq.fig2.noise.out <- array(NA, c(length(rsq.grid),3, ncol = K))

  for (r2.ind in 1:length(rsq.grid)){
    for (i in 1:K){
      Rsq.fig2.out[r2.ind, ,i] <-  matrix(quantile(Rsq.fig2[r2.ind, ,i],probs = c(0.05,0.5,0.95),na.rm = TRUE),nrow = 3, ncol = 1)
      Rsq.fig2.noise.out[r2.ind, ,i] <-  matrix(quantile(Rsq.fig2.noise[r2.ind, ,i],probs = c(0.05,0.5,0.95),na.rm = TRUE),nrow = 3, ncol = 1)
      
      }
    }

### Plot the figures without noise ### 

par(mfrow = c(3,2))
### plot for 1 Factor ### 
plot(rsq.grid, Rsq.fig2.out[ , 1, 1], xlab = expression('True' ~ R^2 ), ylab = expression('adj' ~ R^2), main = 'One factor w/o Noise', type = 'l', col = 'red', xlim = c(0,1), ylim = c(-0.2,1))
lines(rsq.grid, matrix(0,length(rsq.grid),1))
lines(rsq.grid,Rsq.fig2.out[ , 2, 1], lty = 4, col = 'blue')
lines(rsq.grid,Rsq.fig2.out[ , 3, 1], lty = 6, col = 'green')
legend('topleft', c('5%-le','Median','95%-ile'), lty = c(2,4,6), col = c('red', 'blue','green'),cex = 0.6, bg = 'transparent')

plot(rsq.grid, Rsq.fig2.noise.out[ , 1, 1], xlab = expression('True' ~ R^2 ), ylab = expression('adj' ~ R^2), main = 'One factor w/ Noise', type = 'l', col = 'red', xlim = c(0,1), ylim = c(-0.2,1))
lines(rsq.grid, matrix(0,length(rsq.grid),1))
lines(rsq.grid,Rsq.fig2.noise.out[ , 2, 1], lty = 4, col = 'blue')
lines(rsq.grid,Rsq.fig2.noise.out[ , 3, 1], lty = 6, col = 'green')
legend('topleft', c('5%-le','Median','95%-ile'), lty = c(2,4,6), col = c('red', 'blue','green'),cex = 0.6, bg = 'transparent')

### plot for 3 Factors ### 
plot(rsq.grid, Rsq.fig2.out[ , 1, 2], xlab = expression('True' ~ R^2 ), ylab = expression('adj' ~ R^2), main = 'Three factors w/o Noise', type = 'l', col = 'red', xlim = c(0,1), ylim = c(-0.2,1))
lines(rsq.grid, matrix(0,length(rsq.grid),1))
lines(rsq.grid,Rsq.fig2.out[ , 2, 2], lty = 4, col = 'blue')
lines(rsq.grid,Rsq.fig2.out[ , 3, 2], lty = 6, col = 'green')
legend('topleft', c('5%-le','Median','95%-ile'), lty = c(2,4,6), col = c('red', 'blue','green'),cex = 0.6, bg = 'transparent')

plot(rsq.grid, Rsq.fig2.noise.out[ , 1, 2], xlab = expression('True' ~ R^2 ), ylab = expression('adj' ~ R^2), main = 'Three factors w/ Noise', type = 'l', col = 'red', xlim = c(0,1), ylim = c(-0.2,1))
lines(rsq.grid, matrix(0,length(rsq.grid),1))
lines(rsq.grid,Rsq.fig2.noise.out[ , 2, 2], lty = 4, col = 'blue')
lines(rsq.grid,Rsq.fig2.noise.out[ , 3, 2], lty = 6, col = 'green')
legend('topleft', c('5%-le','Median','95%-ile'), lty = c(2,4,6), col = c('red', 'blue','green'),cex = 0.6, bg = 'transparent')


### Plot for 5 Factors ### 
plot(rsq.grid, Rsq.fig2.out[ , 1, 3], xlab = expression('True' ~ R^2 ), ylab = expression('adj' ~ R^2), main = 'Five factors w/o Noise', type = 'l', col = 'red', xlim = c(0,1), ylim = c(-0.2,1))
lines(rsq.grid, matrix(0,length(rsq.grid),1))
lines(rsq.grid,Rsq.fig2.out[ , 2, 3], lty = 4, col = 'blue')
lines(rsq.grid,Rsq.fig2.out[ , 3, 3], lty = 6, col = 'green')
legend('topleft', c('5%-le','Median','95%-ile'), lty = c(2,4,6), col = c('red', 'blue','green'),cex = 0.6, bg = 'transparent')

plot(rsq.grid, Rsq.fig2.noise.out[ , 1, 3], xlab = expression('True' ~ R^2 ), ylab = expression('adj' ~ R^2), main = 'Five factors w/ Noise', type = 'l', col = 'red', xlim = c(0,1), ylim = c(-0.2,1))
lines(rsq.grid, matrix(0,length(rsq.grid),1))
lines(rsq.grid,Rsq.fig2.noise.out[ , 2, 3], lty = 4, col = 'blue')
lines(rsq.grid,Rsq.fig2.noise.out[ , 3, 3], lty = 6, col = 'green')
legend('topleft', c('5%-le','Median','95%-ile'), lty = c(2,4,6), col = c('red', 'blue','green'),cex = 0.6, bg = 'transparent')

