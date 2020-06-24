rm(list = ls())

library('dplyr')
library('lmtest')
library('MASS')
library('gmm')
library('xtable')
library('aTSA')

setwd("C:/Users/mayur/Dropbox/AP Replication Project/Replication Project Code files")

cons.nd <- read.csv('PCND.csv',header = TRUE)
cons.sv <- read.csv('PCESV.csv',header = TRUE)

df.nd <- read.csv('Nondurable deflator.csv',header = TRUE)
df.sv <- read.csv('Services deflator.csv', header = TRUE)



df.full <- cbind(cons.nd, df.nd[,2], cons.sv[,2], df.sv[,2])
names(df.full) <- c('DATE','PCND','DFND','PCSV','DFSV')

df.full$nd.real <- df.full$PCND/df.full$DFND 
df.full$sv.real <- df.full$PCSV/df.full$DFSV

df.full$c.real <- df.full$nd.real + df.full$sv.real 

df.full$c.lag <- lag(df.full$c.real,1) 
df.full$c.g <- log(df.full$c.real/df.full$c.lag)

write.csv(df.full,'cons.csv')





