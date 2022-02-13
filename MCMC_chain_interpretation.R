#### ANALYSIS OF THE CHAIN
load("indexes_Ci_uniform_lambda1_alpha1e4_k2_realData.RData")
load("indexes_Ci_uniform_lambda1.RData")

library(NPflow)
library(mcclust)

df <- list()
for(ii in 1:dim(Ci)[1]){
  df[[ii]] <- Ci[ii,]
}

S=similarityMat(df)

#### PLOT OF CO-CLUSTERING MATRIX

library(ggplot2)
library(ggcorrplot)
library(latex2exp)

x11()
ggcorrplot(S, title = "Co-clustering probability", legend.title = TeX("$p_{ij}$"))

## FOR SIMULATED DATA

x <- 1:100
y <- 1:100
data <- expand.grid(X=x, Y=y)

x11()
ggplot(data, aes(X, Y, fill=as.vector(S))) + 
  geom_tile() + labs(fill=TeX("$p_{ij}$"))+ ggtitle("Co-clustering probability")

## FOR REAL DATA
x <- 1:26
y <- 1:26
data <- expand.grid(X=x, Y=y)

x11()
ggplot(data, aes(X, Y, fill=as.vector(S))) + 
  geom_tile() + labs(fill=TeX("$p_{ij}$"))+ ggtitle("Co-clustering probability")


# heatmap(S)
image(S)

#### CLUSTER ESTIMATION

## SIMULATED DATA

load('Simulated_WP.RData')
load('Functions_WP.RData')


# Choose the model for the simulation:
# simulated data - model 1
data_simul <- data1
time_simul <- time123
rm(list=c('data1','data2','data3','data4','data5'))
rm(list=c('time123','time45'))
rm(list=c('K_2','K_3'))

BL2=cluster_est_binder(data.frame(t(Ci)))

indicii = BL2$c_est

gruppo1 <- which(indicii == 1)
gruppo2 <- which(indicii == 2)

# Plot 
n <- dim(data_simul)[1]

x11()
par(mfrow = c(1,3))

plot(time_simul,data_simul[1,],type = 'l', ylim = c(-3.5,9), lwd = 2, col='black',main = "Simulated data" , ylab='Data',xlab='Time')
for(i in 2:n){
  lines(time_simul,data_simul[i,],type = 'l', col = 'black',lwd = 2)
}
for (i in 80:n){
  lines(time,data[i,],type = 'l', col = 'black', lwd = 2)
}

## SIMULATED DATA
plot(time_simul,data_simul[1,],type = 'l', ylim = c(-3.5,9), col = 'firebrick2', lwd = 2, main ="Main and contaminated processes",ylab='Data',xlab='Time')
for (i in 2:80){
  lines(time_simul,data_simul[i,],type = 'l', col = 'blue',lwd = 2)
}
for (i in 81:100){
  lines(time_simul,data_simul[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
c1 = colMeans(data_simul[1:80,])
c2 = colMeans(data_simul[81:100,])

lines(time_simul,c1,type = 'l', lwd = 3)
lines(time_simul,c2,type = 'l', lwd = 3)

#RESULTS OF BINDER ESTIMATOR

plot(time_simul,data_simul[gruppo2[1],],type = 'l', ylim = c(-3.5,9), col = 'firebrick2', lwd = 2, main = "Estimated after MCMC sampling",ylab='Data',xlab='Time')
for (i in 2:length(gruppo2)){
  lines(time_simul,data_simul[gruppo2[i],],type = 'l', col = 'blue',lwd = 2)
}
for (i in 1:length(gruppo1)){
  lines(time_simul,data_simul[gruppo1[i],],type = 'l', col ='firebrick2',lwd = 2)
}
c1_c = colMeans(data_simul[gruppo1,])
c2_c = colMeans(data_simul[gruppo2,])
lines(time_simul,c1_c,type = 'l', lwd = 3)
lines(time_simul,c2_c,type = 'l', lwd = 3)



##REAL DATA
load("functional_WP.RData")
load('Functions_WP.RData')


#### Application on real data ####
n <- 26
t_points <- 1600
time<-seq(1,t_points)
time_red <- seq(1,1600,3)
f.data_red <- f.data$ausxSL$data[,time_red]
f.Data <- list('data' = f.data_red, 'argvals' = time_red, 'clinical' = data)
rm(time_red)
rm(f.data_red)

load("indexes_Ci_uniform_lambda1_alpha1e4_k2_realData.RData")

BL2=cluster_est_binder(data.frame(t(Ci)))

indicii = BL2$c_est

gruppo1 <- which(indicii == 1)
gruppo2 <- which(indicii == 2)

x11()
par(mfrow = c(1,2))

plot(f.Data$argvals,f.Data$data[1,], ylim = range(f.Data$data) ,type = 'l', lwd = 2, main = "Clinical Data", 
     xlab = 'time', ylab = 'EVOKED POTENTIAL')
for(i in 2:n){
  lines(f.Data$argvals,f.Data$data[i,],type = 'l',lwd = 2)
}

plot(f.Data$argvals,f.Data$data[gruppo1[1],], ylim = range(f.Data$data) ,type = 'l', col = 'firebrick1', lwd = 2, 
     main = "MCMC cluster estimation for k = 2", xlab = 'time', ylab = 'EVOKED POTENTIAL')
for (i in 2:length(gruppo1)){
  lines(f.Data$argvals,f.Data$data[gruppo1[i],],type = 'l', col = 'firebrick1',lwd = 2)
}
for (i in 1:length(gruppo2)){
  lines(f.Data$argvals,f.Data$data[gruppo2[i],],type = 'l', col = 'dodgerblue1',lwd = 2)
}

c1 = colMeans(f.Data$data[gruppo1,])
c1 = colMeans(f.Data$data[gruppo2,])

lines(f.Data$argvals, c1, col='firebrick4', lwd = 3)
lines(f.Data$argvals, c2, col = 'dodgerblue4', lwd = 3)




