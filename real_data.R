#### REAL DATA #### 

#### Setup ####

library(fda.usc)
library(fda)
library(fields)

setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject")
setwd("C:/Users/admin/Documents/R/Project_BS/BayesianProject") #GiuliaR

load("functional_WP.RData")
load('Functions_WP.RData')


#### Application on real data ####
n <- 26
t_points <- 1600
time<-seq(1,t_points)

# Plot of functional mean vs data
X_bar <- colMeans(f.data$ausxSL$data) # functional mean

x11()
plot(t(importMatrix(res, type = c('SL', 'sx'), position = 'au'))[1,], type = "l", ylim = range(f.data$ausxSL$data))
for(i in 1:26){
  lines(t(importMatrix(res, type = c('SL', 'sx'), position = 'au'))[i,], lwd=1)
}
lines(f.data$ausxSL$argvals, X_bar, type = 'l', lwd=3, col = 'firebrick2')


# Covariance matrix for the data: K_hat_R is the covariance matrix computed through the R function cov
K_hat_R <- cov(f.data$ausxSL$data) 

# Plot of the covariance matrix
x11()
image.plot(time,time,K_hat_R,main='R cov function')

# Plot of the covariance matrix K_hat and K_hat_R
x11()
par(mfrow=c(1,2))
image.plot(time,time,K_hat,main='Our cov function')
image.plot(time,time,K_hat_R,main='R cov function')

# Compute eigenvalues and eigenfunctions of the covariance matrix
# eigenvf<-eigen(K_hat)
# lambda<-eigenvf$values
# eigenft<-eigenvf$vectors

# Now the eigenvalues and eigenvectors of the K_hat_R matrix are called as the previous ones
# so that we don't have to rewrite them in the code. The workspace is updated
eigenvf<-eigen(K_hat_R)
values<-eigenvf$values
eigenft<-eigenvf$vectors

# OBS: from the theory we know that only (N-1) eigenvalues are not null where N=number of 
# statistical units. Here we have N=26 and 25 not null eigenvalues.
#lambda[1:26]

# prova dell'approssimazione con una funzione f presa dal dataset.
f_prova <- f.data$ausxSL$data[1,]

x11()
plot(t(importMatrix(res, type = c('SL', 'sx'), position = 'au'))[1,], type = "l", ylim = range(f.data$ausxSL$data[1,]), lwd=2)
lines(f.data$ausxSL$argvals, f_alpha_approx(f_prova,0.1,values,eigenft), type = 'l', lwd=2, col = 'firebrick2')
lines(f.data$ausxSL$argvals, f_alpha_approx(f_prova,1e+4,values,eigenft), type = 'l', lwd=2, col = 'blue')
lines(f.data$ausxSL$argvals, f_alpha_approx(f_prova,1e+5,values,eigenft), type = 'l', lwd=2, col = 'forestgreen')


#### alpha-Mahalanobis distance calculation ####
# Smoothed data
alpha <- 1e+5
f.data_alpha <- matrix(0, nrow = 26, ncol = 1600)
for (i in 1:26){
  f.data_alpha[i,] <- f_alpha_approx(f.data$ausxSL$data[i,],alpha,values,eigenft)
}

x11()
par(mfrow = c(1,2))
plot(f.data$ausxSL$argvals,f.data$ausxSL$data[1,], ylim = range(f.data$ausxSL$data) ,
     type = 'l', lwd = 2, main = "DATA", xlab = 'time', ylab = 'EVOKED POTENTIAL')
for (i in 2:26){
  lines(f.data$ausxSL$argvals, f.data$ausxSL$data[i,], lwd = 2)
}
plot(f.data$ausxSL$argvals, f.data_alpha[1,], ylim = range(f.data$ausxSL$data), 
     type = 'l', lwd = 2, main = 'SMOOTHED DATA', , xlab = 'time', ylab = 'EVOKED POTENTIAL')
for (i in 2:26){
  lines(f.data$ausxSL$argvals, f.data_alpha[i,], lwd = 2)
}


# alpha-Mahalanobis distance matrix 
Mahalanobis_Distance <- matrix(0, nrow = 26, ncol = 26)
for (i in 1:26){
  print(i)
  for (j in 1:26){
    Mahalanobis_Distance[i,j] <- alpha_Mahalanobis(alpha,f.data$ausxSL$data[i,],f.data$ausxSL$data[j,],values,eigenft)
  }
}

x11()
image.plot(1:26,1:26,Mahalanobis_Distance)

# Useless
# # prova
# kernel_estimator(f.data$ausxSL$data,X_bar,n,1,9,t_points) 
# 
# # Covariance matrix for the data
# K_hat <- matrix(0,t_points,t_points)
# for (i in 1:t_points){
#   for (j in 1:t_points){
#     K_hat[i,j] <- kernel_estimator(f.data$ausxSL$data,X_bar,n,i,j,t_points)
#   }
# }
# # sub-matrix (5x5)
# khat_example<-K_hat[1:5,1:5]
# show(khat_example)


#### Reduction of data points ####

time_red <- seq(1,1600,3)
f.data_red <- f.data$ausxSL$data[,time_red]
f.Data <- list('data' = f.data_red, 'argvals' = time_red)
rm(time_red)
rm(f.data_red)

x11()
plot(f.Data$argvals,f.Data$data[1,], ylim = range(f.Data$data), type = 'l', lwd = 2, 
     xlab = 'time', ylab = 'EVOKED PONTENTIAL')
for (i in 2:26){
  lines(f.Data$argvals,f.Data$data[i,], lwd = 2)
}

K_hat_t <- cov(f.Data$data)
eigen_t <- eigen(K_hat)

max(eigenvf$values)
max(eigen_t$values)


#### Transformation of data ####

f.data_t <- list('data' = f.data$ausxSL$data, 'argvals' = f.data$ausxSL$argvals)
scale <- max(max(f.data_t$data),abs(min(f.data_t$data)))
f.data_t$data <- f.data_t$data/scale

x11()
plot(f.data_t$argvals,f.data_t$data[1,], ylim = c(-1.05,1.05), type = 'l', lwd = 2, 
     xlab = 'time', ylab = 'EVOKED PONTENTIAL')
for (i in 2:26){
  lines(f.data_t$argvals,f.data_t$data[i,], lwd = 2)
}

K_hat_t <- cov(f.data_t$data)
eigen_t <- eigen(K_hat_R)
lambda_t <- eigenvf$values
max(lambda_t)
eigenf_t <- eigenvf$vectors

max(eigen_t$values) == max(eigenvf$values)

#### Uniform prior with fixed covariance ####

##### k = 2 ######

f.data_clust_2 <- fda_clustering_mahalanobis_general(n_clust = 2, alpha = alpha,
                                                     cov_matrix = cov(f.Data$data),
                                                     cov_type = 'fixed', toll = 1e-2, 
                                                     data = f.Data$data)
c_opt_2 <- f.data_clust_2$label
show(c_opt_2)
show(f.data_clust_2$loss)

c1 <- f.data_clust_2$centroids[1,]
c2 <- f.data_clust_2$centroids[2,]

data1 <- f.Data$data[which(c_opt_2=='1'),]
rownames(data1) <- c(1:dim(data1)[1])
data2 <- f.Data$data[which(c_opt_2=='2'),]
rownames(data2) <- c(1:dim(data2)[1])

###### Plot #####
x11()
par(mfrow = c(1,2))
plot(f.Data$argvals,f.Data$data[1,], ylim = range(f.Data$data) ,type = 'l', lwd = 2, main = "Data", 
     xlab = 'time', ylab = 'EVOKED POTENTIAL')
for(i in 2:n){
  lines(f.Data$argvals,f.Data$data[i,],type = 'l',lwd = 2)
}

plot(f.Data$argvals,data1[1,], ylim = range(f.Data$data) ,type = 'l', col = 'firebrick2', lwd = 2, 
     main = "Uniform (fixed cov) k = 2", xlab = 'time', ylab = 'EVOKED POTENTIAL')
for (i in 2:dim(data1)[1]){
  lines(f.Data$argvals,data1[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in 1:dim(data2)[1]){
  lines(f.Data$argvals,data2[i,],type = 'l', col = 'dodgerblue3',lwd = 2)
}
lines(f.Data$argvals, c1, lwd = 3)
lines(f.Data$argvals, c2, lwd = 3)


##### k = 3 ######
f.data_clust_3 <- fda_clustering_mahalanobis_general(n_clust = 3, alpha = alpha,
                                                     cov_matrix = cov(f.Data$data),
                                                     cov_type = 'fixed', toll = 1e-2, 
                                                     data = f.Data$data)
c_opt_3 <- f.data_clust_3$label
show(c_opt_3)
show(f.data_clust_3$loss)

c1 <- f.data_clust_3$centroids[1,]
c2 <- f.data_clust_3$centroids[2,]
c3 <- f.data_clust_3$centroids[3,]

data1 <- f.Data$data[which(c_opt_3=='1'),]
data2 <- f.Data$data[which(c_opt_3=='2'),]
data3 <- f.Data$data[which(c_opt_3=='3'),]

###### Plot #####
x11()
par(mfrow = c(1,2))
plot(f.Data$argvals,f.Data$data[1,], ylim = range(f.Data$data) ,type = 'l', lwd = 2, main = "Data",
     xlab = 'time', ylab = 'EVOKED POTENTIAL')
for(i in 2:n){
  lines(f.Data$argvals,f.Data$data[i,],type = 'l',lwd = 2)
}

plot(f.Data$argvals,data1[1,], ylim = range(f.Data$data) ,type = 'l', col = 'dodgerblue3', lwd = 2, 
     main = "Uniform (fixed cov) k = 3", xlab = 'time', ylab = 'EVOKED POTENTIAL')
# plot(f.Data$argvals,data1[1,], ylim = range(f.Data$data) ,type = 'l', col = 'gold', lwd = 2, main = "Uniform (fixed cov) k = 3")
for (i in 2:dim(data1)[1]){
  lines(f.Data$argvals,data1[i,],type = 'l', col = 'dodgerblue3',lwd = 2)
}
for (i in 1:dim(data2)[1]){
  lines(f.Data$argvals,data2[i,],type = 'l', col = 'forestgreen',lwd = 2)
}
for (i in 1:dim(data3)[1]){
  lines(f.Data$argvals,data3[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
lines(f.Data$argvals, c1, lwd = 3)
lines(f.Data$argvals, c2, lwd = 3)
lines(f.Data$argvals, c3, lwd = 3)


##### Cluster merging #####
library(ggplot2)
library(dplyr) # pipe (%>%)
library(tidyr) # gather()
library(tibble) # add_column()
library(gridExtra) # grid.arrange()
require(gtools) # combinations()








#### Uniform prior with covariance updating within clusters ####

##### k = 2 ######

f.data_clust_2up <- fda_clustering_mahalanobis_general(n_clust = 2, alpha = alpha,
                                                     cov_matrix = cov(f.Data$data),
                                                     cov_type = 'updated', toll = 1e-1, 
                                                     data = f.Data$data)
c_opt_2up <- f.data_clust_2up$label
show(c_opt_2up)
show(f.data_clust_2up$loss)

c1 <- f.data_clust_2up$centroids[1,]
c2 <- f.data_clust_2up$centroids[2,]

data1 <- f.Data$data[which(c_opt_2up=='1'),]
data2 <- f.Data$data[which(c_opt_2up=='2'),]

###### Plot #####
x11()
par(mfrow = c(1,2))
plot(f.Data$argvals,f.Data$data[1,], ylim = range(f.Data$data) ,type = 'l', lwd = 2, main = "Data", 
     xlab = 'time', ylab = 'EVOKED POTENTIAL')
for(i in 2:n){
  lines(f.Data$argvals,f.Data$data[i,],type = 'l',lwd = 2)
}

plot(f.Data$argvals,data1[1,], ylim = range(f.Data$data) ,type = 'l', col = 'firebrick2', lwd = 2, 
     main = "Uniform (updated cov) k = 2", xlab = 'time', ylab = 'EVOKED POTENTIAL')
for (i in 2:dim(data1)[1]){
  lines(f.Data$argvals,data1[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in 1:dim(data2)[1]){
  lines(f.Data$argvals,data2[i,],type = 'l', col = 'dodgerblue3',lwd = 2)
}
lines(f.Data$argvals, c1, lwd = 3)
lines(f.Data$argvals, c2, lwd = 3)

##### k = 3 ######
f.data_clust_3up <- fda_clustering_mahalanobis_general(n_clust = 3, alpha = alpha,
                                                       cov_matrix = cov(f.Data$data),
                                                       cov_type = 'updated', toll = 1e-1, 
                                                       data = f.Data$data)
c_opt_3 <- f.data_clust_3up$label
show(c_opt_3)
show(f.data_clust_3up$loss)

c1 <- f.data_clust_3up$centroids[1,]
c2 <- f.data_clust_3up$centroids[2,]
c3 <- f.data_clust_3up$centroids[3,]

data1 <- f.Data$data[which(c_opt_3=='1'),]
data2 <- f.Data$data[which(c_opt_3=='2'),]
data3 <- f.Data$data[which(c_opt_3=='3'),]

###### Plot #####
x11()
par(mfrow = c(1,2))
plot(f.Data$argvals,f.Data$data[1,], ylim = range(f.Data$data) ,type = 'l', lwd = 2, main = "Data",
     xlab = 'time', ylab = 'EVOKED POTENTIAL')
for(i in 2:n){
  lines(f.Data$argvals,f.Data$data[i,],type = 'l',lwd = 2)
}

plot(f.Data$argvals,data1, ylim = range(f.Data$data) ,type = 'l', col = 'dodgerblue3', lwd = 3, 
     main = "Uniform (updated cov) k = 3", xlab = 'time', ylab = 'EVOKED POTENTIAL')
# plot(f.Data$argvals,data1[1,], ylim = range(f.Data$data) ,type = 'l', col = 'dodgerblue3', lwd = 2, 
#      main = "Uniform (updated cov) k = 3", xlab = 'time', ylab = 'EVOKED POTENTIAL')
# for (i in 2:dim(data1)[1]){
#   lines(f.Data$argvals,data1[i,],type = 'l', col = 'dodgerblue3',lwd = 2)
# }
for (i in 1:dim(data2)[1]){
  lines(f.Data$argvals,data2[i,],type = 'l', col = 'forestgreen',lwd = 2)
}
for (i in 1:dim(data3)[1]){
  lines(f.Data$argvals,data3[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
# lines(f.Data$argvals, c1, lwd = 3)
lines(f.Data$argvals, c2, lwd = 3)
lines(f.Data$argvals, c3, lwd = 3)








#### Pitman-Yor EPPF Prior ####

alpha <- alpha
sigma <- 0.25
theta <- 3.66
lambda <- 0.75

##### k = 2 ####
f.data_clust_py_2 <- fda_clustering_pitmanyor(n_clust = 2, alpha, sigma, theta,
                                       lambda, cov_matrix = cov(f.data$ausxSL$data), 
                                       toll = 1e-10, data = f.data$ausxSL$data)

c_opt_py_2 <- f.data_clust_py_2$label
show(c_opt_py_2)
show(f.data_clust_py_2$posterior)

c1 <- f.data_clust_py_2$centroids[1,]
c2 <- f.data_clust_py_2$centroids[2,]

data1 <- f.data$ausxSL$data[which(c_opt_py_2=='1'),]
data2 <- f.data$ausxSL$data[which(c_opt_py_2=='2'),]


###### Plot #####

x11()
par(mfrow = c(1,2))
plot(time,f.data$ausxSL$data[1,], ylim = range(f.data$ausxSL$data) ,type = 'l', lwd = 2, main = "Data")
for(i in 2:n){
  lines(time,f.data$ausxSL$data[i,],type = 'l',lwd = 2)
}

plot(time,data1[1,], ylim = range(f.data$ausxSL$data) ,type = 'l', col = 'gold', lwd = 2, main = "Pitman-Yor k = 2")
for (i in 2:dim(data1)[1]){
  lines(time,data1[i,],type = 'l', col = 'gold',lwd = 2)
}
for (i in 1:dim(data2)[1]){
  lines(time,data2[i,],type = 'l', col = 'forestgreen',lwd = 2)
}



##### k = 3 #####

f.data_clust_py_3 <- fda_clustering_pitmanyor(n_clust = 3, alpha, sigma, theta,
                                              lambda, cov_matrix = cov(f.data$ausxSL$data), 
                                              toll = 1e-10, data = f.data$ausxSL$data)
c_opt_py_3 <- f.data_clust_py_3$label
show(c_opt_py_3)
show(f.data_clust_py_3$posterior)



c1 <- f.data_clust_py_3$centroids[1,]
c2 <- f.data_clust_py_3$centroids[2,]
c3 <- f.data_clust_py_3$centroids[3,]

data1 <- f.data$ausxSL$data[which(c_opt_py_3=='1'),]
data2 <- f.data$ausxSL$data[which(c_opt_py_3=='2'),]
data3 <- f.data$ausxSL$data[which(c_opt_py_3=='3'),]

###### Plot #####
x11()
par(mfrow = c(1,2))
plot(time,f.data$ausxSL$data[1,], ylim = range(f.data$ausxSL$data) ,type = 'l', lwd = 2, main = "Data")
for(i in 2:n){
  lines(time,f.data$ausxSL$data[i,],type = 'l',lwd = 2)
}

plot(time,data1[1,], ylim = range(f.data$ausxSL$data) ,type = 'l', col = 'gold', lwd = 2, main = "Pitman-Yor k = 3")
for (i in 2:dim(data1)[1]){
  lines(time,data1[i,],type = 'l', col = 'gold',lwd = 2)
}
# for (i in 1:dim(data2)[1]){
#   lines(time,data2[i,],type = 'l', col = 'forestgreen',lwd = 2)
# }
lines(time,data2,type = 'l', col = 'forestgreen',lwd = 2)
for (i in 1:dim(data3)[1]){
  lines(time,data3[i,],type = 'l', col = 'firebrick3',lwd = 2)
}





























#### CLUSTERING FUNCTION UPDATING COVARIANCE WITHIN CLUSTERS ####
f.data.clust <- fda_clustering_mahalanobis_updated(n_clust = 2, alpha = 10000,
                                                   cov_matrix = cov(f.data$ausxSL$data),
                                                   toll = 1e-2, data = f.data$ausxSL$data)

c_opt <- f.data.clust$label
show(c_opt)
show(f.data.clust$loss)

c1 <- f.data.clust$centroids[1,]
c2 <- f.data.clust$centroids[2,]
#c3 <- f.data.clust$centroids[3,]
#c4 <- clust$centroids[4,]


data1 <- f.data$ausxSL$data[which(c_opt=='1'),]
data2 <- f.data$ausxSL$data[which(c_opt=='2'),]
#data3 <- f.data$ausxSL$data[which(c_opt=='3'),]

x11()
par(mfrow = c(1,2))
plot(time,f.data$ausxSL$data[1,], ylim = range(f.data$ausxSL$data) ,type = 'l', lwd = 2, main = "Data")
for(i in 2:n){
  lines(time,f.data$ausxSL$data[i,],type = 'l',lwd = 2)
}

plot(time,data1[1,], ylim = range(f.data$ausxSL$data) ,type = 'l', col = 'gold', lwd = 2, main = "Clustered data")
for (i in 2:dim(data1)[1]){
  lines(time,data1[i,],type = 'l', col = 'gold',lwd = 2)
}
for (i in 1:dim(data2)[1]){
  lines(time,data2[i,],type = 'l', col = 'forestgreen',lwd = 2)
}
for (i in 1:dim(data3)[1]){
  lines(time,data3[i,],type = 'l', col = 'gold',lwd = 2)
}


##### Save Workspace ####
setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject")
save.image("functional_WP.RData")

save.image("~/R/Project_BS/BayesianProject/functional_WP.RData") #GiuliaR


