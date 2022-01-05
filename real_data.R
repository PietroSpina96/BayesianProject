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
plot(t(importMatrix(res, type = c('SL', 'sx'), position = 'au'))[1,], type = "l", ylim = c(-250, 250))
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
lambda<-eigenvf$values
eigenft<-eigenvf$vectors

# OBS: from the theory we know that only (N-1) eigenvalues are not null where N=number of 
# statistical units. Here we have N=26 and 25 not null eigenvalues.
#lambda[1:26]

# prova dell'approssimazione con una funzione f presa dal dataset.
f_prova <- f.data$ausxSL$data[1,]

x11()
plot(t(importMatrix(res, type = c('SL', 'sx'), position = 'au'))[1,], type = "l", ylim = c(-250, 250), lwd=2)
lines(f.data$ausxSL$argvals, f_alpha_approx(f_prova,0.1,lambda,eigenft), type = 'l', lwd=2, col = 'firebrick2')
lines(f.data$ausxSL$argvals, f_alpha_approx(f_prova,0.01,lambda,eigenft), type = 'l', lwd=2, col = 'blue')
lines(f.data$ausxSL$argvals, f_alpha_approx(f_prova,0.0001,lambda,eigenft), type = 'l', lwd=2, col = 'forestgreen')


#### alpha-Mahalanobis distance calculation ####
# Smoothed data
alpha <- 1e+4
f.data_alpha <- matrix(0, nrow = 26, ncol = 1600)
for (i in 1:26){
  f.data_alpha[i,] <- f_alpha_approx(f.data$ausxSL$data[i,],alpha,lambda,eigenft)
}

# alpha-Mahalanobis distance matrix 
Mahalanobis_Distance <- matrix(0, nrow = 26, ncol = 26)
for (i in 1:26){
  #print(i)
  for (j in 1:26){
    Mahalanobis_Distance[i,j] <- alpha_Mahalanobis(alpha,f.data$ausxSL$data[i,],f.data$ausxSL$data[j,],lambda,eigenft)
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


#### ORIGINAL CLUSTEING FUNCTIONS ####
f.data.clust <- fda_clustering_mahalanobis(n_clust = 2, alpha = 10000,
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


