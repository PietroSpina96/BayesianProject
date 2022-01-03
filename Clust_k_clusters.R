#### CLUSTERING ON SIMULATED DATA: Cmap ####

setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject") #Pietro
#setwd("C:/Users/admin/Documents/R/Project_BS/BayesianProject") #GiuliaR

load('Simulated_WP.RData')


# Choose the model for the simulation:
# simulated data - model 1
data <- data1
eig <- eig_1
time <- time123
rm(list=c('data1','data2','data3','data4'))
rm(list=c('eig_1','eig_2','eig_3'))
rm(list=c('time123','time4'))

# simulated data - model 2
data <- data2
eig <- eig_2
time <- time123
rm(list=c('data1','data2','data3','data4'))
rm(list=c('eig_1','eig_2','eig_3'))
rm(list=c('time123','time4'))

# simulated data - model 3
data <- data3
eig <- eig_3
time <- time123
rm(list=c('data1','data2','data3','data4'))
rm(list=c('eig_1','eig_2','eig_3'))
rm(list=c('time123','time4'))

# simulated data - model 4. 
# The covariance matrix of data4 is cov_4 and K_4_1 and K_4_2 are the covariance matrices of the two clusters
# No eig has been computed 
data <- data4
time <- time4
rm(list=c('data1','data2','data3','data4'))
rm(list=c('eig_1','eig_2','eig_3'))
rm(list=c('time123','time4'))



#### LOSS FUNCTION ####
# centroids is a matrix with nrow = n_clust and ncol = t_points
# eig is the list of eigenvalues and eigenvectors

gibbs_loss <- function(n_clust, centroids, label , eig, data){
  n <- dim(data)[1]
  res = rep(0,n_clust)
  sum_partial <- 0
  
  for (k in 1:n_clust){
    for (i in 1:n){
      if (label[i] == k){
        sum_partial = alpha_Mahalanobis(alpha,data[i,],centroids[k,],eig$values,eig$vectors)
        res[k] = res[k] + sum_partial
      }
    }
  }
  
  tot = sum(res)
  return(tot)
}


#### ORIGINAL CLUSTERING FUNCTION ####
# The function works with n_clust=k
# alpha is the smoothing parameter of the data
# toll is the tolerance for the while loop
# cov_matrix is the covariance matrix of the data

fda_clustering_mahalanobis <- function(n_clust, alpha, cov_matrix, toll,data){
  
  t_points <- dim(data)[2]
  n <- dim(data)[1]
  
  # index of each centroid randomly defined through sampling
  y0 <- rep(0,n_clust)
  vect_sample <- 1:n
  
  y0[1] <- sample(vect_sample,1)
  
  for (k in 2:n_clust) {
    value <- y0[k-1]
    
    for (i in 1:length(vect_sample)){
      if (vect_sample[i] == value)
        t = i
    }
    
    vect_sample <- vect_sample[-t]
    y0[k] <- sample(vect_sample,1)
  }
  
  # vector of labels
  c_lab <- rep(0,n)
  
  # eigenvalues and eigenfunctions for the alpha-mahalanobis function
  eig <- eigen(cov_matrix)
  values <- eig$values 
  vectors <- eig$vectors
  
  Mahalanobis_Distance <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n){
    for (j in 1:n){
      Mahalanobis_Distance[i,j] <- alpha_Mahalanobis(alpha,data[i,],data[j,],values,vectors)
    }
  }
  
  # i-th unit belongs to cluster_k if the distance(centroid_k,i-th unit) is the smallest one
  Maha_dis <- matrix(0,nrow=n, ncol=n_clust)
  for (i in 1:n){
    for (k in 1:n_clust) {
      Maha_dis[i,k] <- Mahalanobis_Distance[i,y0[k]]
    }
    index <- which.min(Maha_dis[i,])
    c_lab[i] <- index
  }
  
  # define the matrix of the centroids (random centroids)
  centroids_random <- matrix(0,nrow = n_clust,ncol = t_points)
  for (k in 1:n_clust){
    centroids_random[k,] <- data[y0[k],]
  }
  
  loss_value1 <- gibbs_loss(n_clust = n_clust, centroids = centroids_random, label = c_lab, eig = eig, data = data)
  
  # update each centroid as the mean of the clusters data
  centroids_mean<-matrix(0,nrow = n_clust, ncol = t_points)
  for (k in 1:n_clust){
    centroids_mean[k,] <- colMeans(data[which(c_lab==k),])
  }
  
  loss_value2 <- gibbs_loss(n_clust = n_clust, centroids = centroids_mean, label = c_lab, eig = eig, data = data)
  
  while(abs(loss_value1 - loss_value2) >= toll){
    
    c_lab <- rep(0,n)
    
    loss_value1 <- loss_value2
    
    Maha_dis_k <- matrix(0,nrow=n, ncol=n_clust)
    for (i in 1:n){
      for (k in 1:n_clust) {
        Maha_dis_k[i,k] <- alpha_Mahalanobis(alpha,centroids_mean[k,],data[i,],eig$values,eig$vectors)
      }
      index <- which.min(Maha_dis_k[i,])
      c_lab[i] <- index
    }
    
    
    for (k in 1:n_clust){
      centroids_mean[k,] <- colMeans(data[which(c_lab==k),])
    }
    
    loss_value2 <- gibbs_loss(n_clust = n_clust, centroids = centroids_mean, label = c_lab, eig = eig, data = data)
  }
  
  return(list("label" = c_lab, "centroids" = centroids_mean, "loss" = loss_value2))
  
}


##### Application on the simulated data ####
n <- dim(data)[1]
t_points <- dim(data)[2]

k <- 2
alpha <- 0.1

clust <- fda_clustering_mahalanobis(n_clust = k, alpha = alpha, cov_matrix = K_1, toll = 1e-2,  data = data)
c_opt <- clust$label
show(c_opt)  #label switching 

c1 <- clust$centroids[1,]
c2 <- clust$centroids[2,]
c3 <- clust$centroids[3,]
#c4 <- clust$centroids[4,]


# Theoretical optimal plot vs clustering plot
data1 <- data[which(c_opt=='1'),]
data2 <- data[which(c_opt=='2'),]
data3 <- data[which(c_opt=='3'),]
#data4 <- data[which(c_opt=='4'),]

# Plot model 1
x11()
par(mfrow = c(1,2))
plot(time,data[1,],type = 'l', ylim = c(-3.5,9), lwd = 2, main = "Main and contaminated processes")
for(i in 2:(n-c)){
  lines(time,data[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in (n-c+1):n){
  lines(time,data[i,],type = 'l', col = 'blue', lwd = 2)
}

plot(time,data1[1,],type = 'l', ylim = c(-3.5,9), col = 'firebrick2', lwd = 2, main = "Clustered data")
for (i in 2:dim(data1)[1]){
  lines(time,data1[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in 1:dim(data2)[1]){
  lines(time,data2[i,],type = 'l', col = 'blue',lwd = 2)
}
# for (i in 1:dim(data3)[1]){
#  lines(time,data3[i,],type = 'l', col = 'forestgreen',lwd = 2)
# }
# for (i in 1:dim(data4)[1]){
#  lines(time,data4[i,],type = 'l', col = 'lightcyan',lwd = 2)
# }

lines(time,c1,type = 'l', lwd = 3)
lines(time,c2,type = 'l', lwd = 3)
# lines(time,c3,type = 'l', lwd = 3)
# lines(time,c4,type = 'l', lwd = 3)

rm(data1)
rm(data2)
# rm(data3)
# rm(data4)



# ##### Theoretical optimal plot vs clustering plot SMOOTHED ####
# data1 <- f.data_alpha_sim[which(c_opt=='1'),]
# data2 <- f.data_alpha_sim[which(c_opt=='2'),]
# 
# x11()
# par(mfrow = c(1,2))
# 
# plot(time,f.data_alpha_sim[1,],type = 'l', ylim = c(-3,7.5), col = 'firebrick2', lwd = 2, main = "Smooth processes")
# for(i in 2:(n-c)){
#   lines(time,f.data_alpha_sim[i,],type = 'l', col = 'firebrick2',lwd = 2)
# }
# for (i in (n-c+1):n){
#   lines(time,f.data_alpha_sim[i,],type = 'l', col = 'blue', lwd = 2)
# }
# 
# plot(time,data1[1,],type = 'l', ylim = c(-3,7.5), col = 'firebrick2', lwd = 2, main = "Clustered smoothed data")
# for (i in 2:dim(data1)[1]){
#   lines(time,data1[i,],type = 'l', col = 'firebrick2',lwd = 2)
# }
# for (i in 1:dim(data2)[1]){
#   lines(time,data2[i,],type = 'l', col = 'blue',lwd = 2)
# }
# lines(time,f_alpha_approx(c1,alpha,eig$values,eig$vectors),type = 'l', lwd = 4)
# lines(time,f_alpha_approx(c2,alpha,eig$values,eig$vectors),type = 'l', lwd = 4)
# 
# rm(data1)
# rm(data2)


#### WARNING CLUSTERING FUNCTION ####
fda_clustering_mahalanobis_warning <- function(n_clust, alpha, cov_matrix, toll,data){
  
  t_points <- dim(data)[2]
  n <- dim(data)[1]
  
  # index of each centroid randomly defined through sampling
  y0 <- rep(0,n_clust)
  vect_sample <- 1:n
  
  y0[1] <- sample(vect_sample,1)
  
  for (k in 2:n_clust) {
    value <- y0[k-1]
    
    for (i in 1:length(vect_sample)){
      if (vect_sample[i] == value)
        t = i
    }
    
    vect_sample <- vect_sample[-t]
    y0[k] <- sample(vect_sample,1)
  }
  
  # vector of labels
  c_lab <- rep(0,n)
  
  # eigenvalues and eigenfunctions for the alpha-mahalanobis function
  eig <- eigen(cov_matrix)
  values <- eig$values 
  vectors <- eig$vectors
  
  Mahalanobis_Distance <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n){
    for (j in 1:n){
      Mahalanobis_Distance[i,j] <- alpha_Mahalanobis(alpha,data[i,],data[j,],values,vectors)
    }
  }
  
  # i-th unit belongs to cluster_k if the distance(centroid_k,i-th unit) is the smallest one
  Maha_dis <- matrix(0,nrow=n, ncol=n_clust)
  for (i in 1:n){
    for (k in 1:n_clust) {
      Maha_dis[i,k] <- Mahalanobis_Distance[i,y0[k]]
    }
    index <- which.min(Maha_dis[i,])
    c_lab[i] <- index
  }
  
  # define the matrix of the centroids (random centroids)
  centroids_random <- matrix(0,nrow = n_clust,ncol = t_points)
  for (k in 1:n_clust){
    centroids_random[k,] <- data[y0[k],]
  }
  
  loss_value1 <- gibbs_loss(n_clust = n_clust, centroids = centroids_random, label = c_lab, eig = eig, data = data)
  
  # update each centroid as the mean of the clusters data
  centroids_mean<-matrix(0,nrow = n_clust, ncol = t_points)
  for (k in 1:n_clust){
    centroids_mean[k,] <- colMeans(data[which(c_lab==k),])
  }
  
  loss_value2 <- gibbs_loss(n_clust = n_clust, centroids = centroids_mean, label = c_lab, eig = eig, data = data)
  
  while(abs(loss_value1 - loss_value2) >= toll){
    
    c_lab <- rep(0,n)
    
    loss_value1 <- loss_value2
    
    Maha_dis_k <- matrix(0,nrow=n, ncol=n_clust)
    for (i in 1:n){
      for (k in 1:n_clust) {
        Maha_dis_k[i,k] <- alpha_Mahalanobis(alpha,centroids_mean[k,],data[i,],eig$values,eig$vectors)
      }
      index <- which.min(Maha_dis_k[i,])
      c_lab[i] <- index
    }
    
    
    for (k in 1:n_clust){
      centroids_mean[k,] <- colMeans(data[which(c_lab==k),])
    }
    
    loss_value2 <- gibbs_loss(n_clust = n_clust, centroids = centroids_mean, label = c_lab, eig = eig, data = data)
  }
  
  flag <- 0
  for (k in 1:(n_clust-1)){
    for (j in (k+1):n_clust){
      
      diff_centroids <- centroids_mean[k,] - centroids_mean[j,]
      dis_centroids <- norm(as.matrix(diff_centroids),type = 'i')
      
      eps <- # funzione per eps
      if (dis_centroids <= eps )
        flag <- flag + 1
    }
    
  }
  
  if (flag > 0)
    print ("WARNING: the distance between some centroids is very small. Use a smaller K ")
  
  return(list("label" = c_lab, "centroids" = centroids_mean, "loss" = loss_value2))
  
}

# Application on the data
n <- dim(data)[1]
t_points <- dim(data)[2]

k <- 3
alpha <- 0.1

clust <- fda_clustering_mahalanobis_warning(n_clust = k, alpha = alpha, cov_matrix = K_1, toll = 1e-2,  data = data)
c_opt <- clust$label
show(c_opt)  #label switching 

c1 <- clust$centroids[1,]
c2 <- clust$centroids[2,]
c3 <- clust$centroids[3,]


# Theoretical optimal plot vs clustering plot
data1 <- data[which(c_opt=='1'),]
data2 <- data[which(c_opt=='2'),]
data3 <- data[which(c_opt=='3'),]

# Plot model 1
x11()
par(mfrow = c(1,2))
plot(time,data[1,],type = 'l', ylim = c(-3.5,9), lwd = 2, main = "Main and contaminated processes")
for(i in 2:(n-c)){
  lines(time,data[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in (n-c+1):n){
  lines(time,data[i,],type = 'l', col = 'blue', lwd = 2)
}

plot(time,data1[1,],type = 'l', ylim = c(-3.5,9), col = 'firebrick2', lwd = 2, main = "Clustered data")
for (i in 2:dim(data1)[1]){
  lines(time,data1[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in 1:dim(data2)[1]){
  lines(time,data2[i,],type = 'l', col = 'blue',lwd = 2)
}
for (i in 1:dim(data3)[1]){
 lines(time,data3[i,],type = 'l', col = 'forestgreen',lwd = 2)
}

lines(time,c1,type = 'l', lwd = 3)
lines(time,c2,type = 'l', lwd = 3)
lines(time,c3,type = 'l', lwd = 3)

rm(data1)
rm(data2)
rm(data3)



#### CLUSTERING FUNCTION MERGING USING LOSS MINIMIZATION ####
# The parameter eig corresponds to the output of the eigen function (list of eigenvalues and eigenvectors)
# The function works with n_clust=k
# alpha is the smoothing parameter
# toll is the tolerance for the while loop

fda_clustering_mahalanobis_merge <- function(n_clust, alpha, eig, toll,data){
   
   n <- dim(data)[1] 
   t_points <- dim(data)[2]
   
   # index of each centroid randomly defined through sampling
   y0 <- rep(0,n_clust)
   vect_sample <- 1:n
   
   y0[1] <- sample(vect_sample,1)
   
   for (k in 2:n_clust) {
     value <- y0[k-1]
     
     for (i in 1:length(vect_sample)){
       if (vect_sample[i] == value)
         t = i
     }
     
     vect_sample <- vect_sample[-t]
     y0[k] <- sample(vect_sample,1)
   }
   
   # matrix of labels
   c_lab_mat <- matrix(0, nrow = n_clust, ncol = n)
   
   # loss vector
   loss_vec <- rep(0,n_clust)
   
   # eigenvalues and eigenfunctions for the alpha-mahalanobis function
   values <- eig$values
   vectors <- eig$vectors
   
   Mahalanobis_Distance <- matrix(0, nrow = n, ncol = n)
   for (i in 1:n){
     for (j in 1:n){
       Mahalanobis_Distance[i,j] <- alpha_Mahalanobis(alpha,data[i,],data[j,],values,vectors)
     }
   }
   
   for (l in 1:n_clust){
     
     # vector of labels
     c_lab <- rep(0,n)
     
     # starting centroids
     y00 <- rep(0,l)
     for (x in 1:l){
       y00[x] <- y0[x]
     }
    
     # i-th unit belongs to cluster_k if the distance(centroid_k,i-th unit) is the smallest one
     Maha_dis <- matrix(0,nrow=n, ncol=l)
     for (i in 1:n){
       for (k in 1:l) {
         Maha_dis[i,k] <- Mahalanobis_Distance[i,y00[k]]
       }
       index <-which.min(Maha_dis[i,])
       c_lab[i] <- index
     }
     
     # define the matrix of the centroids (random centroids)
     centroids_random <- matrix(0,nrow = l,ncol = t_points )
     for (k in 1:l){
       centroids_random[k,] <- data[y00[k],]
     }
     
     loss_value1 <- gibbs_loss(n_clust = l, centroids = centroids_random, label = c_lab, data = data)
     
     # update each centroid as the mean of the clusters data
     centroids_mean<-matrix(0,nrow = l, ncol = t_points )
     for (k in 1:l){
       centroids_mean[k,] <- colMeans(data[which(c_lab==k),])
     }
      
     loss_value2 <- gibbs_loss(n_clust = l, centroids = centroids_mean, label = c_lab, data = data)
      
       Maha_dis_k <- matrix(0,nrow=n, ncol=l)
       for (i in 1:n ){
         for (k in 1:l) {
           Maha_dis_k[i,k] <- alpha_Mahalanobis(alpha,centroids_mean[k,],data[i,],values,vectors)
         }
         index2 <- which.min(Maha_dis_k[i,])
         c_lab[i] <- index2
       }
       
       loss_value1 <- loss_value2
       
       for (k in 1:l){
         centroids_mean[k,] <- colMeans(data[which(c_lab==k),])
       }
       
       loss_value2 <- gibbs_loss(n_clust = l, centroids = centroids_mean, label = c_lab, data = data)
       
       loss_vec[l] <- loss_value2
       c_lab_mat[l,] <- c_lab
       
       if (l == 1){
         centr_mat <- centroids_mean
       } else {centr_mat <- rbind(centr_mat,centroids_mean)}
     } 
     
     loss <- min(loss_vec)
     index_loss <- match(loss,loss_vec)
     
     if (index_loss == 1){
       a <- 1
     } else {a <- sum(1:(index_loss-1)) + 1}
     
     if (index_loss == 1){
       b <- 1
     } else {b <- sum(1:index_loss)}
     
     return(list("label" = c_lab_mat[index_loss,], "centroids" = centr_mat[a:b,], "loss" = loss, "loss_vector" = loss_vec))
     
   } 

# Application on simulated data
n <- dim(data)[1]
t_points <- dim(data)[2]

k <- 3
alpha <- 0.1

clust <- fda_clustering_mahalanobis_merge(n_clust = k, alpha = alpha, eig = eig, toll = 1e-2,  data = data)
c_opt <- clust$label
loss_vector <- clust$loss_vector
show(c_opt)
show(loss_vector)

c1 <- clust$centroids[1,]
c2 <- clust$centroids[2,]
c3 <- clust$centroids[3,]

data1 <- data[which(c_opt=='1'),]
data2 <- data[which(c_opt=='2'),]
data3 <- data[which(c_opt=='3'),]


x11()
par(mfrow = c(1,2))
plot(time,data[1,],type = 'l', ylim = c(-3.5,9), lwd = 2, main = "Main and contaminated processes")
for(i in 2:(n-c)){
  lines(time,data[i,],type = 'l', lwd = 2)
}
for (i in (n-c+1):n){
  lines(time,data[i,],type = 'l', lwd = 2)
}

plot(time,data1[1,],type = 'l', ylim = c(-3.5,9), col = 'firebrick2', lwd = 2, main = "Clustered data")
for (i in 2:dim(data1)[1]){
  lines(time,data1[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in 1:dim(data2)[1]){
  lines(time,data2[i,],type = 'l', col = 'forestgreen',lwd = 2)
}
for (i in 1:dim(data3)[1]){
  lines(time,data3[i,],type = 'l', col = 'blue',lwd = 2)
}

lines(time,c1,type = 'l', lwd = 3)
lines(time,c2,type = 'l', lwd = 3)
lines(time,c3,type = 'l', lwd = 3)

rm(data1)
rm(data2)
rm(data3)


#### CLUSTERING FUNCTION WITH MERGING OF CLUSTERS ####
# eps is the tolerance for the union of clusters (eps = 0 -> original function)
#     automatically set based on mean centroids distances
# also outputs the comparison plot, but sometimes applies wrong colors
#TODO: fix messed-up colors in some random cases
library(ggplot2)
library(dplyr) # pipe (%>%)
library(tidyr) # gather()
library(tibble) # add_column()
library(gridExtra) # grid.arrange()
require(gtools) # combinations()
# library(viridis) # color palette

fda_clustering_mahalanobis_union <- function(n_clust, alpha, eig, toll, eps=0, data){
  max_clust <- n_clust
  print(sprintf(" ** CLUSTERING: up to k=%d clusters ** ",max_clust))
  n <- dim(data)[1]
  t_points <- dim(data)[2]
  # index of each centroid randomly defined through sampling
  y0 <- rep(0,max_clust)
  vect_sample <- 1:n
  y0[1] <- sample(vect_sample,1)
  
  for (k in 2:max_clust) {
    value <- y0[k-1]
    for (i in 1:length(vect_sample))
      if (vect_sample[i] == value)
        t = i
      vect_sample <- vect_sample[-t]
      y0[k] <- sample(vect_sample,1)
  }
  
  # vector of labels
  print("Calculating a-Mahalanobis distances...")
  c_lab <- rep(0,n)
  #   and eigenfunctions for the alpha-mahalanobis function
  values <- eig$values
  vectors <- eig$vectors
  Mahalanobis_Distance <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n)
    for (j in 1:n)
      Mahalanobis_Distance[i,j] <- alpha_Mahalanobis(alpha,data[i,],data[j,],values,vectors)
  
  # i-th unit belongs to cluster_k if the distance(centroid_k,i-th unit) is the smallest one
  print("Assigning clusters...")
  Maha_dis <- matrix(0,nrow=n, ncol=max_clust)
  for (i in 1:n){
    for (k in 1:max_clust) 
      Maha_dis[i,k] <- Mahalanobis_Distance[i,y0[k]]
    index <- which.min(Maha_dis[i,])
    c_lab[i] <- index
  }
  
  # define the matrix of the centroids (random centroids)
  print("Compute centroids...")
  centroids_random <- matrix(0,nrow = max_clust,ncol = t_points)
  for (k in 1:max_clust)
    centroids_random[k,] <- data[y0[k],]
  loss_value1 <- gibbs_loss(n_clust = max_clust, centroids = centroids_random, 
                            label = c_lab, data = data)
  # update each centroid as the mean of the clusters data
  centroids_mean<-matrix(0,nrow = max_clust, ncol = t_points)
  for (k in 1:n_clust)
    centroids_mean[k,] <- colMeans(data[which(c_lab==k),])
  loss_value2 <- gibbs_loss(n_clust = max_clust, centroids = centroids_mean, 
                            label = c_lab, data = data)
  
  while(abs(loss_value1 - loss_value2) >= toll){
    c_lab <- rep(0,n)
    Maha_dis_k <- matrix(0,nrow=n, ncol=max_clust)
    for (i in 1:n){
      for (k in 1:max_clust) 
        Maha_dis_k[i,k] <- alpha_Mahalanobis(alpha,centroids_mean[k,],data[i,],values,vectors)
      index <- which.min(Maha_dis_k[i,])
      c_lab[i] <- index
    }
    loss_value1 <- loss_value2
    for (k in 1:n_clust)
      centroids_mean[k,] <- colMeans(data[which(c_lab==k),])
  }
  
  # -> PLOT: prepare data to plot (df)
  df <- centroids_mean %>% t() %>% as.data.frame() %>% 
    add_column(x=1:dim(centroids_mean)[2]) %>% ##############TODO: x=time
    gather(group, y, -x)
  cluster_colors <- rainbow(max_clust) # generate colors
  # -> PLOT: current clusters (now just saved, actually plotted later with clustered colors)
  theplot.before <- ggplot(df, aes(x, y, color = group)) + geom_line(size=1)
  # theplot.before
  
  ## compute 'dists' distances matrix (to set smart tolerance + check small distances)
  dists0=matrix(0,max_clust,max_clust) #distances (with original k=max_clust)
  for (k in 1:(max_clust-1)){
    for (j in (k+1):max_clust){
      diff_centroids <- centroids_mean[k,] - centroids_mean[j,]
      dis_centroids <- norm(as.matrix(diff_centroids),type = 'i')
      print(sprintf(" - Clusters %d and %d - centroid distance = %.2f",k,j,dis_centroids))
      dists0[j,k]<-dists0[k,j]<-dis_centroids
    }
  }
  # DO: union of close clusters
  eps=median(dists0) # set epsilon as the median of initial distances (seems to work quite well)
  flag <- if(eps) 0 else 1 # 0 = check closeness (enter the while)
  do.union <- !flag
  clusts <- 1:max_clust # keep track of the mergings, eg. (1,2,3,4,5) -> (1,2,1,1,5)
  print("=> Original clusters vector:")
  print(clusts)
  while (do.union & !flag & n_clust>1){ # check closeness
    print(sprintf("Searching clusters to merge with tolerance eps=%.2f...",eps))
    did.a.merge=FALSE
    # dists=matrix(0,n_clust,n_clust) #distances (updated every while-loop)
    clusts.unique <- clusts %>% unique()
    clusts.combs <- combinations( clusts.unique%>%length(), 2, clusts.unique )
    # for (k in 1:(n_clust-1))
    #   for (j in (k+1):n_clust){
    # for (k in clusts.combs[,1])
    for (i in 1:(clusts.unique%>%length())){
      k<-clusts.combs[i,1]
      j<-clusts.combs[i,2]
      # diff_centroids <- centroids_mean[k,] - centroids_mean[j,]
      # diff_centroids <- centroids_mean[k,] - centroids_mean[j,]
      # d<-dists[j,k]<-dists[k,j]<-norm(as.matrix(diff_centroids),type = 'i')
      # d<-norm(as.matrix(diff_centroids),type = 'i')
      d <- dists0[j,k]
      if (d < eps){ # then merge cluster j into k
        c_lab[ which(c_lab==j) ] <- k
        cluster_colors[j] <- cluster_colors[k]
        clusts[j] <- k
        print(sprintf(" -> MERGED: cluster %d into %d. (centroid distance = %.2f)",j,k,d))
        did.a.merge=TRUE
      }
    }
    if(did.a.merge){
      print("=> Updated clusters vector:")
      print(clusts)
    }
    
    # convert the label vector to observe the levels and then reconvert to return a vector
    c_lab <- as.factor(c_lab)
    labels <- levels(c_lab)
    k_new <- length(labels)
    c_lab <- as.numeric(c_lab)
    
    centroids_mean_post <- matrix(0, nrow=k_new , ncol=t_points )
    for (k in 1:k_new)
      centroids_mean_post[k,] <- colMeans(data[which(c_lab==k),])
    
    if (k_new == n_clust) {
      flag_matr <- matrix(0, nrow = k_new, ncol = t_points)
      
      for (k in 1:k_new)
        for (i in 1:t_points) 
          if (centroids_mean_post[k,i] == centroids_mean[k,i])
            flag_matr[k,i] <- 1
          
          if (sum(flag_matr) == k_new*(t_points))
            flag <- 1
    }
    n_clust <- k_new
    centroids_mean <- centroids_mean_post
  }
  
  # PLOTs
  if(do.union){ # PLOT before+after (unmerged + merged)
    print(sprintf("Merging done. Obtained %d clusters.",n_clust))
    # -> PLOT: prepare data to plot (df)
    df <- centroids_mean_post %>% t() %>% as.data.frame() %>% 
      add_column(x=1:dim(centroids_mean_post)[2]) %>%
      gather(group, y, -x)
    # merged clusters plot
    theplot.after <- ggplot(df, aes(x, y, color = group)) + geom_line(size=1)
    # assign colors
    theplot.before = theplot.before + scale_color_manual(values=cluster_colors)
    theplot.after  = theplot.after  + scale_color_manual(values=cluster_colors %>% unique())
    # show
    grid.arrange(theplot.before, theplot.after, nrow = 1)
  } else { #PLOT before only (didn't merge anything)
    print(sprintf("Done.",n_clust))
    theplot.before
  }
  
  return(list("label" = c_lab,
              "centroids" = centroids_mean,
              "loss" = loss_value2
              # "plot" = grid.arrange(theplot.before, theplot.after, nrow=1)
  ))
}

k <- 5
clust <- fda_clustering_mahalanobis_union(n_clust = k, alpha = alpha, eig = eig,
                                          toll = 1e-10, data = data)


#### CLUSTERING FUNCTION UPDATING COVARIANCE WITHIN CLUSTERS ####
# There is only one while loop (flag_1 and not flag_2)
fda_clustering_mahalanobis_updated <- function(n_clust, alpha, cov_matrix, toll,data){
  
  t_points <- dim(data)[2]
  n <- dim(data)[1]
  
  # vector of labels
  c_lab <- rep(0,n)
  
  # eigenvalues and eigenfunctions for the alpha-mahalanobis function
  eig <- eigen(cov_matrix)
  values <- eig$values 
  vectors <- eig$vectors
  
  # index of each centroid randomly defined through sampling
  y0 <- rep(0,n_clust)
  vect_sample <- 1:n
  
  # while cycle checks that there are no single unit clusters in the initial step
  flag_1 <- 1
  while (flag_1 != 0) {
    
    flag_1 <- 0
    
    y0[1] <- sample(vect_sample,1)
    
    for (k in 2:n_clust) {
      value <- y0[k-1]
      
      for (i in 1:length(vect_sample)){
        if (vect_sample[i] == value)
          t = i
      }
      
      vect_sample <- vect_sample[-t]
      y0[k] <- sample(vect_sample,1)
    }
    
    
    Mahalanobis_Distance <- matrix(0, nrow = n, ncol = n)
    for (i in 1:n){
      for (j in 1:n){
        Mahalanobis_Distance[i,j] <- alpha_Mahalanobis(alpha,data[i,],data[j,]
                                                       ,values,vectors)
      }
    }
    
    # i-th unit belongs to cluster_k if the distance(centroid_k,i-th unit) is the smallest one
    Maha_dis <- matrix(0,nrow=n, ncol=n_clust)
    for (i in 1:n){
      for (k in 1:n_clust) {
        Maha_dis[i,k] <- Mahalanobis_Distance[i,y0[k]]
      }
      index <- which.min(Maha_dis[i,])
      c_lab[i] <- index
    }
    
    # define the matrix of the centroids (random centroids)
    centroids_random <- matrix(0,nrow = n_clust,ncol = t_points)
    for (k in 1:n_clust){
      centroids_random[k,] <- data[y0[k],]
    }
    
    for (k in 1:n_clust){
      if (sum(c_lab == k) == 1)
        flag_1 <- flag_1 + 1   # flag gets updated if there are single unit clusters
    }
    
  }
  
  loss_value1 <- gibbs_loss(n_clust = n_clust, centroids = centroids_random, 
                            label = c_lab, eig = eig,data = data)
  
  # update each centroid as the mean of the clusters data
  centroids_mean<-matrix(0,nrow = n_clust, ncol = t_points)
  
  for (k in 1:n_clust)
    centroids_mean[k,] <- colMeans(data[which(c_lab == k),])
  
  
  loss_value2 <- gibbs_loss(n_clust = n_clust, centroids = centroids_mean, 
                            label = c_lab, eig = eig, data = data)
  
  
  # Create the vector and the matrix that will contain eigenvalues/eigenvectors of a single cluster
  values_k <- rep(0,t_points)
  vector_k <- matrix (0, nrow = t_points, ncol = t_points)
  
  values_matrix <- matrix (0, nrow = t_points, ncol = n_clust)
  vector_matrix <- matrix (0, nrow = (n_clust*t_points), ncol = t_points )
  
  
  while(abs(loss_value1 - loss_value2) >= toll){
    
    loss_value1 <- loss_value2
    
    for (k in 1:n_clust){
      data_k <- data[which(c_lab == k),]
      cov_k <- cov(data_k)
      eig_k <- eigen(cov_k)
      values_matrix[,k] <- abs(eig_k$values)
      vector_matrix[((k-1)*t_points + 1):(k*t_points),] <- eig_k$vectors
    }
    
    
    Maha_dis_k <- matrix(0,nrow=n, ncol=n_clust)
    for (i in 1:n){
      for (k in 1:n_clust) {
        values_k <- values_matrix[,k]
        vector_k <- vector_matrix[((k-1)*t_points + 1):(k*t_points),]
        Maha_dis_k[i,k] <- alpha_Mahalanobis(alpha,centroids_mean[k,],data[i,],values_k,vector_k)
      }
      index <- which.min(Maha_dis_k[i,])
      c_lab[i] <- index
    }
    
    for (k in 1:n_clust){
      if (sum(c_lab == k) == 1) {
        centroids_mean[k,] <- data[which(c_lab == k),]
      }
      else 
        centroids_mean[k,] <- colMeans(data[which(c_lab == k),])
    }
    
    loss_value2 <- gibbs_loss_k(n_clust = n_clust, centroids = centroids_mean, 
                                label = c_lab, values_matrix, vector_matrix, data = data)
    
  }
  
  return(list("label" = c_lab, "centroids" = centroids_mean, "loss" = loss_value2,
              "K" = n_clust))
  
}


gibbs_loss <- function(n_clust, centroids, label , eig, data){
  n <- dim(data)[1]
  res = rep(0,n_clust)
  sum_partial <- 0
  
  for (k in 1:n_clust){
    for (i in 1:n){
      if (label[i] == k){
        sum_partial = alpha_Mahalanobis(alpha,data[i,],centroids[k,],eig$values,eig$vectors)
        res[k] = res[k] + sum_partial
      }
    }
  }
  
  tot = sum(res)
  return(tot)
}

# Modified for eigenvalues/eigenvectros related to the clusters. Both functions are called in the code
gibbs_loss_k <- function(n_clust, centroids, label , values_matrix, vector_matrix, data){
  t_points <- dim(data)[2]
  n <- dim(data)[1]
  
  res = rep(0,n_clust)
  sum_partial <- 0
  
  values_k <- rep(0,t_points)
  vector_k <- matrix(0, t_points, t_points)
  
  for (k in 1:n_clust){
    values_k <- values_matrix[,k]
    vector_k <- vector_matrix[((k-1)*t_points + 1):(k*t_points),]
    for (i in 1:n){
      if (label[i] == k){
        sum_partial = alpha_Mahalanobis(alpha,data[i,],centroids[k,],values_k,vector_k)
        res[k] = res[k] + sum_partial
      }
    }
  }
  
  tot = sum(res)
  return(tot)
}

# Application on the data
n <- dim(data)[1]
t_points <- dim(data)[2]

k <- 3
alpha <- 0

clust <- fda_clustering_mahalanobis_updated(n_clust = k, alpha = alpha, cov_matrix = cov(data),
                                            toll = 1e-2,  data = data)
c_opt <- clust$label
loss_min <- clust$loss
show(c_opt)  #label switching 
show(loss_min)

# Theoretical optimal plot vs clustering plot
c1 <- clust$centroids[1,]
c2 <- clust$centroids[2,]
c3 <- clust$centroids[3,]

data1 <- data[which(c_opt=='1'),]
data2 <- data[which(c_opt=='2'),]
data3 <- data[which(c_opt=='3'),]




#### FUNCTIONS THAT COULD BE USED ####

# # Check
# diff_centroids <- c1 - c2
# diff_centroids <- c3 - c4
# diff_centroids <- c1 - c3
# dis_centroids <- norm(as.matrix(diff_centroids),type='i')
# show(dis_centroids)


# flag <- 1
# while (flag <= n_clust) {
#   
#   data_k <- data[which(c_lab==flag),]
#   eig <- eigen(cov(data_k))
#   eig_dynamic <- paste0("eig_", flag)
#   assign(eig_dynamic, eig, .GlobalEnv)
#   flag <- flag + 1
#   
# }


####  Old version of fda_updated with two while loops
# fda_clustering_mahalanobis_updated <- function(n_clust, alpha, cov_matrix, toll,data){
#   
#   t_points <- dim(data)[2]
#   n <- dim(data)[1]
#   
#   # vector of labels
#   c_lab <- rep(0,n)
#   
#   # eigenvalues and eigenfunctions for the alpha-mahalanobis function
#   eig <- eigen(cov_matrix)
#   values <- eig$values 
#   vectors <- eig$vectors
#   
#   # index of each centroid randomly defined through sampling
#   y0 <- rep(0,n_clust)
#   vect_sample <- 1:n
#   
#   # while cycle checks that there are no single unit clusters in the initial step
#   flag_1 <- 1
#   while (flag_1 != 0) {
#     
#     flag_1 <- 0
#     
#     y0[1] <- sample(vect_sample,1)
#     
#     for (k in 2:n_clust) {
#       value <- y0[k-1]
#       
#       for (i in 1:length(vect_sample)){
#         if (vect_sample[i] == value)
#           t = i
#       }
#       
#       vect_sample <- vect_sample[-t]
#       y0[k] <- sample(vect_sample,1)
#     }
#     
#     
#     Mahalanobis_Distance <- matrix(0, nrow = n, ncol = n)
#     for (i in 1:n){
#       for (j in 1:n){
#         Mahalanobis_Distance[i,j] <- alpha_Mahalanobis(alpha,data[i,],data[j,]
#                                                        ,values,vectors)
#       }
#     }
#     
#     # i-th unit belongs to cluster_k if the distance(centroid_k,i-th unit) is the smallest one
#     Maha_dis <- matrix(0,nrow=n, ncol=n_clust)
#     for (i in 1:n){
#       for (k in 1:n_clust) {
#         Maha_dis[i,k] <- Mahalanobis_Distance[i,y0[k]]
#       }
#       index <- which.min(Maha_dis[i,])
#       c_lab[i] <- index
#     }
#     
#     # define the matrix of the centroids (random centroids)
#     centroids_random <- matrix(0,nrow = n_clust,ncol = t_points)
#     for (k in 1:n_clust){
#       centroids_random[k,] <- data[y0[k],]
#     }
#     
#     for (k in 1:n_clust){
#       if (sum(c_lab == k) == 1)
#         flag_1 <- flag_1 + 1   # flag gets updated if there are single unit clusters
#     }
#     
#   }
#   
#   loss_value1 <- gibbs_loss(n_clust = n_clust, centroids = centroids_random, 
#                             label = c_lab, eig = eig,data = data)
#   
#   # update each centroid as the mean of the clusters data
#   centroids_mean<-matrix(0,nrow = n_clust, ncol = t_points)
#   
#   for (k in 1:n_clust)
#     centroids_mean[k,] <- colMeans(data[which(c_lab == k),])
#   
#   
#   loss_value2 <- gibbs_loss(n_clust = n_clust, centroids = centroids_mean, 
#                             label = c_lab, eig = eig, data = data)
#   
#   
#   # Create the vector and the matrix that will contain eigenvalues/eigenvectors of a single cluster
#   values_k <- rep(0,t_points)
#   vector_k <- matrix (0, nrow = t_points, ncol = t_points)
#   
#   # just to initialize
#   n_clust_new <- n_clust
#   
#   while(abs(loss_value1 - loss_value2) >= toll){
#     
#     loss_value1 <- loss_value2
#     
#     # while cycle to check there are no single unit clusters in the intermediate steps. 
#     flag_2 <- 1
#     while (flag_2 != 0) {
#       
#       flag_2 <- 0
#       
#       values_matrix <- matrix (0, nrow = t_points, ncol = n_clust_new)
#       vector_matrix <- matrix (0, nrow = (n_clust_new*t_points), ncol = t_points )
#       
#       for (k in 1:n_clust_new){
#         data_k <- data[which(c_lab == k),]
#         cov_k <- cov(data_k)
#         eig_k <- eigen(cov_k)
#         values_matrix[,k] <- abs(eig_k$values)
#         vector_matrix[((k-1)*t_points + 1):(k*t_points),] <- eig_k$vectors
#       }
#       
#       c_lab <- rep(0,n)
#       
#       Maha_dis_k <- matrix(0,nrow=n, ncol=n_clust_new)
#       for (i in 1:n){
#         for (k in 1:n_clust_new) {
#           values_k <- values_matrix[,k]
#           vector_k <- vector_matrix[((k-1)*t_points + 1):(k*t_points),]
#           Maha_dis_k[i,k] <- alpha_Mahalanobis(alpha,centroids_mean[k,],data[i,],values_k,vector_k)
#         }
#         index <- which.min(Maha_dis_k[i,])
#         c_lab[i] <- index
#       }
#       
#       # If the number of clusters decreases:
#       clust_new <- levels(factor(c_lab))
#       n_clust_new <- length(clust_new)
#       
#       centroids_mean<-matrix(0,nrow = n_clust_new, ncol = t_points)
#       for (k in 1:n_clust_new){
#         if (sum(c_lab == k) == 1) {
#           centroids_mean[k,] <- data[which(c_lab == k),]
#         }
#         else 
#           centroids_mean[k,] <- colMeans(data[which(c_lab == k),])
#       }
#       
#       for (k in 1:n_clust_new){
#         if (sum(c_lab == k) == 1) {
#           flag_2 <- flag_2 + 1 # flag gets updated if there are single unit clusters
#         }
#       }
#       
#     }
#     
#     loss_value2 <- gibbs_loss_k(n_clust = n_clust_new, centroids = centroids_mean, 
#                                 label = c_lab, values_matrix, vector_matrix, data = data)
#     
#   }
#   
#   return(list("label" = c_lab, "centroids" = centroids_mean, "loss" = loss_value2,
#               "K" = n_clust_new))
#   
# }

# Plots
# Simulated  data plot for model 1,2,3
# x11()
# plot(time,data[1,],type = 'l', ylim = c(-3.5,7.5), col = 'firebrick2', lwd = 2)
# for(i in 2:(n-c)){
#   lines(time,data[i,],type = 'l', col = 'firebrick2',lwd = 2)
# }
# for (i in (n-c+1):n){
#   lines(time,data[i,],type = 'l', col = 'blue', lwd = 2)
# }
# title('Simulated data')
# 
# 
