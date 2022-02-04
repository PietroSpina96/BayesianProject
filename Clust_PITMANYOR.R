#### PITMAN-YOR MODEL ####

#### SETUP ####
setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject") #Pietro
# setwd("C:/Users/admin/Documents/R/Project_BS/BayesianProject") #GiuliaR

load('Simulated_WP.RData')
load('Functions_WP.RData')


# Choose the model for the simulation:
# simulated data - model 1
data <- data1
time <- time123
rm(list=c('data1','data2','data3','data4','data5'))
rm(list=c('time123','time45'))
rm(list=c('K_2','K_3'))

# simulated data - model 2
data <- data2
time <- time123
rm(list=c('data1','data2','data3','data4','data5'))
rm(list=c('time123','time45'))
rm(list=c('K_1','K_3'))

# simulated data - model 3
data <- data3
time <- time123
rm(list=c('data1','data2','data3','data4','data5'))
rm(list=c('time123','time45'))
rm(list=c('K_2','K_1'))

# simulated data - model 4
# The covariance matrix of data4 is cov_4 and K_4_1 and K_4_2 are the covariance matrices of the two clusters
data <- data4
time <- time45
rm(list=c('data1','data2','data3','data4','data5'))
rm(list=c('time123','time45'))
rm(list=c('K_2','K_3','K_1'))

# simulated data - model 5
# The covariance matrix of data5 is cov_5 and K_5_1 and K_5_2 are the covariance matrices of the two clusters
data <- data5
time <- time45
rm(list=c('data1','data2','data3','data4','data5'))
rm(list=c('time123','time45'))
rm(list=c('K_2','K_3','K_1'))

c <- n_cluster2     #number of items in the second cluster for model 1,2,3
n1 <- n_cluster1    #number of items in the first cluster for model 4,5
rm(list=c('n_cluster1','n_cluster2'))


#### Clustering ####
##### Function ####
# This section will be moved in the Functions.R script later on

# # To run from inside
# cov_matrix <- K_1
# n_clust <- 3
# alpha <- 0.1
# sigma <- 0.25
# theta <- 3.66
# lambda <- 0.75
# toll <- 1e-4

fda_clustering_pitmanyor <- function(n_clust, alpha, sigma, theta, lambda, cov_matrix , toll, data){
  
  t_points <- dim(data)[2]
  n <- dim(data)[1]
  
  # index of each centroid randomly defined through sampling
  y0 <- sample(1:n,n_clust,replace = FALSE)
  
  # vector of labels
  c_lab <- rep(0,n)
  
  # covariance matrix must have positive eigenvalues
  delta <- 1e-10
  diag(cov_matrix) <- diag(cov_matrix) + delta
  
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
  
  for (k in 1:n_clust){
    if (is.null(dim(data[which(c_lab == k),])[1]) == TRUE){
      print("Null dimension")
      n_hat <- sample(1:n,1)
      c_lab[n_hat] <- k
    }
  }
  
  # define the matrix of the centroids (random centroids)
  centroids_random <- matrix(0,nrow = n_clust,ncol = t_points)
  for (k in 1:n_clust){
    centroids_random[k,] <- data[y0[k],]
  }
  
  #loss_value1 <- gibbs_loss(n_clust = n_clust, centroids = centroids_random, label = c_lab, eig = eig, data = data)
  loss_value1 <- gibbs_loss_general(n_clust = n_clust, centroids = centroids_random,
                                    label = c_lab, eig = eig, values_matrix = 0, 
                                    vector_matrix = 0, eig_type = 'fixed', data = data)
  
  # Calculation of Pitman-Yor posterior value with random centroids
  post_value1 <- posterior_pitmanyor(n_clust = n_clust, sigma = sigma, theta = theta, 
                                     lambda = lambda, label = c_lab, loss = loss_value1, 
                                     data = data)$posterior
  clust1 <- posterior_pitmanyor(n_clust = n_clust, sigma = sigma, theta = theta, 
                                     lambda = lambda, label = c_lab, loss = loss_value1, 
                                     data = data)$cluster_size
  
  # update each centroid as the mean of the clusters data
  centroids_mean<-matrix(0,nrow = n_clust, ncol = t_points)
  for (k in 1:n_clust){
    if (dim(data[which(c_lab == k),])[1] > 1)
      centroids_mean[k,] <- colMeans(data[which(c_lab == k),])
    else
      centroids_mean[k,] <- data[which(c_lab == k),]
  }
  
  #loss_value2 <- gibbs_loss(n_clust = n_clust, centroids = centroids_mean, label = c_lab, eig = eig, data = data)
  loss_value2 <- gibbs_loss_general(n_clust = n_clust, centroids = centroids_mean,
                                    label = c_lab, eig = eig, values_matrix = 0, 
                                    vector_matrix = 0, eig_type = 'fixed', data = data)
  
  # Calculation of Pitman-Yor posterior value with updated centroids
  post_value2 <- posterior_pitmanyor(n_clust = n_clust, sigma = sigma, theta = theta, 
                                     lambda = lambda, label = c_lab, loss = loss_value2, 
                                     data = data)$posterior
  clust2 <- posterior_pitmanyor(n_clust = n_clust, sigma = sigma, theta = theta, 
                                     lambda = lambda, label = c_lab, loss = loss_value2, 
                                     data = data)$cluster_size
  
  # it <- 0
  # N_it <- 50
  # post_trend <- numeric(N_it)
  #while( abs(log(post_value2) - log(post_value1)) >= toll ){
  while( post_value2 > post_value1 ){
    
    # it <- it + 1
    # print(it)
    
    c_lab <- rep(0,n)
    
    post_value1 <- post_value2
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
      if (is.null(dim(data[which(c_lab == k),])[1]) == TRUE){
        print("Null dimension")
        n_hat <- sample(1:n,1)
        c_lab[n_hat] <- k
      }
    }
    
    for (k in 1:n_clust){
      if (dim(data[which(c_lab == k),])[1] > 1)
        centroids_mean[k,] <- colMeans(data[which(c_lab == k),])
      else
        centroids_mean[k,] <- data[which(c_lab == k),]
    }
    
    #loss_value2 <- gibbs_loss(n_clust = n_clust, centroids = centroids_mean, label = c_lab, eig = eig, data = data)
    loss_value2 <- gibbs_loss_general(n_clust = n_clust, centroids = centroids_mean,
                                      label = c_lab, eig = eig, values_matrix = 0, 
                                      vector_matrix = 0, eig_type = 'fixed', data = data)
    
    # Calculation of Pitman-Yor posterior value with updated centroids
    post_value2 <- posterior_pitmanyor(n_clust = n_clust, sigma = sigma, theta = theta, 
                                       lambda = lambda, label = c_lab, loss = loss_value2, 
                                       data = data)$posterior
    clust2 <- posterior_pitmanyor(n_clust = n_clust, sigma = sigma, theta = theta, 
                                       lambda = lambda, label = c_lab, loss = loss_value2, 
                                       data = data)$cluster_size
    
  }
  
  return(list("label" = c_lab, "centroids" = centroids_mean, 
              "loss" = loss_value1, "posterior" = post_value1, "clusters_dim" = clust1))
  
} 

# Comment: it is not an error to return clust1,post_value1 and loss_value1 due to the while condition.

##### Application on data 1 ####
# Iterative process
alpha <- 0.1
sigma <- 0.25
theta <- 3.66
lambda <- 0.75
n_clust <- 3
cov_matrix <- K_1

# Fix the number of simulations
nsimul <- 15
c_post<-matrix(0, nrow=nsimul, ncol=dim(data)[1])
post_value <- rep(0,nsimul)
post_dim <- matrix(0, nrow=nsimul, ncol = n_clust)

for (j in 1:nsimul){
  print(j)
  posterior1 <- fda_clustering_pitmanyor(n_clust, alpha, sigma, theta, lambda, cov_matrix , toll, data)
  c_post[j,] <- posterior1$label
  post_value[j] <- posterior1$posterior
  post_dim[j,] <- posterior1$clusters_dim
}
#show(post_dim)
#show(post_value)

best_index <-which.max(post_value)
best_posterior <- post_value[best_index]
best_c_opt <-c_post[best_index,]

show(best_posterior)
show(best_c_opt)




# k = 1
clust_py_1 <- fda_clustering_pitmanyor(n_clust = 1, alpha, sigma, theta,
                                       lambda, cov_matrix = K_1, toll = 1e-10, data = data)
# k = 2
clust_py_2 <- fda_clustering_pitmanyor(n_clust = 2, alpha, sigma, theta,
                                       lambda, cov_matrix = K_1, toll = 1e-10, data = data)
# k = 3

clust_py_3.1 <- fda_clustering_pitmanyor(n_clust = 3, alpha, sigma, theta,
                                         lambda, cov_matrix = K_1, toll = 1e-10, data = data)
clust_py_3.2 <- fda_clustering_pitmanyor(n_clust = 3, alpha, sigma, theta,
                                         lambda, cov_matrix = K_1, toll = 1e-10, data = data)

clust_py_3 <- fda_clustering_pitmanyor(n_clust = 3, alpha, sigma, theta,
                                       lambda, cov_matrix = K_1, data = data)

# checking posterior values vs loss values
clust_py_1$loss
clust_py_2$loss
clust_py_3$loss
clust_py_3.1$loss
clust_py_3.2$loss

clust_py_1$posterior
clust_py_2$posterior
clust_py_3$posterior
clust_py_3.1$posterior
clust_py_3.2$posterior

clust_py_3.1$label
clust_py_3.2$label


clust_py_3.1$vector
clust_py_3.2$vector


# Check the labels
c_opt <- clust_py_3$label
show(c_opt)  #label switching 


# Plot the tred of the posterior
library(pracma)
post_trend <- clust_py_3.2$vector
iterations <- 1:length(post_trend)
#length(post_trend)

x11()
par(mfrow=c(1,2))
plot(iterations,post_trend,type='l',col='blue')
semilogy(iterations,post_trend,type = 'l',col='blue')
#plot(iterations,log10(post_trend),type='l',col='blue')


###### Iterations trial #####

it <- 10000
post_vec <- numeric(it)
label_mat <- matrix(0, nrow = it, ncol = 100)

for (i in 1:it){
  clust <- fda_clustering_pitmanyor(n_clust = 3, alpha, sigma, theta,
                                    lambda, cov_matrix = K_1,
                                    toll = 1e-10, data = data)
  post_it[i] <- clust$posterior
  label_mat[i,] <- clust$label
  print(i)
}

max(post_it)

x11()
plot(1:it,post_it,type = 'l')
points(1:it,post_it,col = 'firebrick3',pch = 16)

#### GENERAL STUFF ####
################################################################################
# #### Posterior on Data 1 with fixed covariance ###############################
# # For now this is just a test to see if the new function works
# 
# n <- dim(data)[1]
# n_clust <- k <- 2
# alpha <- 0.1
# lambda <- 0.75
# 
# # Itialize eigen object
# eig <- eigen(K_1)
# 
# # Values initialized by considering example at page 29 of the paper
# sigma <- 0.25
# theta <- 3.66
# 
# # Initialize cetroid matrix (bad clusters)
# c <- 10
# c_lab_bad <- c(rep(1,n-c), rep(2,c))
# centroids_mean <- matrix(0, nrow = 2, ncol = 200)
# centroids_mean[1,] <- colMeans(data[1:(n-c),])
# centroids_mean[2,] <- colMeans(data[(n-c+1):n,])
# 
# post_bad <- posterior_pitmanyor(sigma = sigma, theta = theta, label = c_lab_bad, data = data)
# 
# # Initialize cetroid matrix (good clusters)
# c <- 20
# c_lab_good <- c(rep(1,n-c), rep(2,c))
# centroids_mean <- matrix(0, nrow = 2, ncol = 200)
# centroids_mean[1,] <- colMeans(data[1:(n-c),])
# centroids_mean[2,] <- colMeans(data[(n-c+1):n,])
# 
# post_good <- posterior_pitmanyor(sigma = sigma, theta = theta, label = c_lab_good, data = data)
# 
# # Check 
# post_bad$posterior
# post_good$posterior
# post_good$posterior > post_bad$posterior # as expected the optimal partition has higher posterior
# 
# # Label switching check
# c <- 20
# c_lab_good_2 <- c(rep(2,n-c), rep(1,c))
# centroids_mean <- matrix(0, nrow = 2, ncol = 200)
# centroids_mean[2,] <- colMeans(data[1:(n-c),])
# centroids_mean[1,] <- colMeans(data[(n-c+1):n,])
# 
# post_good_2 <- posterior_pitmanyor(sigma = sigma, theta = theta, label = c_lab_good_2, data = data)
# 
# post_good$posterior
# post_good_2$posterior
# post_good$posterior == post_good_2$posterior

###### Errors ####

# 1) Cluster di una dimensione - errore è nella funzione di clustering
# Error in colMeans(data[which(c_lab == k), ]) : 
#   'x' dev'essere un array con almeno due dimensioni

# 2) posterior diventa 0
# log(post_value2) : NaNs produced

# 3) cluster vuoto (anche se non sono sicuro) - errore è nella funzione della posterior
# Error in clust_obs[k] <- nk : replacement has length zero

# # To prevent: Error 2)
# if(post_value2 == 0){
#   post_value2 <- post_value1
#   print('ERROR: NULL POSTERIOR!')
#   break
# }
# 
# # Save posterior trend
# # post_trend[it] <- post_value2

# Alternative return

# return(list("label" = c_lab, "centroids" = centroids_mean, 
#             "loss" = loss_value2, "posterior" = post_value2,
#             "vector" = post_trend))


--------------------------------------------------------------------------------------------------------------------------
  # # Theoretical optimal plot vs clustering plot
  # c1 <- clust_py_2$centroids[1,]
  # c2 <- clust_py_2$centroids[2,]
  # # c3 <- clust_py_2$centroids[3,]
  # 
  # data1 <- data[which(c_opt=='1'),]
  # data2 <- data[which(c_opt=='2'),]
  # # data3 <- data[which(c_opt=='3'),]
  # 
  # 
  # # Plot 
# x11()
# par(mfrow = c(1,2))
# plot(time,data[1,],type = 'l', ylim = c(-3.5,9), lwd = 2, col='black',main = "Main and contaminated processes")
# for(i in 2:(n1)){
#   lines(time,data[i,],type = 'l', col = 'black',lwd = 2)
# }
# for (i in (n1+1):n){
#   lines(time,data[i,],type = 'l', col = 'black', lwd = 2)
# }
# # legend(x=0.75,y=9.5,ncol=1,box.lwd=1,legend=c('Process 1','Process 2'),fill=c('blue','firebrick2'),x.intersp=0.3,
# #text.col=c('blue','firebrick2'))
# 
# plot(time,data1[1,],type = 'l', ylim = c(-3.5,9), col = 'firebrick2', lwd = 2, main = "Clustered data")
# for (i in 2:dim(data1)[1]){
#   lines(time,data1[i,],type = 'l', col = 'firebrick2',lwd = 2)
# }
# for (i in 1:dim(data2)[1]){
#   lines(time,data2[i,],type = 'l', col = 'blue',lwd = 2)
# }
# # for (i in 1:dim(data3)[1]){
# #  lines(time,data3[i,],type = 'l', col = 'forestgreen',lwd = 2)
# # }
# 
# lines(time,c1,type = 'l', lwd = 3)
# lines(time,c2,type = 'l', lwd = 3)
# # lines(time,c3,type = 'l', lwd = 3)
# 
# #legend(x=0.75,y=9.5,ncol=1,box.lwd=1,legend=c('Process 1','Process 2','Centroids'),fill=c('blue','firebrick2','black'),x.intersp=0.3,
# #text.col=c('blue','firebrick2','black'))
# legend(x=0.6,y=9.5,ncol=1,box.lwd=1,legend=c('Main process','Contaminated process','Centroids'),fill=c('firebrick2','blue','black'),x.intersp=0.3,
#        text.col=c('firebrick2','blue','black'))
# 
# 
# rm(data1)
# rm(data2)
# # rm(data3)