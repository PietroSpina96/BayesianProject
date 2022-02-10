#### CLUSTERING ON SIMULATED DATA: Cmap with PITMAN-YOR ####
# FUNZIONE DA SISTEMARE: non va benissimo. La partizione ottimale che restituisce ha una sola osservazione nel secondo cluster

#### SETUP ####
setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject") #Pietro
setwd("C:/Users/admin/Documents/R/Project_BS/BayesianProject") #GiuliaR

load('Simulated_WP.RData')
load('Functions_WP.RData')


# Choose the model for the simulation:
# simulated data - model 4
# The covariance matrix of data4 is cov_4 and K_4_1 and K_4_2 are the covariance matrices of the two clusters
data <- data4
time <- time45
rm(list=c('data1','data2','data3','data4','data5'))
rm(list=c('time123','time45'))
rm(list=c('K_2','K_3','K_1'))

c <- n_cluster2     #number of items in the second cluster for model 1,2,3
n1 <- n_cluster1    #number of items in the first cluster for model 4,5
rm(list=c('n_cluster1','n_cluster2'))
##################################################################################

posterior_pitmanyor <- function(n_clust, sigma, theta , lambda, label, loss, data){
  
  # Create vector counting the number of observations in each cluster
  cluster_size <- rep(0,n_clust)
  for (k in 1:n_clust){
    cluster_size[k] <- sum(label == k)
  }
  
  post_vec <- rep(0,n_clust)
  for (k in 1:n_clust){
    post_vec[k] <- prod(theta + k*sigma, gamma(cluster_size[k] - sigma))
  }
  
  post <- prod(post_vec,1/(theta + n_clust*sigma),exp(-lambda*loss))
  
  return(list("posterior" = post, "cluster_size" = cluster_size))
  
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

gibbs_loss_updated <- function(n_clust, centroids, label , values_matrix, vector_matrix, data){
  t_points <- dim(data)[2]
  n <- dim(data)[1]
  cluster_size <- rep(0,n_clust)
  res = rep(0,n_clust)
  
  for (k in 1:n_clust){
    cluster_size[k] <- sum(label == k)
  }
  
  values_k <- rep(0,t_points)
  vector_k <- matrix(0, t_points, t_points)
  
  for (k in 1:n_clust){
    if(cluster_size[k] == 0){
      res[k] <- 0
    }
    
    else{
      values_k <- values_matrix[,k]
      vector_k <- vector_matrix[((k-1)*t_points + 1):(k*t_points),]
      for (i in 1:n)
        if (label[i] == k)
          res[k] = res[k] +  alpha_Mahalanobis(alpha,data[i,],centroids[k,],values_k,vector_k)
    }
  }
  
  return(sum(res))
}

fda_clustering_pitmanyor_updated <- function(n_clust, alpha, sigma, theta, lambda, cov_matrix , toll, data){
  
  t_points <- dim(data)[2]
  n <- dim(data)[1]
  
  # vector of labels
  c_lab <- rep(0,n)
  
  # Dimension of clusters
  c_size0 <- c_size1 <- c_size2 <- rep(0,n_clust)
  
  # Create the vector and the matrix that will contain eigenvalues/eigenvectors of a single cluster
  {values_k <- rep(0,t_points)
    vector_k <- matrix (0, nrow = t_points, ncol = t_points)
    
    values_matrix <- matrix (0, nrow = t_points, ncol = n_clust)
    vector_matrix <- matrix (0, nrow = (n_clust*t_points), ncol = t_points )}
  
  # covariance matrix must have positive eigenvalues
  delta <- 1e-10
  diag(cov_matrix) <- diag(cov_matrix) + delta
  
  # eigenvalues and eigenfunctions for the alpha-mahalanobis function
  {eig <- eigen(cov_matrix)
    values <- eig$values 
    vectors <- eig$vectors}
  
  # centroids sampling 
  y0 <- sample(1:n,n_clust,replace = FALSE)
  
  Mahalanobis_Distance <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n){
    for (j in 1:n)
      Mahalanobis_Distance[i,j] <- alpha_Mahalanobis(alpha,data[i,],data[j,]
                                                     ,values,vectors)
  }
  
  # i-th unit belongs to cluster_k if the distance(centroid_k,i-th unit) is the smallest one
  Maha_dis <- matrix(0,nrow=n, ncol=n_clust)
  for (i in 1:n){
    for (k in 1:n_clust)
      Maha_dis[i,k] <- Mahalanobis_Distance[i,y0[k]]
    index <- which.min(Maha_dis[i,])
    c_lab[i] <- index
  }
  
  # Checking for empty clusters
  for (k in 1:n_clust){
    c_size0[k] <- sum(c_lab == k)
  }
  
  for (k in 1:n_clust){
    if (c_size0[k] == 0){
      print("Null dimension")
      n_hat <- sample(1:n,1)
      c_lab[n_hat] <- k
    }
  }
  c_size0 <- as.numeric(table(c_lab)) 
  
  # define the matrix of the centroids (random centroids)
  centroids_random <- matrix(0,nrow = n_clust,ncol = t_points)
  for (k in 1:n_clust){
    centroids_random[k,] <- data[y0[k],]
  }
  
  
  # loss_value1 <- gibbs_loss(n_clust = n_clust, centroids = centroids_random,label = c_lab, eig = eig,data = data)
  loss_value1 <- gibbs_loss_general(n_clust = n_clust, centroids = centroids_random,
                                    label = c_lab, eig = eig, values_matrix = values_matrix, 
                                    vector_matrix = vector_matrix, eig_type = 'fixed', data = data)
  
  # Calculation of Pitman-Yor posterior value with random centroids
  posterior_value1 <-  posterior_pitmanyor(n_clust = n_clust, sigma = sigma, theta = theta, 
                                           lambda = lambda, label = c_lab, loss = loss_value1, data = data)
  post_value1 <- posterior_value1$posterior
  # c_size1 <- posterior_value1$cluster_size
  # c_size1 == c_size0
  
  
  # update each centroid as the mean of the clusters data
  centroids_mean<-matrix(0,nrow = n_clust, ncol = t_points)
  
  for (k in 1:n_clust){
    if (c_size0[k] > 1)
      centroids_mean[k,] <- colMeans(data[which(c_lab == k),])
    else
      centroids_mean[k,] <- data[which(c_lab == k),]
  }
  
  
  #loss_value2 <- gibbs_loss(n_clust = n_clust, centroids = centroids_mean,label = c_lab, eig = eig, data = data)
  loss_value2 <- gibbs_loss_general(n_clust = n_clust, centroids = centroids_mean,
                                    label = c_lab, eig = eig, values_matrix = values_matrix, 
                                    vector_matrix = vector_matrix, eig_type = 'fixed', data = data)
  
  # Calculation of Pitman-Yor posterior value with updated centroids
  posterior_value2 <-  posterior_pitmanyor(n_clust = n_clust, sigma = sigma, theta = theta, 
                                           lambda = lambda, label = c_lab, loss = loss_value2, 
                                           data = data)
  post_value2 <- posterior_value2$posterior
  # c_size2 <- posterior_value2$cluster_size
  # c_size2 == c_size1 == c_size0
  
  for (k in 1:n_clust){
    writeLines(sprintf("%d size cluster %d, outside while",c_size2[k],k))
  }
  
  # Compute eigenvalues and eigenfunctions of each cluster data
  for (k in 1:n_clust){
    data_k <- data[which(c_lab == k),]
    
    if (c_size0[k] == 1){
      cov_k <- cov(data) + diag(rep(delta,t_points))
    }
    
    else{
      cov_k <- cov(data_k) + diag(rep(delta,t_points)) #diagonal correction
    }
    
    eig_k <- eigen(cov_k)
    
    values_matrix[,k] <- abs(eig_k$values)
    vector_matrix[((k-1)*t_points + 1):(k*t_points),] <- eig_k$vectors
  }
  
  # Count how many iterations we will perform in the cycle
  iterations <- 0
  trend_loss <- trend_post <- rep(0,50)
  c_trend <- matrix(0,nrow=50,ncol=n)
  
  #while(post_value2 > post_value1 & is.na(post_value2) == FALSE & is.na(post_value1) == FALSE){
  
  while(abs(loss_value1 - loss_value2) >= toll  && iterations < 50){
    #while(iterations < 50){
    iterations <- iterations + 1
    
    post_value1 <- post_value2
    loss_value1 <- loss_value2
    c_size1 <- c_size2
    
    c_lab <- rep(0,n)
    c_size2 <- rep(0,n_clust)
    
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
    
    # Checking for empty clusters
    for (k in 1:n_clust){
      c_size2[k] <- sum(c_lab == k)
    }
    
    for (k in 1:n_clust){
      if (c_size2[k] > 1)
        centroids_mean[k,] <- colMeans(data[which(c_lab == k),])
      if (c_size2[k] == 1)
        centroids_mean[k,] <- data[which(c_lab == k),]
      else
        centroids_mean[k,] <- colMeans(data)
    }
    
    for (k in 1:n_clust){
      if(c_size2[k] == 0){
        data_k <- vector(mode = 'numeric',length=t_points)
        cov_k <- cov(data) + diag(rep(delta,t_points))
      }
      
      else{
        data_k <- data[which(c_lab == k),]
        
        if (c_size2[k] == 1){
          data_k <- as.matrix(data_k)
          cov_k <- cov(data) + diag(rep(delta,t_points))
        }
        
        if(c_size2[k] > 1){
          cov_k <- cov(data_k) + diag(rep(delta,t_points)) #diagonal correction
        }
        
        eig_k <- eigen(cov_k)
        
        values_matrix[,k] <- abs(eig_k$values)
        vector_matrix[((k-1)*t_points + 1):(k*t_points),] <- eig_k$vectors
      }
    }
    
    #loss_value2 <- gibbs_loss_updated(n_clust = n_clust, centroids = centroids_mean,label = c_lab, values_matrix, vector_matrix, data = data)
    loss_value2 <- gibbs_loss_general(n_clust = n_clust, centroids = centroids_mean,
                                      label = c_lab, eig = eig, values_matrix = values_matrix, 
                                      vector_matrix = vector_matrix, eig_type = 'updated', data = data)
    
    # Calculation of Pitman-Yor posterior value with updated centroids
    posterior_value2 <- posterior_pitmanyor(n_clust = n_clust, sigma = sigma, theta = theta, 
                                            lambda = lambda, label = c_lab, loss = loss_value2, data = data)
    post_value2 <- posterior_value2$posterior
    # c_size2 <- posterior_value2$cluster_size
    
    trend_loss[iterations] <- loss_value2
    trend_post[iterations] <- post_value2
    c_trend[iterations,] <- c_size1
    
  }
  
  return(list("label" = c_lab, "centroids" = centroids_mean, 
              "loss" = loss_value1, "posterior" = post_value1, "clusters_dim" = c_size1,"trend_post"=trend_post,"trend_loss"=trend_loss,"c_trend"=c_trend))
}

# Iterative process
alpha <- 0.1
sigma <- 0.25
theta <- 3.66 #3.66
lambda <- 0.75 #0.75
n_clust <- 2
toll <- 1e-2
cov_matrix <- cov(data)

# Fix the number of simulations
nsimul <-20
c_post<-matrix(0, nrow=nsimul, ncol=dim(data)[1])
post_value <- rep(0,nsimul)
post_dim <- matrix(0, nrow=nsimul, ncol = n_clust)

for (j in 1:nsimul){
  print(j)
  posterior1 <- fda_clustering_pitmanyor_updated(n_clust, alpha, sigma, theta, lambda, cov_matrix , toll, data)
  c_post[j,] <- posterior1$label
  post_value[j] <- posterior1$posterior
  post_dim[j,] <- posterior1$clusters_dim
}

# trend_loss <-posterior1$trend_loss
# trend_post <-posterior1$trend_post
# trend_c <- posterior1$c_trend
# 
# library(pracma)
# x11()
# par(mfrow=c(1,2))
# semilogy((1:50),trend_loss,col='blue')
# semilogy((1:50),trend_post,col='red')

best_index <-which.max(post_value)
best_posterior <- post_value[best_index]
best_c_opt <-c_post[best_index,]

show(best_posterior)
show(best_c_opt)
show(best_index)

# Iterative process
alpha <- 0.1
sigma <- 0.25
theta <- 3.66
lambda <- 0.75
n_clust <- 2
toll <- 1e-2
cov_matrix <- cov(data)

# Fix the number of simulations
nsimul <-100
c_post<-matrix(0, nrow=nsimul, ncol=dim(data)[1])
post_value <- rep(0,nsimul)
post_dim <- matrix(0, nrow=nsimul, ncol = n_clust)

for (j in 1:nsimul){
  print(j)
  posterior1 <- fda_clustering_pitmanyor_updated(n_clust, alpha, sigma, theta, lambda, cov_matrix , toll, data)
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


