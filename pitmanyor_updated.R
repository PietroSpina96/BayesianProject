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
  clust_obs <- as.numeric(table(label))
  
  post_vec <- rep(0,n_clust)
  for (k in 1:n_clust){
    post_vec[k] <- prod(theta + k*sigma, gamma(clust_obs[k] - sigma))
  }
  
  post <- prod(post_vec,1/(theta + n_clust*sigma),exp(-lambda*loss))
  
  return(list("posterior" = post, "cluster_size" = clust_obs))
  
}

fda_clustering_pitmanyor_updated <- function(n_clust, alpha, sigma, theta, lambda, cov_matrix , toll, data){
  
  t_points <- dim(data)[2]
  n <- dim(data)[1]
  
  # vector of labels
  c_lab <- rep(0,n)
  
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
  c_size0 <- as.numeric(table(c_lab))
  
  while (length(c_size0) < n_clust){
    for (k in 1:n_clust){
      print("Null dimension")
      n_hat <- sample(1:n,1)
      c_lab[n_hat] <- k
    }
    c_size0 <-  as.numeric(table(c_lab))
  }
  
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
  c_size1 <- posterior_value1$cluster_size
  
  
  # update each centroid as the mean of the clusters data
  centroids_mean<-matrix(0,nrow = n_clust, ncol = t_points)
  
  for (k in 1:n_clust){
    if (c_size1[k] > 1)
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
  c_size2 <- posterior_value2$cluster_size
  
  for (k in 1:n_clust){
    writeLines(sprintf("%d size cluster %d, outside while",c_size2[k],k))
  }
  
  # Compute eigenvalues and eigenfunctions of each cluster data
  for (k in 1:n_clust){
    data_k <- data[which(c_lab == k),]
    
    if (c_size1[k] == 1){
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
  
  #while(post_value2 > post_value1 & is.na(post_value2) == FALSE & is.na(post_value1) == FALSE){
  
  while(abs(loss_value1 - loss_value2) >= toll & is.na(post_value2) == FALSE & is.na(post_value1) == FALSE  && iterations < 50){
    iterations <- iterations + 1
    
    post_value1 <- post_value2
    loss_value1 <- loss_value2
    
    c_lab <- rep(0,n)
    
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
    
    c_size1 <- as.numeric(table(c_lab))
    
    while (length(c_size1) < n_clust){
      for (k in 1:n_clust){
        print("Null dimension")
        n_hat <- sample(1:n,1)
        c_lab[n_hat] <- k
      }
      c_size1 <-  as.numeric(table(c_lab))
    }
    
    for (k in 1:n_clust){
      if (c_size1[k] > 1)
        centroids_mean[k,] <- colMeans(data[which(c_lab == k),])
      else
        centroids_mean[k,] <- data[which(c_lab == k),]
    }
    
    for (k in 1:n_clust){
      data_k <- data[which(c_lab == k),]
      
      if (c_size1[k] == 1){
        data_k <- as.matrix(data_k)
        cov_k <- cov(data) + diag(rep(delta,t_points))
      }
      
      else{
        cov_k <- cov(data_k) + diag(rep(delta,t_points)) #diagonal correction
      }
      
      eig_k <- eigen(cov_k)
      
      values_matrix[,k] <- abs(eig_k$values)
      vector_matrix[((k-1)*t_points + 1):(k*t_points),] <- eig_k$vectors
    }
    
    #loss_value2 <- gibbs_loss_updated(n_clust = n_clust, centroids = centroids_mean,label = c_lab, values_matrix, vector_matrix, data = data)
    loss_value2 <- gibbs_loss_general(n_clust = n_clust, centroids = centroids_random,
                                      label = c_lab, eig = eig, values_matrix = values_matrix, 
                                      vector_matrix = vector_matrix, eig_type = 'updated', data = data)
    
    # Calculation of Pitman-Yor posterior value with updated centroids
    posterior_value2 <- posterior_pitmanyor(n_clust = n_clust, sigma = sigma, theta = theta, 
                               lambda = lambda, label = c_lab, loss = loss_value2, data = data)
    post_value2 <- posterior_value2$posterior
    c_size2 <- posterior_value2$cluster_size
    
    # writeLines(sprintf("#%d. Loss value: %.2f  /  diff: %.2f",iterations,loss_value2,abs(loss_value2-loss_value1)))
  }
  
  # if (iterations == 50)
  #   writeLines('WARNING: oscillating behaviour. Try again: you will be luckier next time')
  # else
  #   writeLines('Jackpot: Optimal partition found!')
  
  return(list("label" = c_lab, "centroids" = centroids_mean, 
              "loss" = loss_value2, "posterior" = post_value2, "clusters_dim" = c_size2))
}

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







































#### APPLICATION ON THE SIMULATED DATA ####

##### Original clustering functions: fda_clustering_mahalanobis ####
n <- dim(data)[1]
t_points <- dim(data)[2]

k <- 2
alpha <- 0.1

clust <- fda_clustering_mahalanobis_general(n_clust = k, alpha = alpha,
                                            cov_matrix = K_1, cov_type = 'fixed',
                                            toll = 1e-2,  data = data)
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

# Plot 
x11()
par(mfrow = c(1,2))
plot(time,data[1,],type = 'l', ylim = c(-3.5,9), lwd = 2, col='black',main = "Main and contaminated processes")
for(i in 2:(n1)){
  lines(time,data[i,],type = 'l', col = 'black',lwd = 2)
}
for (i in (n1+1):n){
  lines(time,data[i,],type = 'l', col = 'black', lwd = 2)
}
# legend(x=0.75,y=9.5,ncol=1,box.lwd=1,legend=c('Process 1','Process 2'),fill=c('blue','firebrick2'),x.intersp=0.3,
#text.col=c('blue','firebrick2'))

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

#legend(x=0.75,y=9.5,ncol=1,box.lwd=1,legend=c('Process 1','Process 2','Centroids'),fill=c('blue','firebrick2','black'),x.intersp=0.3,
#text.col=c('blue','firebrick2','black'))
legend(x=0.6,y=9.5,ncol=1,box.lwd=1,legend=c('Main process','Contaminated process','Centroids'),fill=c('firebrick2','blue','black'),x.intersp=0.3,
       text.col=c('firebrick2','blue','black'))


rm(data1)
rm(data2)
# rm(data3)
# rm(data4)