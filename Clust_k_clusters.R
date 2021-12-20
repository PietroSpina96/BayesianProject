########################################################################################
########################################################################################
#### CLUSTERING ON SIMULATED DATA: Cmap ####


#setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject")
setwd("C:/Users/admin/Documents/R/Project_BS/BayesianProject") #GiuliaR

load('Simulated_WP.RData')
data<-data1
rm(data1)


#### Loss function ####
#centroids è una matrice avente n_clust righe e ogni riga è un centroide
gibbs_loss <- function(n_clust, centroids, label ,data){
  res = rep(0,n_clust)
  
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


#### Clustering function ####
# The parameter eig corresponds to the output of the eigen function (list of eigenvalues and eigenvectors)
# The function works with n_clust=k
# alpha is the smoothing parameter
# toll is the tolerance for the while loop

fda_clustering_mahalanobis <- function(n_clust, alpha, eig, toll,data){
  
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
      Maha_dis[i,k] <- Mahalanobis_Distance[y0[k],i]
    }
    index <- which.min(Maha_dis[i,])
    c_lab[i] <- index
  }
  
  # define the matrix of the centroids (random centroids)
  centroids_random <- matrix(0,nrow = n_clust,ncol = dim(data)[2])
  for (k in 1:n_clust)
    centroids_random[k,] <- data[y0[k],]
  
  loss_value1 <- gibbs_loss(n_clust = n_clust, centroids = centroids_random, label = c_lab, data = data)
  
  # update each centroid as the mean of the clusters data
  centroids_mean<-matrix(0,nrow = n_clust, ncol = dim(data)[2])
  for (k in 1:n_clust)
    centroids_mean[k,] <- colMeans(data[which(c_lab=='k'),])
  
  loss_value2 <- gibbs_loss(n_clust = n_clust, centroids = centroids_mean, label = c_lab, data = data)
  
  while(abs(loss_value1 - loss_value2) >= toll){
    c_lab <- rep(0,n)
    
    Maha_dis_k <- matrix(0,nrow=n, ncol=n_clust)
    for (i in 1:n){
      for (k in 1:n_clust) {
        Maha_dis_k[i,k] <- alpha_Mahalanobis(alpha,centroids_mean[k,],data[i,],values,vectors)
      }
      index <- which.min(Maha_dis_k[i,])
      c_lab[i] <- index
    }
    
    loss_value1 <- loss_value2
    
    for (k in 1:n_clust){
      centroids_mean[k,] <- colMeans(data[which(c_lab=='k'),])
    }
      
    loss_value2 <- gibbs_loss(n_clust = n_clust, centroids = centroids_mean, label = c_lab, data = data)
  } 
  
  return(list("label" = c_lab, "centroids" = centroids_mean))
  
}

# Application on the simulated data
k <- 2
clust <- fda_clustering_mahalanobis(n_clust = k, alpha = alpha, eig = eig, toll = 1e-6, data = data)
c_opt <- clust$label
c1 <- clust$centroids[1,]
c2 <- clust$centroids[2,]
show(c_opt)  #label switching 


# Theoretical optimal plot vs clustering plot
data1 <- data[which(c_opt=='1'),]
data2 <- data[which(c_opt=='2'),]

x11()
par(mfrow = c(1,2))

plot(time,data[1,],type = 'l', ylim = c(-2,7.5), col = 'firebrick2', lwd = 2, main = "Main and contaminated processes")
for(i in 2:(n-c)){
  lines(time,data[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in (n-c+1):n){
  lines(time,data[i,],type = 'l', col = 'blue', lwd = 2)
}

plot(time,data1[1,],type = 'l', ylim = c(-2,7.5), col = 'firebrick2', lwd = 2, main = "Clustered data")
for (i in 2:dim(data1)[1]){
  lines(time,data1[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in 1:dim(data2)[1]){
  lines(time,data2[i,],type = 'l', col = 'blue',lwd = 2)
}
lines(time,c1,type = 'l', lwd = 4)
lines(time,c2,type = 'l', lwd = 4)

rm(data1)
rm(data2)


# Theoretical optimal plot vs clustering plot SMOOTHED
data1 <- f.data_alpha_sim[which(c_opt=='1'),]
data2 <- f.data_alpha_sim[which(c_opt=='2'),]

x11()
par(mfrow = c(1,2))

plot(time,f.data_alpha_sim[1,],type = 'l', ylim = c(-2,7.5), col = 'firebrick2', lwd = 2, main = "Smooth processes")
for(i in 2:(n-c)){
  lines(time,f.data_alpha_sim[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in (n-c+1):n){
  lines(time,f.data_alpha_sim[i,],type = 'l', col = 'blue', lwd = 2)
}

plot(time,data1[1,],type = 'l', ylim = c(-2,7.5), col = 'firebrick2', lwd = 2, main = "Clustered smoothed data")
for (i in 2:dim(data1)[1]){
  lines(time,data1[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in 1:dim(data2)[1]){
  lines(time,data2[i,],type = 'l', col = 'blue',lwd = 2)
}
lines(time,f_alpha_approx(c1,alpha,eig$values,eig$vectors),type = 'l', lwd = 4)
lines(time,f_alpha_approx(c2,alpha,eig$values,eig$vectors),type = 'l', lwd = 4)

rm(data1)
rm(data2)





