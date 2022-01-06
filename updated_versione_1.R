# Questa funzione è analoga ad updated ma con qualche cambio nel while loop. 
# Non è più presente il while iniziale che controlla i centroidi
# La matrice data da cov(data) non è definita positiva (gli ultimi autovalori sono infinitesimi ma negativi) e dunque per renderla def pos
# ho aggiunto alla diagonale principale un delta = 1e-10 che risulta maggiore degli autovalori più piccoli e sistema le cose.

# Inoltre ho messo anche fuori dal while il calcolo degli autovalori aggiornati al cluster ma solo per il primo loss_value nel while.

# Si runna dall'interno per vedere come la loss vada a caso


setwd("C:/Users/admin/Documents/R/Project_BS/BayesianProject") #GiuliaR

load('Simulated_WP.RData')
load('Functions_WP.RData')

data <- data4
time <- time4
rm(list=c('data1','data2','data3','data4'))
rm(list=c('eig_1','eig_2','eig_3'))
rm(list=c('time123','time4'))

data <- data1
eig <- eig_1
time <- time123
rm(list=c('data1','data2','data3','data4'))
rm(list=c('eig_1','eig_2','eig_3'))
rm(list=c('time123','time4'))

n <- dim(data)[1]
t_points <- dim(data)[2]

k <- 2
alpha <- 0.1

n_clust <- k
cov_matrix <- cov(data)

x11()
plot(time,data[1,],type = 'l', ylim = c(-10,10), col = 'blue', lwd = 2)
for(i in 2:40){
  lines(time,data[i,],type = 'l', col = 'blue',lwd = 2)
}
for (i in (40 + 1):n){
  lines(time,data[i,],type = 'l', col = 'firebrick2', lwd = 2)
}
title('Simulated data - model 4')

# sample <- sample(1:dim(data)[1])
# data <- data[sample,]
# cov_matrix <- cov(data)

fda_clustering_mahalanobis_updated <- function(n_clust, alpha, cov_matrix, toll,data){
  
  t_points <- dim(data)[2]
  n <- dim(data)[1]
  
  # vector of labels
  c_lab <- rep(0,n)
  
  # eigenvalues and eigenfunctions for the alpha-mahalanobis function
  delta <-  1e-10
  for (l in 1:t_points){
    cov_matrix[l,l] <- cov_matrix[l,l] + delta}
  eig <- eigen(cov_matrix)
  values <- eig$values 
  vectors <- eig$vectors
  
  x11()
  image.plot(time,time,cov_matrix)
  
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

  loss_value1 <- gibbs_loss(n_clust = n_clust, centroids = centroids_random, 
                            label = c_lab, eig = eig, data = data)
  
  # update each centroid as the mean of the clusters data
  centroids_mean<-matrix(0,nrow = n_clust, ncol = t_points)
  
  for (k in 1:n_clust){
    centroids_mean[k,] <- colMeans(data[which(c_lab == k),])
  }
  
  
  loss_value2 <- gibbs_loss(n_clust = n_clust, centroids = centroids_mean, 
                            label = c_lab, eig = eig, data = data)
  
  
  # Create the vector and the matrix that will contain eigenvalues/eigenvectors of a single cluster
  values_k <- rep(0,t_points)
  vector_k <- matrix (0, nrow = t_points, ncol = t_points)
  
  values_matrix <- matrix (0, nrow = t_points, ncol = n_clust)
  vector_matrix <- matrix (0, nrow = (n_clust*t_points), ncol = t_points )
  
  for (k in 1:n_clust){
    data_k <- data[which(c_lab == k),]
    cov_k <- cov(data_k)
    
    for (l in 1:t_points)
      cov_k[l,l] <- cov_k[l,l] + delta
    
    eig_k <- eigen(cov_k)
    values_matrix[,k] <- eig_k$values
    vector_matrix[((k-1)*t_points + 1):(k*t_points),] <- eig_k$vectors
  }
  
  while(abs(loss_value1 - loss_value2) >= toll){
    
    loss_value1 <- loss_value2

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
        centroids_mean[k,] <- colMeans(data[which(c_lab == k),])
    }
    
    for (k in 1:n_clust){
      data_k <- data[which(c_lab == k),]
      cov_k <- cov(data_k)
      
      for (l in 1:t_points)
        cov_k[l,l] <- cov_k[l,l] + delta
      
      eig_k <- eigen(cov_k)
      values_matrix[,k] <- eig_k$values
      vector_matrix[((k-1)*t_points + 1):(k*t_points),] <- eig_k$vectors
    }
    
    loss_value2 <- gibbs_loss_k(n_clust = n_clust, centroids = centroids_mean, 
                                label = c_lab, values_matrix, vector_matrix, data = data)
    
  }
  
  return(list("label" = c_lab, "centroids" = centroids_mean, "loss" = loss_value2,
              "K" = n_clust))
  
}
c_opt <- c_lab
loss_min <- loss_value2
show(c_opt)  #label switching 
show(loss_min)

# Theoretical optimal plot vs clustering plot
c1 <- centroids_mean[1,]
c2 <- centroids_mean[2,]
#c3 <- clust$centroids[3,]

data1 <- data[which(c_opt=='1'),]
data2 <- data[which(c_opt=='2'),]
data3 <- data[which(c_opt=='3'),]

x11()
par(mfrow = c(1,2))
plot(time,data[1,],type = 'l', ylim = c(-3.5,9), lwd = 2, main = "Main and contaminated processes")
for(i in 2:(80)){
  lines(time,data[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in (80+1):n){
  lines(time,data[i,],type = 'l', col = 'blue', lwd = 2)
}

plot(time,data1[1,],type = 'l', ylim = c(-3.5,9), col = 'firebrick2', lwd = 2, main = "Clustered data")
for (i in 2:dim(data1)[1]){
  lines(time,data1[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in 1:dim(data2)[1]){
  lines(time,data2[i,],type = 'l', col = 'blue',lwd = 2)
}
lines(time,c1,type = 'l', lwd = 3)
lines(time,c2,type = 'l', lwd = 3)
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
