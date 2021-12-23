#### CLUSTERING ON SIMULATED DATA: Cmap ####

setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject") #Pietro
#setwd("C:/Users/admin/Documents/R/Project_BS/BayesianProject") #GiuliaR

load('Simulated_WP.RData')


# Choose the model for the simulation:
# simulated data - model 1
data <- data1
eig <- eig_1
rm(data1)
rm(data2)
rm(data3)

# simulated data - model 2
data <- data2
eig <- eig_2
rm(data1)
rm(data2)
rm(data3)

# simulated data - model 3
data <- data3
eig <- eig_3
rm(data1)
rm(data2)
rm(data3)


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
# eps is the tolerance for the max-norm for the 'collapse' version. 
 
fda_clustering_mahalanobis <- function(n_clust, alpha, eig, toll, eps, data){
  
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
  
  #   and eigenfunctions for the alpha-mahalanobis function
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
  centroids_random <- matrix(0,nrow = n_clust,ncol = dim(data)[2])
  for (k in 1:n_clust){
    centroids_random[k,] <- data[y0[k],]
  }
   
  loss_value1 <- gibbs_loss(n_clust = n_clust, centroids = centroids_random, label = c_lab, data = data)
  
  # update each centroid as the mean of the clusters data
  centroids_mean<-matrix(0,nrow = n_clust, ncol = dim(data)[2])
  for (k in 1:n_clust){
    centroids_mean[k,] <- colMeans(data[which(c_lab==k),])
  }
  
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
      centroids_mean[k,] <- colMeans(data[which(c_lab==k),])
    }
    
   } 
  
  # union of similar clusters 
  flag <- 0
  while (flag == 0) {
    
    for (k in 1:(n_clust-1)){
      for (j in (k+1):n_clust){
        
        diff_centroids <- centroids_mean[k,] - centroids_mean[j,]
        dis_centroids <- norm(as.matrix(diff_centroids),type = 'i')
        
        if (dis_centroids <= eps){
          indexes <- which(c_lab==j)
          c_lab[indexes] <- k
        }
      }
    }
    
    # convert the label vector to observe the levels and then reconvert to return a vector
    c_lab <- as.factor(c_lab)
    labels <- levels(c_lab)
    k_new <- length(labels)
    c_lab<-as.numeric(c_lab)
 
    centroids_mean_post <- matrix(0, nrow=k_new , ncol=dim(data)[2] )
    for (k in 1:k_new) {
      centroids_mean_post[k,] <- colMeans(data[which(c_lab==k),])
    }
    
    if (k_new == n_clust) {
      flag_matr <- matrix(0, nrow = k_new, ncol = dim(data)[2])
      
      for (k in 1:k_new){
        for (i in 1:dim(data)[2]) {
          if (centroids_mean_post[k,i] == centroids_mean[k,i])
            flag_matr[k,i] <- 1
        }
      }
      if (sum(flag_matr) == k_new*(dim(data)[2]))
        flag <-1
    }
    n_clust <- k_new
    centroids_mean <- centroids_mean_post
    
  }
  
  return(list("label" = c_lab, "centroids" = centroids_mean))
  
}

##### Application on the simulated data ####

# Simulated  data plot 
x11()
plot(time,data[1,],type = 'l', ylim = c(-3.5,7.5), col = 'firebrick2', lwd = 2)
for(i in 2:(n-c)){
  lines(time,data[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in (n-c+1):n){
  lines(time,data[i,],type = 'l', col = 'blue', lwd = 2)
}
title('Simulated data')

k <- 3
#eps <- 1       # model 1
eps <- 1.5         # model 2
#eps <- 1.3        # model 3

clust <- fda_clustering_mahalanobis(n_clust = k, alpha = alpha, eig = eig, toll = 1e-6, eps = eps , data = data)
c_opt <- clust$label
show(c_opt)  #label switching 

c1 <- clust$centroids[1,]
c2 <- clust$centroids[2,]
c3 <- clust$centroids[3,]
#c4 <- clust$centroids[4,]


# # Check
# diff_centroids <- c1 - c2
# diff_centroids <- c3 - c4
# diff_centroids <- c1 - c3
# dis_centroids <- norm(as.matrix(diff_centroids),type='i')
# show(dis_centroids)

# Theoretical optimal plot vs clustering plot
data1 <- data[which(c_opt=='1'),]
data2 <- data[which(c_opt=='2'),]
data3 <- data[which(c_opt=='3'),]
#data4 <- data[which(c_opt=='4'),]

x11()
par(mfrow = c(1,2))
plot(time,data[1,],type = 'l', ylim = c(-3.5,9), col = 'firebrick2', lwd = 2, main = "Main and contaminated processes")
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

lines(time,c1,type = 'l', lwd = 3)
lines(time,c2,type = 'l', lwd = 3)
# lines(time,c3,type = 'l', lwd = 3)


rm(data1)
rm(data2)
# rm(data3)



##### Theoretical optimal plot vs clustering plot SMOOTHED ####
data1 <- f.data_alpha_sim[which(c_opt=='1'),]
data2 <- f.data_alpha_sim[which(c_opt=='2'),]

x11()
par(mfrow = c(1,2))

plot(time,f.data_alpha_sim[1,],type = 'l', ylim = c(-3,7.5), col = 'firebrick2', lwd = 2, main = "Smooth processes")
for(i in 2:(n-c)){
  lines(time,f.data_alpha_sim[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in (n-c+1):n){
  lines(time,f.data_alpha_sim[i,],type = 'l', col = 'blue', lwd = 2)
}

plot(time,data1[1,],type = 'l', ylim = c(-3,7.5), col = 'firebrick2', lwd = 2, main = "Clustered smoothed data")
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




#### K_CLUST FUNCTION MERGING USING LOSS MINIMIZATION ####
fda_clustering_mahalanobis_merge <- function(n_clust, alpha, eig, toll,data){
   
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
       index <-  (Maha_dis[i,])
       c_lab[i] <- index
     }
     
     # define the matrix of the centroids (random centroids)
     centroids_random <- matrix(0,nrow = l,ncol = dim(data)[2])
     for (k in 1:l){
       centroids_random[k,] <- data[y00[k],]
     }
     
     loss_value1 <- gibbs_loss(n_clust = l, centroids = centroids_random, label = c_lab, data = data)
     
     # update each centroid as the mean of the clusters data
     centroids_mean<-matrix(0,nrow = l, ncol = dim(data)[2])
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

 
k <- 3

clust <- fda_clustering_mahalanobis_merge(n_clust = k, alpha = alpha, eig = eig, toll = 1e-2,  data = data)
c_opt <- clust$label
c1 <- clust$centroids[1,]
c2 <- clust$centroids[2,]
c3 <- clust$centroids[3,]
show(c_opt)

data1 <- data[which(c_opt=='1'),]
data2 <- data[which(c_opt=='2'),]
data3 <- data[which(c_opt=='3'),]
#data4 <- data[which(c_opt=='4'),]

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
  lines(time,data2[i,],type = 'l', col = 'blue',lwd = 2)
}


lines(time,c1,type = 'l', lwd = 3)
lines(time,c2,type = 'l', lwd = 3)
lines(time,c3,type = 'l', lwd = 3)


rm(data1)
rm(data2)
rm(data3)

#### ORIGINAL CLUSTERING FUNCTION ####

fda_clustering_mahalanobis_original <- function(n_clust, alpha, eig, toll,data){
  
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
      Maha_dis[i,k] <- Mahalanobis_Distance[i,y0[k]]
    }
    index <- which.min(Maha_dis[i,])
    c_lab[i] <- index
  }
  
  # define the matrix of the centroids (random centroids)
  centroids_random <- matrix(0,nrow = n_clust,ncol = dim(data)[2])
  for (k in 1:n_clust){
    centroids_random[k,] <- data[y0[k],]
  }
  
  loss_value1 <- gibbs_loss(n_clust = n_clust, centroids = centroids_random, label = c_lab, data = data)

  # update each centroid as the mean of the clusters data
  centroids_mean<-matrix(0,nrow = n_clust, ncol = dim(data)[2])
  for (k in 1:n_clust){
    centroids_mean[k,] <- colMeans(data[which(c_lab==k),])
  }
  
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
      centroids_mean[k,] <- colMeans(data[which(c_lab==k),])
    }
    
    loss_value2 <- gibbs_loss(n_clust = n_clust, centroids = centroids_mean, label = c_lab, data = data)
  }
  
  return(list("label" = c_lab, "centroids" = centroids_mean, "loss" = loss_value1))
  
}

k <- 2
fda_clustering_mahalanobis_original(n_clust = k, alpha = alpha, eig = eig, toll = 1e-2,  data = data)$loss

k <- 5
fda_clustering_mahalanobis_original(n_clust = k, alpha = alpha, eig = eig, toll = 1e-4,  data = data)$loss
