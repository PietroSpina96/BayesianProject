#### FUNCTIONS ####

#### SETUP ####

library(ggplot2)
library(dplyr) # pipe (%>%)
library(tidyr) # gather()
library(tibble) # add_column()
library(gridExtra) # grid.arrange()
require(gtools) # combinations()
# library(viridis) # color palette

setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject") #Pietro
#setwd("C:/Users/admin/Documents/R/Project_BS/BayesianProject") #GiuliaR
# setwd("C:\Users\imthe\OneDrive - Politecnico di Milano\6.1 Bayesian\Proj\GitBayesianProject") #Alessio

load('Functions_WP.RData')


##### Basic Functions ####

#scalar product
scalar_prod<- function (f1,f2) {
   # f sono vettori colonna
   t(f1) %*% f2
}

# # equivalently. Not loaded in the workspace
# scalar_prod_norm<- function (f1,f2) {
#    # f sono vettori colonna
#    f1<-as.matrix(f1)
#    f2<-as.matrix(f2)
#    res<- (norm(f1+f2,'2')^2 - norm(f1-f2,'2')^2)/4
#    return(res)
# }


# alpha-mahalanobis distance
# INPUT: 
## alpha: smoothing parameter
## lambda: eigenvalues 
## eigenft: eigenfunctions (of the covariance matrix)

# OUTPUT: alpha-mahalanobis distance (squared norm) between functions f1,f2

alpha_Mahalanobis <- function(alpha,f1,f2,lambda,eigenft) {
   # es. f1 = data[i,]     or    f1 = data[j,]
   #     f2 = data[j,]           f2 = centroids[k,]
   # ( the correction '+sum(log(lambda))', if needed, is better added outside at the end )
   sum( lambda/(lambda+alpha)^2 * scalar_prod(f1-f2,eigenft)^2 )
}



#alpha_approximation of the function f with f_alpha
#same input as above
f_alpha_approx <-function(f,alpha,lambda,eigenft){
   t_points <- length(f)

   coeff<-prod<-res<-rep(0,t_points)
   approx<-matrix(0,t_points,t_points)

   for (j in 1:t_points) {

      coeff[j]<-lambda[j]/(lambda[j]+alpha)
      prod[j]<-scalar_prod(f,eigenft[,j])
      approx[,j]<- as.numeric(coeff[j]*prod[j])*eigenft[,j]
      res<-res+approx[,j]

   }
   return(res)
}

# norm and inner product wrt sample covariance function K
# same input as above
norm2_K <- function (f,lambda,eigenft){
   t_points <- length(f)
   norm_vect<-prod<-rep(0,t_points)
   
   for (j in 1:t_points){
      prod[j]<-(scalar_prod(f,eigenft[,j]))^2
      norm_vect[j]<-prod[j]/lambda[j]
   }
   
   res<-sum(norm_vect)
   return(res)
}

inner_product_K<- function(f,g,lambda,eigenft) {
   t_points <- length(f)
   prod_vect<-prod_f<-prod_g<-rep(0,t_points)
   
   for (j in 1:t_points){
      prod_f[j]<-(scalar_prod(f,eigenft[,j]))
      prod_g[j]<-(scalar_prod(g,eigenft[,j]))
      
      norm_vect[j]<-prod_f[j]*prod_g[j]/lambda[j]
   }
   
   res<-sum(norm_vect)
   return(res)
}


# Covariance function estimator K_hat
# It has been replaced by the R_covariance function: there is a slightly difference in the elements
kernel_estimator <- function(process,med,n,s,t,points){
  est <- rep(0,points)
  for (i in 1:n) {
    est[i] <- (process[i,s] - med[s])*(process[i,t] - med[t])
  }
  estimator = (1/n)*sum(est)

  return(estimator)
}
 
# compute distances between centroids (used in cluster merging)
centroids_dists <- function(centroids_mean,normtype='i'){
   writeLines("Calculating clusters' centroids distances...")
   max_clust <- dim(centroids_mean)[1]
   dists=matrix(0,max_clust,max_clust) #distances (with original k=max_clust)
   for (k in 1:(max_clust-1)){
      for (j in (k+1):max_clust){
         diff_centroids <- centroids_mean[k,] - centroids_mean[j,]
         dis_centroids <- norm(as.matrix(diff_centroids),type = normtype)
         writeLines(sprintf(" - Clusters %d and %d - centroid distance = %.2f",k,j,dis_centroids))
         dists[j,k]<-dists[k,j]<-dis_centroids
      }
   }
   return(dists)
}

##### Model Functions ####
###### Gibbs loss ####
# INPUT:
## n_clust:   number of clusters fixed a priori
## centroids: matrix of the centroids (random or mean) with nrow=n_clust 
## label:     vector of labels
## eig:       list of eigenvalues and eigenfunctions
## data:      dataset

# OUTPUT: value of the gibbs loss
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
 
# Modified for eigenvalues/eigenvectors related to the clusters. Both functions are called in the code
# INPUT: 
## values_matrix: matrix that contains the eigenvalues of all clusters
## vector_matrix: matrix that contains the eigenvectors of all clusters

# OUTPUT: value of the gibbs loss
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
 
###### Clustering Function (Original) ####
# INPUT:
## n_clust:    number of clusters
## alpha:      smoothing parameter of the data
## toll:       tolerance for the while loop
## cov_matrix: covariance matrix of the data

# OUTPUT: list(clusters_vector_label, centroids, loss_value)
fda_clustering_mahalanobis <- function(n_clust, alpha, cov_matrix, toll,data){
   
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


###### Clustering Function + Updating Covs ####
# INPUT: same input of fda_clustering_mahalanobis

fda_clustering_mahalanobis_updated <- function(n_clust, alpha, cov_matrix, toll,data){
   
   t_points <- dim(data)[2]
   n <- dim(data)[1]
   
   # vector of labels
   c_lab <- rep(0,n)
   
   # covariance matrix must have positive eigenvalues
   delta <- 1e-10
   diag(cov_matrix) <- diag(cov_matrix) + delta
   
   # eigenvalues and eigenfunctions for the alpha-mahalanobis function
   eig <- eigen(cov_matrix)
   values <- eig$values 
   vectors <- eig$vectors
   
   # while cycle checks that there are no single unit clusters in the initial step
   flag_1 <- 1
   while (flag_1 != 0) {
      
      flag_1 <- 0
      
      # centroids sampling 
      y0 <- sample(1:n,n_clust,replace = FALSE)
      
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
         
         for (l in 1:t_points)
            cov_k[l,l] <- cov_k[l,l] + delta
         
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
 


######  Merging clusters ####
# looks at centroids_mean and detects if some clusters are too close and should be merged.
# eps = tolerance for the union of clusters (eps = 0 -> no merging),
#       automatically set based on mean centroids distances
## INPUT: centroids_mean   = clusters
##        eps              = merging tolerance (*optional)
## OUTPUT: centroids_mean (reduced)
#TODO: fix messed-up cluster tracking in some random cases

clusters_union <- function(clust, eps=0){
   writeLines("  **  MERGING CLOSE CLUSTERS  **  ")
   centroids_mean <- clust$centroids
   c_lab <- clust$label
   max_clust <- dim(centroids_mean)[1]
   t_points <- dim(centroids_mean)[2]
   if(max_clust<3){
      writeLines(sprintf("Only %d clusters, fine like that.\nExiting.",max_clust))
      return(clust)
   }
   ## compute 'dists' distances matrix (to set smart tolerance + check small distances)
   dists <- centroids_dists(centroids_mean,'2')
   
   # DO: union of close clusters
   clusts <- 1:max_clust # keep track of the mergings, eg. (1,2,3,4,5) -> (1,2,1,1,5)
   writeLines(" => Initial clusters vector:")
   print(clusts)
   
   # if(missing(eps)) eps <- median(dists) # set merging tolerance
   if(missing(eps)) eps <- (median(dists)+mean(dists))/2 # set merging tolerance
   # if(missing(eps)) eps <- 1.678*sd(dists[which(dists>0)]) # set merging tolerance
   do.merge <- TRUE
   n_clust <- max_clust
   while (do.merge){ # check closeness
      writeLines(sprintf("\nSearching clusters to merge with tolerance eps=%.2f...",eps))
      did.a.merge <- FALSE
      clusts.unique <- clusts %>% unique() # eg. (1,2,1,1,5) -> (1,2,5)
      clusts.combs <- combinations( clusts.unique%>%length(), 2, clusts.unique ) # creates a n_clust-by-2 matrix with all combinations of cluster couples
      # for (k in 1:(n_clust-1))
      #   for (j in (k+1):n_clust){
      # for (k in clusts.combs[,1])
      for (i in 1:(clusts.combs%>%dim())[1]){
         k<-clusts.combs[i,1]
         j<-clusts.combs[i,2]
         # diff_centroids <- centroids_mean[k,] - centroids_mean[j,]
         # diff_centroids <- centroids_mean[k,] - centroids_mean[j,]
         # d<-dists[j,k]<-dists[k,j]<-norm(as.matrix(diff_centroids),type = 'i')
         # d<-norm(as.matrix(diff_centroids),type = 'i')
         d <- dists[j,k]
         if (d < eps){ # then merge cluster j into k
            did.a.merge <- TRUE
            c_lab[ which(c_lab==j) ] <- clusts[k]
            clusts[j] <- clusts[k]
            writeLines(sprintf(" -> MERGED: cluster %d into %d. (centroid distance = %.2f)",j,k,d))
         }
      }
      if(did.a.merge){
         writeLines(" => Updated clusters vector:")
         print(clusts)
      }
      
      # how many labels left
      k_new <- c_lab %>% unique() %>% length()
      
      # compute new centroids
      centroids_mean_post <- matrix(0, nrow=k_new , ncol=t_points)
      for (k in 1:k_new)
         centroids_mean_post[k,] <- colMeans(data[which(c_lab==(clusts%>%unique())[k]),])
      
      if (k_new == n_clust) {
         flag_matr <- matrix(0, nrow = k_new, ncol = t_points)
         
         for (k in 1:k_new)
            for (i in 1:t_points) 
               if (centroids_mean_post[k,i] == centroids_mean[k,i])
                  flag_matr[k,i] <- 1
               
         if (sum(flag_matr) == k_new*(t_points))
            do.merge <- FALSE
      }
      n_clust <- k_new
      centroids_mean <- centroids_mean_post
   }
   if(n_clust<max_clust)
      writeLines(sprintf("==> Merging done. Obtained %d clusters (starting from %d).",n_clust,max_clust))
   else
      writeLines("No further merging needed, clusters are good!\nDone.")
   
   return(list("label" = c_lab,
               "centroids" = centroids_mean,
               #"loss" = loss_value2, # maybe again meglio una funzione a parte? o la ricalcolo qua?
               "K" = n_clust,
               "v" = clusts) # only used to manage colors in plot before/after merging
          )
}


###### PLOT clusters ####
# plots (x,y)=(time,centroids_mean) with colors
clusters_plot <- function(time, clust, cols){
   ## INPUT: time              = x (vector: n)
   ##        clust$centroids   = y (matrix: k by n)
   ##        cols              = cluster hex colors (vector: n)(*optional)
   ## OUTPUT: the plot object
   
   if(missing(cols)) cols <- rainbow(dim(clust$centroids)[1])
   
   df <- clust$centroids %>% t() %>% as.data.frame() %>% 
      add_column(x=time) %>%
      gather(group, y, -x)
   
   theplot <- ggplot(df, aes(x, y, color = group)) + 
      geom_line(size=1) +
      scale_color_manual(values=cols) +
      labs(x="time",y="centroids")
   
   return(theplot)
}

# The following function is not loaded in the workspace
#####  Clustering function merging using loss minimization ####
# The parameter eig corresponds to the output of the eigen function (list of eigenvalues and eigenvectors)
# The function works with n_clust=k
# alpha is the smoothing parameter
# toll is the tolerance for the while loop

# fda_clustering_mahalanobis_merge <- function(n_clust, alpha, cov_matrix, toll,data){
#    
#    n <- dim(data)[1] 
#    t_points <- dim(data)[2]
#    
#    # index of each centroid randomly defined through sampling
#    y0 <- rep(0,n_clust)
#    vect_sample <- 1:n
#    
#    y0[1] <- sample(vect_sample,1)
#    
#    for (k in 2:n_clust) {
#       value <- y0[k-1]
#       
#       for (i in 1:length(vect_sample)){
#          if (vect_sample[i] == value)
#             t = i
#       }
#       
#       vect_sample <- vect_sample[-t]
#       y0[k] <- sample(vect_sample,1)
#    }
#    
#    # matrix of labels
#    c_lab_mat <- matrix(0, nrow = n_clust, ncol = n)
#    
#    # loss vector
#    loss_vec <- rep(0,n_clust)
     
     # covariance matrix must have positive eigenvalues
#    delta <- 1e-10
#    for (l in 1:t_points){
#       cov_matrix[l,l] <- cov_matrix[l,l] + delta
# }
#    
#    # eigenvalues and eigenfunctions for the alpha-mahalanobis function
#    values <- eig$values
#    vectors <- eig$vectors
#    
#    Mahalanobis_Distance <- matrix(0, nrow = n, ncol = n)
#    for (i in 1:n){
#       for (j in 1:n){
#          Mahalanobis_Distance[i,j] <- alpha_Mahalanobis(alpha,data[i,],data[j,],values,vectors)
#       }
#    }
#    
#    for (l in 1:n_clust){
#       
#       # vector of labels
#       c_lab <- rep(0,n)
#       
#       # starting centroids
#       y00 <- rep(0,l)
#       for (x in 1:l){
#          y00[x] <- y0[x]
#       }
#       
#       # i-th unit belongs to cluster_k if the distance(centroid_k,i-th unit) is the smallest one
#       Maha_dis <- matrix(0,nrow=n, ncol=l)
#       for (i in 1:n){
#          for (k in 1:l) {
#             Maha_dis[i,k] <- Mahalanobis_Distance[i,y00[k]]
#          }
#          index <-which.min(Maha_dis[i,])
#          c_lab[i] <- index
#       }
#       
#       # define the matrix of the centroids (random centroids)
#       centroids_random <- matrix(0,nrow = l,ncol = t_points )
#       for (k in 1:l){
#          centroids_random[k,] <- data[y00[k],]
#       }
#       
#       loss_value1 <- gibbs_loss(n_clust = l, centroids = centroids_random, label = c_lab, data = data)
#       
#       # update each centroid as the mean of the clusters data
#       centroids_mean<-matrix(0,nrow = l, ncol = t_points )
#       for (k in 1:l){
#          centroids_mean[k,] <- colMeans(data[which(c_lab==k),])
#       }
#       
#       loss_value2 <- gibbs_loss(n_clust = l, centroids = centroids_mean, label = c_lab, data = data)
#       
#       Maha_dis_k <- matrix(0,nrow=n, ncol=l)
#       for (i in 1:n ){
#          for (k in 1:l) {
#             Maha_dis_k[i,k] <- alpha_Mahalanobis(alpha,centroids_mean[k,],data[i,],values,vectors)
#          }
#          index2 <- which.min(Maha_dis_k[i,])
#          c_lab[i] <- index2
#       }
#       
#       loss_value1 <- loss_value2
#       
#       for (k in 1:l){
#          centroids_mean[k,] <- colMeans(data[which(c_lab==k),])
#       }
#       
#       loss_value2 <- gibbs_loss(n_clust = l, centroids = centroids_mean, label = c_lab, data = data)
#       
#       loss_vec[l] <- loss_value2
#       c_lab_mat[l,] <- c_lab
#       
#       if (l == 1){
#          centr_mat <- centroids_mean
#       } else {centr_mat <- rbind(centr_mat,centroids_mean)}
#    } 
#    
#    loss <- min(loss_vec)
#    index_loss <- match(loss,loss_vec)
#    
#    if (index_loss == 1){
#       a <- 1
#    } else {a <- sum(1:(index_loss-1)) + 1}
#    
#    if (index_loss == 1){
#       b <- 1
#    } else {b <- sum(1:index_loss)}
#    
#    return(list("label" = c_lab_mat[index_loss,], "centroids" = centr_mat[a:b,], "loss" = loss, "loss_vector" = loss_vec))
#    
# } 


##### Save WP ####
 
setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject")
save.image("Functions_WP.RData")

save.image("~/R/Project_BS/BayesianProject/Functions_WP.RData") #GiuliaR
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 