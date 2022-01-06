# Questa versione della updated definisce un vettore di label in modo casuale, così da avere immediatamente i centroidi
# e poter calcolare la prima loss tramite autovalori e autovettori del cluster di appartenenza. è leggermente più veloce 
# della versione 1. La seconda loss viene inizializzata a caso come 100 volte il valore della prima, solo per entrare nel while. 

# Nel while procedo come al solito: il calcolo dei nuovi autovalori e autovettori viene effettuato dopo aver calcolato i nuovi clusters. 
# è una modifica rispetto alla versione attuale dell'updated (modifica che è presente anche in versione_1 caricata).

# Anche in questo caso è consigliato runnarla dall'interlo per vedere come la loss vada a caso

# VANTAGGIO: non si richiede di inserire la covaianza dei dati. Non serve perchè assegno direttamente a caso le label e calcolo subito 
# autovett e autoval

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
# cov_matrix <- cov(data)


fda_clustering_mahalanobis_updated_2 <- function(n_clust, alpha, toll, data){
  
  t_points <- dim(data)[2]
  n <- dim(data)[1]
  
  # random vector of labels
  n_for <- n
  n_sample <- rep(0,n_clust)
  
  for (k in 1:n_clust){
    if ( k != n_clust ){
      vect <- 1:n_for
      n_sample[k] <- sample(vect,1)
      n_for <- n_for - n_sample[k]
    }
    n_sample[k] <- n - sum(n_sample)
  }
  
  c_lab <-0
  for (k in 1:n_clust){
    c_k <- rep(k,n_sample[k])
    c_lab <- c(c_lab,c_k)
  }
  c_lab <- c_lab[-1]
  c_lab <- c_lab[sample(1:n)]
  #rm(list=c('n_for','n_sample','c_k','vect'))
  
  # update each centroid as the mean of the clusters data
  centroids_mean<-matrix(0,nrow = n_clust, ncol = t_points)
  for (k in 1:n_clust){
    centroids_mean[k,] <- colMeans(data[which(c_lab == k),])
  }
  
  # Create the vector and the matrix that will contain eigenvalues/eigenvectors of a single cluster
  values_k <- rep(0,t_points)
  vector_k <- matrix (0, nrow = t_points, ncol = t_points)
  
  values_matrix <- matrix (0, nrow = t_points, ncol = n_clust)
  vector_matrix <- matrix (0, nrow = (n_clust*t_points), ncol = t_points )
  
  for (k in 1:n_clust){
    data_k <- data[which(c_lab == k),]
    cov_k <- cov(data_k)
    
    delta <- 1e-10
    for (l in 1:t_points)
      cov_k[l,l] <- cov_k[l,l] + delta
    
    eig_k <- eigen(cov_k)
    values_matrix[,k] <- eig_k$values
    vector_matrix[((k-1)*t_points + 1):(k*t_points),] <- eig_k$vectors
  }
  
  loss_value1 <- gibbs_loss_k(n_clust = n_clust, centroids = centroids_mean, 
                              label = c_lab, values_matrix, vector_matrix, data = data)
  loss_value2 <- 100*loss_value1
  
  while (abs(loss_value1 - loss_value2) >= toll){
    
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
loss_min <- loss_value2 #1270.904
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
for(i in 2:(40)){
  lines(time,data[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in (40+1):n){
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
