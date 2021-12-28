
# DA RIGA 184 NUOVA FUNZIONE PER CLUSTERING

#### IMPORTAZIONE DATI ####
data <- read.table('dati.csv', head=TRUE, sep=',', row.names = "id")
data <- data[,-c(5,9,21,22,23,24,25,26)]
data[,2] <- ifelse(data[,2] == 1, "M", "F")
data[,3] <- ifelse(data[,3] == 1, "TR", "VS")
data[,4] <- ifelse(data[,4] == 1, "DX", "SX")
data[,7] <- ifelse(data[,7] == 'FRMN', "SI", "NO")
data[,11] <- ifelse(data[,11] == 'B', 2, 1)
data[,8] <- as.factor(toupper(data[,8]))
data[,9] <- as.factor(toupper(data[,9]))
data[,12] <- as.factor(toupper(data[,12]))
data[,13] <- as.factor(toupper(data[,13]))
data[data[,12]=='NO ',12] <- c('NO')
data[data[,13]=='NO ',13] <- c('NO')
data[,12] <- factor(data[,12])
data[,13] <- factor(data[,13])  

colnames(data) <- c("AGE", "SEX", "EZIOLOGIA", "LATE", "MESIC", "MESIR", "FRMN",
                    "SLDX", "SLSX", "SLSEP", "SLSEPINC", "MLDX", "MLSX", 
                    "MLSEP", "MLSEP2", "GOSE", "LCF", "DRS")
attach(data)
options("scipen"=100, "digits"=4)


res <- list()
latency <- c("SL", "ML")
late <- c("sx", "dx")
w <- 1
for(i in 1:26)
{
  for(j in 1:2)
  {
    for(k in 1:2)
    {
      fileimp <- paste(paste("(", i, ")/", i, sep=''), latency[j], 
                       paste(late[k], ".asc", sep=''), sep = ' ')
      
      restemp <- try(read.table(file=fileimp, skip = 13, dec= ','))
      if('try-error' %in% class(restemp)) next
      else
      {
        res[[w]] <- c(restemp[1:1600,])
        res[[w]]$sogg <- c(i)
        res[[w]]$latency <- latency[j]
        res[[w]]$late <- late[k]
        w <- w+1
      }   
    }
  }
}

#### GENERO MATRICI PER OGNI CONDIZIONE TESTATA AURICOLARE ####
importMatrix <- function(x, type, position)
{
  temp <- matrix(0, nrow=length(res[[1]][[1]]))
  ntmp <- c()
  j <- 1
  for(i in 1:length(x))
  {
    if(res[[i]]$latency == type[1] && res[[i]]$late == type[2])
    {
      if(position == 'au')
      {
        temp <- cbind(temp, as.vector(res[[i]][[1]]))
        ntmp <- c(ntmp, res[[i]]$sogg)
      }
      else
      {
        if(length(res[[i]]) == 11)
        {
          temp <- cbind(temp, as.vector(res[[i]][[3]]))
          ntmp <- c(ntmp, res[[i]]$sogg)
        }
        else
        {
          if(is.numeric(res[[i]][[4]]))
          {
            temp <- cbind(temp, as.vector(res[[i]][[4]]))
            ntmp <- c(ntmp, res[[i]]$sogg)
          }
        }
      }
    }
  }
  colnames(temp) <- as.character(c(0,ntmp))
  temp[,-1]
}


#### PROVA DI GIULIA C ####

setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject")
load("functional_WP.RData")

library(fdakma)

x <- 1:1600
y0 <- t(importMatrix(res, type = c('SL', 'sx'), position = 'au'))

K=2

set.seed(4)
clust <- kma(
  x=x, y0=y0, n.clust = K, 
  warping.method = 'NOalignment', 
  similarity.method = 'd0.pearson',  
  center.method = 'k-means',
  n.out = dim(y0)[2]
)

#cluster indicators
c=clust$labels

#partition
C1=which(clust$labels==1)
C2=which(clust$labels==2)


n=dim(y0)[1]

#centroids
Xk=clust$y0.centers.final


matplot(t(y0), type='l', xlab='x', ylab='orig.func')

matplot(t(Xk), type='l', xlab='x', ylab='centroids')

#clusters without the ith unit
clust_min_i<-function(i,x,y0){
  clust_min_i <- kma(
    x=x, y0=y0[-i,], n.clust = K, 
    warping.method = 'NOalignment', 
    similarity.method = 'd0.pearson',  
    center.method = 'k-means',
    n.out = dim(y0)[2]
  )
  return (clust_min_i)
}

D=matrix(0,n,K)

for (i in 1:n){
  set.seed(4)
  c_mini=clust_min_i(i,x,y0)$labels              #cluster indicators without the ith unit
  C1_mini=which(c_mini==1)
  C2_mini=which(c_mini==2)
  Xk_mini=clust_min_i(i,x,y0)$y0.centers.final   #centroids without the ith unit
  y0_mini=y0[-i,]
  
  #prob ci=1
  d=0
  for (l in C1){
    d=d+dist(rbind(y0[l,],Xk[1,]),method='euclidean')
  }
  for (j in C1_mini){
    d=d-dist(rbind(y0_mini[j,],Xk_mini[1,]),method='euclidean')
  }
  D[i,1]=d
  
  
  #prob ci=2
  d=0
  for (l in C2){
    d=d+dist(rbind(y0[l,],Xk[2,]),method='euclidean')
  }
  for (j in C2_mini){
    d=d-dist(rbind(y0_mini[j,],Xk_mini[2,]),method='euclidean')
  }
  D[i,2]=d
  
}

D
lambda=10^(-4)
P=exp(-lambda*D)



#### CLUSTERING ON REAL DATA TRIAL ####

setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject")
load("functional_WP.RData")
library(fda.usc)
library(fda)
library(fields)

gibbs_loss <- function(n_clust, centroids, eig ,label ,data){
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
  
  loss_value1 <- gibbs_loss(n_clust = n_clust, centroids = centroids_random, 
                            eig = eig,label = c_lab, data = data)
  
  # update each centroid as the mean of the clusters data
  centroids_mean<-matrix(0,nrow = n_clust, ncol = t_points)
  for (k in 1:n_clust){
    if ( is.null(dim(data[which(c_lab==k),])[1]) == TRUE ) {
      centroids_mean[k,] <- data[which(c_lab==k),]
    }
    else centroids_mean[k,] <- colMeans(data[which(c_lab==k),])
  }
  
  loss_value2 <- gibbs_loss(n_clust = n_clust, centroids = centroids_mean, 
                            eig = eig, label = c_lab, data = data)
  
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
      if ( is.null(dim(data[which(c_lab==k),])[1]) == TRUE ){
        centroids_mean[k,] <- data[which(c_lab==k),]
      }
      else centroids_mean[k,] <- colMeans(data[which(c_lab==k),])
    }
    
    loss_value2 <- gibbs_loss(n_clust = n_clust, centroids = centroids_mean, eig = eig
                              ,label = c_lab, data = data)
  }
  
  return(list("label" = c_lab, "centroids" = centroids_mean, "loss" = loss_value2))
  
}

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









