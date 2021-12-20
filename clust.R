
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


################################################

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

###############################################

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


########################################################################################
########################################################################################
#### CLUSTERING ON SIMULATED DATA: Cmap ####

#setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject")
setwd("C:/Users/admin/Documents/R/Project_BS/BayesianProject") #GiuliaR

load('Simulated_WP.RData')
data<-data1
rm(data1)


#### Loss function ####
gibbs_loss <- function(n_clust, centroids, label ,data){
  res = rep(0,n_clust)
  
  for (k in 1:n_clust){
    for (i in 1:n){
      if (label[i] == k){
        sum_partial = alpha_Mahalanobis(alpha,data[i,],centroids[[k]],eig$values,eig$vectors)
        res[k] = res[k] + sum_partial
      }
    }
  }
  
  tot = sum(res)
  return(tot)
}


#### Clustering function ####
# The parameter eig corresponds to the output of the eigen function (list of eigenvalues and eigenvectors)
# The function works with n_clust=2 and not generic k ( for the moment)
# alpha is the smoothing parameter
# toll is the tolerance for the while loop

fda_clustering_mahalanobis <- function(n_clust, alpha, eig, toll,data){
  
  n <- dim(data)[1]
  
  # index of each centroid randomly defined through sampling
  y0.1 <- sample(1:n,1)
  y0.2 <- sample(1:n,1)
  
  while (y0.2 == y0.1) {
    y0.2 <- sample(1:n,1)
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
  
  # i-th unit belongs to cluster1 if the distance(centroid1, i-th unit) is less than the distance(centroid2,i-th unit)
  for (i in 1:n){
    if (Mahalanobis_Distance[y0.1,i] < Mahalanobis_Distance[y0.2,i]){
      c_lab[i] <- 1
    }
    if (Mahalanobis_Distance[y0.1,i] >= Mahalanobis_Distance[y0.2,i])
      c_lab[i] <- 2
  }
  
  loss_value1 <- gibbs_loss(n_clust = n_clust, centroids = list(data[y0.1,], data[y0.2,]), label = c_lab, data = data)

  # update each centroid as the mean of the clusters data
  centroid1 <- colMeans(data[which(c_lab=='1'),])
  centroid2 <- colMeans(data[which(c_lab=='2'),])
  
  loss_value2 <- gibbs_loss(n_clust = n_clust, centroids = list(centroid1, centroid2), label = c_lab, data = data)
  
  while(abs(loss_value1 - loss_value2) >= toll){
    c_lab <- rep(0,n)
    
    for (i in 1:n){
      if (alpha_Mahalanobis(alpha,centroid1,data[i,],values,vectors) < alpha_Mahalanobis(alpha,centroid2,data[i,],values,vectors)){
        c_lab[i] <- 1
      }
      else c_lab[i] <- 2
    }
    loss_value1 <- loss_value2
    
    centroid1 <- colMeans(data[which(c_lab=='1'),])
    centroid2 <- colMeans(data[which(c_lab=='2'),])
    
    loss_value2 <- gibbs_loss(n_clust = n_clust, centroids = list(centroid1, centroid2), label = c_lab, data = data)
  }
  
  return(list("label" = c_lab, "centroid1" = centroid1, "centroid2" = centroid2))
  
}

# Application on the simulated data
k <- 2
clust <- fda_clustering_mahalanobis(n_clust = k, alpha = alpha, eig = eig, toll = 1e-6, data = data)
c_opt <- clust$label
c1 <- clust$centroid1
c2 <- clust$centroid2
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





