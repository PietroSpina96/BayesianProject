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







#### CLUSTERING ON SIMULATED DATA ####

setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject")
load('Simulated_WP.RData')
data<-data1
rm(data1)



#### Clustering function ####
y0.1 <- sample(1:n,1)
y0.2 <- sample(1:n,1)
while (y0.2 == y0.1) {
  y0.2 <- sample(1:n,1)
}

c_lab <- rep(0,n)

for (i in 1:n){
  if (Mahalanobis_Distance[y0.1,i] < Mahalanobis_Distance[y0.2,i]){
    c_lab[i] <- 1
  }
  if (Mahalanobis_Distance[y0.1,i] > Mahalanobis_Distance[y0.2,i])
    c_lab[i] <- 2
}



#### Loss function ####


gibbs_loss <- function(n_clust, centroid, label ,data){
  
}







