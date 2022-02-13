#### REAL DATA #### 

#### Setup ####

library(fda.usc)
library(fda)
library(fields)

setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject")
setwd("C:/Users/admin/Documents/R/Project_BS/BayesianProject") #GiuliaR

load("functional_WP.RData")
load('Functions_WP.RData')


#### Application on real data ####
n <- 26
t_points <- 1600
time<-seq(1,t_points)

# Plot of functional mean vs data
X_bar <- colMeans(f.data$ausxSL$data) # functional mean

x11()
plot(t(importMatrix(res, type = c('SL', 'sx'), position = 'au'))[1,], type = "l", ylim = range(f.data$ausxSL$data))
for(i in 1:26){
  lines(t(importMatrix(res, type = c('SL', 'sx'), position = 'au'))[i,], lwd=1)
}
lines(f.data$ausxSL$argvals, X_bar, type = 'l', lwd=3, col = 'firebrick2')


# Covariance matrix for the data: K_hat_R is the covariance matrix computed through the R function cov
K_hat_R <- cov(f.data$ausxSL$data) 

# Plot of the covariance matrix
x11()
image.plot(time,time,K_hat_R,main='R cov function')

# Plot of the covariance matrix K_hat and K_hat_R
x11()
par(mfrow=c(1,2))
image.plot(time,time,K_hat,main='Our cov function')
image.plot(time,time,K_hat_R,main='R cov function')

# Compute eigenvalues and eigenfunctions of the covariance matrix
# eigenvf<-eigen(K_hat)
# lambda<-eigenvf$values
# eigenft<-eigenvf$vectors

# Now the eigenvalues and eigenvectors of the K_hat_R matrix are called as the previous ones
# so that we don't have to rewrite them in the code. The workspace is updated
eigenvf<-eigen(K_hat_R)
values<-eigenvf$values
eigenft<-eigenvf$vectors

# OBS: from the theory we know that only (N-1) eigenvalues are not null where N=number of 
# statistical units. Here we have N=26 and 25 not null eigenvalues.
#lambda[1:26]

# Useless
# # prova
# kernel_estimator(f.data$ausxSL$data,X_bar,n,1,9,t_points) 
# 
# # Covariance matrix for the data
# K_hat <- matrix(0,t_points,t_points)
# for (i in 1:t_points){
#   for (j in 1:t_points){
#     K_hat[i,j] <- kernel_estimator(f.data$ausxSL$data,X_bar,n,i,j,t_points)
#   }
# }
# # sub-matrix (5x5)
# khat_example<-K_hat[1:5,1:5]
# show(khat_example)


#### Reduction of data points ####

time_red <- seq(1,1600,3)
f.data_red <- f.data$ausxSL$data[,time_red]
f.Data <- list('data' = f.data_red, 'argvals' = time_red, 'clinical' = data)
rm(time_red)
rm(f.data_red)

x11()
plot(f.Data$argvals,f.Data$data[1,], ylim = range(f.Data$data), type = 'l', lwd = 2, 
     xlab = 'time', ylab = 'EVOKED PONTENTIAL')
for (i in 2:26){
  lines(f.Data$argvals,f.Data$data[i,], lwd = 2)
}

K_hat_t <- cov(f.Data$data)
eigen_t <- eigen(K_hat_t)

max(eigenvf$values)
max(eigen_t$values)


##### alpha-Mahalanobis distance calculation ####
# Setting alpha
f_prova <- f.Data$data[10,]

x11()
plot(f.Data$argvals, f_prova, type = "l", ylim = range(f_prova), lwd=2)
lines(f.Data$argvals, f_alpha_approx(f_prova,1e+3,eigen_t$values,eigen_t$vectors), 
      type = 'l', lwd=2, col = 'firebrick2')
lines(f.Data$argvals, f_alpha_approx(f_prova,1e+4,eigen_t$values,eigen_t$vectors), 
      type = 'l', lwd=2, col = 'blue')
lines(f.Data$argvals, f_alpha_approx(f_prova,1e+5,eigen_t$values,eigen_t$vectors), 
      type = 'l', lwd=2, col = 'forestgreen')



# Smoothed data
alpha <- 1e+4
f.data_alpha <- matrix(0, nrow = 26, ncol = 534)
for (i in 1:26){
  f.data_alpha[i,] <- f_alpha_approx(f.Data$data[i,],alpha,eigen_t$values,eigen_t$vectors)
}

x11()
par(mfrow = c(1,2))
plot(f.Data$argvals,f.Data$data[1,], ylim = range(f.Data$data) ,
     type = 'l', lwd = 2, main = "DATA", xlab = 'time', ylab = 'EVOKED POTENTIAL')
for (i in 2:26){
  lines(f.Data$argvals, f.Data$data[i,], lwd = 2)
}
plot(f.Data$argvals, f.data_alpha[1,], ylim = range(f.Data$data), 
     type = 'l', lwd = 2, main = 'SMOOTHED DATA', , xlab = 'time', ylab = 'EVOKED POTENTIAL')
for (i in 2:26){
  lines(f.Data$argvals, f.data_alpha[i,], lwd = 2)
}


# alpha-Mahalanobis distance matrix 
Mahalanobis_Distance <- matrix(0, nrow = 26, ncol = 26)
for (i in 1:26){
  print(i)
  for (j in 1:26){
    Mahalanobis_Distance[i,j] <- alpha_Mahalanobis(alpha,f.Data$data[i,],
                                                   f.Data$data[j,],eigen_t$values,
                                                   eigen_t$vectors)
  }
}

x11()
image.plot(1:26,1:26,Mahalanobis_Distance)

#### Uniform prior with fixed covariance ####

##### k = 2 ######

f.data_clust_2 <- fda_clustering_mahalanobis_general(n_clust = 2, alpha = alpha,
                                                     cov_matrix = cov(f.Data$data),
                                                     cov_type = 'fixed', toll = 1e-2, 
                                                     data = f.Data$data)
c_opt_2 <- f.data_clust_2$label
show(c_opt_2)
show(f.data_clust_2$loss)

c1 <- f.data_clust_2$centroids[1,]
c2 <- f.data_clust_2$centroids[2,]

data1 <- f.Data$data[which(c_opt_2=='1'),]
rownames(data1) <- c(1:dim(data1)[1])
data2 <- f.Data$data[which(c_opt_2=='2'),]
rownames(data2) <- c(1:dim(data2)[1])

###### Plot #####
x11()
par(mfrow = c(1,2))
plot(f.Data$argvals,f.Data$data[1,], ylim = range(f.Data$data) ,type = 'l', lwd = 2, main = "Data", 
     xlab = 'time', ylab = 'EVOKED POTENTIAL')
for(i in 2:n){
  lines(f.Data$argvals,f.Data$data[i,],type = 'l',lwd = 2)
}

plot(f.Data$argvals,data1[1,], ylim = range(f.Data$data) ,type = 'l', col = 'firebrick3', lwd = 2, 
     main = "Uniform (fixed cov) k = 2", xlab = 'time', ylab = 'EVOKED POTENTIAL')
for (i in 2:dim(data1)[1]){
  lines(f.Data$argvals,data1[i,],type = 'l', col = 'firebrick3',lwd = 2)
}
for (i in 1:dim(data2)[1]){
  lines(f.Data$argvals,data2[i,],type = 'l', col = 'dodgerblue3',lwd = 2)
}
lines(f.Data$argvals, c1, lwd = 3)
lines(f.Data$argvals, c2, lwd = 3)


##### k = 3 ######
f.data_clust_3 <- fda_clustering_mahalanobis_general(n_clust = 3, alpha = alpha,
                                                     cov_matrix = cov(f.Data$data),
                                                     cov_type = 'fixed', toll = 1e-2, 
                                                     data = f.Data$data)
c_opt_3 <- f.data_clust_3$label
show(c_opt_3)
show(f.data_clust_3$loss)

c1 <- f.data_clust_3$centroids[1,]
c2 <- f.data_clust_3$centroids[2,]
c3 <- f.data_clust_3$centroids[3,]

data1 <- f.Data$data[which(c_opt_3=='1'),]
data2 <- f.Data$data[which(c_opt_3=='2'),]
data3 <- f.Data$data[which(c_opt_3=='3'),]

###### Plot #####
x11()
par(mfrow = c(1,2))
plot(f.Data$argvals,f.Data$data[1,], ylim = range(f.Data$data) ,type = 'l', lwd = 2, main = "Data",
     xlab = 'time', ylab = 'EVOKED POTENTIAL')
for(i in 2:n){
  lines(f.Data$argvals,f.Data$data[i,],type = 'l',lwd = 2)
}

plot(f.Data$argvals,data1[1,], ylim = range(f.Data$data) ,type = 'l', col = 'dodgerblue3', lwd = 2, 
     main = "Uniform (fixed cov) k = 3", xlab = 'time', ylab = 'EVOKED POTENTIAL')
# plot(f.Data$argvals,data1[1,], ylim = range(f.Data$data) ,type = 'l', col = 'gold', lwd = 2, main = "Uniform (fixed cov) k = 3")
for (i in 2:dim(data1)[1]){
  lines(f.Data$argvals,data1[i,],type = 'l', col = 'dodgerblue3',lwd = 2)
}
for (i in 1:dim(data2)[1]){
  lines(f.Data$argvals,data2[i,],type = 'l', col = 'forestgreen',lwd = 2)
}
for (i in 1:dim(data3)[1]){
  lines(f.Data$argvals,data3[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
lines(f.Data$argvals, c1, lwd = 3)
lines(f.Data$argvals, c2, lwd = 3)
lines(f.Data$argvals, c3, lwd = 3)


##### Cluster merging #####
library(ggplot2)
library(dplyr) # pipe (%>%)
library(tidyr) # gather()
library(tibble) # add_column()
library(gridExtra) # grid.arrange()
require(gtools) # combinations()

# PLOT data + centroids
plot1<-clusters_plot(f.Data$argvals, f.Data$data, c_opt_3) + labs(title="Clustered data") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot2<-clusters_plot(f.Data$argvals, f.data_clust_3) + labs(title="Clusters centroids") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# grid.arrange(plot1,plot2,nrow=1)
#
# DO merge
f.data_clust_3merge <- clusters_union(f.data_clust_3, data=f.Data$data)
# create colors for plot1 knowing the successive merging
kk <- f.data_clust_3merge$v %>% length
colors1 <- rainbow(kk)
for(i in 1:kk)
  colors1[i] <- colors1[f.data_clust_3merge$v[i]]
colors2 <- colors1 %>% unique
# show both plots
x11(width=10,height=6)
plot1.m <- clusters_plot(f.Data$argvals, f.data_clust_3, colors1) + labs(title="Before merging") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot2.m <- clusters_plot(f.Data$argvals, f.data_clust_3merge, colors2) + labs(title="After merging") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
grid.arrange(plot1,plot2, plot1.m,plot2.m, nrow=2)
#






#### Uniform prior with covariance updating within clusters ####

##### k = 2 ######

f.data_clust_2up <- fda_clustering_mahalanobis_general(n_clust = 2, alpha = alpha,
                                                     cov_matrix = cov(f.Data$data),
                                                     cov_type = 'updated', toll = 1e-1, 
                                                     data = f.Data$data)
c_opt_2up <- f.data_clust_2up$label
show(c_opt_2up)
show(f.data_clust_2up$loss)

c1 <- f.data_clust_2up$centroids[1,]
c2 <- f.data_clust_2up$centroids[2,]

data1 <- f.Data$data[which(c_opt_2up=='1'),]
data2 <- f.Data$data[which(c_opt_2up=='2'),]

###### Plot #####
x11()
par(mfrow = c(1,2))
plot(f.Data$argvals,f.Data$data[1,], ylim = range(f.Data$data) ,type = 'l', lwd = 2, main = "Data", 
     xlab = 'time', ylab = 'EVOKED POTENTIAL')
for(i in 2:n){
  lines(f.Data$argvals,f.Data$data[i,],type = 'l',lwd = 2)
}

plot(f.Data$argvals,data1[1,], ylim = range(f.Data$data) ,type = 'l', col = 'firebrick2', lwd = 2, 
     main = "Uniform (updated cov) k = 2", xlab = 'time', ylab = 'EVOKED POTENTIAL')
for (i in 2:dim(data1)[1]){
  lines(f.Data$argvals,data1[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in 1:dim(data2)[1]){
  lines(f.Data$argvals,data2[i,],type = 'l', col = 'dodgerblue3',lwd = 2)
}
lines(f.Data$argvals, c1, lwd = 3)
lines(f.Data$argvals, c2, lwd = 3)

##### k = 3 ######
f.data_clust_3up <- fda_clustering_mahalanobis_general(n_clust = 3, alpha = alpha,
                                                       cov_matrix = cov(f.Data$data),
                                                       cov_type = 'updated', toll = 1e-1, 
                                                       data = f.Data$data)
c_opt_3up <- f.data_clust_3up$label
show(c_opt_3up)
show(f.data_clust_3up$loss)

c1 <- f.data_clust_3up$centroids[1,]
c2 <- f.data_clust_3up$centroids[2,]
c3 <- f.data_clust_3up$centroids[3,]

data1 <- f.Data$data[which(c_opt_3up=='1'),]
data2 <- f.Data$data[which(c_opt_3up=='2'),]
data3 <- f.Data$data[which(c_opt_3up=='3'),]

###### Plot #####
x11()
par(mfrow = c(1,2))
plot(f.Data$argvals,f.Data$data[1,], ylim = range(f.Data$data) ,type = 'l', lwd = 2, main = "Data",
     xlab = 'time', ylab = 'EVOKED POTENTIAL')
for(i in 2:n){
  lines(f.Data$argvals,f.Data$data[i,],type = 'l',lwd = 2)
}

plot(f.Data$argvals,data1, ylim = range(f.Data$data) ,type = 'l', col = 'dodgerblue3', lwd = 2, 
     main = "Uniform (updated cov) k = 3", xlab = 'time', ylab = 'EVOKED POTENTIAL')
# plot(f.Data$argvals,data1[1,], ylim = range(f.Data$data) ,type = 'l', col = 'dodgerblue3', lwd = 2, 
#      main = "Uniform (updated cov) k = 3", xlab = 'time', ylab = 'EVOKED POTENTIAL')
# for (i in 2:dim(data1)[1]){
#   lines(f.Data$argvals,data1[i,],type = 'l', col = 'dodgerblue3',lwd = 2)
# }
for (i in 1:dim(data2)[1]){
  lines(f.Data$argvals,data2[i,],type = 'l', col = 'forestgreen',lwd = 2)
}
for (i in 1:dim(data3)[1]){
  lines(f.Data$argvals,data3[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
lines(f.Data$argvals, c1, lwd = 5, col = 'dodgerblue3')
lines(f.Data$argvals, c2, lwd = 3)
lines(f.Data$argvals, c3, lwd = 3)

##### Cluster merging #####
library(ggplot2)
library(dplyr) # pipe (%>%)
library(tidyr) # gather()
library(tibble) # add_column()
library(gridExtra) # grid.arrange()
require(gtools) # combinations()

# PLOT data + centroids
plot1<-clusters_plot(f.Data$argvals, f.Data$data, c_opt_3up) + labs(title="Clustered data") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot2<-clusters_plot(f.Data$argvals, f.data_clust_3up) + labs(title="Clusters centroids") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# grid.arrange(plot1,plot2,nrow=1)
#
# DO merge
f.data_clust_3upmerge <- clusters_union(f.data_clust_3up, data=f.Data$argvals)
# create colors for plot1 knowing the successive merging
colors1 <- rainbow(k)
for(i in 1:k)
  colors1[i] <- colors1[f.data_clust_3upmerge$v[i]]
colors2 <- colors1 %>% unique
# show both plots
x11(width=10,height=6)
plot1.m <- clusters_plot(f.Data$argvals, f.data_clust_3up, colors1) + labs(title="Before merging") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot2.m <- clusters_plot(f.Data$argvals, f.data_clust_3upmerge, colors2) + labs(title="After merging") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
grid.arrange(plot1,plot2, plot1.m,plot2.m, nrow=2)
#






#### Pitman-Yor EPPF Prior ####

alpha <- alpha
sigma <- 0.75
theta <- 3.66
lambda <- 1

f.data_clust_py <- fda_clustering_pitmanyor_overall(n_clust = 5, nsimul = 10, alpha = alpha, 
                                                    sigma = sigma, theta = theta, lambda = lambda,
                                                    cov_matrix = cov(f.Data$data), data = f.Data$data)

f.data_clust_py$posterior_all_k
f.data_clust_py$posterior
f.data_clust_py$clusters_number
f.data_clust_py$loss
c_opt_py <- f.data_clust_py$labels
show(c_opt_py)







#### Clinical relevancy #####

##### GOSE ####

gose1 <- f.Data$data[which(f.Data$clinical$GOSE == '1'),]
gose2 <- f.Data$data[which(f.Data$clinical$GOSE == '2'),]
index_gose1 <- which(f.Data$clinical$GOSE == 1)

###### k = 2 Uniform prior with fixed covariance structure ####

data1 <- f.Data$data[which(c_opt_2=='1'),]
data2 <- f.Data$data[which(c_opt_2=='2'),]

show(c_opt_2)
show(f.Data$clinical$GOSE)
show(c_opt_2[index_gose1]) # brutto
show(c_opt_2[-index_gose1])

x11()
par(mfrow = c(1,2))
plot(f.Data$argvals, gose2[1,], ylim = range(f.Data$data) ,type = 'l', lwd = 2, main = "GOSE", 
     xlab = 'time', ylab = 'EVOKED POTENTIAL', col = 'grey')
for(i in 2:dim(gose2)[1]){
  lines(f.Data$argvals, gose2[i,], type = 'l', lwd = 2, col = 'grey')
}
for(i in 1:dim(gose1)[1]){
  lines(f.Data$argvals, gose1[i,], type = 'l', lwd = 2, col = 'firebrick2')
}

plot(f.Data$argvals,data1[1,], ylim = range(f.Data$data) ,type = 'l', col = 'firebrick2', lwd = 2, 
     main = "Uniform (fixed cov) k = 2", xlab = 'time', ylab = 'EVOKED POTENTIAL')
for (i in 2:dim(data1)[1]){
  lines(f.Data$argvals,data1[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in 1:dim(data2)[1]){
  lines(f.Data$argvals,data2[i,],type = 'l', col = 'grey',lwd = 2)
}


###### k = 3 Uniform prior with fixed covariance structure #####

data1 <- f.Data$data[which(c_opt_3=='1'),]
data2 <- f.Data$data[which(c_opt_3=='2'),]
data3 <- f.Data$data[which(c_opt_3=='3'),]

show(c_opt_3)
show(f.Data$clinical$GOSE)
show(c_opt_3[index_gose1]) # bruttissimo 
show(c_opt_3[-index_gose1])

x11()
par(mfrow = c(1,2))
plot(f.Data$argvals, gose2[1,], ylim = range(f.Data$data) ,type = 'l', lwd = 2, main = "GOSE", 
     xlab = 'time', ylab = 'EVOKED POTENTIAL', col = 'grey')
for(i in 2:dim(gose2)[1]){
  lines(f.Data$argvals, gose2[i,], type = 'l', lwd = 2, col = 'grey')
}
for(i in 1:dim(gose1)[1]){
  lines(f.Data$argvals, gose1[i,], type = 'l', lwd = 2, col = 'firebrick2')
}

plot(f.Data$argvals,data1[1,], ylim = range(f.Data$data) ,type = 'l', col = 'grey', lwd = 2, 
     main = "Uniform (fixed cov) k = 3", xlab = 'time', ylab = 'EVOKED POTENTIAL')
for (i in 2:dim(data1)[1]){
  lines(f.Data$argvals,data1[i,],type = 'l', col = 'grey',lwd = 2)
}
for (i in 1:dim(data2)[1]){
  lines(f.Data$argvals,data2[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in 1:dim(data3)[1]){
  lines(f.Data$argvals,data3[i,],type = 'l', col = 'orange2',lwd = 2)
}


###### k = 2 Uniform prior with updating covariance within clusters #####

data1 <- f.Data$data[which(c_opt_2up=='1'),]
data2 <- f.Data$data[which(c_opt_2up=='2'),]

show(c_opt_2up)
show(f.Data$clinical$GOSE)
show(c_opt_2up[index_gose1]) # brutto
show(c_opt_2up[-index_gose1]) 

x11()
par(mfrow = c(1,2))
plot(f.Data$argvals, gose2[1,], ylim = range(f.Data$data) ,type = 'l', lwd = 2, main = "GOSE", 
     xlab = 'time', ylab = 'EVOKED POTENTIAL', col = 'grey')
for(i in 2:dim(gose2)[1]){
  lines(f.Data$argvals, gose2[i,], type = 'l', lwd = 2, col = 'grey')
}
for(i in 1:dim(gose1)[1]){
  lines(f.Data$argvals, gose1[i,], type = 'l', lwd = 2, col = 'firebrick2')
}

plot(f.Data$argvals,data2[1,], ylim = range(f.Data$data) ,type = 'l', col = 'firebrick2', lwd = 2, 
     main = "Uniform (updated cov) k = 2", xlab = 'time', ylab = 'EVOKED POTENTIAL')
for (i in 2:dim(data2)[1]){
  lines(f.Data$argvals,data2[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in 1:dim(data1)[1]){
  lines(f.Data$argvals,data1[i,],type = 'l', col = 'grey',lwd = 2)
}


###### k = 3 Uniform prior with updating covariance within clusters ####

data1 <- f.Data$data[which(c_opt_3up=='1'),]
data2 <- f.Data$data[which(c_opt_3up=='2'),]
data3 <- f.Data$data[which(c_opt_3up=='3'),]

show(c_opt_3up)
show(f.Data$clinical$GOSE)
show(c_opt_3up[index_gose1]) # bruttissimo 
show(c_opt_3up[-index_gose1])

x11()
par(mfrow = c(1,2))
plot(f.Data$argvals, gose2[1,], ylim = range(f.Data$data) ,type = 'l', lwd = 2, main = "GOSE", 
     xlab = 'time', ylab = 'EVOKED POTENTIAL', col = 'grey')
for(i in 2:dim(gose2)[1]){
  lines(f.Data$argvals, gose2[i,], type = 'l', lwd = 2, col = 'grey')
}
for(i in 1:dim(gose1)[1]){
  lines(f.Data$argvals, gose1[i,], type = 'l', lwd = 2, col = 'firebrick2')
}

plot(f.Data$argvals,data1, ylim = range(f.Data$data) ,type = 'l', col = 'grey', lwd = 2, 
     main = "Uniform (updated cov) k = 3", xlab = 'time', ylab = 'EVOKED POTENTIAL')
for (i in 1:dim(data2)[1]){
  lines(f.Data$argvals,data2[i,],type = 'l', col = 'orange2',lwd = 2)
}
for (i in 1:dim(data3)[1]){
  lines(f.Data$argvals,data3[i,],type = 'l', col = 'firebrick2',lwd = 2)
}


###### Pitman-Yor EPPF prior #####



##### LCF ####

lcf1 <- f.Data$data[which(f.Data$clinical$LCF == '1'),]
lcf2 <- f.Data$data[which(f.Data$clinical$LCF == '2'),]
index_lcf1 <- which(f.Data$clinical$LCF == 1)

###### k = 2 Uniform prior with fixed covariance structure ####

data1 <- f.Data$data[which(c_opt_2=='1'),]
data2 <- f.Data$data[which(c_opt_2=='2'),]

show(c_opt_2)
show(f.Data$clinical$LCF)
show(c_opt_2[index_lcf1]) # brutto
show(c_opt_2[-index_lcf1])

x11()
par(mfrow = c(1,2))
plot(f.Data$argvals, lcf2[1,], ylim = range(f.Data$data) ,type = 'l', lwd = 2, main = "LCF", 
     xlab = 'time', ylab = 'EVOKED POTENTIAL', col = 'grey')
for(i in 2:dim(lcf2)[1]){
  lines(f.Data$argvals, lcf2[i,], type = 'l', lwd = 2, col = 'grey')
}
for(i in 1:dim(lcf1)[1]){
  lines(f.Data$argvals, lcf1[i,], type = 'l', lwd = 2, col = 'firebrick2')
}

plot(f.Data$argvals,data1[1,], ylim = range(f.Data$data) ,type = 'l', col = 'firebrick2', lwd = 2, 
     main = "Uniform (fixed cov) k = 2", xlab = 'time', ylab = 'EVOKED POTENTIAL')
for (i in 2:dim(data1)[1]){
  lines(f.Data$argvals,data1[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in 1:dim(data2)[1]){
  lines(f.Data$argvals,data2[i,],type = 'l', col = 'grey',lwd = 2)
}


###### k = 3 Uniform prior with fixed covariance structure #####

data1 <- f.Data$data[which(c_opt_3=='1'),]
data2 <- f.Data$data[which(c_opt_3=='2'),]
data3 <- f.Data$data[which(c_opt_3=='3'),]

show(c_opt_3)
show(f.Data$clinical$LCF)
show(c_opt_3[index_lcf1]) # cluster 2 and 3 seem to show a trend but remember that 
                          # cluster 1 has only two observations
show(c_opt_3[-index_lcf1])

x11()
par(mfrow = c(1,2))
plot(f.Data$argvals, lcf2[1,], ylim = range(f.Data$data) ,type = 'l', lwd = 2, main = "LCF", 
     xlab = 'time', ylab = 'EVOKED POTENTIAL', col = 'grey')
for(i in 2:dim(lcf2)[1]){
  lines(f.Data$argvals, lcf2[i,], type = 'l', lwd = 2, col = 'grey')
}
for(i in 1:dim(lcf1)[1]){
  lines(f.Data$argvals, lcf1[i,], type = 'l', lwd = 2, col = 'firebrick2')
}

plot(f.Data$argvals,data1[1,], ylim = range(f.Data$data) ,type = 'l', col = 'grey', lwd = 2, 
     main = "Uniform (fixed cov) k = 3", xlab = 'time', ylab = 'EVOKED POTENTIAL')
for (i in 2:dim(data1)[1]){
  lines(f.Data$argvals,data1[i,],type = 'l', col = 'grey',lwd = 2)
}
for (i in 1:dim(data2)[1]){
  lines(f.Data$argvals,data2[i,],type = 'l', col = 'orange2',lwd = 2)
}
for (i in 1:dim(data3)[1]){
  lines(f.Data$argvals,data3[i,],type = 'l', col = 'firebrick2',lwd = 2)
}


###### k = 2 Uniform prior with updating covariance within clusters #####

data1 <- f.Data$data[which(c_opt_2up=='1'),]
data2 <- f.Data$data[which(c_opt_2up=='2'),]

show(c_opt_2up)
show(f.Data$clinical$GOSE)
show(c_opt_2up[index_lcf1]) # brutto 
show(c_opt_2up[-index_lcf1])

x11()
par(mfrow = c(1,2))
plot(f.Data$argvals, lcf2[1,], ylim = range(f.Data$data) ,type = 'l', lwd = 2, main = "LCF", 
     xlab = 'time', ylab = 'EVOKED POTENTIAL', col = 'grey')
for(i in 2:dim(lcf2)[1]){
  lines(f.Data$argvals, lcf2[i,], type = 'l', lwd = 2, col = 'grey')
}
for(i in 1:dim(lcf1)[1]){
  lines(f.Data$argvals, lcf1[i,], type = 'l', lwd = 2, col = 'firebrick2')
}

plot(f.Data$argvals,data1[1,], ylim = range(f.Data$data) ,type = 'l', col = 'firebrick2', lwd = 2, 
     main = "Uniform (updated cov) k = 2", xlab = 'time', ylab = 'EVOKED POTENTIAL')
for (i in 2:dim(data1)[1]){
  lines(f.Data$argvals,data1[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in 1:dim(data2)[1]){
  lines(f.Data$argvals,data2[i,],type = 'l', col = 'grey',lwd = 2)
}


###### k = 3 Uniform prior with updating covariance within clusters ####

data1 <- f.Data$data[which(c_opt_3=='1'),]
data2 <- f.Data$data[which(c_opt_3=='2'),]
data3 <- f.Data$data[which(c_opt_3=='3'),]

show(c_opt_3up)
show(f.Data$clinical$LCF)
show(c_opt_3up[index_lcf1]) # brutto
show(c_opt_3up[-index_lcf1])

x11()
par(mfrow = c(1,2))
plot(f.Data$argvals, lcf2[1,], ylim = range(f.Data$data) ,type = 'l', lwd = 2, main = "LCF", 
     xlab = 'time', ylab = 'EVOKED POTENTIAL', col = 'grey')
for(i in 2:dim(lcf2)[1]){
  lines(f.Data$argvals, lcf2[i,], type = 'l', lwd = 2, col = 'grey')
}
for(i in 1:dim(lcf1)[1]){
  lines(f.Data$argvals, lcf1[i,], type = 'l', lwd = 2, col = 'firebrick2')
}

plot(f.Data$argvals, data1[1,], ylim = range(f.Data$data) ,type = 'l', col = 'grey', lwd = 2, 
     main = "Uniform (updated cov) k = 3", xlab = 'time', ylab = 'EVOKED POTENTIAL')
for (i in 1:dim(data2)[1]){
  lines(f.Data$argvals,data2[i,],type = 'l', col = 'orange2',lwd = 2)
}
for (i in 1:dim(data3)[1]){
  lines(f.Data$argvals,data3[i,],type = 'l', col = 'firebrick2',lwd = 2)
}


###### Pitman-Yor EPPF prior #####



##### DRS #########
drs1 <- f.Data$data[which(f.Data$clinical$DRS == '1'),]
drs2 <- f.Data$data[which(f.Data$clinical$DRS == '2'),]
index_drs1 <- which(f.Data$clinical$DRS == 1)

###### k = 2 Uniform prior with fixed covariance structure ####

data1 <- f.Data$data[which(c_opt_2=='1'),]
data2 <- f.Data$data[which(c_opt_2=='2'),]

show(c_opt_2)
show(f.Data$clinical$DRS)
show(c_opt_2[index_drs1]) # brutto
show(c_opt_2[-index_drs1])

x11()
par(mfrow = c(1,2))
plot(f.Data$argvals, drs2[1,], ylim = range(f.Data$data) ,type = 'l', lwd = 2, main = "DSR", 
     xlab = 'time', ylab = 'EVOKED POTENTIAL', col = 'grey')
for(i in 2:dim(drs2)[1]){
  lines(f.Data$argvals, drs2[i,], type = 'l', lwd = 2, col = 'grey')
}
for(i in 1:dim(drs1)[1]){
  lines(f.Data$argvals, drs1[i,], type = 'l', lwd = 2, col = 'firebrick2')
}

plot(f.Data$argvals,data1[1,], ylim = range(f.Data$data) ,type = 'l', col = 'firebrick2', lwd = 2, 
     main = "Uniform (fixed cov) k = 2", xlab = 'time', ylab = 'EVOKED POTENTIAL')
for (i in 2:dim(data1)[1]){
  lines(f.Data$argvals,data1[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in 1:dim(data2)[1]){
  lines(f.Data$argvals,data2[i,],type = 'l', col = 'grey',lwd = 2)
}


###### k = 3 Uniform prior with fixed covariance structure #####

data1 <- f.Data$data[which(c_opt_3=='1'),]
data2 <- f.Data$data[which(c_opt_3=='2'),]
data3 <- f.Data$data[which(c_opt_3=='3'),]

show(c_opt_3)
show(f.Data$clinical$DRS)
show(c_opt_3[index_drs1]) # bruttissimo 
show(c_opt_3[-index_drs1])

x11()
par(mfrow = c(1,2))
plot(f.Data$argvals, drs2[1,], ylim = range(f.Data$data) ,type = 'l', lwd = 2, main = "DSR", 
     xlab = 'time', ylab = 'EVOKED POTENTIAL', col = 'grey')
for(i in 2:dim(drs2)[1]){
  lines(f.Data$argvals, drs2[i,], type = 'l', lwd = 2, col = 'grey')
}
for(i in 1:dim(drs1)[1]){
  lines(f.Data$argvals, drs1[i,], type = 'l', lwd = 2, col = 'firebrick2')
}

plot(f.Data$argvals,data1[1,], ylim = range(f.Data$data) ,type = 'l', col = 'grey', lwd = 2, 
     main = "Uniform (fixed cov) k = 3", xlab = 'time', ylab = 'EVOKED POTENTIAL')
for (i in 2:dim(data1)[1]){
  lines(f.Data$argvals,data1[i,],type = 'l', col = 'grey',lwd = 2)
}
for (i in 1:dim(data2)[1]){
  lines(f.Data$argvals,data2[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in 1:dim(data3)[1]){
  lines(f.Data$argvals,data3[i,],type = 'l', col = 'orange2',lwd = 2)
}


###### k = 2 Uniform prior with updating covariance within clusters #####

data1 <- f.Data$data[which(c_opt_2up=='1'),]
data2 <- f.Data$data[which(c_opt_2up=='2'),]

show(c_opt_2up)
show(f.Data$clinical$DRS)
show(c_opt_2up[index_drs1]) # brutto
show(c_opt_2up[-index_drs1]) 

x11()
par(mfrow = c(1,2))
plot(f.Data$argvals, drs2[1,], ylim = range(f.Data$data) ,type = 'l', lwd = 2, main = "DSR", 
     xlab = 'time', ylab = 'EVOKED POTENTIAL', col = 'grey')
for(i in 2:dim(drs2)[1]){
  lines(f.Data$argvals, drs2[i,], type = 'l', lwd = 2, col = 'grey')
}
for(i in 1:dim(drs1)[1]){
  lines(f.Data$argvals, drs1[i,], type = 'l', lwd = 2, col = 'firebrick2')
}

plot(f.Data$argvals,data2[1,], ylim = range(f.Data$data) ,type = 'l', col = 'firebrick2', lwd = 2, 
     main = "Uniform (updated cov) k = 2", xlab = 'time', ylab = 'EVOKED POTENTIAL')
for (i in 2:dim(data2)[1]){
  lines(f.Data$argvals,data2[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in 1:dim(data1)[1]){
  lines(f.Data$argvals,data1[i,],type = 'l', col = 'grey',lwd = 2)
}


###### k = 3 Uniform prior with updating covariance within clusters ####

data1 <- f.Data$data[which(c_opt_3up=='1'),]
data2 <- f.Data$data[which(c_opt_3up=='2'),]
data3 <- f.Data$data[which(c_opt_3up=='3'),]

show(c_opt_3up)
show(f.Data$clinical$DRS)
show(c_opt_3up[index_drs1]) # bruttissimo 
show(c_opt_3up[-index_drs1])

x11()
par(mfrow = c(1,2))
plot(f.Data$argvals, drs2[1,], ylim = range(f.Data$data) ,type = 'l', lwd = 2, main = "DSR", 
     xlab = 'time', ylab = 'EVOKED POTENTIAL', col = 'grey')
for(i in 2:dim(drs2)[1]){
  lines(f.Data$argvals, drs2[i,], type = 'l', lwd = 2, col = 'grey')
}
for(i in 1:dim(drs1)[1]){
  lines(f.Data$argvals, drs1[i,], type = 'l', lwd = 2, col = 'firebrick2')
}

plot(f.Data$argvals,data1, ylim = range(f.Data$data) ,type = 'l', col = 'grey', lwd = 2, 
     main = "Uniform (updated cov) k = 3", xlab = 'time', ylab = 'EVOKED POTENTIAL')
for (i in 1:dim(data2)[1]){
  lines(f.Data$argvals,data2[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in 1:dim(data3)[1]){
  lines(f.Data$argvals,data3[i,],type = 'l', col = 'orange2',lwd = 2)
}


###### Pitman-Yor EPPF prior #####






#### Save Workspace ####
setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject")
save.image("functional_WP.RData")

save.image("~/R/Project_BS/BayesianProject/functional_WP.RData") #GiuliaR


