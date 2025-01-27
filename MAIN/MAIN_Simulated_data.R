#### UNIFORM PRIOR MODEL ####

#### SETUP ####
setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject") #Pietro
# setwd("C:/Users/admin/Documents/R/Project_BS/BayesianProject") #GiuliaR


load('Workspaces/Simulated_WP.RData')
load('Workspaces/Functions_WP.RData')


# Choose the model for the simulation:
# simulated data - model 1
data <- data1
time <- time123
rm(list=c('data1','data2','data3','data4','data5'))
rm(list=c('time123','time45'))
rm(list=c('K_2','K_3'))

# simulated data - model 2
data <- data2
time <- time123
rm(list=c('data1','data2','data3','data4','data5'))
rm(list=c('time123','time45'))
rm(list=c('K_1','K_3'))

# simulated data - model 3
data <- data3
time <- time123
rm(list=c('data1','data2','data3','data4','data5'))
rm(list=c('time123','time45'))
rm(list=c('K_2','K_1'))

# simulated data - model 4
# The covariance matrix of data4 is cov_4 and K_4_1 and K_4_2 are the covariance matrices of the two clusters
data <- data4
time <- time45
rm(list=c('data1','data2','data3','data4','data5'))
rm(list=c('time123','time45'))
rm(list=c('K_2','K_3','K_1'))

# simulated data - model 5
# The covariance matrix of data5 is cov_5 and K_5_1 and K_5_2 are the covariance matrices of the two clusters
data <- data5
time <- time45
rm(list=c('data1','data2','data3','data4','data5'))
rm(list=c('time123','time45'))
rm(list=c('K_2','K_3','K_1'))

c <- n_cluster2     #number of items in the second cluster for model 1,2,3
n1 <- n_cluster1    #number of items in the first cluster for model 4,5
rm(list=c('n_cluster1','n_cluster2'))

##### Fixed covariance structure ####
n <- dim(data)[1]
t_points <- dim(data)[2]

k <- 2
alpha <- 0.1

clust <- fda_clustering_mahalanobis_general(n_clust = k, alpha = alpha,
                                            cov_matrix = cov(data), cov_type = 'updated',
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

###### Plot #####
x11()
par(mfrow = c(1,2))
plot(time,data[1,],type = 'l', ylim = c(-3.5,9), lwd = 2, col='black',main = "Main and contaminated processes", ylab='Data',xlab='Time')
for(i in 2:(n-c)){
  lines(time,data[i,],type = 'l', col = 'black',lwd = 2)
}
for (i in (n-c+1):n){
  lines(time,data[i,],type = 'l', col = 'black', lwd = 2)
}
#legend(x=0.75,y=9.5,ncol=1,box.lwd=1,legend=c('Process 1','Process 2'),fill=c('blue','firebrick2'),x.intersp=0.3,
       #text.col=c('blue','firebrick2'))

plot(time,data1[1,],type = 'l', ylim = c(-3.5,9), col = 'firebrick2', lwd = 2, main = "Clustered data",ylab='Data',xlab='Time')
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
legend(x=0.6,y=9.5,ncol=1,box.lwd=1,legend=c('Main process','Contaminated process','Centroids'),fill=c('blue','firebrick2','black'),x.intersp=0.3,
       text.col=c('blue','firebrick2','black'))

x11()
par(mfrow = c(1,2))
plot(time,data[1,],type = 'l', ylim = c(-3.5,9), lwd = 2, col='black',main = "Main and contaminated processes", ylab='Data',xlab='Time')
for(i in 2:n1){
  lines(time,data[i,],type = 'l', col = 'black',lwd = 2)
}
for (i in (n1+1):n){
  lines(time,data[i,],type = 'l', col = 'black', lwd = 2)
}
#legend(x=0.75,y=9.5,ncol=1,box.lwd=1,legend=c('Process 1','Process 2'),fill=c('blue','firebrick2'),x.intersp=0.3,
#text.col=c('blue','firebrick2'))

plot(time,data1[1,],type = 'l', ylim = c(-3.5,9), col = 'firebrick2', lwd = 2, main = "Clustered data",ylab='Data',xlab='Time')
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

legend(x=0.75,y=9.5,ncol=1,box.lwd=1,legend=c('Process 1','Process 2','Centroids'),fill=c('firebrick2','blue','black'),x.intersp=0.3,
       text.col=c('firebrick2','blue','black'))
#legend(x=0.65,y=9.5,ncol=1,box.lwd=1,legend=c('Main process','Contaminated process','Centroids'),fill=c('firebrick2','blue','black'),x.intersp=0.3,
#text.col=c('firebrick2','blue','black'))



rm(data1)
rm(data2)
# rm(data3)
# rm(data4)


###### Cluster Merging (a posteriori) ####
# fda_clustering_mahalanobis + clusters_union + clusters_plot
n <- dim(data)[1]
t_points <- dim(data)[2]
k <- 3
alpha <- 0.1
clust1 <- fda_clustering_mahalanobis_general(n_clust = k, alpha = alpha,
                                             cov_matrix = cov(data), cov_type = 'fixed',
                                             toll = 1e-2,  data = data)
c_opt <- clust1$label # labels
c1 <- clust1$centroids[1,]
c2 <- clust1$centroids[2,]
# PLOT data + centroids
plot1<-clusters_plot(time, data, c_opt) + labs(title="Clustered data") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot2<-clusters_plot(time, clust1) + labs(title="Clusters centroids") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# grid.arrange(plot1,plot2,nrow=1)
#
# DO merge
clust2 <- clusters_union(clust1,data=data)
# create colors for plot1 knowing the successive merging
colors1 <- rainbow(k)
for(i in 1:k)
  colors1[i] <- colors1[clust2$v[i]]
colors2 <- colors1 %>% unique
# show both plots
x11(width=10,height=6)
plot1.m <- clusters_plot(time, clust1, colors1) + labs(title="Before merging") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot2.m <- clusters_plot(time, clust2, colors2) + labs(title="After merging") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
grid.arrange(plot1,plot2, plot1.m,plot2.m, nrow=2)
#
 
##### Updating covariance within clusters ####
n <- dim(data)[1]
t_points <- dim(data)[2]

k <- 2
alpha <- 0.1

clust <- fda_clustering_mahalanobis_updated(n_clust = k, alpha = alpha, cov_matrix = cov(data),
                                            toll = 1e-2,  data = data)
c_opt <- clust$label
loss_min <- clust$loss
show(c_opt)  #label switching 
show(loss_min)

# Theoretical optimal plot vs clustering plot
c1 <- clust$centroids[1,]
c2 <- clust$centroids[2,]
c3 <- clust$centroids[3,]

data1 <- data[which(c_opt=='1'),]
data2 <- data[which(c_opt=='2'),]
data3 <- data[which(c_opt=='3'),]



#### PITMAN-YOR MODEL ####

#### SETUP ####
setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject") #Pietro
# setwd("C:/Users/admin/Documents/R/Project_BS/BayesianProject") #GiuliaR

load('Simulated_WP.RData')
load('Functions_WP.RData')


# Choose the model for the simulation:
# simulated data - model 1
data <- data1
time <- time123
rm(list=c('data1','data2','data3','data4','data5'))
rm(list=c('time123','time45'))
rm(list=c('K_2','K_3'))

# simulated data - model 2
data <- data2
time <- time123
rm(list=c('data1','data2','data3','data4','data5'))
rm(list=c('time123','time45'))
rm(list=c('K_1','K_3'))

# simulated data - model 3
data <- data3
time <- time123
rm(list=c('data1','data2','data3','data4','data5'))
rm(list=c('time123','time45'))
rm(list=c('K_2','K_1'))

n2 <- n_cluster2     #number of items in the second cluster for model 1,2,3
n1 <- n_cluster1    #number of items in the first cluster for model 4,5
rm(list=c('n_cluster1','n_cluster2'))


##### Application on simulated data ####
# Parameters
alpha <- 0.1
sigma <- 0.50 #0.25
theta <- -0.4 #3.66
lambda <- 0.75 #0.75

# To execute the algorithm for a specific cluster number
posterior_k <- fda_clustering_pitmanyor_kfixed(n_clust = 2 ,nsimul = 20,
                                               alpha = alpha, sigma = sigma, theta = theta,
                                               lambda = lambda,cov_matrix = K_1, data = data)

post_value_k <- posterior_k$posterior
post_labels_k <- posterior_k$labels
centroids <- posterior_k$centroids

show(post_value_k)
show(post_labels_k)


# To execute the algorithm for more clusters. Fix the clusters and simulations number
posterior_overall <- fda_clustering_pitmanyor_overall(n_clust = 4 ,nsimul = 10,
                                                      alpha = alpha, sigma = sigma, theta = theta,
                                                      lambda = lambda,cov_matrix = K_1, data = data)

post_value_allk <- posterior_overall$posterior_all_k

post_value_overall <- posterior_overall$posterior
post_labels <- posterior_overall$labels
number_clusters<- posterior_overall$clusters_number
centroids <- posterior_overall$centroids

show(post_value_allk)
show(post_value_overall)
show(post_labels)
show(number_clusters)


##### Plot ####
c1 <- centroids[1,]
c2 <- centroids[2,]
c3 <- centroids[3,]
#c4 <- centroids[4,]


# Theoretical optimal plot vs clustering plot
data1 <- data[which(post_labels=='1'),]
data2 <- data[which(post_labels=='2'),]
data3 <- data[which(post_labels=='3'),]
#data4 <- data[which(post_labels=='4'),]

# Plot for simulated data 1,2,3 
n <- dim(data)[1]
x11()
par(mfrow = c(1,2))
plot(time,data[1,],type = 'l', ylim = c(-3.5,9), lwd = 2, col='black',main = "Main and contaminated processes")
for(i in 2:n2){
  lines(time,data[i,],type = 'l', col = 'black',lwd = 2)
}
for (i in (n2+1):n){
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





