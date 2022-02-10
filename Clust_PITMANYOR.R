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

