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


#### Posterior on Data 1 with fixed covariance #################################################
# For now this is just a test to see if the new function works

n <- dim(data)[1]
n_clust <- k <- 2
alpha <- 0.1
lambda <- 0.75

# Itialize eigen object
eig <- eigen(K_1)

# Values initialized by considering example at page 29 of the paper
sigma <- 0.25
theta <- 3.66

# Initialize cetroid matrix (bad clusters)
c <- 10
c_lab_bad <- c(rep(1,n-c), rep(2,c))
centroids_mean <- matrix(0, nrow = 2, ncol = 200)
centroids_mean[1,] <- colMeans(data[1:(n-c),])
centroids_mean[2,] <- colMeans(data[(n-c+1):n,])

post_bad <- posterior_pitman_yor(sigma = sigma, theta = theta, label = c_lab_bad, data = data)

# Initialize cetroid matrix (good clusters)
c <- 20
c_lab_good <- c(rep(1,n-c), rep(2,c))
centroids_mean <- matrix(0, nrow = 2, ncol = 200)
centroids_mean[1,] <- colMeans(data[1:(n-c),])
centroids_mean[2,] <- colMeans(data[(n-c+1):n,])

post_good <- posterior_pitman_yor(sigma = sigma, theta = theta, label = c_lab_good, data = data)

# Check 
post_bad$posterior
post_good$posterior
post_good$posterior > post_bad$posterior # as expected the optimal partition has higher posterior

# Label switching check
c <- 20
c_lab_good_2 <- c(rep(2,n-c), rep(1,c))
centroids_mean <- matrix(0, nrow = 2, ncol = 200)
centroids_mean[2,] <- colMeans(data[1:(n-c),])
centroids_mean[1,] <- colMeans(data[(n-c+1):n,])

post_good_2 <- posterior_pitman_yor(sigma = sigma, theta = theta, label = c_lab_good_2, data = data)

post_good$posterior
post_good_2$posterior
post_good$posterior == post_good_2$posterior


