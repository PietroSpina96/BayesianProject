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


#### Posterior calculation on Model 1 with optimal c_lab and fixed structure ####
# For now this is just a test to see if the new function works

n <- dim(data)[1]
c_lab <- c(rep(1,n-c), rep(2,c))
n_clust <- k <- 2
alpha <- 0.1
lambda <- 0.25

# Initialize cetroid matrix
centroids_mean <- matrix(0, nrow = 2, ncol = 200)
centroids_mean[1,] <- colMeans(data[1:(n-c),])
centroids_mean[2,] <- colMeans(data[(n-c+1):n,])

# Itialize eigen object
eig <- eigen(K_1)

# Values initialized by considering example at page 29 of the paper
sigma <- 0.25
theta <- 3.66

posterior_pitman_yor(sigma = sigma, theta = theta, label = c_lab, data = data)

##### Checking function stuff #####
label <- c_lab
clust_obs <- rep(0,n_clust)
for (k in 1:n_clust){
  nk <- dim(data[which(label == k),])[1]
  clust_obs[k] <- nk
}

post_vec <- rep(0,n_clust)
for (k in 1:n_clust){
  post_vec[k] <- prod(theta + k*sigma, gamma(clust_obs[k] - sigma)) 
  # Con lambda = 1 questo l'ultimo valore è piccolissimo
}
post_vec

post <- prod(post_vec,1/(theta + n_clust*sigma), exp(-lambda*gibbs_loss(n_clust = n_clust , 
                                                                        centroids = centroids_mean, 
                                                                        label = label, eig = eig,data = data)))
post
