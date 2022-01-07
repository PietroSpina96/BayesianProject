#### Setup ####

library(fda.usc)
library(fda)
library(fields)
library(roahd)
# library(rstan)

# Commentate i set delle directory non vostre
#setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject")
setwd("C:/Users/admin/Documents/R/Project_BS/BayesianProject") #GiuliaR

load("Simulated_WP.RData")
# load('Functions_WP.RData') 

###################### DATA SIMULATION - MODEL 1 ###############################

n <- 100
c <- 20
t_points <- 200

data <- matrix(0, nrow = n, ncol = t_points)

cov_M1 <- function(s,t){
  K <- 0.3 * exp(-abs(s - t)/0.3) 
  return(K)
}

# Covariance Matrix
time <- seq(0,1,1/(t_points-1))
K_1 <- matrix(0, nrow = t_points, ncol = t_points)
for (i in 1:t_points){
  for (j in 1:t_points){
    K_1[i,j] <- cov_M1(time[i],time[j])
  }
}

x11()
image.plot(time,time,K_1,main='Cov matrix mod-1')

# Create simulated data
set.seed(1234)
m <- rep(0,t_points)
random_process_1 <- generate_gauss_fdata(n, m, Cov = K_1)

main_proc_1 <- function(t,points){
  X <- rep(0,points)
  X[t] <- 30 * t * (1-t)^(3/2)
}

cont_proc_1 <- function(t,points){
  X <- rep(0,points)
  X[t] <- 30 * t^(3/2) * (1-t)
}

for (i in 1:(n-c)){
  for (j in 1:t_points){
    data[i,j] <- main_proc_1(time[j],t_points)
  }
}

for (i in (n-c+1):n){
  for (j in 1:t_points){
    data[i,j] <- cont_proc_1(time[j],t_points)
  }
}

data1 <- matrix(0, nrow = n, ncol = t_points)
for (i in 1:n){
  for (j in 1:t_points){
    data1[i,j] <- data[i,j] + random_process_1[i,j]
  }
}


# Simulated data plot
x11()
plot(time,data1[1,],type = 'l', ylim = c(-2,7.5), col = 'firebrick2', lwd = 2)
for(i in 2:(n-c)){
  lines(time,data1[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in (n-c+1):n){
  lines(time,data1[i,],type = 'l', col = 'blue', lwd = 2)
}
title('Simulated data - model 1')


# Smoothed data
eig_1 <- eigen(K_1)
values_1 <- eig_1$values
vectors_1 <- eig_1$vectors


alpha <- 0.1
f.data_alpha_sim_1 <- matrix(0, nrow = n, ncol = t_points)
for (i in 1:n){
  f.data_alpha_sim_1[i,] <- f_alpha_approx(data1[i,],alpha,values_1,vectors_1)
}

x11()
plot(time, data1[1,], type = 'l', lwd = 2,xlab='Time',ylab='Values',main='Comparison of alpha for mod-1: best alpha=0.1')
lines(time, f_alpha_approx(data1[1,],1,values_1,vectors_1), type = 'l', lwd = 2, col = 'firebrick2')
lines(time, f_alpha_approx(data1[1,],0.1,values_1,vectors_1), type = 'l', lwd = 2, col = 'blue')
lines(time, f_alpha_approx(data1[1,],0.01,values_1,vectors_1), type = 'l', lwd = 2, col = 'forestgreen')
legend("topright",ncol=1,box.lwd=1,legend=c('alpha = 1','alpha = 0.1','alpha = 0.01'),fill=c('firebrick2','blue','forestgreen'),x.intersp=0.3,
       text.col=c('firebrick2','blue','forestgreen'))


# alpha-Mahalanobis distance matrix 
Mahalanobis_Distance_1 <- matrix(0, nrow = n, ncol = n)
for (i in 1:n){
  for (j in 1:n){
    Mahalanobis_Distance_1[i,j] <- alpha_Mahalanobis(alpha,data1[i,],data1[j,],values_1,vectors_1)
  }
}

x11()
image.plot(1:n,1:n,Mahalanobis_Distance_1)

x11()
plot(time,f.data_alpha_sim_1[1,],type = 'l', ylim = c(-2,7.5), col = 'firebrick2', lwd = 2)
for(i in 2:(n-c)){
  lines(time,f.data_alpha_sim_1[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in (n-c+1):n){
  lines(time,f.data_alpha_sim_1[i,],type = 'l', col = 'blue', lwd = 2)
}

time123 <- time

#################### DATA SIMULATION - MODEL 2 #################################

n <- 100
c <- 20
t_points <- 200

data <- matrix(0, nrow = n, ncol = t_points)

cov_M2 <- function(s,t){
  K <- exp(-abs(s - t)) 
  return(K)
}

# Covariance Matrix
time <- seq(0,1,1/(t_points-1))
K_2 <- matrix(0, nrow = t_points, ncol = t_points)
for (i in 1:t_points){
  for (j in 1:t_points){
    K_2[i,j] <- cov_M2(time[i],time[j])
  }
}

# Create simulated data
set.seed(154)
m <- rep(0,t_points)   # mean
random_process_2 <- generate_gauss_fdata(n, m, Cov = K_2)
mu_2 <- runif (1, min = 0.25, max = 0.75)
u_2 <- rbinom(1, 1, prob = 0.5)

main_proc_2 <- function(t,points){
  X <- rep(0,points)
  X[t] <- 4*t
}

cont_proc_2 <- function(t,mu,u,points){
  X <- rep(0,points)
  X[t] <- 4*t + 1.8*(-1)^(u) + ((0.02*pi)^(-0.5))*exp(-((t-mu)^2)/0.02)
}

for (i in 1:(n-c)){
  for (j in 1:t_points){
    data[i,j] <- main_proc_2(time[j],t_points)
  }
}

for (i in (n-c+1):n){
  for (j in 1:t_points){
    data[i,j] <- cont_proc_2(time[j],mu = mu_2,u = u_2,t_points)
  }
}

data2 <- matrix(0, nrow = n, ncol = t_points)
for (i in 1:n){
  for (j in 1:t_points){
    data2[i,j] <- data[i,j] + random_process_2[i,j]
  }
}

# Simulated data plot
x11()
plot(time,data2[1,],type = 'l', ylim = c(-3.5,9), col = 'firebrick2', lwd = 2)
for(i in 2:(n-c)){
  lines(time,data2[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in (n-c+1):n){
  lines(time,data2[i,],type = 'l', col = 'blue', lwd = 2)
}
title('Simulated data - model 2')


# Smoothed data
eig_2 <- eigen(K_2)
values_2 <- eig_2$values
vectors_2 <- eig_2$vectors

alpha <- 0.1
f.data_alpha_sim_2<- matrix(0, nrow = n, ncol = t_points)
for (i in 1:n){
  f.data_alpha_sim_2[i,] <- f_alpha_approx(data2[i,],alpha,values_2,vectors_2)
}

x11()
plot(time, data2[1,], type = 'l', lwd = 2)
lines(time, f_alpha_approx(data2[1,],1,values_2,vectors_2), type = 'l', lwd = 2, col = 'firebrick2')
lines(time, f_alpha_approx(data2[1,],0.1,values_2,vectors_2), type = 'l', lwd = 2, col = 'blue')
lines(time, f_alpha_approx(data2[1,],0.01,values_2,vectors_2), type = 'l', lwd = 2, col = 'forestgreen')
title ('Curves comparison for alpha: best alpha=0.1')

# alpha-Mahalanobis distance matrix 
Mahalanobis_Distance_2 <- matrix(0, nrow = n, ncol = n)
for (i in 1:n){
  for (j in 1:n){
    Mahalanobis_Distance_2[i,j] <- alpha_Mahalanobis(alpha,data2[i,],data2[j,],values_2,vectors_2)
  }
}

x11()
image.plot(1:n,1:n,Mahalanobis_Distance_2)

x11()
plot(time,f.data_alpha_sim_2[1,],type = 'l', ylim = c(-3.5,9), col = 'firebrick2', lwd = 2)
for(i in 2:(n-c)){
  lines(time,f.data_alpha_sim_2[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in (n-c+1):n){
  lines(time,f.data_alpha_sim_2[i,],type = 'l', col = 'blue', lwd = 2)
}
title('Smoothed simulated data - model 2')


#################### DATA SIMULATION - MODEL 3 #################################

n <- 100
c <- 20
t_points <- 200

data <- matrix(0, nrow = n, ncol = t_points)

cov_M3 <- function(s,t){
  K <- exp(-abs(s - t)) 
  return(K)
}

# Covariance Matrix
time <- seq(0,1,1/(t_points-1))
K_3 <- matrix(0, nrow = t_points, ncol = t_points)
for (i in 1:t_points){
  for (j in 1:t_points){
    K_3[i,j] <- cov_M3(time[i],time[j])
  }
}

# Create simulated data
set.seed(123)
m <- rep(0,t_points)   # mean
random_process_3 <- generate_gauss_fdata(n, m, Cov = K_3)
mu_3 <- runif (1, min = 0.25, max = 0.75)

main_proc_3 <- function(t,points){
  X <- rep(0,points)
  X[t] <- 4*t
}

cont_proc_3 <- function(t,mu,points){
  X <- rep(0,points)
  X[t] <- 4*t + 2*sin(4*(t + mu)*pi)
}

for (i in 1:(n-c)){
  for (j in 1:t_points){
    data[i,j] <- main_proc_3(time[j],t_points)
  }
}

for (i in (n-c+1):n){
  for (j in 1:t_points){
    data[i,j] <- cont_proc_3(time[j],mu = mu_3,t_points)
  }
}

data3 <- matrix(0, nrow = n, ncol = t_points)
for (i in 1:n){
  for (j in 1:t_points){
    data3[i,j] <- data[i,j] + random_process_3[i,j]
  }
}

# Simulated data plot
x11()
plot(time,data3[1,],type = 'l', ylim = c(-3,7.5), col = 'firebrick2', lwd = 2)
for(i in 2:(n-c)){
  lines(time,data3[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in (n-c+1):n){
  lines(time,data3[i,],type = 'l', col = 'blue', lwd = 2)
}
title('Simulated data - model 3')


# Smoothed data
eig_3 <- eigen(K_3)
values_3 <- eig_3$values
vectors_3 <- eig_3$vectors

alpha <- 0.1
f.data_alpha_sim_3<- matrix(0, nrow = n, ncol = t_points)
for (i in 1:n){
  f.data_alpha_sim_3[i,] <- f_alpha_approx(data3[i,],alpha,values_3,vectors_3)
}

x11()
plot(time, data3[1,], type = 'l', lwd = 2)
lines(time, f_alpha_approx(data3[1,],1,values_3,vectors_3), type = 'l', lwd = 2, col = 'firebrick2')
lines(time, f_alpha_approx(data3[1,],0.1,values_3,vectors_3), type = 'l', lwd = 2, col = 'blue')
lines(time, f_alpha_approx(data3[1,],0.01,values_3,vectors_3), type = 'l', lwd = 2, col = 'forestgreen')
title ('Curves comparison for alpha: best alpha=0.1')

# alpha-Mahalanobis distance matrix 
Mahalanobis_Distance_3 <- matrix(0, nrow = n, ncol = n)
for (i in 1:n){
  for (j in 1:n){
    Mahalanobis_Distance_3[i,j] <- alpha_Mahalanobis(alpha,data3[i,],data3[j,],values_3,vectors_3)
  }
}

x11()
image.plot(1:n,1:n,Mahalanobis_Distance_3)

x11()
plot(time,f.data_alpha_sim_3[1,],type = 'l', ylim = c(-3,7.5), col = 'firebrick2', lwd = 2)
for(i in 2:(n-c)){
  lines(time,f.data_alpha_sim_3[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in (n-c+1):n){
  lines(time,f.data_alpha_sim_3[i,],type = 'l', col = 'blue', lwd = 2)
}
title('Smoothed simulated data - model 3')

###################### DATA SIMULATION - MODEL 4 ###############################
n <- 100
n1 <- 40

t_points <- 200
time <- seq(0,1,(1)/(t_points - 1))

cov_M4_1 <- function(s,t){
  K <- 0.3 * exp(-abs(s - t)/0.3) 
  return(K)
}

cov_M4_2 <- function(s,t){
  K <- 1.5*exp(-abs(s - t)/3) 
  return(K)
}

K_4_1 <- matrix(0, nrow = t_points, ncol = t_points)
for (i in 1:t_points){
  for (j in 1:t_points){
    K_4_1[i,j] <- cov_M4_1(time[i],time[j])
  }
}

K_4_2 <- matrix(0, nrow = t_points, ncol = t_points)
for (i in 1:t_points){
  for (j in 1:t_points){
    K_4_2[i,j] <- cov_M4_2(time[i],time[j])
  }
}

x11()
par(mfrow=c(1,2))
image.plot(time,time,K_4_1,main='Cov matrix mod-4 cluster-1')
image.plot(time,time,K_4_2,main='Cov matrix mod-4 cluster-2')

set.seed(123564)
m <- rep(0,t_points)
random_process_cluster_1 <- generate_gauss_fdata(n1, m, Cov = K_4_1)
random_process_cluster_2 <- generate_gauss_fdata(n - n1, m, Cov = K_4_2)
random_process_cluster <- rbind(random_process_cluster_1,random_process_cluster_2)

cluster_1 <- function(t,points){
  X <- rep(0,points)
  X[t] <- sin(t)
}

cluster_2 <- function(t,points){
  X <- rep(0,points)
  X[t] <- 30 * t * (1-t)^(3/2)
}

data4 <- matrix(0, nrow = n, ncol = t_points)
for (i in 1:n1){
  for (j in 1:t_points){
    data4[i,j] <- cluster_1(time[j],t_points) + random_process_cluster[i,j]
  }
}

for (i in (n1 + 1):n){
  for (j in 1:t_points){
    data4[i,j] <- cluster_2(time[j],t_points) + random_process_cluster[i,j]
  }
}


# Simulated data plot
x11()
plot(time,data4[1,],type = 'l', ylim = c(-5,10), col = 'blue', lwd = 2,xlab='Time',ylab='Values')
for(i in 2:n1){
  lines(time,data4[i,],type = 'l', col = 'blue',lwd = 2)
}
for (i in (n1 + 1):n){
  lines(time,data4[i,],type = 'l', col = 'firebrick2', lwd = 2)
}
title('Simulated data - model 4')

# Covariance matrix of the data
cov_4 <- cov(data4)

x11()
image.plot(time,time,cov_4,main='Covariance matrix')

# Smoothed data for the second cluster
eig_4_2 <- eigen(K_4_2)
values_4_2 <- eig_4_2$values
vectors_4_2 <- eig_4_2$vectors

alpha <- 0.1
x11()
plot(time, data4[n1+1,], type = 'l', lwd = 2)
lines(time, f_alpha_approx(data4[n1+1,],1,values_4_2,vectors_4_2), type = 'l', lwd = 2, col = 'firebrick2')
lines(time, f_alpha_approx(data4[n1+1,],0.1,values_4_2,vectors_4_2), type = 'l', lwd = 2, col = 'blue')
lines(time, f_alpha_approx(data4[n1+1,],0.01,values_4_2,vectors_4_2), type = 'l', lwd = 2, col = 'forestgreen')
title('Best value for alpha is 0.1 for the second cluster')

# Smoothed data for the first cluster
eig_4_1 <- eigen(K_4_1)
values_4_1 <- eig_4_1$values
vectors_4_1 <- eig_4_1$vectors

x11()
plot(time, data4[1,], type = 'l', lwd = 2)
lines(time, f_alpha_approx(data4[1,],1,values_4_1,vectors_4_1), type = 'l', lwd = 2, col = 'firebrick2')
lines(time, f_alpha_approx(data4[1,],0.1,values_4_1,vectors_4_1), type = 'l', lwd = 2, col = 'blue')
lines(time, f_alpha_approx(data4[1,],0.01,values_4_1,vectors_4_1), type = 'l', lwd = 2, col = 'forestgreen')
title('Best value for alpha is 0.1 for the first cluster')


# alpha-Mahalanobis distance matrix 
Mahalanobis_Distance_4 <- matrix(0, nrow = n, ncol = n)
for (i in 1:n1){
  for (j in 1:n1){
    Mahalanobis_Distance_4[i,j] <- alpha_Mahalanobis(alpha,data4[i,],data4[j,],values_4_1,vectors_4_1)
  }
}
for (i in (n1+1):n){
  for (j in (n1+1):n){
    Mahalanobis_Distance_4[i,j] <- alpha_Mahalanobis(alpha,data4[i,],data4[j,],values_4_2,vectors_4_2)
  }
}

x11()
image.plot(1:n,1:n,Mahalanobis_Distance_4)

time4 <- time

f.data_alpha_sim_4 <- matrix(0, nrow = n, ncol = t_points)
for (i in 1:n1){
  f.data_alpha_sim_4[i,] <- f_alpha_approx(data4[i,],alpha,values_4_1,vectors_4_1)
}
for (i in (n1+1):n){
  f.data_alpha_sim_4[i,] <- f_alpha_approx(data4[i,],alpha,values_4_2,vectors_4_2)
}
x11()
plot(time,f.data_alpha_sim_4[1,],type = 'l', ylim = range(data4) ,col = 'firebrick2', lwd = 2)
for(i in 2:n1){
  lines(time,f.data_alpha_sim_4[i,],type = 'l', col = 'firebrick2',lwd = 2)
}
for (i in (n1+1):n){
  lines(time,f.data_alpha_sim_4[i,],type = 'l', col = 'blue', lwd = 2)
}
title('Smoothed data')




##### FINAL WORKSPACE #### 
# Remove useless variables and functions
rm(list=c('data','random_process_1','random_process_2','random_process_3'))
rm(list=c('random_process_cluster_1','random_process_cluster_2','random_process_cluster'))
rm(list=c('mu_2','mu_3','u_2','i','j','m'))
rm(list = c('n','t_points','time','delta'))
rm(list=c('values_1','values_2','values_3','vectors_1','vectors_2','vectors_3'))
rm(list=c('values_4_1','values_4_2','vectors_4_1','vectors_4_2'))

rm(list=c('cont_proc_1','cont_proc_2','cont_proc_3'))
rm(list=c('main_proc_1','main_proc_2','main_proc_3'))
rm(list=c('cov_M1','cov_M2','cov_M3','cov_M4_1','cov_M4_2'))
rm(list=c('cluster_1','cluster_2'))

rm(list=c('clusters_plot','clusters_union','fda_clustering_mahalanobis'))
rm(list=c('fda_clustering_mahalanobis_union','fda_clustering_mahalanobis_updated','gibbs_loss','gibbs_loss_k'))

# Save Workspace
#setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject")
#save.image("Simulated_WP.RData")

save.image("~/R/Project_BS/BayesianProject/Simulated_WP.RData") #GiuliaR
