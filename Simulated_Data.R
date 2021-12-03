#### Setup ####

library(fda.usc)
library(fda)
library(fields)
library(roahd)

# Commentate i set delle directory non vostre
setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject")
#setwd("C:/Users/admin/Documents/R/Project_BS/potenziali_evocati")


#### DATA SIMULATION - MODEL 1 ####

n <- 100
c <- 15
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

# Create simulated data

m <- rep(0,t_points)
random_process <- generate_gauss_fdata(n, m, Cov = K_1)

main_proc <- function(t,points){
  X <- rep(0,points)
  X[t] <- 30 * t * (1-t)^(3/2)
}

cont_proc <- function(t,points){
  X <- rep(0,points)
  X[t] <- 30 * t^(3/2) * (1-t)
}

for (i in 1:(n-c)){
  for (j in 1:t_points){
    data[i,j] <- main_proc(time[j],t_points)
  }
}

for (i in (n-c+1):n){
  for (j in 1:t_points){
    data[i,j] <- cont_proc(time[j],t_points)
  }
}

data1 <- matrix(0, nrow = n, ncol = t_points)
for (i in 1:n){
  for (j in 1:t_points){
    data1[i,j] <- data[i,j] + random_process[i,j]
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


#### alpha-Mahalanobis distance calculation ####

# Smoothed data

eig <- eigen(K_1)
values <- eig$values
vectors <- eig$vectors

alpha <- 1e+4
f.data_alpha_sim <- matrix(0, nrow = n, ncol = t_points)
for (i in 1:n){
  f.data_alpha_sim[i,] <- f_alpha_approx(data1[i,],alpha,values,vectors)
}

x11()
plot(time, data1[1,], type = 'l', lwd = 2)
lines(time, f_alpha_approx(data1[1,],1e+4,values,vectors), type = 'l', lwd = 2, col = 'firebrick2')
lines(time, f_alpha_approx(data1[1,],1e+6,values,vectors), type = 'l', lwd = 2, col = 'blue')
lines(time, f_alpha_approx(data1[1,],1e+8,values,vectors), type = 'l', lwd = 2, 'forestgreen')

# alpha-Mahalanobis distance matrix 

Mahalanobis_Distance <- matrix(0, nrow = n, ncol = n)
for (i in 1:n){
  for (j in 1:n){
    Mahalanobis_Distance[i,j] <- alpha_Mahalanobis(alpha,data1[i,],data1[j,],values,vectors)
  }
}

x11()
image.plot(1:n,1:n,Mahalanobis_Distance)

