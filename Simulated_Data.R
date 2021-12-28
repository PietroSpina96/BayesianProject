#### Setup ####

library(fda.usc)
library(fda)
library(fields)
library(roahd)

# Commentate i set delle directory non vostre
#setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject")
setwd("C:/Users/admin/Documents/R/Project_BS/BayesianProject") #GiuliaR

load("Simulated_WP.RData")


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


#### alpha-Mahalanobis distance calculation ####
# functions

#alpha_approximation
f_alpha_approx <-function(f,alpha,lambda,eigenft){
  t_points <- length(f)
  
  coeff<-prod<-res<-rep(0,t_points)
  approx<-matrix(0,t_points,t_points)
  
  for (j in 1:t_points) {
    
    coeff[j]<-lambda[j]/(lambda[j]+alpha)
    prod[j]<-scalar_prod(f,eigenft[,j])
    approx[,j]<- as.numeric(coeff[j]*prod[j])*eigenft[,j]
    res<-res+approx[,j]
    
  }
  return(res)
}

#alpha-mahalanobis distance
alpha_Mahalanobis <- function(alpha,f1,f2,lambda, eigenft) {
  t_points <- length(f1)
  dis<-coeff<-prod<-rep(0,t_points)
  
  for (j in 1:t_points){
    
    coeff[j]<-lambda[j]/(lambda[j]+alpha)^2
    prod[j]<-(scalar_prod(f1-f2,eigenft[,j]))^2
    dis[j]<-coeff[j]*prod[j]
    
  }
  res<-sum(dis)
  return(res)
}

#scalar product
scalar_prod<- function (f1,f2) {
  # f sono vettori colonna
  res<- t(f1)%*%f2
  return(res)
}


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
plot(time, data1[1,], type = 'l', lwd = 2)
lines(time, f_alpha_approx(data1[1,],1,values_1,vectors_1), type = 'l', lwd = 2, col = 'firebrick2')
lines(time, f_alpha_approx(data1[1,],0.1,values_1,vectors_1), type = 'l', lwd = 2, col = 'blue')
lines(time, f_alpha_approx(data1[1,],0.01,values_1,vectors_1), type = 'l', lwd = 2, col = 'forestgreen')

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

save.image("~/R/Project_BS/BayesianProject/Simulated_WP.RData") #GiuliaR


#################### DATA SIMULATION - MODEL 4 #################################
# The following model is taken from the fda library 

data4 <- t(CanadianWeather$dailyAv[,,1])
n <- dim(data4)[1]
t_points <- 365
time<- 1:365

# Plot of the data
x11()
matplot(t(data4),type='l',xlab='Day',ylab='Temperature')

x11()
plot(time,data4[1,],ylim=c(-40,30))
for (i in 2:dim(data4)[1])
  lines(time,data4[i,])

# Covariance matrix and plot
K_4 <- cov(data4)

x11()
image.plot(time,time,K_4,main='Covariance matrix')

# Eigenvalues and eigenfunctions of the covariance matrix
eig_4 <- eigen(K_4)
values_4 <- eig_4$values
vectors_4 <- eig_4$vectors

x11()
plot(time, data4[1,], type = 'l', lwd = 2)
lines(time, f_alpha_approx(data4[1,],1,values_4,vectors_4), type = 'l', lwd = 2, col = 'firebrick2')
lines(time, f_alpha_approx(data4[1,],10,values_4,vectors_4), type = 'l', lwd = 2, col = 'blue')
lines(time, f_alpha_approx(data4[1,],100,values_4,vectors_4), type = 'l', lwd = 2, col = 'forestgreen')
title ('Curves comparison for alpha: best alpha=10')


alpha <- 10
f.data_alpha_sim_4<- matrix(0, nrow = n, ncol = t_points)
for (i in 1:n){
  f.data_alpha_sim_4[i,] <- f_alpha_approx(data4[i,],alpha,values_4,vectors_4)
}

# alpha-Mahalanobis distance matrix 
Mahalanobis_Distance_4 <- matrix(0, nrow = n, ncol = n)
for (i in 1:n){
  for (j in 1:n){
    Mahalanobis_Distance_4[i,j] <- alpha_Mahalanobis(alpha,data4[i,],data4[j,],values_4,vectors_4)
  }
}

x11()
image.plot(1:n,1:n,Mahalanobis_Distance_4)

x11()
plot(time,f.data_alpha_sim_4[1,],type = 'l', ylim = c(-40,30), col = 'firebrick2', lwd = 2)
for(i in 2:n)
  lines(time,f.data_alpha_sim_4[i,],type = 'l', col = 'firebrick2',lwd = 2)
title('Smoothed simulated data - model 4')



################ DATA SIMULATION - MODEL 5 #####################################
# The following model is taken from Secchi's TDE (19.07.2019). There are three clusters

data5 <- read.table('data_model5.txt',header=TRUE)
data5 <- data5[,1:365]
data5 <- as.matrix(data5)
t_points <- 365
time <- 1:365
n <- dim(data5)[1]

x11()
matplot(t(data5),type='l',main='Data5',xlab='time',ylab='Values',ylim=range(data5))

x11()
plot(time,data5[1,],ylim=c(9,28),type='l')
for (i in 2:n)
  lines(time,data5[i,])

# Covariance matrix and plot
K_5 <- cov(data5)

x11()
image.plot(time,time,K_5,main='Covariance matrix')

# Eigenvalues and eigenfunctions of the covariance matrix
eig_5 <- eigen(K_5)
values_5 <- eig_5$values
vectors_5 <- eig_5$vectors

x11()
plot(time, data5[1,], type = 'l', lwd = 2)
lines(time, f_alpha_approx(data5[1,],1,values_5,vectors_5), type = 'l', lwd = 2, col = 'firebrick2')
lines(time, f_alpha_approx(data5[1,],0.1,values_5,vectors_5), type = 'l', lwd = 2, col = 'blue')
lines(time, f_alpha_approx(data5[1,],0.01,values_5,vectors_5), type = 'l', lwd = 2, col = 'forestgreen')
title ('Curves comparison for alpha: best alpha=10')


alpha <- 10
f.data_alpha_sim_5<- matrix(0, nrow = n, ncol = t_points)
for (i in 1:n){
  f.data_alpha_sim_5[i,] <- f_alpha_approx(data5[i,],alpha,values_5,vectors_5)
}

# alpha-Mahalanobis distance matrix 
Mahalanobis_Distance_5 <- matrix(0, nrow = n, ncol = n)
for (i in 1:n){
  for (j in 1:n){
    Mahalanobis_Distance_5[i,j] <- alpha_Mahalanobis(alpha,data5[i,],data5[j,],values_5,vectors_5)
  }
}

x11()
image.plot(1:n,1:n,Mahalanobis_Distance_5)

x11()
plot(time,f.data_alpha_sim_5[1,],type = 'l', ylim = c(-40,30), col = 'firebrick2', lwd = 2)
for(i in 2:n)
  lines(time,f.data_alpha_sim_5[i,],type = 'l', col = 'firebrick2',lwd = 2)
title('Smoothed simulated data - model 5')


################ DATA SIMULATION - MODEL 6 #####################################
library(fdakma)
library(fda)
library(kma)

data(kma.data)
data6 <- kma.data$y0
# y1_6 <- kma.data$y1

time <- as.vector(kma.data$x)
t_points <- 200
n <- dim(data6)[1]

x11()
matplot(t(data6),type='l')

K_6 <- cov(data6)

x11()
image.plot(time,time,K_6,main='Covariance matrix')

# Eigenvalues and eigenfunctions of the covariance matrix
eig_6 <- eigen(K_6)
values_6 <- eig_6$values
vectors_6 <- eig_6$vectors

alpha <- 0
f.data_alpha_sim_6<- matrix(0, nrow = n, ncol = t_points)
for (i in 1:n){
  f.data_alpha_sim_6[i,] <- f_alpha_approx(data6[i,],alpha,values_6,vectors_6)
}

x11()
plot(time, data6[1,], type = 'l', lwd = 2)
lines(time, f_alpha_approx(data6[1,],0,values_6,vectors_6), type = 'l', lwd = 2, col = 'firebrick2')
lines(time, f_alpha_approx(data6[1,],0.1,values_6,vectors_6), type = 'l', lwd = 2, col = 'blue')
lines(time, f_alpha_approx(data6[1,],0.01,values_6,vectors_6), type = 'l', lwd = 2, col = 'forestgreen')
title ('Curves comparison for alpha: best alpha=0')

# alpha-Mahalanobis distance matrix 
Mahalanobis_Distance_6 <- matrix(0, nrow = n, ncol = n)
for (i in 1:n){
  for (j in 1:n){
    Mahalanobis_Distance_6[i,j] <- alpha_Mahalanobis(alpha,data6[i,],data6[j,],values_6,vectors_6)
  }
}

x11()
image.plot(1:n,1:n,Mahalanobis_Distance_6)



# Remove useless variables
rm(list=c('data','kma.data','random_process_1','random_process_2','random_process_3'))
rm(list=c('mu_2','mu_3','u_2'))

##### Save Workspace ####
#setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject")
#save.image("Simulated_WP.RData")

save.image("~/R/Project_BS/BayesianProject/Simulated_WP.RData") #GiuliaR
