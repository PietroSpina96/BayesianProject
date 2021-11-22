#### Setup ####

library(fda.usc)
library(fda)
library(fields)

# Commentate i set delle directory non vostre
setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati")
setwd('C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/potenziali_evocati')
#setwd("C:/Users/admin/Documents/R/Project_BS/potenziali_evocati")

load('functional_WP.RData')

#### Kernel estimator ####

n <- 26
t_points <- 1600
X_bar <- colMeans(f.data$ausxSL$data) # functional mean

# Plot of functional mean vs data
x11()
plot(t(importMatrix(res, type = c('SL', 'sx'), position = 'au'))[1,], type = "l", ylim = c(-250, 250))
for(i in 1:26){
  lines(t(importMatrix(res, type = c('SL', 'sx'), position = 'au'))[i,], lwd=1)
}
lines(f.data$ausxSL$argvals, X_bar, type = 'l', lwd=3, col = 'firebrick2')


# Covariance function estimator K_hat
kernel_estimator <- function(process,med,n,s,t,points){
  est <- rep(0,points)
  for (i in 1:n) {
    est[i] <- (process[i,s] - med[s])*(process[i,t] - med[t])
  }
  estimator = (1/n)*sum(est)
  
  return(estimator)
}

kernel_estimator(f.data$ausxSL$data,X_bar,n,1,9,t_points) # prova

K_hat <- matrix(0,t_points,t_points)
for (i in 1:t_points){
  # print(i)
  for (j in 1:t_points){
    K_hat[i,j] <- kernel_estimator(f.data$ausxSL$data,X_bar,n,i,j,t_points)
  }
}

# Plot of the covariance matrix
time<-seq(1,t_points)
x11()
image.plot(time,time,K_hat)

# sub-matrix (5x5)
khat_example<-K_hat[1:5,1:5]
show(khat_example)

#### alpha-Mahalanobis distance and other usefull functions ####
# functions
scalar_prod<- function (f1,f2) {
  # f sono vettori colonna
  res<- t(f1)%*%f2
  return(res)
}


# equivalentemennte
#scalar_prod_norm<- function (f1,f2) {
# f sono vettori colonna
#f1<-as.matrix(f1)
#f2<-as.matrix(f2)
#res<- (norm(f1+f2,'2')^2 - norm(f1-f2,'2')^2)/4
#return(res)
#}


# Compute eigenvalues and eigenfunctions of the covariance matrix
eigenvf<-eigen(K_hat)
lambda<-eigenvf$values
eigenft<-eigenvf$vectors


# OBS: from the theory we know that only (N-1) eigenvalues are not null where N=number of 
# statistical units. Here we have N=26 and 25 not null eigenvalues.
lambda[1:26]

alpha_Mahalanobis <- function(alpha,f1,f2,lambda, eigenft) {
  #calcolo la distanza al quadrato
  
  dis<-coeff<-prod<-rep(0,t_points)
  for (j in 1:t_points){
    
    coeff[j]<-lambda[j]/(lambda[j]+alpha)^2
    
    prod[j]<-(scalar_prod(f1-f2,eigenft[,j]))^2
    
    dis[j]<-coeff[j]*prod[j]
  }
  
  res<-sum(dis)
  return(res)
}

# approximation of the function f with f_alpha
f_alpha_approx <-function(f,alpha,lambda,eigenft){
  coeff<-prod<-rep(0,t_points)
  approx<-res<-matrix(0,t_points,t_points)
  
  for (j in 1:t_points) {
    
    coeff[j]<-lambda[j]/(lambda[j]+alpha)
    prod[j]<-scalar_prod(f,eigenft[,j])
    
    approx[,j]<- (coeff[j]*prod[j])*eigenft[,j]
    
    res[,j]<-res[,j] + approx[,j]
    
  }
  return(res[,t_points])
}

# norm and inner product wrt sample covariance function K
norm2_K <- function (f,lambda,eigenft){
  norm_vect<-prod<-rep(0,t_points)
  
  for (j in 1:t_points){
    prod[j]<-(scalar_prod(f,eigenft[,j]))^2
    norm_vect[j]<-prod[j]/lambda[j]
  }
  res<-sum(norm_vect)
  return(res)
}

inner_product_K<- function(f,g,lambda,eigenft) {
  prod_vect<-prod_f<-prod_g<-rep(0,t_points)
  
  for (j in 1:t_points){
    prod_f[j]<-(scalar_prod(f,eigenft[,j]))
    prod_g[j]<-(scalar_prod(g,eigenft[,j]))
    
    norm_vect[j]<-prod_f[j]*prod_g[j]/lambda[j]
  }
  res<-sum(norm_vect)
  return(res)
}



# prova dell'approssimazione con una funzione f presa dal dataset.
# problem: the eigenfunctions have very small values and this leads the approximation of f to be very small and far from the true one. 
x11()
plot(t(importMatrix(res, type = c('SL', 'sx'), position = 'au'))[1,], type = "l", ylim = c(-250, 250))
lines(f.data$ausxSL$argvals, X_bar, type = 'l', lwd=3, col = 'firebrick2')
lines(f.data$ausxSL$argvals, f_prova_alpha, type = 'l', lwd=3, col = 'blue')

f_prova<-X_bar
f_prova_alpha<-f_alpha_approx(f_prova,alpha=1e-8,lambda,eigenft)


# Example of sum for the first two eigenfunctions
coeff_1<-lambda[1]/(lambda[1]+ 1e-8)
prod_1<-scalar_prod(f_prova,eigenft[,1])
approx_1<- (coeff_1*prod_1)*eigenft[,1]

coeff_2<-lambda[2]/(lambda[2]+ 1e-8)
prod_2<-scalar_prod(f_prova,eigenft[,2])
approx_2<- (coeff_2*prod_2)*eigenft[,2]

res<-approx_2 + approx_1


##### Save Workspace ####
setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/BayesianProject")
save.image("functional_WP.RData")
