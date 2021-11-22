library(fda.usc)

setwd("C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati")

#### Kernel estimator ####

n <- 26
X_bar <- colMeans(f.data$ausxSL$data) # functional mean

# Plot of functional mean vs data
x11()
plot(t(importMatrix(res, type = c('SL', 'sx'), position = 'au'))[1,], type = "l", ylim = c(-250, 250))
for(i in 1:26){
  lines(t(importMatrix(res, type = c('SL', 'sx'), position = 'au'))[i,], lwd=1)
}
lines(f.data$ausxSL$argvals, X_bar, type = 'l', lwd=3, col = 'firebrick2')


# Covariance function estimator K_hat
kernel_estimator <- function(process,med,n,s,t){
  est <- rep(0,1600)
  for (i in 1:n) {
    est[i] <- (process[i,s] - med[s])*(process[i,t] - med[t])
  }
  estimator = (1/n)*sum(est)
  
  return(estimator)
}

kernel_estimator(f.data$ausxSL$data,X_bar,n,1,9) # prova

K_hat <- matrix(0,1600,1600)
for (i in 1:1600){
  print(i)
  for (j in 1:1600){
    K_hat[i,j] <- kernel_estimator(f.data$ausxSL$data,X_bar,n,i,j)
  }
}

#### alpha-Mahalanobis Distance ####

alpha_Mahalanobis <- function(alpha,f1,f2,K) {
  
}