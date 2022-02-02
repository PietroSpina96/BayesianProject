
discrepancy_within <- function(x, centroid, alpha, eig){
  #compute square discrepancy
  res <- 0
  n <- dim(x)[1]
  for (i in 1:n){
    sum_partial = alpha_Mahalanobis(alpha,x[i,],centroid,eig$values,eig$vectors)
    res = res + sum_partial
  }
  return(res)
}


prob_i_k <- function(obs_index, clust_index, Ci, data, lambda, alpha, eig, verbose=FALSE){
  ## compute fullconditional term for uniform prior
  
  idx_k = which(Ci == clust_index) #indexes of elements of cluster k
  if(!(obs_index %in% idx_k)){#check if our 'i' is already in idx_k, if not, add it
    idx_k = c(obs_index,idx_k)
    obs_idx = 1
  } else {
    obs_idx = match(obs_index, idx_k)#if it is in clust k, get position of the first(and only) match 
  }
  
  x_k = data[idx_k,]
  centroid = colMeans(x_k)
  
  
  start_time <- Sys.time()
  disc_in = discrepancy_within(x_k, centroid, alpha, eig)
  end_time <- Sys.time()
  
  if(verbose){
    print(paste("Computing full discrepancy in ",end_time - start_time,"s"))
  }
  
  
  x_k_no_i = x_k[-obs_idx,]
  centroid_no_i = colMeans(x_k_no_i)
  
  start_time <- Sys.time()
  disc_no_i <- discrepancy_within(x_k_no_i, centroid_no_i, alpha, eig)
  end_time <- Sys.time()
  
  if(verbose){
    print(paste("Computing discrepancy without obs ",obs_index," in ",end_time - start_time,"s"))
  }
  
  
  return(exp(-lambda*(disc_in - disc_no_i)))
}


prob_i_k_PY <- function(obs_index, clust_index, Ci, data, sig, lambda, alpha, eig, verbose=FALSE){
  ## compute fullconditional term for PY-EPPF prior
  
  
  idx_k = which(Ci == clust_index) #indexes of elements of cluster k
  if(!(obs_index %in% idx_k)){#check if our 'i' is already in idx_k, if not, add it
    idx_k = c(obs_index,idx_k)
    obs_idx = 1
  } else {
    obs_idx = match(obs_index, idx_k)#if it is in clust k, get position of the first(and only) match 
  }
  
  x_k = data[idx_k,]
  centroid = colMeans(x_k)
  
  
  start_time <- Sys.time()
  disc_in = discrepancy_within(x_k, centroid, alpha, eig)
  end_time <- Sys.time()
  
  if(verbose){
    print(paste("Computing full discrepancy in ",end_time - start_time,"s"))
  }
  
  
  x_k_no_i = x_k[-obs_idx,]
  centroid_no_i = colMeans(x_k_no_i)
  
  start_time <- Sys.time()
  disc_no_i <- discrepancy_within(x_k_no_i, centroid_no_i, alpha, eig)
  end_time <- Sys.time()
  
  if(verbose){
    print(paste("Computing discrepancy without obs ",obs_index," in ",end_time - start_time,"s"))
  }
  
  
  return((length(idx_k)-1-sig)*exp(-lambda*(disc_in - disc_no_i)))
}


update_k_i <- function(obs_idx, Ci, n_clust, data, sig, lambda, alpha, eig){
  # Sample a value for the cluster, compatible with the given expression of the full conditional
  probs = rep(0,n_clust)
  for(ii in 1:n_clust)
    probs[ii] = prob_i_k_PY(obs_idx, ii, Ci, data, sig, lambda, alpha, eig,verbose=F)
    
  k = sample(n_clust, 1, replace = TRUE, prob = probs)
  return(k)
}




#### MH for our distro

gibbs_sampler_PY <- function(N, N_burnin, x0, data, lambda, alpha, eig, sig, verbose = F){
  #compute MH samples of the posterior for cluter indexes Ci
  #INPUTS
  # N = number of samples
  # N_burnin = number of initial elements to discard
  # x0 = initial labels
  
  k = 2 #FIXED HERE FOR NOW
  n <- dim(data)[1]
  
  samplee = matrix(0, nrow = N, ncol = n)
  posz = 0
  iterMH = 1
  
  
  if(verbose){
    print("## Begin Gibbs Iterations ##")
  }
  
  X_temp = x0
  
  while(posz < N){
    
    gibbs_time_init <- Sys.time()
    
    # for(t in 1:n){ #update one component at the time
    #   #sample the value according to the full conditional
    #   init_time <- Sys.time()
    #   X_temp[t] = update_k_i(t, X_temp, k, data, sig, lambda, alpha, eig)
    #   finish_time <- Sys.time()
    #   print(paste("Sampled indicator for obs",t," in ",finish_time - init_time,"s"))
    # }
    
    for(t in 1:n)#update one component at the time
      X_temp[t] = update_k_i(t, X_temp, k, data, sig, lambda, alpha, eig)
    
    if(iterMH > N_burnin){
        posz = iterMH - N_burnin
        samplee[posz,] = X_temp
      }
    
    
    gibbs_time_finish <-Sys.time()
    
    if(verbose){
      print(paste("Gibbs iteration N = ",iterMH," in ",gibbs_time_finish-gibbs_time_init,"s"))
    }
    iterMH = iterMH + 1
  }
  return(samplee)
}


Ci <- gibbs_sampler_PY(11000, 1000, c_opt, data, .1, alpha, eigen(K_1), 0.25, verbose = T)

# likhood_i_k(79, 1, c_opt, data, .1, alpha, eig, F)


avg_ci <- colMeans(Ci)
colore = avg_ci + rep(2, dim(data)[1])

plot(c_opt - avg_ci)


x11()
#par(bg = 'gold')
plot(time,data[1,],type = 'l', ylim = c(-3.5,9), lwd = 2, main = "Main and contaminated processes")
for(i in 2:n){
  lines(time,data[i,],type = 'l', col = colore[i],lwd = 2)
}
