
discrepancy_within <- function(x, centroid, alpha, eig){
  #compute square discrepancy
  res <- 0
  n <- dim(x)[1]
  if(!is.null(n)){
    for (i in 1:n){
      sum_partial = alpha_Mahalanobis(alpha,x[i,],centroid,eig$values,eig$vectors)
      res = res + sum_partial
    }
  }else{
    res = alpha_Mahalanobis(alpha,x,centroid,eig$values,eig$vectors)
  }
  return(res)
}


prob_i_k <- function(obs_index, clust_index, Ci, data, lambda, alpha, verbose=FALSE){
  ## compute fullconditional term for PY-EPPF prior
  
  
  idx_k = which(Ci == clust_index) #indexes of elements of cluster k
  if(!(obs_index %in% idx_k)){#check if our 'i' is already in idx_k, if not, add it
    idx_k = c(obs_index,idx_k)
    obs_idx = 1
  } else {
    obs_idx = match(obs_index, idx_k)#if it is in clust k, get position of the first(and only) match 
  }
  
  x_k = data[idx_k,]
  
  #update covariance eigs
  eig = eigen(cov(x_k)) #not efficient here but the code is easier to write and more readable
  
  if(length(idx_k)>1) #i.e. at least 2 observations
    centroid = colMeans(x_k)
  else
    return(-1) # if it's the only observation it can't be reallocated
  
  
  start_time <- Sys.time()
  disc_in = discrepancy_within(x_k, centroid, alpha, eig)
  end_time <- Sys.time()
  
  if(verbose){
    print(paste("Computing full discrepancy in ",end_time - start_time,"s"))
  }
  
  x_k_no_i = x_k[-obs_idx,]
  
  if(length(idx_k)==2){ # i.e. exactly 2 obs in the cluster
    disc_no_i <- 0
    if(verbose){
      print("Since length(idx_k)==2 computation of discrepancy without obs is not necessary! ")
    }
  }else{
    centroid_no_i = colMeans(x_k_no_i)
    start_time <- Sys.time()
    disc_no_i <- discrepancy_within(x_k_no_i, centroid_no_i, alpha, eig)
    end_time <- Sys.time()
    
    if(verbose){
      print(paste("Computing discrepancy without obs ",obs_index," in ",end_time - start_time,"s"))
    }
  }
  return(exp(-lambda*(disc_in - disc_no_i)))
}



update_k_i <- function(obs_idx, Ci, n_clust, data, lambda, alpha, verbose=F){
  # Sample a value for the cluster, compatible with the given expression of the full conditional
  probs = rep(0,n_clust)
  for(ii in 1:n_clust){
    val = prob_i_k(obs_idx, ii, Ci, data, lambda, alpha, verbose=F)
    if(val == -1){
      if(verbose)
        print("Atomic cluster detected!")
      return(ii) #if it's the only obs in that cluster it is not reallocated
    }else{
      probs[ii] = val
    }
  }
  
  k = sample(n_clust, 1, replace = TRUE, prob = probs)
  return(k)
}




#### MH for our distro

gibbs_sampler <- function(N, N_burnin, x0, data, lambda, alpha, k, verbose = F){
  #compute MH samples of the posterior for cluter indexes Ci
  #INPUTS
  # N = number of samples
  # N_burnin = number of initial elements of the chain to be discarded
  # x0 = initial labels
  # data = functional data [number of observations * number of points of a single observation ]
  # lambda = weigth of the likelihood wrt the prior
  # alpha = smoothing factor
  # k = number of cluster
  
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
    
    for(t in 1:n)#update one component at the time
      X_temp[t] = update_k_i(t, X_temp, k, data, lambda, alpha, verbose)
    
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

### TEST 
c_opt_test = rep(2,100)
for(ii in c(95,96,97,98,99,100))
  c_opt_test[ii]=1
c_opt_test

#test run
Ci <- gibbs_sampler(10, 1, c_opt, data, 1, alpha, 2, verbose = T)

#actual run
Ci <- gibbs_sampler(8000, 800, c_opt, data, 1, alpha, 2, verbose = T)


save(Ci, file="indexes_Ci_uniform_lambda1_updated.RData")

save.image("chain_Ci_data1_uniform.RData")

