#### Gibbs sampler WITHOUT update of covariance matrix eigenvalues/eigenvectors



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


prob_i_k_unif_fixed <- function(obs_index, clust_index, Ci, data, lambda, alpha, eig, verbose=FALSE){
  ## compute fullconditional term for PY-EPPF prior
  
  
  idx_k = which(Ci == clust_index) #indexes of elements of cluster k
  if(!(obs_index %in% idx_k)){#check if our 'i' is already in idx_k, if not, add it
    idx_k = c(obs_index,idx_k)
    obs_idx = 1
  } else {
    obs_idx = match(obs_index, idx_k)#if it is in clust k, get position of the first(and only) match 
  }
  
  x_k = data[idx_k,]
 
  
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



update_k_i_unif_fixed <- function(obs_idx, Ci, n_clust, data, lambda, alpha, eig, verbose=F){
  # Sample a value for the cluster, compatible with the given expression of the full conditional
  probs = rep(0,n_clust)
  for(ii in 1:n_clust){
    val = prob_i_k_unif_fixed(obs_idx, ii, Ci, data, lambda, alpha, eig, verbose=F)
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
gibbs_sampler_unif_fixed <- function(N, N_burnin, x0, data, lambda, alpha, k, verbose = F){
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
  
  eig = eigen(cov(data))
  
  while(posz < N){
    
    gibbs_time_init <- Sys.time()
    
    for(t in 1:n)#update one component at the time
      X_temp[t] = update_k_i_unif_fixed(t, X_temp, k, data, lambda, alpha, eig, verbose)
    
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

# DO RUN
k <- 2
alpha <- 1e4
lambda <- 1

#TEST RUN ON CLINICAL DATA
Ci <- gibbs_sampler_unif_fixed(2, 1, c_opt_2, f.Data$data, lambda, alpha, k, verbose = T)
#REAL RUN ON CLINICAL DATA
Ci <- gibbs_sampler_unif_fixed(10000, 1000, c_opt_2, f.Data$data, lambda, alpha, k, verbose = T)

#SAVE RESULTS
save(Ci, file=paste("indexes_Ci_uniform_lambda",lambda,"_alpha",alpha,"_k",k,"_realData.RData", sep=""))

# SHOW misclassification probabilities
plot_misclassification_probabilities(f.Data$argvals, f.Data$data, alpha, Ci)
