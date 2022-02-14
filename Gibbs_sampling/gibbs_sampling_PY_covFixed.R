
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


prob_i_k_PY <- function(obs_index, clust_index, Ci, data, sig, lambda, alpha, theta, eig, verbose=FALSE){
  ## compute fullconditional term for PY-EPPF prior
  
  N <- dim(data)[1] #total numer of observations
  
  if(!(clust_index%in%Ci)){ # check if it is a new cluster
    # print(paste("posterior value for clsuter",clust_index,"=",(clust_index*sig+theta)/(N + theta)))
    return((clust_index*sig+theta)/(N + theta))
  }else{#if it is not a new cluster
    idx_k = which(Ci == clust_index) #indexes of elements of cluster k
    if(!(obs_index %in% idx_k)){#check if our 'i' is already in idx_k, if not, add it
      idx_k = c(obs_index,idx_k)
      obs_idx = 1 #location of the observation into consideration
    } else {
      obs_idx = match(obs_index, idx_k)#if it is in clust k, get position of the first(and only) match 
    }
    
   
    
    x_k = data[idx_k,]
    if(length(idx_k)>1) #i.e. at least 2
      centroid = colMeans(x_k)
    else{
      # centroid = x_k
      # return(((clust_index-1)*sig+theta)/(N + theta))
      return(1-sig)
      
      if(verbose){
        print("Cluster had only 1 observation")
      }
    }
      
      start_time <- Sys.time()
      disc_in = discrepancy_within(x_k, centroid, alpha, eig)
      end_time <- Sys.time()
    
    
    if(verbose){
      print(paste("Computing full discrepancy in ",end_time - start_time,"s"))
    }
    # print(paste("numerosità x_k = ",length(idx_k)))
    # print(paste("osservazione rimossa = ",obs_idx))
    # print(paste(" numerosità X_k_no_i = ", dim(x_k_no_i)[1]))
    
      x_k_no_i = x_k[-obs_idx,]
      if(length(idx_k) > 2){
        #i.e. if x_k_no_i has 1 row
        centroid_no_i = colMeans(x_k_no_i)
        start_time <- Sys.time()
        disc_no_i <- discrepancy_within(x_k_no_i, centroid_no_i, alpha, eig)
        end_time <- Sys.time()
        if(verbose){
          print(paste("Computing discrepancy without obs ",obs_index," in ",end_time - start_time,"s"))
        }
      }else{
        disc_no_i <- 0
        if(verbose){
          print(paste("Computation of discrepancy without obs not necessary, idx_k length:",length(idx_k)))
        }
      }
      
      return((length(idx_k)-sig)*exp(-lambda*(disc_in - disc_no_i)))
    
  }
}

calcolaIndici <- function(Ci){
  v <- unique(Ci)
  flag = 0
  t = 1
  to_add = max(v)+1
  while(t < max(v) &  flag==0){
    if(!(t%in%v)){
      to_add = t
      flag = 1
    }
    t = t+1
  }
  return(c(v,to_add))
}


update_k_i_PY <- function(obs_idx, Ci, data, sig, lambda, alpha, theta, eig, verbose){
  # Sample a value for the cluster, compatible with the given expression of the full conditional
  indici = calcolaIndici(Ci)
  probs = rep(0,length(indici))
  for(ii in indici) #for all the cluster up to now but possibly a new one
    probs[ii] =prob_i_k_PY(obs_idx, ii, Ci, data, sig, lambda, alpha, theta, eig,verbose=F) #FOR PY
  
  idx0 = 0
  if(verbose){
    print("probs provvisorio:")
    print(probs)
  }
  
  for(tt in 1:length(probs)){
    if(probs[tt]<1e-6 | is.na(probs[tt]))
      idx0 = c(idx0,tt)
  }
  if(tail(idx0,n=1)!=0){
    idx0 = idx0[-1]
    if(verbose){
      print(paste("removing observation#:",idx0))
    }
    probs = probs[-idx0]
  }
  
  
  
  if(verbose){
    print(paste("Values of probs obs#",obs_idx))
    print(probs)
    print("indici")
    print(indici)
    #print(paste("Probabilità di allocazione elem #",obs_idx," in #",length(indici)," clusters: ",probs))
  }
  
  k = sample(indici, 1, replace = TRUE, prob = probs)
  return(c(k,length(indici)))
}



gibbs_sampler_PY <- function(N, N_burnin, x0, data, lambda, alpha, eig, sig, theta, verbose = F){
  #compute MH samples of the posterior for cluter indexes Ci
  #INPUTS
  # N = number of samples
  # N_burnin = number of initial elements to discard
  # x0 = initial labels
  
  #k = 2 #FIXED HERE FOR NOW
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
    max_n_clust = 1
    
    for(t in 1:n){#update one component at the time
       v_out = update_k_i_PY(t, X_temp, data, sig, lambda, alpha, theta, eig, verbose=F)
       X_temp[t] = v_out[1]
       max_n_clust = max(max_n_clust,v_out[2])
    
    }
      
    
    if(iterMH > N_burnin){
      posz = iterMH - N_burnin
      samplee[posz,] = X_temp
    }
    
    
    gibbs_time_finish <-Sys.time()
    
    if(verbose){
      print(paste("Gibbs iteration N = ",iterMH," in ",gibbs_time_finish-gibbs_time_init,"s"))
      print(paste("Actual number of cluster:",max_n_clust))
    }
    iterMH = iterMH + 1
  }
  return(samplee)
}


#### TEST SIMULATED DATA

Ci <- gibbs_sampler_PY(3, 0, c_opt_test, data, .7, alpha, eigen(K_1), 0.25, .5, verbose = T)

#### TEST REAL DATA

Ci <- gibbs_sampler_PY(10, 0, c_opt_3up, f.Data$data, .3, alpha, eigen(cov(f.Data$data)), 0.25, .1, verbose = T)


#FULL RUN SIMULATED DATA
Ci <- gibbs_sampler_PY(3, 0, c_opt_test, data, .7, alpha, eigen(K_1), 0.25, .5, verbose = T)


#FULL RUN REAL DATA
Ci <- gibbs_sampler_PY(10000, 1000, c_opt_3up, f.Data$data, .4, alpha, eigen(cov(f.Data$data)), 0.75, .1, verbose = T)


save(Ci, file="indexes_Ci_PY_TEST.RData")



























#### OLD CODE

prob_i_k <- function(obs_index, clust_index, Ci, data, sig, lambda, alpha, eig, verbose=FALSE){
  ## compute fullconditional term for PY-EPPF prior
  
  
  idx_k = which(Ci == clust_index) #indexes of elements of cluster k
  if(!(obs_index %in% idx_k)){#check if our 'i' is already in idx_k, if not, add it
    idx_k = c(obs_index,idx_k)
    obs_idx = 1
  } else {
    obs_idx = match(obs_index, idx_k)#if it is in clust k, get position of the first(and only) match 
  }
  
  x_k = data[idx_k,]
  if(length(idx_k)>1)
    centroid = colMeans(x_k)
  else
    centroid = x_k
  
  
  start_time <- Sys.time()
  disc_in = discrepancy_within(x_k, centroid, alpha, eig)
  end_time <- Sys.time()
  
  if(verbose){
    print(paste("Computing full discrepancy in ",end_time - start_time,"s"))
  }
  # print(paste("numerosità x_k = ",length(idx_k)))
  # print(paste("osservazione rimossa = ",obs_idx))
  # print(paste(" numerosità X_k_no_i = ", dim(x_k_no_i)[1]))
  if(length(idx_k)>1){
    x_k_no_i = x_k[-obs_idx,]
    if(length(idx_k) > 2){
      #i.e. if x_k_no_i has 1 row
      centroid_no_i = colMeans(x_k_no_i)
    }else{
      centroid_no_i = x_k_no_i
    }
    start_time <- Sys.time()
    disc_no_i <- discrepancy_within(x_k_no_i, centroid_no_i, alpha, eig)
    end_time <- Sys.time()
    
    if(verbose){
      print(paste("Computing discrepancy without obs ",obs_index," in ",end_time - start_time,"s"))
    }
    return(exp(-lambda*(disc_in - disc_no_i)))
  }else{
    disc_no_i <- 0
    
    if(verbose){
      print("Cluster had only 1 observation, disc_no_i set to zero")
    }
    
    return(exp(-lambda*(disc_in - disc_no_i))) # approx to be checked with mario
  }
}


