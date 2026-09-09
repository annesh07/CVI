#' CVI implementation for one set of initial parameters
#'
#' @param config List of inputs that are generated if not user-provided
#' @param X the data matrix
#' @param N samples of X
#' @param D dimensions of X
#' @param T0 variational clusters
#' @param prior_shape_alpha shape parameter of Gamma prior for the DP
#' concentration parameter alpha. Default is 0.001
#' @param prior_rate_alpha rate parameter of Gamma prior for the DP
#' concentration parameter alpha. Default is 0.001
#' @param post_shape_alpha initial value for posterior update of shape parameter
#' for alpha. Default is 0.001
#' @param post_rate_alpha initial value for posterior update of ratee parameter
#' for alpha. Default is 0.001
#' @param prior_mean_eta mean vector of MVN prior for the DP mean parameters.
#' Default is zero vector
#' @param post_mean_eta initial value of posterior update for the DP mean
#' parameter
#' @param fixed_variance covariance matrix of the data is considered known (fixed)
#' or unknown.
#' @param covariance_type covariance matrix is considered diagonal or full.
#' @param cluster_specific_covariance covariance matrix is specific to a cluster
#' allocation or it is same over all cluster choices.
#' @param variance_prior_type For unknown and full covariance matrix, choice of
#' matrix prior is either Inverse-Wishart ('IW') or Cholesky-decomposed
#' ('decomposed'). For unknown, full and cluster-specific covariance matrix,
#' choice of matrix prior is either Inverse-Wishart ('IW'), element-wise Gamma
#' and Laplace distributed ('sparse') or element-wise Gamma and Normal
#' distributed ('off-diagonal normal')
#' @param maxit Maximum number of iterations for variational updates
#' @param varargs List of case specific parameters
#'
#' @return a `list` with the following elements:
#'   - `alpha`: posterior DP concentration parameter
#'   - `Cluster number`: number of clusters from posterior probability allocation matrix
#'   - `Cluster Proportion`: cluster proportions from posterior probability allocation matrix
#'   - `log Probability matrix`: log of posterior probability allocation matrix
#'   - `ELBO`: Optimisation of the ELBO function
#'   - `Iterations`: Number of iterations required for convergence
#'   
#' @export
run_single <- function(config, X, N, D, T0, prior_shape_alpha, prior_rate_alpha,
                       post_shape_alpha, post_rate_alpha, 
                       prior_mean_eta, post_mean_eta,
                       fixed_variance, covariance_type, 
                       cluster_specific_covariance,
                       variance_prior_type, maxit, varargs){
  params <- list()
  inverts <- list()
  
  params$N <- N 
  params$D <- D
  params$T0 <- T0
  params$prior_mean_eta <- prior_mean_eta
  params$prior_shape_alpha <- prior_shape_alpha
  params$prior_rate_alpha <- prior_rate_alpha
  params$post_shape_alpha <- post_shape_alpha
  params$post_rate_alpha <- post_rate_alpha
  params$post_mean_eta <- post_mean_eta
  
  #log-probability matrix based on input
  params$log_prob_matrix <- config$log_prob_matrix
  params$P <- t(apply(exp(params$log_prob_matrix), 1, function(x){x/sum(x)}))
  RP <- Rfast::colsums(params$P)
  
  #updating the parameter list based on the conditions
  if(covariance_type == "diagonal") {
    
    if(fixed_variance) {
      params$prior_precision_scalar_eta <- varargs$prior_precision_scalar_eta
      params$post_precision_scalar_eta <- varargs$post_precision_scalar_eta
      params$cov_data <- varargs$cov_data
      params_check(params, fixed_variance, covariance_type,
                   cluster_specific_covariance,
                   variance_prior_type)
      
      C00 <- diag(D)/varargs$prior_precision_scalar_eta #covariance of DP mean parameters
      inverts[["inv_C0"]] <- Rfast::spdinv(varargs$cov_data)
      inverts[["inv_C00"]] <- Rfast::spdinv(C00)
      
    } else {
      params$prior_shape_scalar_cov <- varargs$prior_shape_scalar_cov
      params$prior_rate_scalar_cov <- varargs$prior_rate_scalar_cov
      params$post_shape_scalar_cov <- varargs$post_shape_scalar_cov
      params$post_rate_scalar_cov <- varargs$post_rate_scalar_cov
      params$post_precision_scalar_eta <- varargs$post_precision_scalar_eta
      params$prior_precision_scalar_eta <- varargs$prior_precision_scalar_eta
      params_check(params, fixed_variance, covariance_type,
                   cluster_specific_covariance,
                   variance_prior_type)
      
      C00 <- diag(D)/varargs$prior_precision_scalar_eta #covariance of DP mean parameters
      inverts[["inv_C00"]] <- Rfast::spdinv(C00)
      
    }
    
  } else if(covariance_type == "full") {
    
    if(fixed_variance) {
      params$post_cov_eta <- varargs$post_cov_eta
      params$cov_data <- varargs$cov_data
      params$prior_cov_eta <- varargs$prior_cov_eta
      params_check(params, fixed_variance, covariance_type,
                   cluster_specific_covariance,
                   variance_prior_type)
      
      inverts[["inv_C0"]] <- Rfast::spdinv(varargs$cov_data)
      inverts[["inv_C00"]] <- Rfast::spdinv(varargs$prior_cov_eta)
      
      
    } else {
      if(!cluster_specific_covariance) {
        if(variance_prior_type == "IW"){
          params$prior_df_cov <- varargs$prior_df_cov
          params$prior_scale_cov <- varargs$prior_scale_cov
          params$post_df_cov <- varargs$post_df_cov
          params$post_scale_cov <- varargs$post_scale_cov
          params$post_cov_eta <- varargs$post_cov_eta
          params$prior_cov_eta <- varargs$prior_cov_eta
          params_check(params, fixed_variance, covariance_type,
                       cluster_specific_covariance,
                       variance_prior_type)
          
          inverts[["inv_V0"]] <- Rfast::spdinv(varargs$prior_scale_cov)
          inverts[["inv_C00"]] <- Rfast::spdinv(varargs$prior_cov_eta)
          
          
        } else if (variance_prior_type == "decomposed"){
          params$prior_shape_diag_decomp <- varargs$prior_shape_diag_decomp
          params$prior_rate_diag_decomp <- varargs$prior_rate_diag_decomp
          params$prior_mean_offdiag_decomp <- varargs$prior_mean_offdiag_decomp
          params$prior_var_offdiag_decomp <- varargs$prior_var_offdiag_decomp
          params$post_shape_diag_decomp <- varargs$post_shape_diag_decomp
          params$post_rate_diag_decomp <- varargs$post_rate_diag_decomp
          params$post_mean_offdiag_decomp <- varargs$post_mean_offdiag_decomp
          params$post_var_offdiag_decomp <- varargs$post_var_offdiag_decomp
          params$post_cov_eta <- varargs$post_cov_eta
          params$prior_cov_eta <- varargs$prior_cov_eta
          params_check(params, fixed_variance, covariance_type,
                       cluster_specific_covariance,
                       variance_prior_type)
          
          inverts[["inv_C00"]] <- Rfast::spdinv(varargs$prior_cov_eta)
          
          
        } else {
          stop("'variance_prior_type' can only be either 'IW' or 'decomposed'
               when 'cluster_specific_covariance' is FALSE")
        }
        
      }else{
        if(variance_prior_type == "IW"){
          params$prior_df_cs_cov <- varargs$prior_df_cs_cov
          params$prior_scale_cs_cov <- varargs$prior_scale_cs_cov
          params$post_df_cs_cov <- varargs$post_df_cs_cov
          params$post_scale_cs_cov <- varargs$post_scale_cs_cov
          params$scaling_cov_eta <- varargs$scaling_cov_eta
          params_check(params, fixed_variance, covariance_type,
                       cluster_specific_covariance,
                       variance_prior_type)
          
          
        } else if (variance_prior_type == "sparse"){
          params$prior_shape_d_cs_cov <- config$prior_shape_d_cs_cov
          params$prior_rate_d_cs_cov <- config$prior_rate_d_cs_cov
          params$prior_var_offd_cs_cov <- varargs$prior_var_offd_cs_cov
          params$post_shape_d_cs_cov <- varargs$post_shape_d_cs_cov
          params$post_rate_d_cs_cov <- varargs$post_rate_d_cs_cov
          params$post_var_offd_cs_cov <- varargs$post_var_offd_cs_cov
          params$scaling_cov_eta <- varargs$scaling_cov_eta
          
          params_check(params, fixed_variance, covariance_type,
                       cluster_specific_covariance,
                       variance_prior_type)
          
          
        } else if (variance_prior_type == "off-diagonal normal"){
          params$prior_shape_d_cs_cov <- varargs$prior_shape_d_cs_cov
          params$prior_rate_d_cs_cov <- varargs$prior_rate_d_cs_cov
          params$prior_var_offd_cs_cov <- varargs$prior_var_offd_cs_cov
          params$post_shape_d_cs_cov <- varargs$post_shape_d_cs_cov
          params$post_rate_d_cs_cov <- varargs$post_rate_d_cs_cov
          params$post_mean_offd_cs_cov <- varargs$post_mean_offd_cs_cov
          params$scaling_cov_eta <- varargs$scaling_cov_eta
          params_check(params, fixed_variance, covariance_type,
                       cluster_specific_covariance,
                       variance_prior_type)
          
          
        } else {
          stop("'variance_prior_type' can only be either 'IW' or 'decomposed'
               when 'cluster_specific_covariance' is TRUE")
        }
      }
      
    }
  } else {
    stop("covariance_type can only be either 'diagonal' or 'full'.")
  }
  
  #store the output of ELBO function for every iteration of updates
  elbo_values <- vector("list", maxit + 1)
  elbo_values[[1]] <- ELBO_function(fixed_variance, covariance_type,
                                    cluster_specific_covariance,
                                    variance_prior_type, X, inverts, params)
  
  for (m in 1:maxit){
    
    updated_params <- CVI_update_function(fixed_variance, covariance_type,
                                          cluster_specific_covariance,
                                          variance_prior_type,
                                          X, inverts, params)
    
    params <- updated_params
    elbo_values[[m+1]] <- ELBO_function(fixed_variance, covariance_type,
                                        cluster_specific_covariance,
                                        variance_prior_type, X, inverts, params)
    if (abs(sum(elbo_values[[m]]) - sum(elbo_values[[m + 1]])) < 0.000001 ){
      break
    }
    
    if (m > 3) {
      if (abs(sum(elbo_values[[m+1]]) - sum(elbo_values[[m-2]])) < 1e-8 ||
          abs(sum(elbo_values[[m+1]]) - sum(elbo_values[[m-1]])) < 1e-8) {
        
        message("Stopping: ELBO oscillating or stuck.")
        break
      }
    }
    message("outer loop: ", m,"\n", elbo_values[[m + 1]], '\n', sep="")
  }
  elbo_values <- elbo_values[seq_len(m + 1)]
  W1 <- params$post_shape_alpha
  W2 <- params$post_rate_alpha
  Plog <- params$log_prob_matrix
  
  alpha0 <- W1/W2 #posterior concentration parameter
  clustering <- apply(Plog, MARGIN = 1, FUN=which.max)
  clust <- table(clustering) #clusters with proportions
  clustnum <- length(unique(clustering)) #number of clusters
  
  post_distribution = list()
  L1 <- params$post_mean_eta
  
  if(covariance_type == "diagonal") {
    
    if(fixed_variance) {
      L2 <- params$post_precision_scalar_eta
      
      post_distribution[["Mean"]] = L1/c(L2)
      
    } else {
      G1 <- params$post_shape_scalar_cov
      G2 <- params$post_rate_scalar_cov
      L2 <- params$post_precision_scalar_eta <- varargs$post_precision_scalar_eta
      
      post_distribution[["Mean"]] = L1/c(L2)
      post_distribution[["Precision"]] = (G1/G2)*diag(D)
      #G2/G1 instead of G1/G2 because prior on the precision diagonal scalar
    }
    
  } else if(covariance_type == "full") {
    
    if(fixed_variance) {
      L2 <- params$post_cov_eta
      L21 <- matrix(0, nrow = T0, ncol = D)
      for (i in 1:T0){
        L21[i,] = mat_mult(L1[i,, drop = FALSE], L2[,,i])
      }
      
      post_distribution[["Mean"]] = L21
      
    } else {
      if(!cluster_specific_covariance) {
        if(variance_prior_type == "IW"){
          nu <- params$post_df_cov
          V <- params$post_scale_cov
          L2 <- params$post_cov_eta
          
          L21 <- matrix(0, nrow = T0, ncol = D)
          for (i in 1:T0){
            L21[i,] = mat_mult(L1[i,, drop = FALSE], L2[,,i])
          }
          
          post_distribution[["Mean"]] = L21
          post_distribution[["Precision"]] = nu*V
          
        } else if (variance_prior_type == "decomposed"){
          a1 <- params$post_shape_diag_decomp
          b1 <- params$post_rate_diag_decomp
          mu1 <- params$post_mean_offdiag_decomp
          c1 <- params$post_var_offdiag_decomp
          L2 <- params$post_cov_eta
          
          mean_lower <- matrix(0, nrow = D, ncol = D) #mean matrix of the decomposed
          mean_lower[lower.tri(mean_lower, diag = FALSE)] <- mu1
          sigma_lower <- matrix(0, nrow = D, ncol = D) #var matrix of the decomposed
          sigma_lower[lower.tri(sigma_lower, diag = FALSE)] <- c1
          mean_L <- mean_lower + diag(sqrt(1/b1)*sqrt(pi)/beta(a1,0.5))
          diag(sigma_lower) <- (1/b1)*(a1 - (sqrt(pi)/beta(a1,0.5))^2)
          #expected inverse of C0; covariance matrix of data
          inv_C0 <- mat_mult(mean_L, t(mean_L)) +
            diag(Rfast::rowsums(sigma_lower))
          
          L21 <- matrix(0, nrow = T0, ncol = D)
          for (i in 1:T0){
            L21[i,] = mat_mult(L1[i,, drop = FALSE], L2[,,i])
          }
          
          post_distribution[["Mean"]] = L21
          post_distribution[["Precision"]] = inv_C0
          
        } else {
          stop("'variance_prior_type' can only be either 'IW' or 'decomposed'
               when 'cluster_specific_covariance' is FALSE")
        }
        
      }else{
        if(variance_prior_type == "IW"){
          nu1 <- params$post_df_cs_cov
          V1 <- params$post_scale_cs_cov
          
          V1_inv <- array(apply(V1, 3, function(x){spdinv(x)}), dim = dim(V1))
          #expectation of inverse of data covariance matrix
          inv_C0 <- sweep_3D(V1_inv, nu1, c(D, D, T0))
          
          post_distribution[["Mean"]] = L1
          post_distribution[["Precision_cs"]] = inv_C0
          
        } else if (variance_prior_type == "sparse"){
          a1 <- params$post_shape_d_cs_cov
          B1 <- params$post_rate_d_cs_cov
          
          #expectation of inverse of C0, data covariance matrix
          inv_C0 <- array(0, c(D, D, T0))
          for (i in 1:T0){
            inv_C0[,,i] <- temp <- Rfast::Diag.matrix(D, a1[1,i]/B1[i,])
          }
          
          post_distribution[["Mean"]] = L1
          post_distribution[["Precision_cs"]] = inv_C0
          post_distribution[["Precision_cs_a1"]] = a1
          post_distribution[["Precision_cs_B1"]] = B1
          
        } else if (variance_prior_type == "off-diagonal normal"){
          a1 <- params$post_shape_d_cs_cov
          B1 <- params$post_rate_d_cs_cov
          C1 <- params$post_mean_offd_cs_cov
          
          #expectation of inverse of data covariance matrix
          inv_C0 <- array(0, c(D, D, T0))
          for (i in 1:T0){
            inv_C0[,,i] <- temp <- Rfast::Diag.fill(C1[,,i], a1[1,i]/B1[i,])}
          
          post_distribution[["Mean"]] = L1
          post_distribution[["Precision_cs"]] = inv_C0
          
        } else {
          stop("'variance_prior_type' can only be either 'IW' or 'decomposed'
               when 'cluster_specific_covariance' is TRUE")
        }
      }
      
    }
  } else {
    stop("covariance_type can only be either 'diagonal' or 'full'.")
  }
  
  
  posterior <- c(list("alpha" = alpha0, "Cluster number" = clustnum,
                      "Cluster Proportion" = clust,
                      "log Probability matrix" = Plog), post_distribution)
  logBayes <- as.list(elbo_values[length(elbo_values)])$e_data - 
    as.list(elbo_values[1])$e_data
  optimisation <- list("ELBO" = elbo_values,
                       "Iterations" = (length(elbo_values)-1),
                       "logBF" = logBayes)
  
  output <-  list("posterior" = posterior, "optimisation" = optimisation)
  class(output) <- "CVIoutput_partial"
  
  return(output)
}
