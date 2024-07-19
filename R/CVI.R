#' Title Collapsed Variational Inference (CVI)
#'
#' @param N number of observed data
#' @param D dimension of each observed data
#' @param T0 total clusters fixed for the variational distribution
#' @param s1 shape parameter for the prior distribution of alpha
#' @param s2 rate parameter for the prior distribution of alpha
#' @param L20 precision parameter for the prior distribution of eta_i's
#' @param X the observed data, N x D matrix
#' @param W1 shape parameter for the posterior distribution of alpha
#' @param W2 rate parameter for the posterior distribution of alpha
#' @param L1 1st parameter for the posterior distribution of eta_i's
#' @param L2 Precision parameter for the posterior distribution of eta_i's
#' @param Plog Logarithm of Latent allocation probability matrix
#' @param maxit Maximum number of iterations the updates are run
#'
#' @return the posterior value of alpha, the number of clusters from the
#' latent allocation probability, the proportions of each cluster and the
#' latent allocation probability matrix
#' @export
#'
#' @examples CVI(N = 100, D = 2, T0 = 10, s1 = 0.01, s2 = 0.01,
#'             L20 = 0.01,
#'             X = matrix(c(rep(0, 50), rep(10, 50), rep(0, 50), rep(10, 50)),
#'                 nrow = 100, ncol = 2),
#'             W1 = 0.01, W2 = 0.01,
#'             L1 = matrix(0.01, nrow = 10, ncol = 2),
#'             L2 = matrix(0.01, nrow = 10, ncol = 1),
#'             Plog = matrix(-3, nrow = 100, ncol = 10), maxit = 1000)
CVI <- function(N, D, T0, s1, s2, L20, X, W1, W2, L1, L2, Plog, maxit){
  C0 <- diag(D)
  #the mean vector of the parameters eta_i
  Mu0 <- matrix(c(rep(0,D)), nrow=1)
  #the covariance matrix of the parameters eta_i
  C00 <- diag(D)/L20
  Mu00 <- Mu0%*%solve(C00)
  #store the output of ELBO function for every iteration of updates
  f <- list()
  f[[1]] <- ELBO(N, D, T0, s1, s2, L20, X, W1, W2, L1, L2, Plog)

  for (m in 1:maxit){
    #updating the latent probability values
    P0 <- exp(Plog)
    #different updates for i = 1, i = {2, ..., T0-1} and i = T0
    L21 <- sweep(L1, 1, L2, "/")
    P230 <- L21 %*% t(X)
    P231 <- diag(-0.5*L21 %*% t(L21))
    for (n in 1:N){
      #update of the n^th vector is done by considering all the vecors except
      #the n^th one
      P1 <- P0[-n,]

      p_eni <- colSums(P1)
      p_vni <- colSums(P1*(1-P1))
      p_enj <- rowSums(apply(P1, 1, f0))
      p_vnj <- rowSums(apply(P1, 1, f1))

      P20 <- log(1 + p_eni) - p_vni/((1 + p_eni)^2) - log(1 + p_eni + p_enj +
        (W1/W2)) + (p_vni + p_vnj + (W1/(W2^2)))/((1 + p_eni + p_enj + (W1/W2))^2)

      P21 <- log((W1/W2) + p_enj) - (p_vnj + (W1/(W2^2)))/(((W1/W2) + p_enj)^2) -
        log(1 + p_eni + p_enj + (W1/W2)) +
        (p_vni + p_vnj + (W1/(W2^2)))/((1 + p_eni + p_enj + (W1/W2))^2)
      P22 <- c(0, cumsum(P21)[1:(T0-1)])

      # P231 <- rep(NA, T0)
      # for (i in 1:T0){
      #   P231[i] <- -0.5*L21[i,, drop=FALSE] %*% t(L21[i,, drop=FALSE])
      # }
      P2 <- P20 + P22 + P230[,n] + P231 - 0.5*D/L2
      #log-sum-exp trick
      p0 <- max(P2)
      Plog[n,] <- P2 - p0 - log(sum(exp(P2 - p0)))
    }
    #labelling the probability matrix so that non-zero cluster allocations
    #are present in the beginning of the matrix
    P00 <- exp(Plog)
    Csum <- colSums(P00)
    index <- which(Csum > 0.000001)
    l0 <- length(index)
    for (l in 1:l0){
      Plog[, c(l, index[l])] <- Plog[, c(index[l], l)]
    }
    #final updated and labelled probability matrix
    Pf0 <- exp(Plog)

    for (i in 1:T0){
      #update of the 1st parameter vector of eta_i
      L1[i,] <- Mu00 + t(Pf0[, i, drop=FALSE])%*%X
      #update of the 2nd parameter value of eta_i
      L2[i, 1] <- L20 + sum(Pf0[, i])
    }

    #update of the shape parameter of alpha
    W1 <- s1 + l0 - 1
    #update of the rate parameter of alpha
    a0 <- l0/log(N)
    a_eni <- colSums(Pf0[,1:l0])
    a_vni <- colSums(Pf0[,1:l0]*(1-Pf0[,1:l0]))
    a_enj <- rowSums(apply(Pf0[,1:l0], 1, f0))
    a_vnj <- rowSums(apply(Pf0[,1:l0], 1, f1))
    W20 <- log(a0 + a_eni[1:(l0 - 1)] + a_enj[1:(l0 - 1)]) -
      0.5*(a_vni[1:(l0 - 1)] + a_vnj[1:(l0 - 1)])/((a0 + a_eni[1:(l0 - 1)]
      + a_enj[1:(l0 - 1)])^2) - log(a0 + a_enj[1:(l0 - 1)]) +
      0.5*a_vnj[1:(l0 - 1)]/((a0 + a_enj[1:(l0 - 1)])^2)
    W21 <- log(a0 + a_eni[l0]) - 0.5*a_vni[l0]/((a0 + a_eni[l0])^2) -
      log(a0)
    W2 <- s2  + sum(W20) + W21

    f[[m+1]] <- ELBO(N, D, T0, s1, s2, L20, X, W1, W2, L1, L2, Plog)
    if (abs(sum(f[[m]]) - sum(f[[m + 1]])) < 0.000001){
      break
    }
    message("outer loop: ", m,"\n", f[[m + 1]], '\n', sep="")
  }

  alpha0 <- W1/W2
  clustering <- apply(Plog, MARGIN = 1, FUN=which.max)
  clust <- table(clustering)
  clustnum <- length(unique(clustering))

  posterior <- list("alpha"=alpha0, "Clusters"=clustnum, "Proportions"=clust, "Clustering" = Plog)
  optimisation <- list("ELBO" = f)

  output <-  list("posterior" = posterior, "optimisation" = optimisation)
  class(output) <- "CVIoutput"

  return(output)
}
