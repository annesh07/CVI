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
#' @param Plog Latent probability values
#' @param maxit Maximum number of iterations the updates are run
#'
#' @return the posterior value of alpha, the number of clusters from latent
#' probability allocation and the proportions of each cluster
#' @export
#'
#' @examples CVI(N = 50, D = 2, T0 = 10, s1 = 0.01, s2 = 0.01, L20 = 0.01,
#' X = matrix(1, nrow = 50, ncol = 2), W1 = 0.01, W2 = 0.01,
#' L1 = matrix(0.01, nrow = 10, ncol = 2), L2 = matrix(0.01, nrow = 10, ncol = 1
#' ), Plog = matrix(-3, nrow = 50, ncol = 10), maxit = 1000)
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
    for (n in 1:N){
      #update of the n^th vector is done by considering all the vecors except
      #the n^th one
      P1 <- P0[-n,]
      P2 <- rep(NA, T0)
      #first element

      #expectation of the indicator variables N_i, where i = 1
      p1qni <- sum(P1[, 1])
      #variance of the corresponding indicator variables
      p1vqni <- sum(P1[, 1]*(1 - P1[, 1]))
      #expectation of the indicator variables N_>i, where i = 1
      p1qnj <- sum(P1[, 2:T0])
      #variance of the corresponding variables
      I <- 1
      p1vqnj <- sum(apply(P1, 1, f0, I=I))
      #update of P[n, 1]
      P2[1] <- log(1 + p1qni) - 0.5*p1vqni/((1 + p1qni)^2) - log(1 + p1qni + p1qnj
        + (W1/W2)) + 0.5*(p1vqni + p1vqnj + (W1/(W2^2)))/((1 + p1qni + p1qnj +
        (W1/W2))^2) + (L1[1,, drop=FALSE]/L2[1,1]) %*% t(X[n,, drop=FALSE]) -
        0.5*((L1[1,, drop=FALSE]/L2[1,1]) %*% (t(L1[1,, drop=FALSE])/L2[1,1]) +
        D/L2[1,1])

      #other elements except the last
      for (i in 2:(T0 - 1)){
        #expectation of the indicator variables N_i
        pqni <- sum(P1[, i])
        #variance of the corresponding indicator variables
        pvqni <- sum(P1[, i]*(1 - P1[, i]))
        #expectation of the indicator variables N_>i
        pqnj <- sum(P1[, (i + 1):T0])
        #variance of the corresponding variables
        I <- i
        pvqnj <- sum(apply(P1, 1, f0, I=I))
        # 1st part of the update of P[n, i]
        P20 <- log(1 + pqni) - 0.5*pvqni/((1 + pqni)^2) - log(1 + pqni + pqnj +
          (W1/W2)) + 0.5*(pvqni + pvqnj + (W1/(W2^2)))/((1 + pqni + pqnj +
          (W1/W2))^2)

        #for i>1, the update includes j terms where j<i
        P21 <- rep(NA, (i - 1))
        for (j in 1:(i - 1)){
          #expectation of N_j
          pjqni <- sum(P1[, j])
          #corresponding variance
          pjvqni <- sum(P1[, j]*(1 - P1[, j]))
          # expectation of N_>j
          pjqnj <- sum(P1[, (j + 1):T0])
          #corresponding variance
          I=j
          pjvqnj <- sum(apply(P1, 1, f0, I=I))
          P21[j] <- log((W1/W2) + pjqnj) - 0.5*(pjvqnj + (W1/(W2^2)))/(((W1/W2) +
            pjqnj)^2) - log(1 + (W1/W2) + pjqni + pjqnj) + 0.5*(pjvqni + pjvqnj +
            (W1/(W2^2)))/((1 + (W1/W2) + pjqni + pjqnj)^2)
        }
        #update of P[n, i]
        P2[i] <- P20 + sum(P21) +
          (L1[i,, drop=FALSE]/L2[i,1]) %*% t(X[n,, drop=FALSE]) -
          0.5*((L1[i,, drop=FALSE]/L2[i,1]) %*% (t(L1[i,, drop=FALSE])/L2[i,1])
               + D/L2[i,1])
      }

      #last term
      Pt0 <- log(1 + sum(P1[, T0])) - 0.5*sum(P1[, T0]*(1 - P1[, T0]))/((1 +
        sum(P1[, T0]))^2) - log(1 + sum(P1[, T0]) + (W1/W2)) +
        0.5*(sum(P1[, T0]*(1 - P1[, T0]))+(W1/(W2^2)))/((1 + sum(P1[, T0]) +
        (W1/W2))^2)

      #similar to previous step, j = T0-1 terms to be computed here as well
      Pt1 <- rep(NA, (T0 - 1))
      for (j in 1:(T0 - 1)){
        ptqni <- sum(P1[, j])
        ptvqni <- sum(P1[, j]*(1 - P1[, j]))
        ptqnj <- sum(P1[,(j+1):T0])
        I <- j
        ptvqnj <- sum(apply(P1, 1, f0, I=I))
        Pt1[j] <- log((W1/W2) + ptqnj) - 0.5*(ptvqnj + (W1/(W2^2)))/(((W1/W2) +
          ptqnj)^2) - log(1 + (W1/W2) + ptqni + ptqnj) +
          0.5*(ptvqni + ptvqnj + (W1/(W2^2)))/((1 + (W1/W2) + ptqni + ptqnj)^2)
      }
      #update of P[n, T0]
      P2[T0] <- Pt0 + sum(Pt1) +
        (L1[T0,, drop=FALSE]/L2[T0,1]) %*% t(X[n,, drop=FALSE]) -
        0.5*((L1[T0,, drop=FALSE]/L2[T0,1]) %*%
        (t(L1[T0,, drop=FALSE])/L2[T0,1]) + D/L2[T0,1])
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
    W20 <- rep(NA, l0)
    for (i in 1:(l0-1)){
      #expectation of N_i
      wqni <- sum(Pf0[, i])
      #corresponding variance
      wvqni <- sum(Pf0[, i]*(1 - P00[, i]))
      #expectation of N_>i
      wqnj <- sum(Pf0[, (i + 1):T0])
      #corresponding variance
      I <- i
      wvqnj <- sum(apply(Pf0, 1, f0, I = I))

      W20[i] <- log(0.0000001 + wqni + wqnj) -
        0.5*(wvqnj + wvqni)/((0.0000001 + wqni + wqnj)^2) -
        log(0.0000001 + wqnj) + 0.5*wvqnj/((0.0000001 + wqnj)^2)
    }
    W20[l0] <- log(0.0000001 + sum(Pf0[, l0])) -
      0.5*sum(Pf0[, l0]*(1 - Pf0[, l0]))/((0.0000001 + sum(Pf0[, l0]))^2) -
      log(1.0000001)
    W2 <- s2 + sum(W20)

    f[[m+1]] <- ELBO(N, D, T0, s1, s2, L20, X, W1, W2, L1, L2, Plog)
    if (abs(sum(f[[m]]) - sum(f[[m + 1]])) < 0.000001){
      break
    }
    cat("outer loop: ", m,"\n", sep="")
    print(f[[m + 1]])
    cat('\n')
  }

  plot(sapply(f, sum)[-1], type="l")
  alpha0 <- W1/W2
  clustering <- apply(Plog, MARGIN = 1, FUN=which.max)
  clust <- table(clustering)
  clustnum <- length(unique(clustering))
  list0 <- list("alpha"=alpha0, "Clusters"=clustnum, "Proportions"=clust)
  return(list0)
}
