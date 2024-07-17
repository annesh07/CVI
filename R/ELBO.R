#' Title Expected lower bound (ELBO)
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
#'
#' @return value of the ELBO function
#' @export
#'
#' @examples ELBO(N = 100, D = 2, T0 = 10, s1 = 0.01, s2 = 0.01, L20 = 0.01,
#' X = matrix(c(rep(0, 50), rep(10, 50), rep(0, 50), rep(10, 50)), nrow = 100,
#' ncol = 2), W1 = 0.01, W2 = 0.01, L1 = matrix(0.01, nrow = 10, ncol = 2),
#' L2 = matrix(0.01, nrow = 10, ncol = 1),
#' Plog = matrix(-3, nrow = 100, ncol = 10))
ELBO <- function(N, D, T0, s1, s2, L20, X, W1, W2, L1, L2, Plog){
  Mu0 <- matrix(c(rep(0,D)), nrow=1)
  C00 <- diag(D)/L20
  Mu00 <- Mu0%*%solve(C00)
  C0 <- diag(D)
  W1 <- W1
  W2 <- W2
  L1 <- matrix(L1, nrow = T0, ncol = D)
  L2 <- matrix(L2, nrow = T0, ncol = 1)
  Plog <- matrix(Plog, nrow = N, ncol = T0)

  P0 <- exp(Plog)
  #the alpha
  e0 <- s1*log(s2) - lgamma(s1) + (s1 - 1)*(digamma(W1)-log(W2)) - s2*(W1/W2)

  #the z's
  eni <- colSums(P0)
  vni <- colSums(P0*(1 - P0))
  enj <- rowSums(apply(P0, 1, f0))
  vnj <- rowSums(apply(P0, 1, f1))
  e10 <- lgamma(1 + eni) + 0.5*trigamma(1 + eni)*vni +
    lgamma((W1/W2) + enj) + 0.5*trigamma((W1/W2) + enj)*((W1/(W2^2)) + vnj) -
    lgamma(1 + (W1/W2)+eni+enj) -
    0.5*trigamma(1 + (W1/W2)+eni+enj)*((W1/(W2^2))+vni+vnj)
  e1 <- T0*(digamma(W1) - log(W2)) + sum(e10)

  #the eta's
  e20 <- rep(NA, T0)
  inv_C00 <- solve(C00)
  for (i in 1:T0){
    e20[i] <- -D/2*log(2*pi) + D*0.5*log(L20) -
      0.5*(((L1[i,, drop=FALSE]/L2[i,1]) %*% inv_C00%*%(t(L1[i,, drop=FALSE])/
      L2[i,1])) + sum(diag(inv_C00/L2[i,1]))) + Mu0 %*% inv_C00 %*%
      (t(L1[i,, drop=FALSE])/L2[i,1]) - 0.5*Mu0 %*% inv_C00 %*% t(Mu0)
  }
  e2 <- sum(e20)

  #the X's
  #e30 <- matrix(NA, nrow = N, ncol = T0)
  inv_C0 <- solve(C0)
  L21 <- sweep(L1, 1, L2, "/")
  e30 <- (-0.5*X %*% inv_C0 %*% t(X)) %*% P0
  e31 <- (L21 %*% inv_C0 %*% t(X)) %*% P0
  e32 <- P0 %*% (- 0.5*L21 %*% inv_C0 %*% t(L21))
  e33 <- sweep(P0, 2, -0.5*L2*sum(diag(inv_C0)), "*")
  e3 <- sum(P0*(-0.5*D*log(2*pi) - 0.5*determinant(C0, logarithm=TRUE)$modulus))
    + sum(e30) + sum(e31) + sum(e32) + sum(e33)

  #the variationa distributions
  e40 <- W1*log(W2) - lgamma(W1) + (W1-1)*(-log(W2) + digamma(W1)) - W1
  e41 <- sum(exp(Plog)*Plog)
  e42 <- sum(D*log(L2)/2) - 0.5*D*T0*(log(2*pi)+1)
  e4 <- e40 + e41 + e42

  return(c("e0"=e0, "e1"=e1, "e2"=e2, "e3"=e3, "me4"=-e4))
}
