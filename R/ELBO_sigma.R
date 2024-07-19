#' Title Expected lower bound (ELBO) for varying sigma^2
#'
#' @param N number of observed data
#' @param D dimension of each observed data
#' @param T0 total clusters fixed for the variational distribution
#' @param s1 shape parameter for the prior distribution of alpha
#' @param s2 rate parameter for the prior distribution of alpha
#' @param L20 precision parameter for the prior distribution of eta_i's
#' @param b1 shape parameter for prior distribution of sigma^2
#' @param b2 rate parameter for prior distribution of sigma^2
#' @param X the observed data, N x D matrix
#' @param W1 shape parameter for the posterior distribution of alpha
#' @param W2 rate parameter for the posterior distribution of alpha
#' @param L1 1st parameter for the posterior distribution of eta_i's
#' @param L2 Precision parameter for the posterior distribution of eta_i's
#' @param G1 shape parameter for the posterior distribution of sigma^2
#' @param G2 rate parameter for the posterior distribution of sigma^2
#' @param Plog Latent probability values
#'
#' @return value of the ELBO function
#' @export
#'
#' @examples ELBO_sigma(N = 100, D = 2, T0 = 10, s1 = 0.01, s2 = 0.01,
#'              L20 = 0.01, b1 = 0.01, b2 = 0.01,
#'              X = matrix(c(rep(0, 50), rep(10, 50), rep(0, 50), rep(10, 50)),
#'                  nrow = 100, ncol = 2),
#'              W1 = 0.01, W2 = 0.01, L1 = matrix(0.01, nrow = 10, ncol = 2),
#'              L2 = matrix(0.01, nrow = 10, ncol = 1), G1 = 0.01, G2 = 0.01,
#'              Plog = matrix(-3, nrow = 100, ncol = 10))
ELBO_sigma <- function(N, D, T0, s1, s2, L20, b1, b2, X, W1, W2, L1, L2, G1, G2,
                       Plog){
  Mu0 <- matrix(c(rep(0,D)), nrow=1)
  C00 <- diag(D)/L20

  G1 <- G1
  G2 <- G2
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
  # e20 <- rep(NA, T0)
  inv_C00 <- solve(C00)
  L21 <- sweep(L1, 1, L2, "/")
  # for (i in 1:T0){
  #   e20[i] <- -D/2*log(2*pi) + D*0.5*log(L20) -
  #     0.5*(((L1[i,, drop=FALSE]/L2[i,1]) %*% inv_C00%*%(t(L1[i,, drop=FALSE])/
  #     L2[i,1])) + sum(diag(inv_C00/L2[i,1]))) + Mu0 %*% inv_C00 %*%
  #     (t(L1[i,, drop=FALSE])/L2[i,1]) - 0.5*Mu0 %*% inv_C00 %*% t(Mu0)
  # }
  e20 <- diag(-0.5*L21 %*% inv_C00 %*% t(L21))
  e21 <- - 0.5*sum(diag(inv_C00))/L2
  e22 <- Mu0 %*% inv_C00 %*% t(L21)
  e2 <-  T0*(-D/2*log(2*pi) + D*0.5*log(L20)- 0.5*Mu0 %*% inv_C00 %*% t(Mu0)) +
    sum(e20) + sum(e21) + sum(e22)
  # e2 <- sum(e20)

  #the X's
  #e30 <- matrix(NA, nrow = N, ncol = T0)
  #L21 <- sweep(L1, 1, L2, "/")
  e30 <- sweep(P0, 1, -0.5*(G1/G2)*diag(X %*% t(X)), "*")
  e31 <- P0*t((G1/G2)*(L21 %*% t(X)))
  e32 <- sweep(P0, 2, -0.5*(G1/G2)*diag(L21 %*% t(L21)), "*")
  e33 <- sweep(P0, 2, -0.5*(G1/G2)*D/L2, "*")
  e3 <- sum(P0*(-0.5*D*log(2*pi) + 0.5*D*(digamma(G1) - log(G2))))
        + sum(e30) + sum(e31) + sum(e32) + sum(e33)

  #the sigma^2
  e4 <- b1*log(b2) - lgamma(b1) + (b1 - 1)*(digamma(G1) - log(G2)) - b2*G1/G2

  #the variationa distributions
  e50 <- W1*log(W2) - lgamma(W1) + (W1-1)*(digamma(W1) - log(W2)) - W1
  e51 <- sum(exp(Plog)*Plog)
  e52 <- sum(D*log(L2)/2) - 0.5*D*T0*(log(2*pi)+1)
  e53 <- G1*log(G2) - lgamma(G1) + (G1-1)*(digamma(G1) - log(G2)) - G1
  e5 <- e50 + e51 + e52 + e53

  return(c("e0"=e0, "e1"=e1, "e2"=e2, "e3"=e3, "e4"=e4, "me5"=-e5))
}
