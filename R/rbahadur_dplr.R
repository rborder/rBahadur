#' Binary random variates with Diagonal Plus Low Rank (dplr) correlations
#'
#' Generate second Bahadur order multivariate Bernoulli random variates with Diagonal Plus Low Rank (dplr) correlation structures.
#'
#' @param n number of observations
#' @param mu vector of means
#' @param U outer product component matrix
#'
#' @details This generates multivariate Bernoulli (MVB) random vectors with mean vector 'mu' and correlation matrix \eqn{C = D + U U^T} where \eqn{D} is a diagonal matrix with values dictated by 'U'. 'mu' must take values in the open unit interval and 'U' must induce a valid second Bahadur order probability distribution. That is, there must exist an MVB probability distribution with first moments 'mu' and standarized central second moments \eqn{C} such that all higher order central moments are zero.
#'
#' @return An \eqn{n}-by-\eqn{m} matrix of binary random variates, where \eqn{m} is the length of 'mu'.
#' @export
#'
#' @examples
rbahadur_dplr <- function(n, mu, U) {

  ## place holder for assert statements

  M <- length(mu)

  k <- matrix(NaN, nrow=n, ncol=M)
  rand_U <- matrix(runif(M*n), nrow=n, ncol=M)

  # initial step
  p <- rep(mu[1],n)
  k[ ,1] <- as.numeric(rand_U[,1] <= p)

  tmp_bool <- (k[ ,1]==0)
  p <- tmp_bool*(1-p) + (!tmp_bool)*p
  Bk0 <- tmp_bool*(1-mu[1]) + (!tmp_bool)*mu[1]
  Bk1 <- tmp_bool*(-1) + (!tmp_bool)*1

  x <- Bk1*U[1]/p
  c <- 1

  # recursive steps
  for (m in 2:(M-1)) {
    p <- mu[m] + x * U[m]
    if (any(p < 0 | p > 1)) {
      stop('Infeasible probability. Fix inputs.')
    }
    k[ ,m] <- (rand_U[ ,m] <= p)

    tmp_bool <- (k[ ,m]==0)
    p <- tmp_bool*(1-p) + (!tmp_bool)*p
    Bk0 <- tmp_bool*(1-mu[m]) + (!tmp_bool)*mu[m]
    Bk1 <- tmp_bool*(-1) + (!tmp_bool)*1

    x <- (x*Bk0 + c*Bk1*U[m])/p
    c <- (Bk0/p)*c
  }
  p <- mu[M] + x*U[M]
  if (any(p < 0 | p > 1)) {
    stop(mu[M],' ',p,' Infeasible probability. Fix inputs.')
  }
  k[ ,M] <-(rand_U[ ,M] <= p)

  return(k)
}


