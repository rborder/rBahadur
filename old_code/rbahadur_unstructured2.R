#' Binary random variates with unstructured correlations
#'
#' Generate second Bahadur order multivariate Bernoulli random variates with unstructured correlations.
#'
#' @param m number of observations
#' @param mu vector of means
#' @param C correlation matrix
#'
#' @details This generates multivariate Bernoulli (MVB) random vectors with mean vector 'mu' and correlation matrix 'C'. 'mu' must take values in the open unit interval and 'C' must induce a valid second Bahadur order probability distribution. That is, there must exist an MVB probability distribution with first moments 'mu' and standardized central second moments 'C' such that all higher order central moments are zero.
#'
#' @return An \eqn{m}-by-\eqn{m} matrix of binary random variates, where \eqn{m} is the length of 'mu'.
#' @export
#'
#' @examples
rbahadur_unstructured2 <- function(n, mu, C) {
  # rmvb_vec_unpatterned <- function(n, mu, C) {

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

  v <- matrix(Bk1/p,ncol=1)
  c <- 1

  # recursive steps
  for (m in 2:(M-1)) {
    p <- drop(mu[m] + v%*%C[1:(m-1),m])
    #         p <- mu[m] + apply(v,1,function(x) x%*% C[1:(m-1),m])
    if (any(p < 0 | p > 1)) {
      stop('Infeasible probability. Fix inputs.')
    }
    k[ ,m] <- (rand_U[ ,m] <= p)

    tmp_bool <- (k[ ,m]==0)
    p <- tmp_bool*(1-p) + (!tmp_bool)*p
    Bk0 <- tmp_bool*(1-mu[m]) + (!tmp_bool)*mu[m]
    Bk1 <- tmp_bool*(-1) + (!tmp_bool)*1
    #         v <- cbind(v*rep(Bk0/p,each=m-1), matrix(c*Bk1/p,nrow=n))
    v <- cbind(v*Bk0/p, + c*Bk1/p)
    c <- (Bk0/p)*c
  }
  #     p <- mu[M] + apply(v,1,function(x) x%*% C[1:(M-1),M])
  p <- drop(mu[M] + v%*%C[1:(M-1),M])
  if (any(p < 0 | p > 1)) {
    stop(mu[M],' ',p,' Infeasible probability. Fix inputs.')
  }
  k[ ,M] <-(rand_U[ ,M] <= p)

  return(k)
}
