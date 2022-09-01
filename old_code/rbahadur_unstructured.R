#' Binary random variates with unstructured correlations
#'
#' Generate second Bahadur order multivariate Bernoulli random variates with unstructured correlations.
#'
#' @param n number of observations
#' @param mu vector of means
#' @param C correlation matrix
#'
#' @details This generates multivariate Bernoulli (MVB) random vectors with mean vector 'mu' and correlation matrix 'C'. 'mu' must take values in the open unit interval and 'C' must induce a valid second Bahadur order probability distribution. That is, there must exist an MVB probability distribution with first moments 'mu' and standardized central second moments 'C' such that all higher order central moments are zero.
#'
#' @return An \eqn{n}-by-\eqn{m} matrix of binary random variates, where \eqn{m} is the length of 'mu'.
#' @export
#'
#' @examples
rb_unstr <- function(n, mu, C) {

  M <- length(mu)
  k <- matrix(NaN, nrow=n, ncol=M)
  rand_U <- matrix(runif(M*n), nrow=n, ncol=M)

  # initial step
  p <- rep(mu[1],n)
  k[ ,1] <- as.numeric(rand_U[,1] <= p)


  Bk0_values_matrix <- rbind(1-mu,mu)
  Bk1_bvec <- c(-1,1)

  tmp_bool <- (k[ ,1]==0)
  p <- tmp_bool*(1-p) + (!tmp_bool)*p

  Bk0_mat <- Bk1_mat <- matrix(as.numeric(NA),
                               nrow=M,
                               ncol=n)
  Bk0 <- Bk0_mat[1,] <- Bk0_values_matrix[cbind(k[ ,1]+1,
                                                rep(1,n))]

  Bk1 <- Bk1_mat[1,] <- Bk1_bvec[k[ ,1]+1]

  v <- matrix(Bk1/p,ncol=1)
  c <- 1

  # recursive steps
  for (m in 2:(M-1)) {
    p <- drop(mu[m] + v%*%C[1:(m-1),m])
    if (any(p < 0 | p > 1)) {
      stop('Infeasible probability. Fix inputs.')
    }
    k[ ,m] <- (rand_U[ ,m] <= p)

    tmp_bool <- (k[ ,m]==0)
    p <- tmp_bool*(1-p) + (!tmp_bool)*p
    Bk0_mat[m,] <- Bk0 <- Bk0_values_matrix[cbind(k[ ,m]+1,
                                                  rep(m,n))]
    Bk1_mat[m,] <- Bk1 <-  Bk1_bvec[k[ ,m]+1]
    v <- cbind(v*Bk0/p, + c*Bk1/p)
    c <- (Bk0/p)*c
  }
  p <- drop(mu[M] + v%*%C[1:(M-1),M])
  if (any(p < 0 | p > 1)) {
    stop(mu[M],' ',p,' Infeasible probability. Fix inputs.')
  }
  k[ ,M] <-(rand_U[ ,M] <= p)

  return(k)
}
