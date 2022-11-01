#' Simulate genotype/phenotype data under equilibrium univariate AM.
#'
#' @param h2_0 generation zero (panmictic) heritability
#' @param r cross-mate phenotypic correlation
#' @param m number of biallelic causal variants
#' @param n sample size
#' @param min_MAF minimum minor allele frequency for causal variants
#'
#' @return A list including the following objects:
#' * `y`: phenotype vector
#' * `g`: heritable component of the phenotype vector
#' * `X`: matrix of diploid genotypes
#' * `AF`: vector of allele frequences
#' * `beta_std`: standardized genetic effects
#' * `beta_raw`: unstandardized genetic effects
#' @export
#'
#' @examples
#' set.seed(1)
#' h2_0 = .5; m = 200; n = 1000; r =.5
#'
#' ## simulate genotype/phenotype data
#' sim_dat <- am_simulate(h2_0, r, m, n)
#' str(sim_dat)
#'
#' ## empirical h2 vs expected equilibrium h2
#' (emp_h2 <- var(sim_dat$g)/var(sim_dat$y))
#' h2_eq(r, h2_0)

am_simulate <- function(h2_0, r, m, n, min_MAF=.1) {
  ## draw standardized diploid allele substitution effects
  beta <- scale(rnorm(m))*sqrt(h2_0 / m)
  ## draw allele frequencies
  AF <- runif(m, min_MAF, 1 - min_MAF)
  ## compute unstandardized effects
  beta_unscaled <- beta/sqrt(2*AF*(1-AF))
  ## generate corresponding haploid quantities
  beta_hap <- rep(beta, each=2)
  AF_hap <- rep(AF, each=2)
  ## compute equilibrium outer product covariance component
  U <- am_covariance_structure(beta, AF, r)
  ## draw multivariate Bernoulli haplotypes
  H <- rb_dplr(n, AF_hap, U)
  ## convert haplotypes to diploid genotypes
  X = (H[,seq(1,2*m,2)]+H[,seq(2,2*m,2)])
  ## compute genetic phenotypes
  g = X %*% beta_unscaled
  ## compute full phenotype
  y = g + rnorm(n, 0, sqrt(1 - h2_0))
  return(list(
    y = y,
    g = g,
    X = X,
    AF = AF,
    beta_std = beta,
    beta_raw = beta_unscaled
    ))
}

