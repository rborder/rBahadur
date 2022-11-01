#' Compute Diagonal plus Low Rank covariance structure under equilibrium assortative mating
#' #'
#' @importFrom stats runif rnorm
#'
#' @param beta vector of standardized diploid allele-substitution effects
#' @param AF vector of allele frequencies
#' @param r cross-mate phenotypic correlation
#'
#' @return Vector 'U' such that $D + U U^T$ corresponds to the expected haploid LD-matrix given the specified genetic architecture (encoded by 'beta' and 'AF') and cross-mate phenotypic correlation 'r'. It is assumed that the total phenotypic variance at generation zero is one.
#' @examples
#' set.seed(1)
#' h2_0 = .5; m = 200; n = 1000; r =.5; min_MAF=.1
#' betas <- rnorm(m,0,sqrt(h2_0/m))
#' afs <- runif(m, min_MAF, 1-min_MAF)
#' output <- am_covariance_structure(betas, afs, r)
#' @export
am_covariance_structure <- function(beta, AF, r) {
  ## obtain haploid substitution effects, variances
  beta_hap <- rep(beta,each=2)
  sd_hap <- rep(sqrt(AF*(1-AF)), each=2)
  ## compute equilibrium variance components
  h2_0 <- sum(beta**2)
  vgeq = vg_eq(r=r, h2_0, h2_0)
  rgeq = rg_eq(r=r, h2_0)
  vtot = vgeq+(1-h2_0)
  ## compute outer product component
  U <- sqrt(vtot/2)/(2*beta_hap*sqrt(r)) * (sqrt(4*beta_hap**2 *r/vtot + (1-rgeq)^2)-(1-rgeq) )* sd_hap
  return(U)
}

