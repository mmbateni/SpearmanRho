ac1ranks <- function(x) {
  # Calculate ranks
  ranks <- rank(x, ties.method = "average")
  
  # Calculate lag-1 autocorrelation of ranks
  ac_ranks <- acf(ranks, plot = FALSE, type = "correlation", lag.max = 1)$acf[2]
  
  # Transform rank autocorrelation to Pearson autocorrelation
  # Standard conversion under bivariate normality: r = 2 * sin(pi/6 * rho_s)
  ac1 <- 2 * sin((pi / 6) * ac_ranks)
  
  return(ac1)
}
