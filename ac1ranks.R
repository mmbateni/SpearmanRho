##This founction calculates lag-1 autocorrelation of Ranks of data and
#transform it to the autocorrelation of Normalized data.
ac1ranks <- function(x){
  
  # Calculate ranks
  ranks <- rank(x, ties.method = "average")
  
  # Calculate autocorrelation of ranks at lag 1
  ac_ranks <- acf(ranks, plot = FALSE, type = "correlation", lag.max = 1)$acf[2]
  
  # Transform to autocorrelation of normalized data
  ac1 <- ac_ranks
  
  return(ac1)
}