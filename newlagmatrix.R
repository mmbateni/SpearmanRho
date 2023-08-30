newlagmatrix <- function(x, nlags, c = 0) {
  # Construction of a matrix of lags (X) and a vector (Y) for use in an autoregression
  #
  # USAGE:
  #     newlagmatrix(x, nlags, c)
  #
  # INPUTS:
  #     x     - The dependent variable (Tx1)
  #     nlags - The number of lags (scalar)
  #     c     - [OPTIONAL] 1 if you want to include a constant; 0 is default
  #
  # OUTPUTS:
  #     y     - (n-p) by 1 vector of contemporaneous values
  #     x     - (n-p) by p (or p+1 if c=1) matrix of lags and possibly constant.
  #
  # COMMENTS:
  #     Original name, 'lagmatrix' conflicts with a Matlab file, so newlagmatrix replaces original 
  
  if (length(c) != 1 || any(c != 0 && any(c != 1))) {
    stop('C must be 1 or 0')
  }
  
  if (length(nlags) != 1 || floor(nlags) != nlags || any(nlags < 0)) {
    stop('NLAGS must be a positive integer')
  }
  
  T <- length(x)
  K <- length(x)
  if (min(T, K) != 1) {
    stop('X must be a vector')
  }
  
  if (T < K) {
    x <- t(x) # Transpose if a column vector is passed
  }
  
  if (nlags > 0) {
    nlags <- nlags + 1
    newX <- c(x, rep(0, nlags))
    lagmatrix <- matrix(rep(newX, nlags), nrow = nlags)
    lagmatrix <- lagmatrix[1:(nrow(lagmatrix) - nlags), ]
    y <- lagmatrix[, 1]
    x <- lagmatrix[, 2:nlags]
    if (c == 1) {
      x <- cbind(1, x)
    }
  } else {
    if (c == 1) {
      y <- x
      x <- matrix(1, nrow = T, ncol = 1)
    } else {
      y <- x
      x <- matrix(nrow = T, ncol = 0)
    }
  }
  
  return(list(y = y, x = x))
}
