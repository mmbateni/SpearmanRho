############################################
#### SPEARMAN'S RHO TEST   ##### ############
############################################
#################
### Performs Spearman's rho test of the null hypothesis of trend
### absence in the vector V,  against the alternative of trend. 
### The result of the test is returned in Td = 1 indicates positive trend
### Td = -1 indicates negative trends, i.e., 
### a rejection of the null hypothesis at the alpha significance level. Td = 0 indicates
### a failure to reject the null hypothesis at the alpha significance level.
############################################
### INPUTS
#V = time series [vector]
#alpha =  significance level of the test [scalar]
############################################
#The significance level of a test is a threshold of probability a agreed
#to before the test is conducted. A typical value of alpha is 0.05. If the p-value of a test is less than alpha,
#the test rejects the null hypothesis. If the p-value is greater than alpha, there is insufficient evidence 
#to reject the null hypothesis. 
############################################
### OUTPUTS
#Td = test result: [1] Positive Trend; [-1] Negative Trend;  [0] Insufficient evidence to reject the null hypothesis
#p_value = p-value of the test
############################################
#The p-value of a test is the probability, under the null hypothesis, of obtaining a value
#of the test statistic as extreme or more extreme than the value computed from
#the sample.
##################################################
#### References 
# Daniels, H. E. (1950), Rank correlation and population models, J. R. Stat.
# Soc. Ser. B, 12, 171– 181.
# Khaliq, M. N., T. B. M. J. Ouarda, P. Gachon, L. Sushama, and A. St-Hilaire
# (2009), Identification of hydrological trends in the presence of serial and 
# cross correlations: A review of selected methods and their application to 
# annual flow regimes of Canadian rivers, J. Hydrol., 368, 117 – 130, 
# doi:10.1016/j.jhydrol.2009.01.035

spearman_rho <- function(v, alpha = 0.05, ac_correction = FALSE){
  
  # Rank the data
  v <- v[order(v)]
  ranks <- rank(v)
  
  # Calculate Spearman's rho 
  d <- ranks - seq_along(v)
  rho <- 1 - 6 * sum(d^2) / (length(v) * (length(v)^2 - 1))
  
  # Apply autocorrelation correction if specified
  if(ac_correction){
    r <- acf(v, plot = FALSE)$acf[2]
    var_s_mod <- (1-r^2)^2 / (length(v) - 1)
  } else {
    var_s_mod <- 1 / (length(v) - 1) 
  }
  
  # Calculate the test statistic
  z <- rho / sqrt(var_s_mod)
  
  # Calculate the p-value
  p_value <- 2 * pnorm(abs(z), lower.tail = FALSE) 
  
  # Determine if the trend is significant
  significant <- abs(z) > qnorm(1 - alpha/2)
  trend <- ifelse(significant, sign(rho), 0)
  
  # Print output
  cat(sprintf("p-value: %.3f\n", p_value))
  if(trend == 1){
    cat(sprintf("Positive trend detected at %.1f%% significance\n", (1 - alpha)*100))
  } else if (trend == -1){
    cat(sprintf("Negative trend detected at %.1f%% significance\n", (1 - alpha)*100))
  } else {
    cat("No significant trend detected\n")
  }
  
  return(list(trend = trend, p_value = p_value))
  
}