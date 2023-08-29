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
### From Matlab Help ###############
#The significance level of a test is a threshold of probability a agreed
#to before the test is conducted. A typical value of alpha is 0.05. If the p-value of a test is less than alpha,
#the test rejects the null hypothesis. If the p-value is greater than alpha, there is insufficient evidence 
#to reject the null hypothesis. 
############################################
### OUTPUTS
#Td = test result: [1] Positive Trend; [-1] Negative Trend;  [0] Insufficient evidence to reject the null hypothesis
#p_value = p-value of the test
############################################
### From Matlab Help ###############
#The p-value of a test is the probability, under the null hypothesis, of obtaining a value
#of the test statistic as extreme or more extreme than the value computed from
#the sample.
##################################################
#### References 
#Daniels, H. E. (1950), Rank correlation and population models, J. R. Stat.
# Soc. Ser. B, 12, 171– 181.
#Khaliq, M. N., T. B. M. J. Ouarda, P. Gachon, L. Sushama, and A. St-Hilaire
#(2009), Identification of hydrological trends in the presence of serial and 
#cross correlations: A review of selected methods and their application to 
#annual flow regimes of Canadian rivers, J. Hydrol., 368, 117 – 130, 
#doi:10.1016/j.jhydrol.2009.01.035
##################################################
spearman_rho_modified <- function(v, alpha = 0.05, ac_correction = FALSE){
  
  # Rank the data
  v <- v[order(v)]
  ranks <- rank(v)
  
  n <- length(v)
  
  # Calculate Spearman's rho
  d <- ranks - seq_along(v)
  rho <- 1 - 6 * sum(d^2) / (n^3 - n)
  
  # Autocorrelation correction
  if(ac_correction){
    
    # Estimate autocorrelation 
    acf_v <- acf(v, plot = FALSE)$acf[2]
    
    # Calculate correction factor
    t <- -n:n
    cf <- function(t, acf_v) acf_v^abs(t) 
    cr <- cf(t, acf_v)
    cr2 <- (2 - 2*cr)^-0.5
    
    # Estimate variance 
    ee <- 0
    for(i in 1:n){
      for(l in 1:n){
        if(l == i) next
        for(u in 1:n){
          if(l == u) next
          e1 <- cr[(i-l)+n+1]
          e2 <- cr[(i-u)+n+1]
          d1 <- cr2[(l-u)+n+1]
          
          for(j in 1:n){
            if(i == j) next
            e3 <- cr[(j-u)+n+1]
            e4 <- cr[(j-l)+n+1]
            d2 <- cr2[(j-i)+n+1]
            
            ee <- ee + asin((e3 - e2 - e4 + e1)*d2*d1) * (i-1) * (l-1)
          }
        }  
      }
    }
    
    var_s_mod <- ee/(2*pi) * 144/(n*(n^2-1))^2
    
  } else {
    var_s_mod <- 1/(n - 1)
  }
  
  # Calculate test statistic
  z <- rho / sqrt(var_s_mod) 
  
  # p-value
  p <- 2 * pnorm(abs(z), lower.tail = FALSE) 
  
  # Determine significant trend
  trend <- ifelse(abs(z) > qnorm(1 - alpha/2), sign(rho), 0)
  
  return(list(trend = trend, p_value = p))
  
}