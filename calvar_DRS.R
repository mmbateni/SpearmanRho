calvar_DRS <- function(n, r){
  
  t <- -n:n
  
  cf <- function(t, r) r^abs(t)
  cr <- cf(t, r)
  
  cr2 <- (2 - 2*cr)^-0.5
  
  ee <- 0
  
  for(i in 1:n){
    for(l in 1:n){
      if(l == i) next
      
      for(u in 1:n){
        if(l == u) next
        
        e1 <- cr[i - l + n + 1]
        e2 <- cr[i - u + n + 1]
        d1 <- cr2[l - u + n + 1]
        
        for(j in 1:n){
          if(i == j) next
          
          e3 <- cr[j - u + n + 1]
          e4 <- cr[j - l + n + 1]
          d2 <- cr2[j - i + n + 1]
          
          e <- (e3 - e2 - e4 + e1) * d2 * d1
          
          if(abs(e) > 1) e <- sign(e)
          
          ee <- ee + asin(e) * (i-1) * (l-1)
        }
      }
    }
  }
  
  v <- ee/(2*pi) * 144/(n*(n^2-1))^2
  
  return(v)
}