init.theta=function(train,cl,k,p)
{ 
  n.train <- nrow(train)
  levels <- unique(cl)
  skew <- matrix(0, k, p)
  for (i in 1:k) {
    skew[i, ] <- apply(train[cl == levels[i], , drop = FALSE], 
                       2, skewness)
    skew[i, ] <- sign(skew[i, ]) * (abs(skew[i, ])^(1/3))
  }
  total.skew <- colSums(skew)
  total.skew <- ifelse(is.na(total.skew), 1, total.skew)
  
  theta=1- abs(total.skew - min(total.skew) - 0.01 )/(max(total.skew) - min(total.skew) )
  return(theta)
}