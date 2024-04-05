
kappa1 <- function(x) {
  k1 <- rep(1/2,length(x))
  tt <- which(abs(x)<=0.0001)
  k1[tt] <- 1/2 + x[tt]/12 - x[tt]^3/720 + x[tt]^5/30240
  tt <- which(abs(x)>0.0001)
  k1[tt] <- 1 - 1/(x[tt]) - 1/(1-exp(x[tt]))
  return(k1)
}

kappa2 <- function(x) {
  k2 <- rep(1/12,length(x))
  tt <- which(abs(x)<=0.015)
  k2[tt] <- 1/12 - x[tt]^2/240 + x[tt]^4/6048
  tt <- which(abs(x)>0.015)
  k2[tt] <- 1/(x[tt])^2 + 1/(2-2*cosh(x[tt]))
  return(k2)
}



fisher.scoring <- function(y,x,initial){
  
  # Initialization
  beta0 <- initial
  beta1 <- beta0
  
  # Main loop
  x <- cbind(rep(1,nrow(x)),x)
  counter <- 0
  repeat{
    counter <- counter+1
    eta <- as.vector(x%*%beta0)
    
    kp <- kappa1(eta)
    kpp <- kappa2(eta)
    kpp[which(kpp<0.01)] <- 0.01
    
    W <- kpp
    Z <- x%*%beta0+(y-kp)/W
    fit <- lm(Z~x-1,weights=W)
    beta1 <- fit$coef
    epsilon <- sqrt(sum((beta0-beta1)^2)/sum(beta0^2))
    print(paste("Epsilon: ", epsilon, sep=""))
    if(epsilon<=0.01) break
    if(counter==100) {print("no convergence"); break}
    beta0 <- beta1
  }
  list(beta.hat=as.vector(beta1),zstat=summary(fit)[["coefficients"]][, "t value"])
}