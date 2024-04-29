
my_k_prime<- function(x)
{ if (abs(x)>=0.0001){
  result<- 1 / (1- exp(-x)) - 1/x}
  else
  {result<-1/2 + x/12 - x^3/720 +x^5/30240}
return(result)}


my_k_second <- function(x, stable=TRUE)
{ if (abs(x) >=0.015){
  result<- 1/x^2 - exp(x)/ ( (1-exp(x))^2 )}
  else 
  {result<- 1/12 - x^2/240 + x^4/6048}
  if (stable==TRUE)
  {if (result<0.01)
    {result<-0.01}}
return (result)}


my_f<-function(X, beta_old, k_second)
{X_beta <- X %*% beta_old


# Apply the function ??''() element-wise to X_beta
k_second_X_beta <- sapply(X_beta, k_second)

# Construct a diagonal matrix with the resulting elements
diag_matrix <- diag(k_second_X_beta)

# Multiply X by the diagonal matrix
result <- t(X) %*% diag_matrix %*% X
return(result)  
}

my_s<-function(X, beta_old, y, k_prime)
  {Ey <-  sapply(X%*%beta_old, k_prime) 
  Ey<-matrix(Ey, nrow = nrow(X), ncol = 1)
  result<- t(X) %*% (y-Ey)
  
  return(result)
  }





weighted_ls_objective <- function(beta, Z, X, W) {
  residuals <- Z - X %*% beta
  return(t(residuals) %*% W %*% residuals)
}


#loss gradient
gradient <- function(beta, Z, X, W) {
  return(-2 * t(X) %*% W %*% (Z - X %*% beta))
}



# Function to compute penalty ###IT is automatically symemtric, so it is strong, double theta to have weak!!!!!!!
penalty <- function(beta, interactions, lambda) {
  total_penalty <- 0
  #print(length(interactions))
  i=0
  for (interaction in interactions) {
    # Compute penalty for current interaction
    #print(interaction)
    penalty_term <- max(abs(beta[interaction[[1]] ]) , sum(abs(beta[interaction[2:length(interaction)] ])) )
    #print(beta[interaction[-1] ])
    # Add lambda times penalty term to total penalty
    total_penalty <- total_penalty + lambda * penalty_term
    #print('total penalty now:')
    #print(total_penalty)
  }
  #print(total_penalty)
  total_penalty<- total_penalty +lambda* sum(abs(beta[41:length(beta)])) ##lambda/2 if weak ##TAKE CARE 41 hardcoded
  return(total_penalty)
}


beta=rep(0,515)
beta[1:40]=1
beta[41]=3
beta[42]=2
#lambda=1
#penalty(beta, list_of_interactions,lambda)

add_gradient_max_abs_l1 <- function(b, indices, gradient, lambda) {
  abs_values <- abs(b[indices])
  max_abs <- max(abs_values[1],sum(abs_values[2:length(abs_values)]) )
  #print(length(abs_values))
  #print(max_abs)
  #print(abs_values[1])
  #print(sum(abs_values[2:length(abs_values)]))
  #print('--------')
  if (abs_values[1] >= max_abs)
  {gradient[indices[1]]<-gradient[indices[1]] + lambda* sign(b[indices[1]])}
  if (sum(abs_values[2:length(abs_values)]) >= max_abs)
  {gradient[indices[2:length(indices)]]<-gradient[indices[2:length(indices)]] + lambda*sign(b[indices[2:length(indices)]])}
  
  return(gradient)
}

penalty_gradient <- function(beta, interactions, lambda) {
  # Initialize gradient vector
  gradient <- rep(0, length(beta))
  
  # Compute gradient for each interaction
  for (interaction in interactions) {
    gradient<-add_gradient_max_abs_l1(beta, interaction, gradient, lambda )
    }
  #print(gradient)
  
  # Compute gradient for the rest of beta (beyond interaction indexes)
  gradient[41:length(beta)] <- gradient[41:length(beta)] + lambda * sign(beta[41:length(beta)]) #lambda instead of l/2 because strong
  
  return(gradient)
}

total_gradient<-function(beta, interactions, lambda, Z,W,X){
  result<-gradient(beta, Z, X, W) + penalty_gradient(beta, interactions, lambda)
  return(result)}

total_loss<-function(beta, Z, X, W, interactions, lambda)
{ result<-weighted_ls_objective(beta, Z, X, W) + penalty(beta, interactions, lambda)
return(result)}

#beta=rep(1e-6,515)
#beta[1:40]=1
#beta[41]=3
#beta[42]=1
#lambda=1
#penalty_gradient(beta, list_of_interactions,lambda)




#### add the gradient of penalty in the algorithm. Consider second derivatives as 0 ?!
Fisher_scoring_alg<-function(beta_init, X, y, interactions, f=my_f, s=my_s, max_iter=1e2, tol=1e-2, k_second=my_k_second, k_prime=my_k_prime,lambda=1, verbose= TRUE)
{ 
  X <- cbind(rep(1,nrow(X)),X)
  beta_new<-beta_init
  for (k in c(1:max_iter)) {
    
  beta_old<-beta_new
  beta_new<-beta_old + solve(f(X, beta_old, k_second)) %*% (s(X,beta_old, y, k_prime)-c(0,penalty_gradient(beta=beta_old[-1], interactions=interactions,lambda=lambda)))
  #print(solve(f(X, beta_old, k_second)) %*% s(X,beta_old, y, k_prime))
  if (verbose == TRUE)
  {if (k %%2==1)
  {print(beta_old[1:5]-beta_new[1:5])
    }}
  
  #print(max(abs(beta_old-beta_new)))
  #print(mean(abs(beta_old-beta_new)))
  print(sqrt(sum((beta_old-beta_new)^2)/sum(beta_old^2)))
  if (sqrt(sum((beta_old-beta_new)^2)/sum(beta_old^2)) <= tol )
    {return (beta_new)}
  } 
print("It did not converge. Last residual is:")
print(abs(beta_new-beta_old))
return(beta_new)
  
}
#####


fisher_scoring <- function(y,X,beta_init, tol, max_iter){
  
  # Initialization
  beta0 <- beta_init
  beta1 <- beta0
  
  # Main loop
  X <- cbind(rep(1,nrow(X)),X)
  counter <- 0
  repeat{
    counter <- counter+1
    eta <- as.vector(X%*%beta0)
    
    kp <- sapply(eta, my_k_prime)
    kpp <- sapply(eta, my_k_second)
    kpp[which(kpp<0.01)] <- 0.01#### why like this? for stability?
    
    W <- kpp
    Z <- X%*%beta0+(y-kp)/W
    fit <- lm(Z~X-1,weights=W)
    beta1 <- fit$coef
    epsilon <- sqrt(sum((beta0-beta1)^2)/sum(beta0^2))
    print(paste("Epsilon: ", epsilon, sep=""))
    if(epsilon<=tol) break
    if(counter==max_iter) {print("no convergence"); break}
    beta0 <- beta1
  }
  list(beta.hat=as.vector(beta1),zstat=summary(fit)[["coefficients"]][, "t value"])
}





fisher_scoring_lasso <- function(y,X,beta_init, tol, max_iter, interactions, lambda){
  
  # Initialization
  beta0 <- beta_init
  beta1 <- beta0
  
  # Main loop
  X <- cbind(rep(1,nrow(X)),X)
  counter <- 0
  repeat{
    counter <- counter+1
    eta <- as.vector(X%*%beta0)
    
    kp <- sapply(eta, my_k_prime)
    kpp <- sapply(eta, my_k_second)
    kpp[which(kpp<0.01)] <- 0.01#### why like this? for stability?
    
    W <- diag(kpp)
    Z <- X%*%beta0+(y-kp)/kpp
    #fit <- lm(Z~X-1,weights=W)
    #beta1 <- fit$coef
    
    beta1 <- optim(par = beta0, fn = total_loss, gr = total_gradient, Z = Z, X = X, W = W, interactions=interactions, lambda=lambda, method ='BFGS',control = list(maxit = 50))$par
    
    epsilon <- sqrt(sum((beta0-beta1)^2)/sum(beta0^2))
    print(paste("Epsilon: ", epsilon, sep=""))
    if(epsilon<=tol) break
    if(counter==max_iter) {print("no convergence"); break}
    beta0 <- beta1
  }
  list(beta.hat=as.vector(beta1))
}






pred<-function(X,beta)
{return(X%*%beta)}

r2 <- function(actual, predicted) {
  # Calculate the mean of the actual values
  mean_actual <- mean(actual)
  
  # Calculate the total sum of squares
  total_sum_squares <- sum((actual - mean_actual)^2)
  
  # Calculate the residual sum of squares
  residual_sum_squares <- sum((actual - predicted)^2)
  
  # Calculate R-squared
  r_squared <- 1 - (residual_sum_squares / total_sum_squares)
  
  return(r_squared)
}



# Set seed for reproducibility
#set.seed(123)

# Number of samples
#n <- 1000

# Generate X matrix with 3 columns
#X <- matrix(runif(n * 3), ncol = 3)


# Generate beta vector with 3 coefficients
#beta <- c(0.5, -0.3, 0.8)

# Generate noise for X and Y
#noise <- rnorm(n, mean = 0, sd = 0.5)

# Calculate X * beta
#Y <- X %*% beta + noise


# Create dataset
#dataset <- data.frame(X1 = X[,1], X2 = X[,2], X3 = X[,3], Y = Y)

# Display the first few rows of the dataset
#head(dataset)


#beta_result<-Fisher_scoring_alg(beta_init=c(0.5,-0.3,0.8), X=X, y=Y, f=my_f, s=my_s, max_iter=1e4, tol=1e-6, k_second=my_k_second, k_prime=my_k_prime, verbose= TRUE)
#beta_result


