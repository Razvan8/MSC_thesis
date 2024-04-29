### FUNCTIONS THAT I USE IN WEAK_HIERNET



## ASSERT FUNCTION
assert <- function(condition, message) {
  if (!condition) stop(message)
}


###USEFUL FUNCTION FOR LASSO LOSS AND WEAK HIER NET

set_Z_values_to_zero <- function(matrix, p) {
  for (i in 1:p) {
    # Calculate the column index corresponding to i,i
    col_index <- (i - 1) * p + i
    
    # Set values to 0 for the selected column
    matrix[, col_index] <- 0
  }
  
  # Return the modified matrix
  return(matrix)
}


get_positions <- function(p, j) {
  # Generate all pairs
  all_pairs <- expand.grid(1:p, 1:p)
  #print(all_pairs)
  
  # Find positions where THE second column is equal to j
  positions <- which( all_pairs[, 2] == j)
  
  # Return the positions
  return(positions)
}

##HIERARCHICAL LASSO LOSS
hierarchical_lasso_loss <- function(X, y, Beta_plus, Beta_minus, beta0, Theta, lambda) {
  # Check if X has the correct dimensions
  assert <- function(condition, message) {
    if (!condition) stop(message)
  }
  
  #DIMENSIONS: X= nxp, Beta= px1, Theta= pxp, Z= nx p^2, with 0 in 0, p, 2p...p2 positions
  
  assert(ncol(X) > 1, "X must have at least 2 columns for main effects and interactions.")
  assert(length(Beta_plus) == ncol(X), "Beta must have length equal to the number of columns in X.")
  assert(length(Beta_minus) == ncol(X), "Beta must have length equal to the number of columns in X.")
  #assert(dim(Theta) == c(ncol(X), ncol(X)), "Theta must have dimensions ncol(X) x ncol(X).")
  Beta0 <- matrix(beta0, nrow = ncol(X), ncol = 1)
  vec_theta =  c(t(Theta)) 
  Z=t(apply(X, 1, function(x) outer(x, x, "*")))
  
  
  
  
  
  # Check if the diagonal elements are zero
  assert(all(diag(Theta) == 0), "Interaction coefficients on the diagonal must be zero.")
  
  
  
  # Calculate the least squares loss
  loss <- sum( (y - Beta0-X %*% (Beta_plus- Beta_minus)- Z%*%vec_theta)^2/2) /  length(y) 
  print(loss)
  
  # Add L1 norm penalties for main effects and interactions
  loss <- loss  + lambda * ( sum(abs(Beta_plus)) +  sum(abs(Beta_minus)) +sum(abs(Theta))/2 )
  print(loss)
  
  # Add small L2 penalty
  epsilon <- 1e-8 * lambda
  loss <- loss + epsilon / 2 * ( sum(Beta_plus^2) + sum(Beta_minus^2) + sum(vec_theta^2) )
  
  
  return(loss)
}



##Example test HIER LASSO LOSS
#X=matrix(c(2,2,1,1), nrow = 2, ncol = 2, byrow = TRUE)
#y=matrix(c(8,4), nrow = 2, ncol = 1)
#Beta_plus=matrix(c(3,3), nrow = 2, ncol = 1)
#Beta_minus=matrix(c(1,1), nrow = 2, ncol = 1)
#Theta=matrix(c(0,1,1,0), nrow=2, ncol = 2)
#beta0=0
#lambda=1
#hierarchical_lasso_loss(X, y, Beta_plus, Beta_minus, beta0, Theta, lambda)


###SOFT THRESHOLDING OPERATOR
Soft_thresholding <- function(c, lambda) {
  assert <- function(condition, message) {
    if (!condition) stop(message)
  }
  
  assert(lambda >= 0, "lambda cannot be negative.")
  
  # Apply soft thresholding component-wise
  result <- sign(c) * pmax(abs(c) - c(lambda), 0)
  
  return(result)
}




## RELU function
RELU <- function(x) {
  if (x>=0){
    return(x)}
  else {
    return (0)
  }
}



### Find knots - HELPING FUNCTION FOR ONE ROW -b) , c) from alg
evaluate_knots <- function(Theta_tilda_j, beta_tilda_plus_j, beta_tilda_minus_j, lambda, j, f, t = 1) {
  p <- length(Theta_tilda_j)
  
  set_elements <- numeric(0)
  
  for (k in 1:p)  { 
    if (k!=j){ ##discard theta_jj
      element <- abs(Theta_tilda_j[k]) / t - lambda / 2
      set_elements <- c(set_elements, element)}
  }
  
  set_elements <- c(set_elements, -beta_tilda_plus_j/t, -beta_tilda_minus_j/t)
  
  selected_elements <- set_elements[set_elements >= 0]
  
  selected_elements<-c(selected_elements,0) #### ADD O TO THE LIST ### IT SHOULD NOT INFLUENCE (TAKE CARE HERE!!!!!!!!!!)
  selected_elements<- sort(selected_elements) ## sort hem increasing
  #print(selected_elements)
  
  result_vector <- sapply(selected_elements, f)## Evaluate f(p) for every p
  #cat("selected elements : ", selected_elements)
  #print("  ")
  #cat("result vector : ", result_vector)
  
  
  return(list("knots" = selected_elements,"evaluated"=result_vector))
}





# Function to find adjacent knots and calculate alpha_hat - HELPING FUNCTION FOR ONEROW
find_adjacent_knots_and_alpha <- function(knots, evaluated) {
  # Find the index where evaluated > 0
  positive_index <- max(which(evaluated > 0))
  #print(which(evaluated > 0))
  #cat("positive index", positive_index)
  
  # Ensure there is a positive index and the next index exists
  if (!is.na(positive_index) && (positive_index + 1) <= length(knots)) {
    # Get corresponding knots
    p1 <- knots[positive_index]
    p2 <- knots[positive_index + 1]
    
    # Calculate alpha_hat
    #alpha_hat <- -evaluated[positive_index] * (evaluated[positive_index + 1] - evaluated[positive_index]) / (p2 - p1) ##initial one
    alpha_hat <- (p1*evaluated[positive_index + 1] - p2*evaluated[positive_index]) / (evaluated[positive_index + 1] - evaluated[positive_index])
    #print(positive_index)
    return(alpha_hat)
  } else {
    # Return an informative message if no such adjacent knots are found
    stop("No adjacent knots found with evaluated > 0 and evaluated < 0 at consecutive positions.")
  }
}




### FINAL RETURN - HELPING FUNCTION FOR ONEROW
final_return<- function(beta_tilda_plus_j, beta_tilda_minus_j, Theta_tilda_j, lambda,j, t, alpha_hat)
  
{ beta_hat_plus_j = RELU(beta_tilda_plus_j +t * alpha_hat -t*lambda ) ###bc -t*lambda has to be in beta tilda !!!
beta_hat_minus_j = RELU(beta_tilda_minus_j +t * alpha_hat -t*lambda) ##
#cat("In final return alpha hat, beta tilda plus j, t, beta hat plus", alpha_hat,beta_tilda_plus_j,t, beta_hat_plus_j)

Theta_hat_j = Soft_thresholding(Theta_tilda_j, t*(lambda/2 + alpha_hat))
#print(Theta_hat_j)
if (j==1){ 
  #print("1")
  Theta_hat_j <- c( 0, Theta_hat_j)#### THIS IS TO ADD ) BACK AT POSITION J
  return (list("beta_hat_plus_j" = beta_hat_plus_j, "beta_hat_minus_j" = beta_hat_minus_j , "Theta_hat_j" = Theta_hat_j)) }
if (j ==length(Theta_tilda_j)+1)
{ #print("len max")
  Theta_hat_j <- c(Theta_hat_j, 0)   #### THIS IS TO ADD ) BACK AT POSITION J
  return (list("beta_hat_plus_j" = beta_hat_plus_j, "beta_hat_minus_j" = beta_hat_minus_j , "Theta_hat_j" = Theta_hat_j)) }
#print("final ret")
Theta_hat_j <- c(Theta_hat_j[1:(j-1)], 0, Theta_hat_j[j:length(Theta_hat_j)])#### THIS IS TO ADD ) BACK AT POSITION J
return (list("beta_hat_plus_j" = beta_hat_plus_j, "beta_hat_minus_j" = beta_hat_minus_j , "Theta_hat_j" = Theta_hat_j))
}





# function for ONEROW
ONEROW <- function(beta_tilda_plus_j, beta_tilda_minus_j, Theta_tilda_j, lambda,j,  t=1) {
  
  f <- function(alpha) {
    #cat(' suma ',sum(abs(Soft_thresholding( Theta_tilda_j, t*(lambda/2+alpha)  ) ) ))
    return( sum(abs(Soft_thresholding( Theta_tilda_j, t*(lambda/2+alpha)  ) ) ) - RELU(beta_tilda_plus_j+t*alpha) - RELU(beta_tilda_minus_j + t* alpha)  )### 1
  }
  
  if ( f(0)<=0)### 1 a)
  {#print('a')
    alpha_hat<-0
    return(final_return(beta_tilda_plus_j = beta_tilda_plus_j, beta_tilda_minus_j = beta_tilda_minus_j, Theta_tilda_j =  Theta_tilda_j[-j],
                        lambda =  lambda, j=j, t=t, alpha_hat = alpha_hat) ) ### 2)
  }   
  
  
  ### 1b) ### 1c)
  knots_evaluated <- evaluate_knots( beta_tilda_plus_j= beta_tilda_plus_j,beta_tilda_minus_j= beta_tilda_minus_j, Theta_tilda_j=Theta_tilda_j,
                                     lambda=lambda,j=j, f=f, t=t)
  ### 1d)
  knots=knots_evaluated$knots
  evaluated=knots_evaluated$evaluated
  #print( evaluated)
  
  zero_indices <- which(evaluated == 0)
  
  if (length(zero_indices) > 0) {
    #print('d')
    alpha_hat <- knots[zero_indices[1]] # i.e. alpha_hat =p s.t. f(p)=0
    cat("alpha hat",alpha_hat)
    return (final_return(beta_tilda_plus_j = beta_tilda_plus_j, beta_tilda_minus_j = beta_tilda_minus_j, Theta_tilda_j =  Theta_tilda_j[-j],
                         lambda =  lambda, j=j, t=t, alpha_hat = alpha_hat) )}   ### 2)
  
  
  alpha_hat <- find_adjacent_knots_and_alpha(knots, evaluated)    ### 1e)
  #print('e')
  
  ### STEP 2 
  return (final_return(beta_tilda_plus_j = beta_tilda_plus_j, beta_tilda_minus_j = beta_tilda_minus_j, Theta_tilda_j =  Theta_tilda_j[-j],
                       lambda =  lambda, j=j, t=t, alpha_hat = alpha_hat)) ### 2)
  
}


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





######################### WEAKHIERNET CLASS #########################################

# Define the WeakHierNet class
WeakHierNetUnscaled <- function(X, Beta_plus_init, Beta_minus_init, Theta_init, y, lambda, t=1, tol=1e-2, max_iter=5000, eps=1e-8, Z=NULL) {
  
  
 
    self = list()
    self$Beta_hat_plus <- Beta_plus_init
    self$Beta_hat_minus <- Beta_minus_init
    self$Theta_hat <- Theta_init

  
  
  
  
  # Public method for fitting the model
  fit <- function( X, Beta_plus_init, Beta_minus_init, Theta_init, y, lambda, t=1, tol=1e-2, max_iter=5000, eps=1e-8, Z=NULL) {
    
    eps<-1e-8*lambda
    p <- ncol(X)
    n<-nrow(X)
    #CONSTRUCT Z
    if (is.null(Z))
    {Z=t(apply(X, 1, function(x) outer(x, x, "*")))
    ####### MAKE ALL COLS 
    # Count the occurrences of -1 in each column
    minus_ones_count <- colSums(Z == -1)
    
    # Identify columns with no occurrences of -1
    no_minus_ones_columns <- which(minus_ones_count == 0)
    
    # Set columns with no -1 to 0
    Z[, no_minus_ones_columns] <- 0
    }
    ####
    
    # Initialize variables
    delta <- 1-t*eps
    Beta_hat_plus <- Beta_plus_init
    Beta_hat_minus <- Beta_minus_init
    Theta_hat <- Theta_init
    vec_Theta_hat <-  c(t(Theta_hat)) 
    r_hat <- matrix(-1, nrow = n, ncol = 1)
    
    
    for (k in 2:max_iter+1) {
      vec_Theta_hat <-  c(t(Theta_hat)) 
      
      r_hat_old <- r_hat
      #print(length(y))
      #print(dim(X))
      #print(length(Beta_hat_minus))
      #print(length(Beta_hat_plus))
      r_hat <-y- X %*% (Beta_hat_plus - Beta_hat_minus) - Z %*% vec_Theta_hat / 2
      r_hat<- -  r_hat/2 ############################### TAKE CARE ! why like this? why/2 ########################
      if (k%%30==3)
      {cat(' Loss rhat',mean(r_hat^2))}
      #cat("Theta", Theta_hat) 
      
      
      
      
      for (j in 1:p) { #Take care dimensions match
        j_cols= get_positions(p,j)

        #cat('mean(r_hat): ', mean(r_hat), " ; ")
        
        ###CHECK HERE!!!
        result_ONEROW <-  ONEROW(delta * Beta_hat_plus[j] - t* t(X[, j]) %*% r_hat , 
                                 delta * Beta_hat_minus[j] + t * t(X[, j]) %*%r_hat ,
                                 delta * Theta_hat[,j] - t* t(Z[, j_cols])%*%r_hat, lambda=lambda ,j=j,  t=t)
        
        Beta_hat_plus[j] <- result_ONEROW$beta_hat_plus_j
        
        Beta_hat_minus[j] <- result_ONEROW$beta_hat_minus_j
        
        Theta_hat[,j] <-  result_ONEROW$Theta_hat_j
      }
      if (k>=3)
      {if (mean((r_hat_old- r_hat)^2) <=tol) #TAKE CARE RHAT IS VECTOR  !!! Sum because already scaled !!!!!!
      {  cat("Converged at iteration ",k)
        self$Beta_hat_plus=Beta_hat_plus
        self$Beta_hat_minus = Beta_hat_minus
        self$Theta_hat=Theta_hat
        return(self) }
      }
      
    }
    self$Beta_hat_plus=Beta_hat_plus
    self$Beta_hat_minus = Beta_hat_minus
    self$Theta_hat=Theta_hat
    cat("It has not converged. The difference between last 2 residuals is:", abs(r_hat[k-1]- r_hat[k-2]))
    return(self) 
  }
  
  
  # method for predicting
  
  predict <- function(self, new_X, Z=NULL) {
    if (!is.numeric(new_X) || !is.matrix(new_X)) {
      stop("Input 'new_X' must be a numeric matrix.")
    }
    
    if (!is.numeric(self$Beta_hat_plus) || !is.numeric(self$Beta_hat_minus)) {
      stop("Coefficients 'Beta_hat_plus' and 'Beta_hat_minus' must be numeric.")
    }
    ## Later add also rescaling
    # Implement prediction code here
    p <- ncol(new_X)
    n<-nrow(new_X)
    #CONSTRUCT Z
    if (is.null(Z))
    {Z=t(apply(X, 1, function(x) outer(x, x, "*")))
      ####### MAKE ALL COLS  0 where same main effect (i.e have 0 of -1)
      # Count the occurrences of -1 in each column
      minus_ones_count <- colSums(Z == -1)
      
      # Identify columns with no occurrences of -1
      no_minus_ones_columns <- which(minus_ones_count == 0)
      
      # Set columns with no -1 to 0
      Z[, no_minus_ones_columns] <- 0
    }

    ####
    y_pred= X %*% (self$Beta_hat_plus - self$Beta_hat_minus) + Z %*% c(t(self$Theta_hat)) / 2
    return(y_pred)
  }
  
  R2_score <- function(self, new_X, y_true, verbose= TRUE, Z=NULL) {
    
    
    y_pred = predict(self,new_X, Z)
    
    if (verbose == TRUE)
    {cat ("r2 score is ", r2(y_true, y_pred))
      plot(y_pred, y_true)}
    
    return(r2(y_true, y_pred))
  }
  

  return(list( fit = fit, predict = predict, R2_score = R2_score, self = self))
}









library(Metrics)

library(hierNet)
library(caret)
library(dplyr)
library(Metrics)








