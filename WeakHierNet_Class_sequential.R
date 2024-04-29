### FUNCTIONS THAT I USE IN WEAK_HIERNET



## ASSERT FUNCTION
assert <- function(condition, message) {
  if (!condition) stop(message)
}


###USEFUL FUNCTION FOR LASSO LOSS AND WEAK HIER NET



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
hierarchical_lasso_seq_loss <- function(X, y,psi, lambda) {
  # Check if X has the correct dimensions
  assert <- function(condition, message) {
    if (!condition) stop(message)
  }
  
  #DIMENSIONS: X= nxp, Beta= px1, Theta= pxp, Z= nx p^2, with 0 in 0, p, 2p...p2 positions
  
  assert(ncol(X) > 1, "X must have at least 2 columns for main effects and interactions.")
  
  
  
  # Calculate the least squares loss
  loss <- sum( (y - X %*% psi)^2/2) 
  print(loss)
  
  # Add L1 norm penalties for main effects and interactions
  loss <- loss  + lambda *  sum(abs(psi)) 
  print(loss)
  
  # Add small L2 penalty
  epsilon <- 1e-8 * lambda
  loss <- loss + epsilon / 2 *  sum(psi^2)
  return(loss)
}



#Example test HIER LASSO LOSS
X=matrix(c(2,2,1,1), nrow = 2, ncol = 2, byrow = TRUE)
print(X)
y=matrix(c(9,4), nrow = 2, ncol = 1)
psi=matrix(c(1,3), nrow = 2, ncol = 1)
lambda=1
hierarchical_lasso_seq_loss(X, y, psi, lambda)


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
evaluate_knots_seq <- function(psi_tilda_jk, lambda, f, t = 1,c=1) {
  p <- length(psi_tilda_jk)
  
  set_elements <- numeric(0)
  
  for (k in 1:p)  { 
    if (abs(psi_tilda_jk[k])>0){ ##discard impossible combinations
      element <- abs(psi_tilda_jk[k]) / t - lambda / 6 #element to add in knot set
      set_elements <- c(set_elements, element)}
  }
  
  selected_elements <- set_elements[set_elements >= 0] ### e >= dar egal e deja verificat?!
  
  #selected_elements<-c(selected_elements,0) #### ADD O TO THE LIST ### IT SHOULD NOT INFLUENCE (TAKE CARE HERE!!!!!!!!!!)
  selected_elements<- sort(selected_elements) ## sort them increasing
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

final_return_seq<- function(psi_hat_jk, lambda, j, t, alpha_hat) # returns psi_jk with all the elements inside (vector with 40 elements)
{ 
psi_hat_jk = Soft_thresholding(psi_hat_jk, t*(lambda/6 + alpha_hat))
return (psi_hat_jk)
}

##RUN TESTS ###############################################################################################################################

get_all_possible_kj<-function(l1=21,l2=14,l3=3,l4=2)
 {range1<-c(1:l1)
  range2<-c((l1+1):(l1+l2))
  range3<-c((l1+l2+1):(l1+l2+l3))
  range4<-c((l1+l2+l3+1):(l1+l2+l3+l4))
  all_possible = c(range1, range2, range3, range4)
  possible_kj=list()

## ij 

   for ( range in list(range1, range2, range3, range4)  ){
    for (i in range){
      for (j in setdiff(all_possible,range) ){
        possible_kj <- c(possible_kj, list(c(i, j)))
      }}}
  return(possible_kj)
}


#get_all_possible_kj(l1=21,l2=14,l3=3,l4=2)




get_possible_positions<- function(l1=21,l2=14,l3=3,l4=2) #get possible positions only once and return a list of tuples of these positions
{
  
  
  list_possible_combinations <- list()
  #### case 1 additive aryle-halyde and (base or ligand)
  for (i in c(1:l1)) { #aditive
    for (j in c((l1+1):(l1+l2) ) ) { #aryl halide
      for (k in c( (l1+l2+1): (l1+l2+l3+l4) ) ) {  #base/ligand
        list_possible_combinations <- append(list_possible_combinations, list(c(i,j,k)))
      }}}
  
  
  ### case 2 additive base ligand
  for (i in c(1:l1)) { #aditive
    for (j in c( (l1+l2+1): (l1+l2+l3) ) ){ #base
      for (k in c( (l1+l2+l3+1): (l1+l2+l3+l4) ) ){ #ligand
        list_possible_combinations<-  append(list_possible_combinations, list(c(i,j,k)))
        
      }}}
  
  
  ### case 3 aryle-halide base ligand
  for (i in c((l1+1):(l1+l2))) { #aryl-halide
    for (j in c( (l1+l2+1): (l1+l2+l3) ) ){ #base
      for (k in c( (l1+l2+l3+1): (l1+l2+l3+l4) ) ){ #ligand
        list_possible_combinations<-  append(list_possible_combinations, list(c(i,j,k)))
        
      }}}
  
  return(list_possible_combinations)
}



my_table <- array(1:64, dim = c(4, 4, 4))

print(dim(xxx.all))

print(colnames(xxx.all)[2134])

psi_value_from_table_position<-function (table,i,j,k)
{return( (table[i,j,k] + table[i,k,j] + table [j,i,k] +table[j,k,i] + table[k,i,j] + table[k,j,i] )/6)}

print(get_possible_positions(l1=3,l2=2,l3=2,l4=2))

table_position_to_vector_index<- function(position_tuple, l1=21,l2=14,l3=2,l4=3) ## takes into account only possible combinations
{
  
  sum_comb_l234<-l2*l3+l2*l4+l3*l4
  l34<-l3+l4
  l12<-l1+l2
  l123<-l1+l2+l3
  
  x<-position_tuple[1]
  y<-position_tuple[2]
  z<-position_tuple[3]
  
  #case 1: additive aryl-halide base/lig
  if (x<l1+1 & y <l1+l2+1 )
    {position_psi<-(x-1)*sum_comb_l234 + (y-l1-1)*l34 + (z-l12) }
  
  #case 3: additive base ligand
  if (x<l1+1 & y<l1+l2+l3+1 & y>l1+l2 )
    {position_psi<- (x-1)*sum_comb_l234 + l2*l34 + (y-l12-1)*l4 + (z-l123)}
  
  
  #case 4: aryl-halide base ligand
  if (x>l1)
  {position_psi<- l1*sum_comb_l234 + (x-1-l1)*l3*l4 + (y-l12-1)*l4 +(z-l123)}
  
  return(position_psi)
  
  }

#table_position_to_vector_index(c(2,23,37))



get_psi_vec<-function(psi,l1=21,l2=14,l3=2,l4=3)
{
  assert(all(dim(psi)==l1+l2+l3+l4), "Dimensions are not ok")
  
  psi_vec<-array(0, dim=c(l1*l2*l3+l1*l2*l4+l1*l3*l4+l2*l3*l4) )
  
  ### CASE 1 have additive, aryl halide and anything else
  for (i in c(1:l1)) { #add
    for (j in c((l1+1):(l1+l2) ) ) { #halide
      for (k in c( (l1+l2+1): (l1+l2+l3+l4) ) ) {  #base\lig
        #cat("i:", i, ", j:", j, ", k:", k ,'\n')
        psi_vec[table_position_to_vector_index(c(i,j,k),l1=l1,l2=2,l3=l3,l4=l4)]<-psi_value_from_table_position(psi,i,j,k)
      }}}
  ### CASE 2 have additive,  base, ligand
  for (i in c(1:l1)) { #add
    for (j in c((l1+l2+1):(l1+l2+l3) ) ) { #base
      for (k in c( (l1+l2+l3+1): (l1+l2+l3+l4) ) ) {  #ligand
        #cat("i:", i, ", j:", j, ", k:", k ,'\n')
        psi_vec[table_position_to_vector_index(c(i,j,k),l1=l1,l2=2,l3=l3,l4=l4)]<-psi_value_from_table_position(psi,i,j,k)
      }}}
  
  ### CASE 3 have halide,  base, ligand
  for (i in c((l1+1):(l1+l2)) ) { #halide
    for (j in c((l1+l2+1):(l1+l2+l3) ) ) { #base
      for (k in c( (l1+l2+l3+1): (l1+l2+l3+l4) ) ) {  #ligand
        #cat("i:", i, ", j:", j, ", k:", k ,'\n')
        psi_vec[table_position_to_vector_index(c(i,j,k),l1=l1,l2=2,l3=l3,l4=l4)]<-psi_value_from_table_position(psi,i,j,k)
      }}}
  
  return(psi_vec)
}


psi<-  array(1, dim = c(9, 9, 9))
# 123 45 67 89
# Set some positions to 1
psi[1, 4, 6] <- 6
psi[1, 4, 7] <- 6
psi[5, 7, 8] <- 6

list(c(1,2,3),c(2,3))

get_psi_vec(psi = psi,l1=3,l2=2,l3=2,l4=2)

set_0s_psi<-function(psi,l1=21,l2=14,l3=2,l4=3)

{range1<-c(1:l1)
 range2<-c((l1+1):(l1+l2))
 range3<-c((l1+l2+1):(l1+l2+l3))
 range4<-c((l1+l2+l3+1):(l1+l2+l3+l4))
  
 ## ij same range
 for ( range in list(range1, range2, range3, range4) ){
     for (i in  range) {
     for (j in range)  {   
    psi[i,j,]<-0
    psi[i, ,j]=0
    psi[,i,j]=0
     }}}
 
 return(psi)
}

set_0s_psi(psi, l1=3,l2=2,l3=2,l4=2)


###############################################################################################################################################




# function for ONEROW
ONEROW_SEQ <- function(psi_tilda_jk, lambda,j,k, theta_bound_jk,  t=1, c=1) {
  
  f <- function(alpha) {
    #cat(' suma ',sum(abs(Soft_thresholding( Theta_tilda_j, t*(lambda/2+alpha)  ) ) ))
    return( sum(abs(Soft_thresholding( psi_tilda_jk, t*(lambda/6+alpha)  ) ) ) - c*abs(theta_bound_jk) )### 1
  }
  if ( f(0)<=0)  ### 1 a)
  {print('a')
    alpha_hat<-0
    return(final_return_seq(psi_hat_jk = psi_tilda_jk, lambda =  lambda, j=j, t=t, alpha_hat = alpha_hat) ) ### 2)
  }   
  
  
  ### 1b) ### 1c)
  knots_evaluated <- evaluate_knots_seq( psi_tilda_jk=psi_tilda_jk, lambda=lambda, f=f, t = t,c=c)
  ### 1d)
  knots=knots_evaluated$knots
  evaluated=knots_evaluated$evaluated
  #print( evaluated)
  
  zero_indices <- which(evaluated == 0)
  
  if (length(zero_indices) > 0) {
    print('d')
    alpha_hat <- knots[zero_indices[1]] # i.e. alpha_hat =p s.t. f(p)=0
    cat("alpha hat",alpha_hat)
    return (final_return_seq(psi_hat_jk = psi_tilda_jk, lambda =  lambda, j=j, t=t, alpha_hat = alpha_hat) )}   ### 2)
  
  
  alpha_hat <- find_adjacent_knots_and_alpha(knots, evaluated)    ### 1e)
  
  print('e')
  cat("alpha hat",alpha_hat)
  
  ### STEP 2 
  return (final_return_seq(psi_hat_jk = psi_tilda_jk, lambda =  lambda, j=j, t=t, alpha_hat = alpha_hat)) ### 2)
  
}

psi_tilda_jk=c(0,1,2,3,5)
lambda=6
j=5
k=7
theta_bound_jk=4
t=1



ONEROW_SEQ(psi_tilda_jk=psi_tilda_jk, lambda=lambda,j=j,k=k, theta_bound_jk=theta_bound_jk,  t=t, c=1)



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




# Create a 40x40x40 array filled with random numbers
my_array <- array(runif(), dim = c(40, 40, 40))

# Calculate the sum of all elements in the array
total_sum <- sum(my_array)

# Print the total sum
print(sum(my_array*my_array))

print(colnames(xxx.all))

######################### WEAKHIERNET CLASS #########################################

# Define the WeakHierNet class
WeakHierNet <- function(X, psi_init, y, theta_bound, lambda, t=1, tol=1e-6, max_iter=5000, eps=1e-8) {
  
  
 
    self = list()
    self$psi_hat <- psi_init
    self$theta_bound <- theta_bound
    

  
  
  # Public method for fitting the model
  fit <- function( X, psi_init, y, lambda, X_mean_init, t=1, tol=1e-2, max_iter=5000, eps=1e-8,l1=21,l2=14,l3=2,l4=3) {
    
    eps<-1e-8*lambda
    p <- ncol(X)
    n<-nrow(X)
    X <- scale(X)#standardize X
    y <- scale(y, center = TRUE, scale = FALSE) #center y # MAYBE ALREADY SCALED ACTUALLY???????????????????????
    
    # Initialize variables
    delta <- 1-t*eps
    psi_hat <- set_0s_psi(psi_init) #matrix form
    r_hat <- matrix(-1, nrow = n, ncol = 1)
    
    
    for (it in 2:max_iter+1) {
      
      vec_psi_hat <- get_psi_vec(psi_hat,l1=l1,l2=l2,l3=l3,l4=l4) ###select only good positions and take as sum /6
      r_hat_old <- r_hat
      r_hat <-y - X %*%   vec_psi_hat 
      r_hat<- -  r_hat ############################### TAKE CARE ! why like this? ?? why rhat/2 or not??? ########################
      if (k%%10==1)
      {cat(' Loss rhat',mean(r_hat^2))} ###/2 to be like in paper ????????????????????????????/ depends how is r_hat
      #cat("Theta", Theta_hat) 
      
      
      
      possible_kj <- get_all_possible_kj(l1=l1, l2=l2,l3=l3,l4=l4)#### ALL possible combinations 
      for (kj in possible_kj) { #possible positions kj
        k<-kj[1]
        j<-kj[2]
        
        
        ###CHECK HERE!!!
        Xkj<-get_cols_Xjk(X[i_cols],k,j)
        psi_hat[k,j] <-  ONEROW_SEQ(delta * psi_hat[k,j,] - t* t(Xkj)%*%r_hat, lambda=lambda ,j=j,  t=t) 
        
      }
      if (it>=3)
      {if (mean((r_hat_old- r_hat)^2) <=tol) #TAKE CARE RHAT IS VECTOR  !!! Sum because already scaled !!!!!!
      {  cat("Converged at iteration ",it)
        self$psi_hat=psi_hat
        return(self) }
      }
      
    }

    self$psi_hat=psi_hat
    cat("It has not converged. The difference between last 2 residuals is:", abs(r_hat[it-1]- r_hat[it-2]))
    return(self) 
  }
  
  
  # method for predicting
  
  predict <- function(self, new_X, mean_scale) {
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
    X <- scale(new_X)#standardize X
    y_pred= X %*%  self$vec_psi_hat + mean_scale  ##predicts already with scale   ###vezi cum faci vec_psi_hat eventual salveaza cu self
    return(y_pred)
  }
  
  R2_score <- function(self, new_X, y_true, verbose= TRUE) { 
    y_pred = predict(self,new_X) 
    
    if (verbose == TRUE)
    {cat ("r2 score is ", r2(y_true, y_pred))
      plot(y_pred, y_true)}
    
    return(r2(y_true, y_pred))
  }
  

  return(list( fit = fit, predict = predict, R2_score = R2_score, self = self))
}

##### How to use the class ######

## First create a basic dataset###


## CREATE DATASET
# Set a seed for reproducibility
set.seed(123)

# Number of rows in the dataset
n <- 200

# Generate predictor variables
col1 <- rnorm(n, mean=1,sd=10)
col2 <- rnorm(n, mean = 2,sd=10)

# Generate noise
noise <- rnorm(n, mean = 0, sd = 0.1)

# Generate response variable
response <- 0 * col1 + 0 * col2 + noise  + 2 *col1*col2

 #Combine predictors and response into a data frame
#synthetic_data <- data.frame(col1, col2, response)
X <- cbind(col1,col2)
y <- response



############ use the class on syntetic dataset ###


Beta_plus_init <- matrix(c(0, 0), ncol=1)
Beta_minus_init <- matrix(c(0, 0), ncol =1)
Theta_init <- matrix(100, nrow = 2, ncol = 2)
lambda <- 100
t<-0.0001
eps=1e-8

# Example usage:

# Create an instance of the WeakHierNet class
myWeakHierNet <- WeakHierNet(X=X, Beta_plus_init= Beta_plus_init  , Beta_minus_init = Beta_minus_init, 
                             Theta_init = Theta_init, y=y, lambda = lambda, t=t, tol=1e-6)

# Fit the model
fitted=myWeakHierNet$fit(X=X, y=y, lambda=lambda, t = t, tol = 1e-6, max_iter = 5000, eps = 1e-8, Beta_plus_init = Beta_plus_init,
                                   Beta_minus_init =  Beta_minus_init, Theta_init =  Theta_init)

# Make predictions
new_X <- X
predictions <- myWeakHierNet$predict(fitted, as.matrix(new_X))

myWeakHierNet$R2_score(self=fitted, new_X= as.matrix(new_X), y_true = y, verbose = TRUE)

fitted


### Use the library

library(Metrics)

library(hierNet)
library(caret)
library(dplyr)
library(Metrics)


lambda=1e2
fit=hierNet(as.matrix(X),y, lam = lambda, diagonal = FALSE, step = 1e-4)
yhat=predict(fit,as.matrix(X))
print(paste("r2- hiernet library:", r2(y, yhat)))
fit$bp
fit$bn
fit


sum(abs(colSums(xxx.all)))

dim(xxx.all)
