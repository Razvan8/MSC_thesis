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
#X=matrix(c(2,2,1,1), nrow = 2, ncol = 2, byrow = TRUE)
#print(X)
#y=matrix(c(9,4), nrow = 2, ncol = 1)
#psi=matrix(c(1,3), nrow = 2, ncol = 1)
#lambda=1
#hierarchical_lasso_seq_loss(X, y, psi, lambda)


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
  
  #cat("psitilda_jk", psi_tilda_jk)
  
  set_elements <- 0##add also 0 to knots, it has positive value in evaluated
  
  for (k in 1:p)  { 
    if (abs(psi_tilda_jk[k])>0){ ##discard impossible combinations
      element <- abs(psi_tilda_jk[k]) / t - lambda / 6 #element to add in knot set
      set_elements <- c(set_elements, element)}
  }
  
  selected_elements <- set_elements[set_elements >= 0] ### cred ca e nevoie si de 0?!?!
  
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
  #print(evaluated)
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
    alpha_hat <- (p1*evaluated[positive_index + 1] - p2*evaluated[positive_index]) / (evaluated[positive_index + 1] - evaluated[positive_index]) #new
    #print(positive_index)
    return(alpha_hat)
  } else {
    # Return an informative message if no such adjacent knots are found
    stop("No adjacent knots found with evaluated > 0 and evaluated < 0 at consecutive positions.")
  }
}

final_return_seq<- function(psi_hat_jk, lambda, t, alpha_hat) # returns psi_jk with all the elements inside (vector with 40 elements)
{ 
psi_hat_jk = Soft_thresholding(psi_hat_jk, t*(lambda/6 + alpha_hat))
return (psi_hat_jk)
}

## FUNCTIONS ##########################################################

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


get_all_possible_kj(l1=1,l2=1,l3=1,l4=1)




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



#my_table <- array(1:64, dim = c(4, 4, 4))



psi_value_from_table_position<-function (table,i,j,k)
{return( (table[i,j,k] + table[i,k,j] + table [j,i,k] +table[j,k,i] + table[k,i,j] + table[k,j,i] )/6)}

#print(get_possible_positions(l1=3,l2=2,l3=2,l4=2))

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
        psi_vec[table_position_to_vector_index(c(i,j,k),l1=l1,l2=l2,l3=l3,l4=l4)]<-psi_value_from_table_position(psi,i,j,k)
      }}}
  ### CASE 2 have additive,  base, ligand
  for (i in c(1:l1)) { #add
    for (j in c((l1+l2+1):(l1+l2+l3) ) ) { #base
      for (k in c( (l1+l2+l3+1): (l1+l2+l3+l4) ) ) {  #ligand
        #cat("i:", i, ", j:", j, ", k:", k ,'\n')
        psi_vec[table_position_to_vector_index(c(i,j,k),l1=l1,l2=l2,l3=l3,l4=l4)]<-psi_value_from_table_position(psi,i,j,k)
      }}}
  
  ### CASE 3 have halide,  base, ligand
  for (i in c((l1+1):(l1+l2)) ) { #halide
    for (j in c((l1+l2+1):(l1+l2+l3) ) ) { #base
      for (k in c( (l1+l2+l3+1): (l1+l2+l3+l4) ) ) {  #ligand
        #cat("i:", i, ", j:", j, ", k:", k ,'\n')
        psi_vec[table_position_to_vector_index(c(i,j,k),l1=l1,l2=l2,l3=l3,l4=l4)]<-psi_value_from_table_position(psi,i,j,k)
      }}}
  
  return(psi_vec)
}


# Create the 3-dimensional array
#array_data <- array(data = c(
  # Values for the first layer
  #0, 0.000000, 0.0000000, 0.000000,
  #0, 0.000000, -3.3490767, 1.899415,
  #0, 2.887088, 0.0000000, 3.100585,
  #0, 1.206102, 0.1402117, 0.000000,
  
  # Values for the second layer
  #0.0000000, 0, 0.8724023, -0.07197379,
  #0.0000000, 0, 0.0000000, 0.00000000,
  #-3.7229584, 0, 0.0000000, 0.50213101,
  #0.7168255, 0, 1.8879085, 0.00000000,
  
  # Values for the third layer
  #0.000000, 3.247535, 0, 0.3154444,
  #2.306129, 0.000000, 0, 1.5593346,
  #0.000000, 0.000000, 0, 0.0000000,
  #-1.263373, -1.752465, 0, 0.0000000,
  
  # Values for the fourth layer
  #0.0000000, 1.4078473, 1.676665, 0,
  #0.3514781, 0.0000000, 1.668121, 0,
  #3.0134272, 0.1127374, 0.000000, 0,
  #0.0000000, 0.0000000, 0.000000, 0
#), dim = c(4, 4, 4))

# Print the array
#print(array_data)


#psi<-  array(1, dim = c(4, 4, 4))
# 123 45 67 89
# Set some positions to 1
#psi[1, 2, 3] <- 6
#psi[1, 3, 2] <- 6
#psi[3, 1, 2] <- 6

#get_psi_vec(psi = array_data,l1=1,l2=1,l3=1,l4=1)











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



get_range<- function(x,l1=21,l2=14,l3=2,l4=3)
{if (x<=l1)
   {return(c(1:l1))}
 if (x<=l1+l2)
    {return(c( (l1+1) : (l1+l2) ))}
 if (x<=l1+l2+l3)
 {return(c( (l1+l2+1) : (l1+l2+l3) ))}
  return(c( (l1+l2+l3+1) : (l1+l2+l3+l4) ))
}



#get_cols_Xkj<-function(X_main,k,j,l1=21,l2=14,l3=2,l4=3) #####De verificat!!!!!!!!!!!!!
#{Xkj<- X_main[,k]*X_main[,j]
#cat("dim Xkj init", dim(Xkj))
#X_cols_kj <- X_main*Xkj
#X_cols_kj<-scale(X_cols_kj)
#range_k<-get_range(k,l1=l1,l2=l2,l3=l3,l4=l4)
#range_j=get_range(j,l1=l1,l2=l2,l3=l3,l4=l4)
#X_cols_kj[,c(range_k,range_j)]<-0
#return(X_cols_kj) #should have 40 columns in my case, all scaled
#}

#get_cols_Xkj(matrix(c(1,2,3,4,5,6,7,8),nrow=32,ncol=8),1,5,2,2,2,2)



get_cols_Xkj<-function(X,k,j,l1=21,l2=14,l3=2,l4=3) #####De verificat!!!!!!!!!!!!! ###assumes data like buchwald-hartwig form
  
{ X_kj<-matrix(0,nrow=nrow(X), ncol = l1+l2+l3+l4)
  range1<-c(1:l1)
  range2<-c((l1+1):(l1+l2))
  range3<-c((l1+l2+1):(l1+l2+l3))
  range4<-c((l1+l2+l3+1):(l1+l2+l3+l4))
  all_possible = c(range1, range2, range3, range4)
  
  range_k<-get_range(k,l1=l1,l2=l2,l3=l3,l4=l4)
  range_j<-get_range(j,l1=l1,l2=l2,l3=l3,l4=l4)
  range_kj<-c(range_k, range_j)
  range_possible<-setdiff(all_possible,range_kj)
  for (i in range_possible){
    X_kj[,i]<-X[,table_position_to_vector_index(sort(c(i,j,k)),l1=l1,l2=l2,l3=l3,l4=l4)]
  }
  return(X_kj)
}

#get_cols_Xkj(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7), nrow=4, ncol =7),3,4, 2,1,1,1)

#print(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7),nrow=4, ncol=7))


###############################################################################################################################################




# function for ONEROW
ONEROW_SEQ <- function(psi_tilda_jk, lambda, theta_bound_jk,  t=1, c=1) {
  
  f <- function(alpha) {
    #print('alpha')
    #print(alpha)
    #cat(' suma ',sum(abs(Soft_thresholding( Theta_tilda_j, t*(lambda/2+alpha)  ) ) ))
    result<- sum(abs(Soft_thresholding( psi_tilda_jk, t*(lambda/6+alpha)  ) ) ) - c*abs(theta_bound_jk) ### 1
    if (abs(result)<1e-7) ###this is new
    {return(0)}## to solve numerical instability
    return(result)
  }
  if ( f(0)<=0)  ### 1 a)
  {#print('a')
    alpha_hat<-0
    return(final_return_seq(psi_hat_jk = psi_tilda_jk, lambda = lambda, t=t, alpha_hat = alpha_hat) ) ### 2)
  }   
  
  
  ### 1b) ### 1c)
  knots_evaluated <- evaluate_knots_seq( psi_tilda_jk=psi_tilda_jk, lambda=lambda, f=f, t = t,c=c)
  ### 1d)
  knots=knots_evaluated$knots
  evaluated=knots_evaluated$evaluated
  #print(knots_evaluated)
  
  zero_indices <- which(evaluated == 0)
  
  if (length(zero_indices) > 0) {
    #print('d')
    alpha_hat <- knots[zero_indices[1]] # i.e. alpha_hat =p s.t. f(p)=0
    #cat("alpha hat",alpha_hat)
    return (final_return_seq(psi_hat_jk = psi_tilda_jk, lambda =lambda, t=t, alpha_hat = alpha_hat) )}   ### 2)
  
  
  alpha_hat <- find_adjacent_knots_and_alpha(knots, evaluated)    ### 1e)
  
  #print('e')
  #cat("alpha hat",alpha_hat)
  
  ### STEP 2 
  return (final_return_seq(psi_hat_jk = psi_tilda_jk, lambda =  lambda, t=t, alpha_hat = alpha_hat)) ### 2)
  
}

#psi_tilda_jk=c(3,2,0,0,0)
#lambda=6
#theta_bound_jk=2
#t=1

#ONEROW_SEQ(psi_tilda_jk=psi_tilda_jk, lambda=lambda, theta_bound_jk=theta_bound_jk,  t=t)



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
WeakHierNet_seq <- function(X, psi_init, y, theta_bound, lambda, t=1, tol=1e-6, max_iter=5000, eps=1e-8,l1=21,l2=14,l3=2,l4=3) {
  
  
 
    self = list()
    self$psi_hat <- psi_init
    self$theta_bound <- theta_bound
    self$vec_psi_hat<-get_psi_vec(psi_init,l1=l1,l2=l2,l3=l3,l4=l4)
    self$means_X<-colMeans(X)
    self$stds_X<-apply(X,2,sd)
    self$mean_y<-mean(y)
    
    

  
  
  # Public method for fitting the model
  fit <- function( X, psi_init, y, lambda,theta_bound, t=1, tol=1e-2, max_iter=5000, eps=1e-8,l1=21,l2=14,l3=2,l4=3) {
    
    eps<-1e-8*lambda
    p <- ncol(X)
    n<-nrow(X)
    X <- scale(X)#standardize X
    y <- scale(y, center = TRUE, scale = FALSE) #center y # MAYBE ALREADY SCALED ACTUALLY???????????????????????
    
    # Initialize variables
    delta <- 1-t*eps
    psi_hat <- set_0s_psi(psi_init,l1=l1,l2=l2,l3=l3,l4=l4) #matrix form
    r_hat <- matrix(-1, nrow = n, ncol = 1)
    
    
    for (it in 2:max_iter+1) {
      
      vec_psi_hat <- get_psi_vec(psi_hat,l1=l1,l2=l2,l3=l3,l4=l4) ###select only good positions and take as sum /6
      r_hat_old <- r_hat
      r_hat <-y - X %*%   vec_psi_hat 
      r_hat<- -  r_hat ############################### TAKE CARE ! why like this? ?? why rhat/2 or not??? ########################
      if (it%%50==1)
      {cat(' Loss rhat',mean(r_hat^2))} ###/2 to be like in paper ????????????????????????????/ depends how is r_hat
      #cat("Theta", Theta_hat) 
      
      
      
      possible_kj <- get_all_possible_kj(l1=l1, l2=l2,l3=l3,l4=l4)#### ALL possible combinations
      #print("possible kj")
      #print(possible_kj)
      for (kj in possible_kj) { #possible positions kj
        k<-kj[1]
        j<-kj[2]
        #cat("k,j",c(k,j))
        
        
        ###CHECK HERE!!!
        Xkj<-get_cols_Xkj(X=X,k=k,j=j,l1=l1,l2=l2,l3=l3,l4=l4)
        psi_hat[,k,j] <-  ONEROW_SEQ(delta * psi_hat[,k,j] - t* t(Xkj)%*%r_hat, lambda=lambda, t=t, theta_bound_jk = theta_bound[j,k]) 
        
      }
      if (it>=3)
      {if (mean((r_hat_old- r_hat)^2) <=tol) #TAKE CARE RHAT IS VECTOR  !!! Sum because already scaled !!!!!!
      {  cat("Converged at iteration ",it)
        self$psi_hat=psi_hat
        self$vec_psi_hat=get_psi_vec(psi_hat,l1=l1,l2=l2,l3=l3,l4=l4)
        return(self) }
      }
      
    }

    self$psi_hat=psi_hat
    self$vec_psi_hat=get_psi_vec(psi_hat,l1=l1,l2=l2,l3=l3,l4=l4)
    print(self$vec_psi_hat)
    cat("It has not converged. The difference between last 2 residuals is:", abs(r_hat[it-1]- r_hat[it-2]))
    return(self) 
  }
  
  
  # method for predicting
  
  predict <- function(self, new_X) {
    
    ## Later add also rescaling
    # Implement prediction code here
    p <- ncol(new_X)
    n<-nrow(new_X)
    X <- (new_X-self$means_X)/self$stds_X ### scale X
    y_pred= X %*%  self$vec_psi_hat + self$mean_y
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


row1<-c(1,1,1,1)
row2<-c(1,1,1,-1)
row3<-c(1,1,-1,1)
row4<-c(1,1,-1,-1)

row5<-c(1,-1,1,1)
row6<-c(1,-1,1,-1)
row7<-c(1,-1,-1,1)
row8<-c(1,-1,-1,-1)

row9<-c(-1,1,1,1)
row10<-c(-1,1,1,-1)
row11<-c(-1,1,-1,1)
row12<-c(-1,1,-1,-1)

row13<-c(-1,-1,1,1)
row14<-c(-1,-1,1,-1)
row15<-c(-1,-1,-1,1)
row16<-c(-1,-1,-1,-1)

#X<-rbind(row1, row2, row3, row4, row5, row6, row7, row8, row9, row10, row11, row12, row13, row14, row15, row16)



###MAKE 2^4 times each matrix

#for (i in 1:5) {
  #X <- rbind(X, X)}





# Create a beta vector of length 32
#beta <- numeric(4)
#for (i in 1:length(beta)) {
#  beta[i] <- i %% 5
#}

#print(beta)
#noise <- rnorm(512, mean = 0, sd = 0.1)
#print(dim(X))
# Generate response variable
#y <- X%*%beta+noise


#library(caret)

# Assuming X and y are your predictor variables and response variable, respectively

# Set the seed for reproducibility
set.seed(123)

# Split the data into training and testing sets
#train_index <- createDataPartition(y, p = 0.8, list = FALSE)
#X_train <- X[train_index, ]
#X_test <- X[-train_index, ]
#y_train <- y[train_index]
#y_test <- y[-train_index]


#print(length(y_train))
#print(dim(X_train))

#lm_model <- lm(y_train ~ ., data = as.data.frame(X_train))

#Predict on X_test
#predictions_test <- predict(lm_model, newdata = as.data.frame(X_test))
#predictions_train<-predict(lm_model, newdata =  as.data.frame(X_train))

#print(r2(y_train, predictions_train))
#print(r2(y_test,predictions_test))







############ use the class on syntetic dataset ###


#lambda <-0
#t<-0.0005
#eps=1e-8
#psi_init<-array(rnorm(4*4*4, mean = 0, sd = 2), dim = c(4,4,4))
#theta_bound<-matrix(2,nrow=4,ncol=4)
#theta_bound[c(1,2,3),c(1,2,3)]<-1






# Example usage:

# Create an instance of the WeakHierNet class
#myWeakHierNet_seq <- WeakHierNet_seq(X=X_train, psi_init=psi_init, y=y_train, theta_bound=theta_bound, lambda=lambda, t=t, tol=1e-8, max_iter=5000, eps=1e-8,
 #                                    l1=1,l2=1,l3=1,l4=1)

# Fit the model
#fitted=myWeakHierNet_seq$fit(X=X_train, psi_init=psi_init, y=y_train, theta_bound=theta_bound, lambda=lambda, t=t, tol=1e-8, max_iter=5000, 
 #                            eps=1e-8,l1=1,l2=1,l3=1,l4=1)
#print(fitted$vec_psi_hat)
#print(beta)
#print(fitted$psi_hat)

# Make predictions
#new_X <- X_test

#print(dim(X_train))
#print(dim(X_test))
#print(length(colMeans(X_train)))
#(X_test-colMeans(X_train)) /apply(X_train, 2, sd)
#print(length(fitted$vec_psi_hat))
#pred_train<-scale(X_train)%*%fitted$vec_psi_hat+mean(y_train)
#pred_test<- ( (X_test-colMeans(X_train))/apply(X_train,2,sd) ) %*%fitted$vec_psi_hat+mean(y_train)
#print(r2(y_train,pred_train))
#print(r2(y_test,pred_test))

#myWeakHierNet_seq$R2_score(self=fitted, new_X= as.matrix(new_X), y_true = y_test, verbose = TRUE)

#fitted


### Use on another synthetic dataset is columns continous values from normal distribution

# Set the number of predictors
#num_predictors <- 68

# Set the beta values
#beta <- c(5, 3, 7, 9, 10, rep(1e-4, num_predictors - 5))

# Generate the predictor variables
#set.seed(123)  # for reproducibility
#predictors <- matrix(rnorm(num_predictors * 1000, mean = 0, sd = 5), ncol = num_predictors)

# Generate the noise
#noise <- rnorm(1000, mean = 0, sd = 0.3)

# Generate the response variable
#response <- predictors %*% beta + noise

# Combine predictors and response into a data frame
#data <- data.frame(response = response, predictors)


# Split predictors and response into training and testing sets
#set.seed(123)  # for reproducibility
#train_index <- sample(1:1000, 700)  # 70% for training
#test_index <- setdiff(1:1000, train_index)

#X_train <- predictors[train_index, ]
#X_test <- predictors[test_index, ]
#y_train <- response[train_index]
#y_test <- response[test_index]

#lm_model <- lm(y_train ~ ., data = as.data.frame(X_train))

# Make predictions on the testing data
#predictions_test <- predict(lm_model, newdata = as.data.frame(X_test))

# Calculate the R-squared value on the testing data
#r2(y_test, predictions_test)



##
#lambda <-3e2
#t<-8e-5
#eps=1e-8
#psi_init<-array(rnorm(11*11*11, mean = 0, sd = 2), dim = c(11,11,11))
#theta_bound<-matrix(1,nrow=11,ncol=11)






# Example usage:

# Create an instance of the WeakHierNet class
#myWeakHierNet_seq <- WeakHierNet_seq(X=X_train, psi_init=psi_init, y=y_train, theta_bound=theta_bound, lambda=lambda, t=t, tol=1e-7, max_iter=5000, eps=1e-8,
 #                                    l1=5,l2=2,l3=2,l4=2)

# Fit the model
#fitted=myWeakHierNet_seq$fit(X=X_train, psi_init=psi_init, y=y_train, theta_bound=theta_bound, lambda=lambda, t=t, tol=1e-7, max_iter=5000, 
#                             eps=1e-8,l1=5,l2=2,l3=2,l4=2)
#print(fitted$vec_psi_hat/fitted$stds_X)
#print(beta)
#print(fitted$psi_hat)

# Make predictions
#new_X <- X_test
#print(dim(X_test))
#pred_train<-scale(X_train)%*%fitted$vec_psi_hat+mean(y_train)
#pred_test<-scale(X_test, center = FALSE)%*%fitted$vec_psi_hat+mean(y_train)
#pred_test2<- ((X_test-colMeans(X_train))/apply(X_train,2,sd) ) %*%fitted$vec_psi_hat+mean(y_train)
#print(r2(y_train,pred_train))
#print(r2(y_test,pred_test2))



#print('USING METHODS FROM MY CLASS')

#print(myWeakHierNet_seq$R2_score(fitted,X_test,y_test))




################ GANDURI ###################

## FA COLS CU MEAN 0 prin tehnica cu -1 in 
## Vf COD iar
### Fa pred si score methods


###TAKE CARE WHICH KNOTS I USE FOR e)
####
#### HOW TO DO SCALE AND MEAN IN PRED AND R2 SCORE
#### VEZI CUM E FECETUL LUI THETA BOUND DUPA SCALE SAU CAND E MIC
##why small numerical instability or variance from one run to another

