swap_increase <- function(x, y) {
  if (x>y)
  {temp <- x
  x <- y
  y <- temp}
  return(list(x = x, y = y))
}

find_index_group <- function(lst, num) {
  if (num <=lst[1] && num>0)
  {return (1)}
  total_sum <- 0
  for (i in seq_along(lst)) {
    total_sum <- total_sum + lst[i]
    if (num > total_sum && num <= total_sum + lst[i + 1]) {
      return(i+1)
    }
  }
  print("returned NA")
  return(NA)  # If no such index is found
}

# Example usage
#my_list <- c(1,1,1)
#my_number <- 3
#index <- find_index_group(my_list, my_number)
#if (!is.na(index)) {
  #cat("Index found:", index, "\n")
#} else {
  #cat("No such index found\n")
#}



interaction_ab<-function(idx_a, idx_b, num_categories_list=c(22,15,3,4), main_lengths=44)
{partial_sums=rep(0,length(num_categories_list))

#Will use to compute idx
for (i in c(1:length(num_categories_list)) )
{partial_sums[i]<-sum(num_categories_list[1:i])}

befores=rep(0,length(num_categories_list))
for (i in c(2:length(num_categories_list)) )
  befores[i]<- befores[i-1]+num_categories_list[i-1]*(sum(num_categories_list)-partial_sums[i-1])


if (idx_a>idx_b)
{ab=swap_increase(idx_a,idx_b)
idx_a=ab$x#smaller one
idx_b=ab$y}#bigger one

group_a<-find_index_group(num_categories_list,idx_a) # i.e 2 from 4 groups
group_b<-find_index_group(num_categories_list, idx_b)

if (group_a==1)
{positions_before_a<-(idx_a-1)*(sum(num_categories_list)-partial_sums[group_a])}
else if (group_a>1){
  positions_before_a<-befores[group_a]+(idx_a-partial_sums[group_a-1]-1)*(sum(num_categories_list)-partial_sums[group_a])
}
positions_from_b<-idx_b-partial_sums[group_a]

if (group_a==group_b){print("There are interactions from same feature, not ok!")}

idx_ab<-positions_before_a +positions_from_b

return (idx_ab+main_lengths)}

#print(interaction_ab(idx_a=25,idx_b=44))

create_all_interactions<- function( X_main,num_categories_list=c(22,15,3,4), main_lengths=44) #maybe optimize and use this only
{Beta_interaction_list<-c()
X_all<-X_main
for (idx_a in c(1:sum(num_categories_list[1:length(num_categories_list)-1])) )
{
  group_a<-find_index_group(num_categories_list,idx_a)
  for (idx_b in c(sum(num_categories_list[1:group_a]) :main_lengths) )
    
  {
    group_b<-find_index_group(num_categories_list,idx_b)
    if (group_b != group_a)
    {
      
      X_all<-cbind(X_all,X_all[,idx_a]*X_all[,idx_b]) 
      idx_ab<-interaction_ab(idx_a, idx_b,num_categories_list, main_lengths)
      Beta_interaction_list<- cbind(Beta_interaction_list,c(idx_a, idx_b, idx_ab))
    }}}
return (list("X"=X_all, "beta_interactions"=Beta_interaction_list))
}




cap <- function(beta, X, y,beta_interactions, lambda=0) {
  result <- (1/dim(X)[1]) * sum((y - X %*% beta)^2)  
  for (i in c(1:ncol(beta_interactions)) )
  { 
    result<-result+lambda*max(abs( beta[beta_interactions[,i]] )) ## add terms like (b1 b2 b12)
    result<-result+ lambda* abs(beta[beta_interactions[3,i] ]) # add all terms like b12
    
  }
  print(result)
  return(result)
}

gradient_max_abs <- function(b) {
  abs_values <- abs(b)
  max_abs <- max(abs_values)
  gradient <- ifelse(abs_values == max_abs, sign(b), 0)
  return(gradient)
}


gradient_cap<- function(beta, X, y,beta_interactions, lambda=0)
{gr<-2/dim(X)[1] * ( t(X)%*%(X%*%beta-y) )
for (i in c(1:ncol(beta_interactions)) )
{ 
  gr[beta_interactions[,i]]<- gr[beta_interactions[,i]]+lambda*gradient_max_abs(beta[beta_interactions[,i]])## add terms like (b1 b2 b12)
  gr[beta_interactions[3,i]]<-  gr[beta_interactions[3,i]]+ lambda* sign(beta[beta_interactions[3,i] ]) # add all terms like b12
}
return(gr)
}
