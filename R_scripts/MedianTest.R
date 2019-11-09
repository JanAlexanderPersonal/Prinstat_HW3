# FUNCTION: calc.median.diff: --------------
# Function: calculate the difference in median between two groups.
#
# param: ind = the indices of the first group
# param: vec = the complete list with data from group 1 and group 2
#
calc.median.diff <- function(ind, vec){
  # make both groups
  group1 <- vec[ind]
  group2 <- vec[!(1:length(vec)) %in% ind]
  # calculate the difference in means
  return(median(group1) - median(group2))
}

# FUNCTION: median.test: --------------
# Function: calculate the p-value of the median difference between
# vector x and vector y, based on the null hypothesis F1(x) = F2(x)
#
# When the number of combinations is sufficiently small to do a full 
# permutation test (limit at 10000), the full permutation test will
# be performed. 
#
# param: x = vector values for group 1
# param: y = vector values group 2
#
median.test <- function(x, y){
  N <- 10000 												# maximum number of permutations
  realization <- median(x) - median(y)			# The difference in median observed for sets x and y
  len_x <- length(x)
  len_y <- length(y)
  vec <- c(x, y)
  len_vec <- len_x + len_y
  if(choose(n = len_vec, len_x) < N){
    #A limited number of combinations: full permutation test
    median.diff <- combn(len_vec, 
                         len_x, 
                         function(ind){
                           return(calc.median.diff(ind, vec))
                         })
  } else {
    #Too many combinations - A sample test will be performed.
    median.diff <- replicate(N, 
                             calc.median.diff(
                               sample(c(1:len_vec),len_x), 
                               vec)
    )
  }
  # The distribution will be symmetrical. The p-value can be calculated by taking
  # the absolute value.
  return(mean(abs(realization) <= abs(median.diff)))
}