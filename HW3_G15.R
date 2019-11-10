# HOMEWORK 3
# Group nr 15
#
#
# Install and load packages
packages <- c('coin','ggplot2', 'reshape2', 'rmutil')
# install.packages(packages)
lapply(packages, library, character.only = TRUE)
#
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
  # calculate the difference in medians
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
  N <- 10000
  realization <- median(x) - median(y)
  len_x <- length(x)
  len_y <- length(y)
  vec <- c(x, y)
  len_vec <- len_x + len_y
  if(choose(n = len_vec, len_x) < N){
    # A limited number of combinations: full permutation test 
    median.diff <- combn(len_vec, 
                         len_x, 
                         function(ind){
                           return(calc.median.diff(ind, vec))
                         })
  } else {
    # Too many combinations - A sample test will be performed.
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

calc.rejection <- function(N, n, delta, p.val = 0.05){
  distributions <- c('rt', 'rt', 'rexp', 'rlogis', 'rnorm', 'runif', 'rlaplace')
  dist_names <- c('t3', 't5', 'exp', 'logis', 'norm', 'unif', 'laplace')
  dist_var <- c(3, 5/3, 1, 2*pi^2/6, 1, 1/12, 2)
  dist_args <- list(list(n, 3), list(n, 5), list(n), list(n), list(n), list(n), list(n))
  for(i in 1:length(distributions)) {
    cat(paste('power calculation for distribution : ',
              distributions[i], '\n'))
    p.med <- p.wmw <- p.t <- c()
    
    # perform N tests on randomly drawn samples Y1 and Y2
    for(j in 1:N) {
      Y1 <- do.call(what = distributions[i], 
                    args = dist_args[[i]])
      Y2 <- do.call(what = distributions[i], 
                    args = dist_args[[i]]) + delta * dist_var[i] / 2
      df <- data.frame(rep(c('A', 'B'), each = n), c(Y1, Y2))
      colnames(df) <- c("group", "Y")
      p.med[j] <- median.test(x = Y1, y = Y2)
      p.wmw[j] <- wilcox.test(Y1, Y2, exact = TRUE)$p.value
      p.t[j] <- pvalue(oneway_test(Y~group, data = df, 
                                   distribution = approximate(nresample = 10000)))
    }
    # calculate the proportion of H0 rejections
    power.med[i] <- mean(p.med < p.val)
    power.wmw[i] <- mean(p.wmw < p.val)
    power.t[i] <- mean(p.t < p.val)
  }
  df_power <- data.frame(Distribution = dist_names,
                         med = rep(x = 0.0, length(distributions)),
                         wmw = rep(x = 0.0, length(distributions)),
                         t = rep(x = 0.0, length(distributions)))
  
  df_power$med <- power.med
  df_power$wmw <- power.wmw
  df_power$t <- power.t
  df_power <- melt(df_power, id.vars = 'Distribution')
  return(df_power)
}



# Simulations
N <- 20        # test simulations
power.med <- power.wmw <- power.t <- c()

for(n in c(20, 40)){
# For simulations under H0 (delta = 0) and Ha (delta = 1)
for(delta in 0:1) {
  if(delta) {
    title <- paste0('Power_', n, '_', N) }
  else{
    title <- paste0('TypeI_error_', n, '_', N) }
    
  df_power <- calc.rejection(N, n, delta)
  
  write.csv(df_power, paste(title, '.csv', sep=''))
  
  p <- ggplot(data = df_power, aes(x = Distribution, y = value, fill = variable)) + 
    geom_bar(position = 'dodge', stat = 'identity') + 
    ggtitle(title) +
    xlab('Distribution') + 
    ylab('Power') +
    coord_flip()
  
  ggsave( paste(title, '.png', sep = ''), 
          plot = p, 
          device = 'png'
          )
}
}