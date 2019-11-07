# HOMEWORK 3
# Group nr 15
#
#
library(coin)
library(ggplot2)
library(reshape2)
library(rmutil)
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
  # calculate the difference in means
  return(median(group1) - median(group2))
}

# FUNCTION: median.test: --------------
# Function: calculate the p-value of the median difference between
# vector x and vector y, based on the null hypothesis F1(x) = F2(x)
#
# When the number of combinations is sufficiently small to do a full 
# permutation test (limit at 15000), the full permutation test will
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
  if(choose(n = len_x+len_y, len_x) < N){
    #cat('A limited number of combinations: full permutation test')
    median.diff <- combn(len_vec, 
                         len_x, 
                         function(ind){
                           return(calc.median.diff(ind, vec))
                         })
  } else {
    #cat('Too many combinations - A sample test will be performed.\n')
    median.diff <- replicate(N, 
                             calc.median.diff(
                               sample(c(1:len_vec),len_x), 
                               vec)
    )
  }
  if(FALSE){
    hist(abs(median.diff), main = paste('p = ', 
                                        mean(abs(realization) <= abs(median.diff))))
    abline(v = abs(realization), col = "red",
           lty = 1, lwd = 2)
  }
  # The distribution will be symmetrical. The p-value can be calculated by taking
  # the absolute value.
  return(mean(abs(realization) <= abs(median.diff)))
}

# FUNCTION: calculate.power: --------------
# Function: calculate the proportion of tests that yield a p-value
# smaller than 0.05 for the median.test, the Wilcoxon exact
# test and the permutation t-test after repetitive sampling.
# 
# param: delta = effect size
# param: d = random distribution function
# dist_arg = arguments or the distribution function
# N_tests = the amount of repetitions
# p_value = significance level
#
calculate.power <- function(delta, d, dist_arg, N_tests = 150, p_value = 0.05){
  cat('power calculation for distribution : ')
  cat(d)
  cat('\n')
  p.perm <- p.wmw <- p.t <- c()
  for(i in 1:N_tests){
    Y1 <- do.call(what = d, 
                  args = dist_arg)
    Y2 <- do.call(what = d, 
                  args = dist_arg) + delta
    df <- data.frame(rep(c('A', 'B'), each = n), c(Y1, Y2))
    colnames(df) <- c("group", "Y")
    p.perm[i] <- median.test(x = Y1, y = Y2)
    p.wmw[i] <- wilcox.test(Y1, Y2, exact = TRUE)$p.value
    p.t[i] <- pvalue(oneway_test(Y~group, data = df, 
                                 distribution = approximate(nresample = 10000)))
  }
  if(FALSE){
    hist(p.perm)
    hist(p.wmw)
  }
  return(c(mean(p.perm < p_value),
           mean(p.wmw < p_value),
           mean(p.t < p_value)
  )
  )
}

# Simulations
n <- 20
distributions <- c('rt', 'rt', 'rexp', 'rlogis', 'rnorm', 'runif', 'rlaplace')
dist_names <- c('t3', 't5', 'exp', 'logis', 'norm', 'unif', 'laplace')
dist_args <- list(list(n, 3), list(n, 5), list(n), list(n), list(n), list(n), list(n))

power.perm <- power.wmw <- power.t <- c()

for(delta in c(0, 0.5, 0.75)){
  for(i in 1:7){
    N_tests <- 1000
    d <- distributions[i]
    res <- calculate.power(delta, d, dist_args[[i]], N_tests = N_tests)
    power.perm[i] <- res[1]
    power.wmw[i] <- res[2]
    power.t[i] <- res[3]
  }
  df_power <- data.frame(Distribution = dist_names,
                         Perm = rep(x = 0.0, length(distributions)),
                         wmw = rep(x = 0.0, length(distributions)),
                         t = rep(x = 0.0, length(distributions)))
  df_power$Perm <- power.perm
  df_power$wmw <- power.wmw
  df_power$t <- power.t
  df_power <- melt(df_power, id.vars = 'Distribution')
  
  write.csv(df_power, paste('delta_', delta*100, 'percent.csv', sep=''))
  
  if(delta != 0){
    title <- paste('Power of perm vs wmw for different distributions for delta = ', delta,
                   '\n(n =', n, ', # simulations =', N_tests, ')')
  }else{
    title <- paste('Proportion of type I error - delta = ', delta,
                   '\n(n =', n, ', # simulations =', N_tests, ')')
  }
  
  plot1 <- ggplot(data = df_power, aes(x = Distribution, y = value, fill = variable)) + 
    geom_bar(position = 'dodge', stat = 'identity') + 
    ggtitle(title) +
    xlab('Distribution') + 
    ylab('Power') + 
    coord_flip()

  ggsave( paste('delta', delta*100, '.png', sep = ''), 
          plot = plot1, 
          device = 'png'
          )
}
