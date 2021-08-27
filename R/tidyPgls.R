tidyPgls <- function(model) {
  # This function restructures the output of a pgls model and additionally
  # returns 95% confidence intervals for easier exporting.
  
  if (class(model) != "pgls") stop("Not an object of class 'pgls' ")
  
  summ <- summary(model)
  collection <- data.frame(cbind(rbind(summ$coefficients[2,]), rbind(summ$coefficients[1,])))
  colnames(collection) <- c("slope_estimate", "slope_std_err", "slope_t", "slope_P", "intercept_estimate", "intercept_std_err", "intercept_t", "intercept_P")
  
  collection$y <- model$varNames[1]
  collection$x <- model$varNames[2]
  collection$fstat <- summ$fstatistic[1]
  collection$adj_r2 <- summ$adj.r.squared
  
  # Move x and y to the first two columns
  col_idx <- match("x", names(collection))
  collection <- collection[, c(col_idx, (1:ncol(collection))[-col_idx])]
  col_idx <- match("y", names(collection))
  collection <- collection[, c(col_idx, (1:ncol(collection))[-col_idx])]  
  
  # Calculate 95% confidence interval for the slope
  # if sample size is > 30, use normal distribution
  # if not, use t distribution
  n <- model$n
  if(n < 30) {
    error <- qt(0.975,df=n-1)*collection$slope_std_err # t distribution
  } else {
    error <- qnorm(0.975)*collection$slope_std_err # normal distribution
  }
  collection$slope_CI95_lower <- collection$slope_estimate-error
  collection$slope_CI95_upper <- collection$slope_estimate+error
  
  return(collection)
  
}

