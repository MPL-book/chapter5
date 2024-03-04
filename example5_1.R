## Simulating data from a stratified Cox model

# write a function to generate data
gen_strat_data <- function(nsample, nstrata, prob_strata, beta){
  
  #beta = c(1, -0.5)
  
  istrata <- sample(1:nstrata, size=nsample, prob = prob_strata, replace = T)
  
  X <- cbind(rbinom(nsample, 1, 0.5), rnorm(nsample, 1, 1))
  
  y_i <- rep(0, nsample)
  
  for(i in 1:nsample){
    neg.log.u = -log(runif(1))
    mu_term = exp(X[i,]%*%beta)
    
    if(istrata[i] == 1){
      y_i[i] = (neg.log.u/(mu_term))^(1/3)
    }else{
      y_i[i] = (3 * (neg.log.u) / (mu_term))
    }
    
  }
  
  t_i = delta_i = rep(0, nsample)
  
  for(i in 1:nsample){
    c_i = rexp(1,1.5)
    delta_i[i] = as.numeric(y_i[i] < c_i)
    t_i[i] = min(c_i, y_i[i])
  }
  
  data = data.frame(c(1:nsample), t_i, delta_i, X, istrata)
  return(data)
  
}

dat.strat <- gen_strat_data(500, 2, c(0.5, 0.5), c(1, -0.5))
table(dat.strat$delta_i)

# log-log survival function
library(survival)
unstrat.fit <- coxph(Surv(t_i, delta_i) ~ X1 + X2 + factor(istrata), data = dat.strat)
summary(unstrat.fit)

km_strat <- survfit(Surv(t_i, delta_i) ~ istrata, data = dat.strat)
plot(km_strat, fun = "cloglog", xlab = "Time using log", ylab = "log-log survival",
       main = "log-log curves by istrat")

# include a time-dependent covariate
unstrat.fit.tvc <- coxph(Surv(t_i, delta_i) ~ X1 + X2 + factor(istrata) + tt(istrata), data = dat.strat,
                          tt = function(x, t,...) t*(x-1))
summary(unstrat.fit.tvc)

# test using schoenfeld residuals
cox.zph(unstrat.fit)
