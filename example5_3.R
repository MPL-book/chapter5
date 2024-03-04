## A simulation example

# Write a function to generate interval censored survival data from two strata
gen_strat_interval_censoring <- function(nsample, nstrata, pi_E, prob_strata, a1, a2){
  
  # Set true coefficients
  beta <- c(1, -0.5)
  
  # Generate strata and covariates
  istrata <- sample(1:nstrata, size=nsample, prob = prob_strata, replace = T)
  X <- cbind(rbinom(nsample, 1, 0.5), rnorm(nsample, 1, 1))
  
  # Generate true event times
  y_i <- rep(0, nsample)
  for(i in 1:nsample){
    neg.log.u <- -log(runif(1))
    mu_term <- exp(X[i,]%*%beta)
    if(istrata[i] == 1){
      y_i[i] <- (neg.log.u/(mu_term))^(1/3)
    }else{
      y_i[i] <- (2 * (neg.log.u) / (mu_term))^(1/2)
    }
  }
  
  # Generate interval censoring times
  TL <- TR <- rep(0, nsample)
  
  for(i in 1:nsample){
    U_E <- runif(1)
    U_L <- runif(1,0,1)
    U_R <- runif(1,U_L,1)
    
    strata_i <- istrata[i]
    time <- y_i[i]
    
    event <- as.numeric(U_E < pi_E[strata_i])
    interval <- as.numeric(a1[strata_i]*U_L <= time & time <= a2[strata_i]*U_R & U_E >= pi_E[strata_i])
    right <- as.numeric(a2[strata_i]*U_R < time & U_E >= pi_E[strata_i])
    left <- as.numeric(time < a1[strata_i]*U_L & U_E >= pi_E[strata_i])
    
    TL[i] <- (time)*event + (a1[strata_i]*U_L)*interval + (a2[strata_i]*U_R)*right + (0)*left
    TR[i] <- (time)*event + (a2[strata_i]*U_R)*interval + (0)*right + (a1[strata_i]*U_L)*left
    
    if(TR[i] == 0){
      TR[i] <- Inf
    }
    if(TL[i] == 0){
      TL[i] <- -Inf
    }
  }
  
  data <- data.frame(c(1:nsample), TL, TR, X, istrata)
  return(data)
}

# Run simulations
save <- matrix(0, nrow = 200, ncol = 82)
save[,1] <- c(1:200)

source("coxph.R")
source("methods.R")
source("functions.R")

for(s in 1:200){
  
  # Generate data
  data <- gen_strat_interval_censoring(200, 2, pi_E = c(0.3, 0.3), c(0.50, 0.50), c(0.4, 0.9), c(1.1, 2))
  data.surv <- Surv(time = data$TL, time2 = data$TR, type = "interval2")
  
  # Save details of data
  save[s,2] <- sum(as.numeric(data.surv[data$istrata == 1,3] == 1))
  save[s,3] <- sum(as.numeric(data.surv[data$istrata == 1,3] == 0))
  save[s,4] <- sum(as.numeric(data.surv[data$istrata == 1,3] == 2))
  save[s,5] <- sum(as.numeric(data.surv[data$istrata == 1,3] == 3))
  
  save[s,6] <- sum(as.numeric(data.surv[data$istrata == 2,3] == 1))
  save[s,7] <- sum(as.numeric(data.surv[data$istrata == 2,3] == 0))
  save[s,8] <- sum(as.numeric(data.surv[data$istrata == 2,3] == 2))
  save[s,9] <- sum(as.numeric(data.surv[data$istrata == 2,3] == 3))
  
  # Fit MPL model
  try1 <- coxph_mpl(data.surv ~ data$X1 + data$X2, istrata = 6, data = data, control = coxph_mpl.control(basis = "msplines"))
  
  save[s,10] <- mpl.fit$coef$Beta[1]
  save[s,11] <- mpl.fit$coef$Beta[2]
  
  save[s,12] <- mpl.fit$se$Beta$se.Eta_M2HM2[1]
  save[s,13] <- mpl.fit$se$Beta$se.Eta_M2HM2[2]
  
  save[s,14:18] <- mpl.fit$s_obj[[1]]$knots_strat$Alpha
  
  save[s,26:31] <- mpl.fit$s_obj[[1]]$M_theta_m1
  save[s,37:42] <- mpl.fit$se$Theta$se.Eta_M2HM2[1:6]
  
  save[s,48:53] <- mpl.fit$s_obj[[2]]$M_theta_m1
  save[s,59:64] <- mpl.fit$se$Theta$se.Eta_M2HM2[7:12]
  
  ## Fit PL model with midpoint imputation
  data$midpoint_time <- 0
  data$midpoint_time[which(data.surv[,3] == 0)] <- data$TL[which(data.surv[,3] == 0)]
  data$midpoint_time[which(data.surv[,3] == 1)] <- data$TL[which(data.surv[,3] == 1)]
  data$midpoint_time[which(data.surv[,3] == 3)] <- data$TL[which(data.surv[,3] == 3)] + (data$TR[which(data.surv[,3] == 3)]  - data$TL[which(data.surv[,3] == 3)])/2
  
  data$midpoint_event <- 0
  data$midpoint_event[which(data.surv[,3] != 0)] <- 1
  
  data.surv.midpoint <- Surv(time = data$midpoint_time, event = data$midpoint_event)
  ph.fit <-  coxph(data.surv.midpoint ~ data$X1 + data$X2 + strata(data$istrata), data = data)
  
  save[s,70] <- ph.fit$coefficients[1]
  save[s,71] <- ph.fit$coefficients[2]
  
  save[s,72] <- sqrt(diag(ph.fit$var))[1]
  save[s,73] <- sqrt(diag(ph.fit$var))[2]
  
}
