## A pseudo melanoma recurrence example

# Write a function to generate data similar to the melanoma data
gen_strat_IC_melanoma <- function(nsample, nstrata, pi_E, prob_strata, a1, a2){
  
  # True parameter values
  beta <- c(0.25, 0.45, 0.55, 0, 0.45, 0.45, 0.25)
  
  # Strata membership
  istrata <- sample(1:nstrata, size=nsample, prob = prob_strata, replace = T)
  
  # Shared covariates
  Body_Site = sample(1:4, size=nsample, prob = c(0.25, 0.25, 0.25, 0.25), replace = T)
  X <- cbind(rbinom(nsample, 1, 0.5), rbinom(nsample, 1, 0.5), 
             rbinom(nsample, 1, 0.5), as.numeric(Body_Site == 1), 
             as.numeric(Body_Site == 2), as.numeric(Body_Site == 3), 
             as.numeric(Body_Site == 4))
  
  # Generate event times from four different distributions
  y_i <- rep(0, nsample)
  
  for(i in 1:nsample){
    neg.log.u <- -log(runif(1))
    mu_term <- exp(X[i,]%*%beta)
    
    if(istrata[i] == 1){
      y_i[i] = 2 * ( (neg.log.u/mu_term)^(1/1.5) )
    }else if(istrata[i] == 2){
      y_i[i] = ( (neg.log.u/(0.5 * mu_term)) )
    }else if(istrata[i] == 3){
      y_i[i] = 1.8 * ( (neg.log.u/mu_term) )
    }else if(istrata[i] == 4){
      y_i[i] = 2 * ( (neg.log.u/mu_term)^(1/1.1) )
    }
  }
  
  # Generate interval censoring times
  TL <- TR <- rep(0, nsample)
  
  for(i in 1:nsample){
    #uniform variables
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

# Generate a sample dataset
mel.data.sim <- gen_strat_IC_melanoma(5000, 4, c(0,0,0,0), c(0.5, 0.1, 0.3, 0.1), 
                                      c(0.8, 0.8, 0.7, 0.7), c(1.2, 1.3, 1.4, 1.5))

# calculate midpoint
mel.data.sim.mp <- mel.data.sim
id1 <- which(is.infinite(mel.data.sim.mp$TL))
id2 <- which(!is.infinite(mel.data.sim.mp$TR))
id3 <- which(!is.infinite(mel.data.sim.mp$TR))
id4 <- which(!is.infinite(mel.data.sim.mp$TR))
id5 <- which(!is.infinite(mel.data.sim.mp$TR))
mel.data.sim.mp$TL[id1] <- 0
mel.data.sim.mp$midpoint <- mel.data.sim.mp$TL
mel.data.sim.mp$midpoint[id2] <- mel.data.sim.mp$midpoint[id3] 
+ (mel.data.sim.mp$TR[id4] - mel.data.sim.mp$TL[id5])/2

# calculate event status
mel.data.sim.mp$event <- as.numeric(!is.infinite(
  mel.data.sim.mp$TR))

# create right censored Surv object
mp.mel.surv <- Surv(mel.data.sim.mp$midpoint, mel.data.sim.mp$event)

# use cox.zph() to check for violation of PH assumption by istrata
(cox.zph(coxph(mp.mel.surv ~ X1 + X2 + X3 + X5 + X6 + X7 + istrata, data = mel.data.sim.mp)))

# refit a stratified Cox model
ph.fit <- coxph(mp.mel.surv ~ X1 + X2 + X3 + X5 + X6 + X7 + strata(istrata), data = mel.data.sim.mp)
summary(ph.fit)

# fit using interval censoring MPL model

# create an interval censored Surv object
IC.data.surv <- Surv(mel.data.sim$TL, mel.data.sim$TR, type = "interval2")

# fit stratified model using MPL method
mpl.ic.fit <- coxph_mpl(IC.data.surv ~ X1 + X2 + X3 + X5 + X6 + X7, istrata = 11, data = mel.data.sim,
                        control = coxph_mpl.control(basis = "msplines", n.knots = c(8,2)))
summary(mpl.ic.fit)



