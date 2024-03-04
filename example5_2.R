# Fitting stratified Cox model using partial likelihood

strat.fit <- coxph(Surv(t_i, delta_i) ~ X1 + X2 + strata(istrata), data = dat.strat)
summary(strat.fit)
cox.zph(strat.fit)

