testme <- function(){
  library(devtools)
  load_all()
  lung2=lung
  lung2$strata_var=sample(c(1,2), size = nrow(lung),prob = c(.5,.5), replace = T) #rep(1, nrow(lung))
  fit_mpl <- coxph_mpl(Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss, istrata = ncol(lung2), data = lung2, max.iter=c(500,7500,1e6)) # istrata = ncol(lung2) selects the last column of lung2, ie lung2$strata_var
  cat("Max diff with expected results:\t")
  cat(max(abs(fit_mpl$Beta - c(0.014574961, -0.510842685, -0.011757287, -0.001759738))))
  cat("\n\nBeta expected:")
  print(round(c(0.014574961, -0.510842685, -0.011757287, -0.001759738),4))
  cat("Beta estimate:")
  print(round(fit_mpl$Beta,4))
  # invisible()
  return(fit_mpl)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
#strata is the number of the col_variable in data used as stratification variable
coxph_mpl=function(formula,data,istrata,subset,na.action,control,...){
  #
  formula=update.formula(formula, paste("~ . +", names(data)[istrata]))
  mc = match.call(expand.dots = FALSE)
  m  = match(c("formula","data","subset","na.action") ,names(mc),0)
  mc = mc[c(1,m)]    
  if (m[1]==0){stop("A formula argument is required")}
  data.name = if(m[2]!=0){deparse(match.call()[[3]])}else{"-"}
  mc[[1]] = as.name("model.frame")
  mc$formula = if(missing(data)) {
    terms(formula)
  } else {
    terms(formula, data=data)
  }
  mf = eval(mc,parent.frame())
  if (any(is.na(mf))) stop("Missing observations in the model variables")
  if (nrow(mf) ==0) stop("No (non-missing) observations")
  mt = attr(mf,"terms")
  # Y
  y0    = model.extract(mf, "response")
  type = attr(y0, "type")
  if(!inherits(y0, "Surv")){stop("Response must be a survival object")}
  ##
  if (attr(y0,which = "type")=="right"){
    left=y0[,1]
    right=rep(NA, nrow(y0))
    icase = which(y0[,2]==1)
    right[icase] = y0[icase,1]
    y0 = Surv(left, right, type="interval2")
  } else if (type!="interval"){
    stop("\nPlease create the survival object using the option type='interval2' in the Surv function.\n")
  }
  ##
  # X
  X0   = model.matrix(mt, mf) #l, contrasts)
  #
  strata=as.factor(X0[,ncol(X0)])
  X0=X0[,-ncol(X0)]
  ##
  X0   = X0[,!apply(X0, 2, function(x) all(x==x[1])), drop=FALSE]

  lev.strata=levels(strata)
  ns = length(lev.strata)
  s_obj = vector("list", ns)
  ## Common to all strata (beta related)
  p           <- ncol(X0)    
  M_beta_p1  <- matrix(0,nrow=p,ncol=1)
  
  for (is in 1:ns){
    iis=which(strata==lev.strata[is])
    y=y0[iis]
    X=X0[iis,]
    ##
    t_i1        = y[,1L]  
    t_i2        = y[,2L] ## CHECK THIS - code not ready for inteval censored?
    n       = length(t_i1)
    ctype   = matrix(NA, nrow=n, ncol=4)
    colnames(ctype) = c("r","e","l","i")
    for(tw in 1:4){ctype[,tw] = y[,3L]==(tw-1)}
    n.ctype     = apply(ctype,2,sum)
    ctypeTF     = n.ctype>0    
    observed    = y[,3L]==1L
    n.obs       = sum(y[,3L]!=0)    
    # control arguments
    extraArgs <- list(...)
    if (length(extraArgs)) {
      controlargs <- names(formals(coxph_mpl.control)) 
      m <- pmatch(names(extraArgs), controlargs, nomatch=0L)
      if (any(m==0L))
        stop(gettextf("Argument(s) %s not matched", names(extraArgs)[m==0L]),
             domain = NA, call. = F)
    }    
    # if (missing(control)) control <- coxph_mpl.control(n.obs, ...)
    control <- coxph_mpl.control(n.obs, ...)
    
    # ties 
    t_i1.obs  = t_i1[observed]    
    ties     = duplicated(t_i1.obs)
    if(any(ties)){
      if(control$ties=="epsilon"){
        if(length(control$seed)>0){
          old <- .Random.seed
          on.exit({.Random.seed <<- old})
          set.seed(control$seed)
        }
        t_i1.obs[ties] = t_i1.obs[ties]+runif(sum(ties),-1e-11,1e-11)
        t_i1[observed] = t_i1.obs
      }else{    
        t_i1.obs = t_i1.obs[!ties]
        n.obs   = length(t_i1.obs)
      }
    }
    #X
    if(is.null(ncol(X))){
      X   = matrix(0,n,1)
      noX = TRUE
    }else{  noX = FALSE}
    ###################################
    ## WARNING: if correction is used, it has to be by strata
      # mean_j      <- apply(X, 2, mean)
      # XC          <- X - rep(mean_j, each=n) 
      XC <- X
      # knot sequence and psi matrices
      knots  <- knots_mpl(control,
                         c(t_i1[ctype[,"i"]],t_i2[ctype[,"i"]],t_i1[ctype[,"e"]]-1e-3,t_i1[ctype[,"e"]]+1e-3,
                           t_i1[ctype[,"r"]],t_i1[ctype[,"l"]]))
      
      ###                                    
      ### Estimation                         
      ###                                    
      
      
      K                        <- control$max.iter
      ##
      s_obj[[is]]$knots_strat       <- knots
      s_obj[[is]]$time        <- t_i1
      s_obj[[is]]$time_obs    <- t_i1.obs
      s_obj[[is]]$s_lambda     <- s_lambda     <- control$smooth
      s_obj[[is]]$s_kappa      <- control$kappa
      s_obj[[is]]$s_t1         <- knots$Alpha[1]
      s_obj[[is]]$s_tn         <- max(knots$Alpha)
      s_obj[[is]]$M_R_mm       <- penalty_mpl(control,knots)
      s_convlimit             <- control$tol    
      s_obj[[is]]$M_X_nop     <- XC[ctype[,2],,drop=F]
      s_obj[[is]]$M_tX_nop    <- t(s_obj[[is]]$M_X_nop)
      s_obj[[is]]$M_psi_nom   <- basis_mpl(t_i1,knots,control$basis,control$order,which=1)[ctype[,2],,drop=F] 
      s_obj[[is]]$M_tpsi_nom  <- t(s_obj[[is]]$M_psi_nom)
      s_obj[[is]]$M_Psi_nom   <- basis_mpl(t_i1,knots,control$basis,control$order,which=2)[ctype[,2],,drop=F] 
      s_obj[[is]]$M_tPsi_nom  <- t(s_obj[[is]]$M_Psi_nom)                                                                 
      s_obj[[is]]$M_X_nrp     <- XC[ctype[,1],,drop=F]
      s_obj[[is]]$M_tX_nrp    <- t(s_obj[[is]]$M_X_nrp)
      s_obj[[is]]$M_Psi_nrm   <- basis_mpl(t_i1,knots,control$basis,control$order,which=2)[ctype[,1],,drop=F]     
      s_obj[[is]]$M_tPsi_nrm  <- t(s_obj[[is]]$M_Psi_nrm)                                                     
      s_obj[[is]]$M_X_nlp     <- XC[ctype[,3],,drop=F]
      s_obj[[is]]$M_tX_nlp    <- t(s_obj[[is]]$M_X_nlp)
      s_obj[[is]]$M_Psi_nlm   <- basis_mpl(t_i1,knots,control$basis,control$order,which=2)[ctype[,3],,drop=F]                                     
      s_obj[[is]]$M_tPsi_nlm  <- t(s_obj[[is]]$M_Psi_nlm)                                         
      s_obj[[is]]$M_X_nip     <- XC[ctype[,4],,drop=F]
      s_obj[[is]]$M_tX_nip    <- t(s_obj[[is]]$M_X_nip)
      s_obj[[is]]$M_Psi1_nim  <- basis_mpl(t_i1,knots,control$basis,control$order,which=2)[ctype[,4],,drop=F]                                     
      s_obj[[is]]$M_Psi2_nim  <- basis_mpl(t_i2,knots,control$basis,control$order,which=2)[ctype[,4],,drop=F]                                     
      s_obj[[is]]$M_tPsi1_nim <- t(s_obj[[is]]$M_Psi1_nim)                                     
      s_obj[[is]]$M_tPsi2_nim <- t(s_obj[[is]]$M_Psi2_nim)                                     
      s_obj[[is]]$m            <- ncol(s_obj[[is]]$M_psi_nom)
      s_obj[[is]]$M_Rstar_ll   <- rbind(matrix(0,p,p+s_obj[[is]]$m),cbind(matrix(0,s_obj[[is]]$m,p),s_obj[[is]]$M_R_mm))    
      
      ###
      ### initialise
      ###
      M_beta_p1_OLD <- M_beta_p1 <- rep(0, p)
      s_obj[[is]]$M_theta_m1 <- matrix(1,nrow=s_obj[[is]]$m,ncol=1)
      s_obj[[is]]$s_df       <- -1    
      
      ## shortcuts
      ##
      s_obj[[is]]$M_mu_no1 <- exp(s_obj[[is]]$M_X_nop%*%M_beta_p1)
      s_obj[[is]]$M_mu_nr1   <- exp(s_obj[[is]]$M_X_nrp%*%M_beta_p1)
      s_obj[[is]]$M_mu_nl1   <- exp(s_obj[[is]]$M_X_nlp%*%M_beta_p1)
      s_obj[[is]]$M_mu_ni1   <- exp(s_obj[[is]]$M_X_nip%*%M_beta_p1)
      ##
      s_obj[[is]]$M_h0_no1    <- s_obj[[is]]$M_psi_nom%*%s_obj[[is]]$M_theta_m1
      s_obj[[is]]$M_H0_no1    <- s_obj[[is]]$M_Psi_nom%*%s_obj[[is]]$M_theta_m1
      s_obj[[is]]$M_H0_nr1    <- s_obj[[is]]$M_Psi_nrm%*%s_obj[[is]]$M_theta_m1
      s_obj[[is]]$M_H0_nl1    <- s_obj[[is]]$M_Psi_nlm%*%s_obj[[is]]$M_theta_m1
      s_obj[[is]]$M_H01_ni1   <- s_obj[[is]]$M_Psi1_nim%*%s_obj[[is]]$M_theta_m1
      s_obj[[is]]$M_H02_ni1   <- s_obj[[is]]$M_Psi2_nim%*%s_obj[[is]]$M_theta_m1
      s_obj[[is]]$Rtheta      <- s_obj[[is]]$M_R_mm%*%s_obj[[is]]$M_theta_m1
      ##
      s_obj[[is]]$thetaRtheta <- thetaRtheta <- t(s_obj[[is]]$M_theta_m1)%*%s_obj[[is]]$Rtheta
      s_obj[[is]]$TwoLRtheta  <- s_obj[[is]]$s_lambda*2*s_obj[[is]]$Rtheta            
      ##
      s_obj[[is]]$M_H_no1  <- s_obj[[is]]$M_H0_no1 * as.numeric(s_obj[[is]]$M_mu_no1)
      ##
      s_obj[[is]]$M_H_nr1  <- s_obj[[is]]$M_H0_nr1 * s_obj[[is]]$M_mu_nr1
      s_obj[[is]]$M_H_nl1  <- s_obj[[is]]$M_H0_nl1 * s_obj[[is]]$M_mu_nl1
      s_obj[[is]]$M_H1_ni1 <- s_obj[[is]]$M_H01_ni1* s_obj[[is]]$M_mu_ni1
      s_obj[[is]]$M_H2_ni1 <- s_obj[[is]]$M_H02_ni1* s_obj[[is]]$M_mu_ni1
      s_obj[[is]]$M_S_no1  <- exp(-s_obj[[is]]$M_H_no1)
      ##
      s_obj[[is]]$M_S_nl1  <- exp(-s_obj[[is]]$M_H_nl1)
      s_obj[[is]]$M_S1_ni1 <- exp(-s_obj[[is]]$M_H1_ni1)
      s_obj[[is]]$M_S2_ni1 <- exp(-s_obj[[is]]$M_H2_ni1)
      # avoid division by 0
      s_obj[[is]]$M_S_nl1[s_obj[[is]]$M_S_nl1==1]   <- 1-control$epsilon[1]
      s_obj[[is]]$M_S1_ni1[s_obj[[is]]$M_S1_ni1==1] <- 1-control$epsilon[1]
      s_obj[[is]]$M_S2_ni1[s_obj[[is]]$M_S2_ni1==1] <- 1-control$epsilon[1]
      ##
      s_obj[[is]]$M_S1mS2_ni1 <- s_obj[[is]]$M_S1_ni1-s_obj[[is]]$M_S2_ni1  
      s_obj[[is]]$M_S1mS2_ni1[s_obj[[is]]$M_S1mS2_ni1<control$epsilon[2]] <- control$epsilon[2]
  }
  ### outer loop
  ctype   = matrix(NA, nrow=nrow(X0), ncol=4)
  colnames(ctype) = c("r","e","l","i")
  for(tw in 1:4){ctype[,tw] = y0[,3L]==(tw-1)}
  
  
  
  full.iter = 0
  thr <- 1-control$epsilon[1]
  thr2 <- control$epsilon[2]
  K         = ifelse(control$max.iter[1]>1,control$max.iter[2],control$max.iter[3])
  for(iter in 1:control$max.iter[1]){
    # cat(iter, "/", control$max.iter[1], "\n", sep='')
    # loglik
    s_lik <- s_lik_OLD <- sum(sapply(s_obj, compute_likelihood))
    ### inner loop
    for(k in 1:K){
      
      ## update betas
      s_omega       <- 1
      M_beta_p1_OLD <- M_beta_p1 
      # s_lik_OLD     <- s_lik 
      M_gradbeta_p1 <- gradBeta(s_obj)
      M_hessbeta_p1 <- hessianBeta(s_obj)      
      M_stepbeta_p1 <- chol2inv(chol(M_hessbeta_p1))%*%M_gradbeta_p1
      #
      M_beta_p1     <- M_beta_p1_OLD+s_omega*M_stepbeta_p1
      s_obj <- update_s_obj_beta(s_obj, M_beta_p1, thr, thr2)
      s_lik <- sum(sapply(s_obj, compute_likelihood))
      ## if likelihood decreases
      if(s_lik<s_lik_OLD){
        i       = 0            
        while(s_lik<s_lik_OLD){
          # update value of omega
          s_omega <- update_omega(s_omega, s_obj)
          #
          M_beta_p1 = M_beta_p1_OLD+s_omega*M_stepbeta_p1
          s_obj <- update_s_obj_beta(s_obj, M_beta_p1, thr, thr2)
          s_lik <- sum(sapply(s_obj, compute_likelihood))
          #
          i = i+1
          if(i>500){warning("max no of iteration reached for beta estimate");break}                
        }
      }
      ## end beta
      ############################
      ## update thetas
      s_omega = 1
      for (is in 1:ns){
        s_obj[[is]]$M_theta_m1_OLD = s_obj[[is]]$M_theta_m1 
      }
        #
      s_obj <- update_s_obj_steptheta(s_obj, thr, thr2, s_nu)
      s_obj <- update_s_obj_theta(s_obj, thr, thr2, s_omega)
      # loglik
      s_lik_OLD      = s_lik
      s_lik <- sum(sapply(s_obj, compute_likelihood))
      ## if likelihood decreases
      if(s_lik<s_lik_OLD){
        i       = 0
        while(s_lik<s_lik_OLD){
          # update value of omega
          s_omega <- update_omega(s_omega, s_obj)
          s_obj <- update_s_obj_theta(s_obj, thr, thr2, s_omega)
          # loglik
          s_lik <- sum(sapply(s_obj, compute_likelihood))
          #
          i = i+1
          if(i>500){break}                
        }    
      }
      # if(all(c(abs(M_beta_p1-M_beta_p1_OLD))< 0.1)){break}# s_convlimit)){break}
      if(all(c(abs(M_beta_p1-M_beta_p1_OLD),sapply(s_obj, function(is)abs(is$M_theta_m1-is$M_theta_m1_OLD)))<s_convlimit)){break}
      if(control$max.iter[1]==1){
        if(any(k==seq(0,control$max.iter[3],control$max.iter[2]))){
            control$epsilon[2] = control$epsilon[2]*10
        }
      }
    }
    # print(M_beta_p1)
    # H matrix 
    # H = HRinv = matrix(0,p+m,p+m)
    
    
    M_hessbeta_components <- lapply(s_obj, function(ss){
      # H matrix 
      H_sep = matrix(0,p+ss$m,p+ss$m)
      H_sep[1:p,1:p] =
        ss$M_tX_nop%*%diag(c(ss$M_H_no1))%*%ss$M_X_nop+
        ss$M_tX_nrp%*%diag(c(ss$M_H_nr1))%*%ss$M_X_nrp+
        ss$M_tX_nlp%*%diag(c((ss$M_S_nl1/(1-ss$M_S_nl1)^2*ss$M_H_nl1^2-ss$M_S_nl1/(1-ss$M_S_nl1)*ss$M_H_nl1)))%*%ss$M_X_nlp+
        ss$M_tX_nip%*%diag(c((ss$M_S1_ni1*ss$M_S2_ni1/ss$M_S1mS2_ni1^2*(ss$M_H2_ni1-ss$M_H1_ni1)^2+(ss$M_S1_ni1*ss$M_H1_ni1-ss$M_S2_ni1*ss$M_H2_ni1)/ss$M_S1mS2_ni1)))%*%ss$M_X_nip
      H_sep[1:p,(p+1):(p+ss$m)] =
        ss$M_tX_nop%*%diag(c(ss$M_mu_no1))%*%ss$M_Psi_nom+
        ss$M_tX_nrp%*%diag(c(ss$M_mu_nr1))%*%ss$M_Psi_nrm+
        ss$M_tX_nlp%*%diag(c((ss$M_S_nl1/(1-ss$M_S_nl1)^2*ss$M_H_nl1-ss$M_S_nl1/(1-ss$M_S_nl1))*ss$M_mu_nl1))%*%ss$M_Psi_nlm+
        ss$M_tX_nip%*%diag(c(ss$M_S1_ni1*ss$M_S2_ni1/ss$M_S1mS2_ni1^2*(ss$M_H2_ni1-ss$M_H1_ni1)*ss$M_mu_ni1))%*%(ss$M_Psi2_nim-ss$M_Psi1_nim)+
        ss$M_tX_nip%*%diag(c(ss$M_S1_ni1/ss$M_S1mS2_ni1*ss$M_mu_ni1))%*%ss$M_Psi1_nim-
        ss$M_tX_nip%*%diag(c(ss$M_S2_ni1/ss$M_S1mS2_ni1*ss$M_mu_ni1))%*%ss$M_Psi2_nim
      H_sep[(p+1):(p+ss$m),1:p] = t(H_sep[1:p,(p+1):(p+ss$m)])
      H_sep[(p+1):(p+ss$m),(p+1):(p+ss$m)] =
        ss$M_tpsi_nom%*%diag(c(1/ss$M_h0_no1^2))%*%ss$M_psi_nom+
        ss$M_tPsi_nlm%*%diag(c(ss$M_S_nl1/(1-ss$M_S_nl1)^2*ss$M_mu_nl1^2))%*%ss$M_Psi_nlm+
      (ss$M_tPsi2_nim-ss$M_tPsi1_nim)%*%diag(c(ss$M_S1_ni1*ss$M_S2_ni1/ss$M_S1mS2_ni1^2*ss$M_mu_ni1^2))%*%(ss$M_Psi2_nim-ss$M_Psi1_nim)
      H_sep
    })
    M_hessbeta <- Reduce(`+`, M_hessbeta_components)
    for (is in 1:ns){
      HRinv = matrix(0, p+s_obj[[is]]$m, p+s_obj[[is]]$m)
      H_sep_smooth = M_hessbeta_components[[is]]
      #
      s_obj[[is]]$s_lambda_old   = s_obj[[is]]$s_lambda
      s_obj[[is]]$s_df_old       = s_obj[[is]]$s_df
      s_obj[[is]]$s_sigma2_old   = 1/(2*s_obj[[is]]$s_lambda_old)
      #pos            = c(if(noX){FALSE}else{rep(TRUE,p)},M_theta_m1>control$min.theta)
      pos            = c(if(noX){FALSE}else{rep(TRUE,p)},(s_obj[[is]]$M_theta_m1>control$min.theta & apply(H_sep_smooth[(p+1):(p+s_obj[[is]]$m),(p+1):(p+s_obj[[is]]$m)],2,sum)>0))
      #MM: (G+Q)^-1
      # temp           = try(chol2inv(chol(H[pos,pos]+(1/s_obj[[is]]$s_sigma2_old)*s_obj[[is]]$M_Rstar_ll[pos,pos])),silent=T)
      temp           = try(solve(H_sep_smooth[pos,pos]+(1/s_obj[[is]]$s_sigma2_old)*s_obj[[is]]$M_Rstar_ll[pos,pos]),silent=T)
      #if(class(temp)!="try-error"&!any(is.infinite(temp))){
        HRinv[pos,pos]=temp
      #}else{
      #  HRinv[pos,pos]=MASS::ginv(H_sep_smooth[pos,pos])
      #}
      
      s_obj[[is]]$s_df           = s_obj[[is]]$m-sum(diag(HRinv%*%s_obj[[is]]$M_Rstar_ll))/s_obj[[is]]$s_sigma2_old
      s_obj[[is]]$s_sigma2       = c(t(s_obj[[is]]$M_theta_m1)%*%s_obj[[is]]$M_R_mm%*%s_obj[[is]]$M_theta_m1/s_obj[[is]]$s_df)
      s_obj[[is]]$s_lambda       = 1/(2*s_obj[[is]]$s_sigma2)
      s_obj[[is]]$TwoLRtheta     = s_obj[[is]]$s_lambda*2*s_obj[[is]]$Rtheta
      s_obj[[is]]$df_conv = (abs(s_obj[[is]]$s_df-s_obj[[is]]$s_df_old)<(control$tol))
    }
    df_conv=T; for (is in 1:ns) df_conv=df_conv * s_obj[[is]]$df_conv; df_conv = as.logical(df_conv)
    full.iter      = full.iter+k
    
    # cat(full.iter, control$max.iter[3], k, control$max.iter[2], df_conv, "\n")
    cat("Ext iters: ", iter, "/", control$max.iter[1], "\t", k, " int iters", "\tllik = ", s_lik, "\n", sep='')
    
    if((full.iter/iter)>(control$max.iter[2]*.975)){control$epsilon[2] = control$epsilon[2]*10}
    if(full.iter>control$max.iter[3]){break}
    if((k<control$max.iter[2]) & df_conv) {break}
  }
  s_lambda        = control$smooth #= s_lambda_old
  #for (ss in s_obj){
  #  s_correction    = c(exp(-ss$mean_j%*%M_beta_p1)) 
  #  ss$M_thetatilde_m1 = ss$M_theta_m1
  #  ss$M_theta_m1      = ss$M_theta_m1*s_correction
  #}
  
  
  ###                                    
  ### Inference                          
  ###                                    
  # # M_corr_ll = cbind(rbind(diag(rep(1,p)),s_correction*matrix(rep(M_thetatilde_m1,p)*rep(-mean_j,each=m),ncol=p)),
  #                   # rbind(matrix(0,ncol=m,nrow=p),diag(rep(s_correction,m))))
  
  m_full = 0
  m_each = c(0)
  pos_full = rep(TRUE, p)
  for(is in 1:ns){
    m_full = m_full + s_obj[[is]]$m
    m_each = c(m_each, s_obj[[is]]$m)
    pos_full = c(pos_full, c(if(noX){FALSE}else{rep(TRUE,p)},(s_obj[[is]]$M_theta_m1>control$min.theta & apply(M_hessbeta_components[[is]][(p+1):(p+s_obj[[is]]$m),(p+1):(p+s_obj[[is]]$m)],2,sum)>0))[-c(1:p)])
  }
  
  H = matrix(0,p+m_full,p+m_full)
  
  for(is in 1:ns){
    H[1:p,1:p] = H[1:p,1:p] + 
      s_obj[[is]]$M_tX_nop%*%diag(c(s_obj[[is]]$M_H_no1))%*%s_obj[[is]]$M_X_nop+
      s_obj[[is]]$M_tX_nrp%*%diag(c(s_obj[[is]]$M_H_nr1))%*%s_obj[[is]]$M_X_nrp+
      s_obj[[is]]$M_tX_nlp%*%diag(c((s_obj[[is]]$M_S_nl1/(1-s_obj[[is]]$M_S_nl1)^2*s_obj[[is]]$M_H_nl1^2-s_obj[[is]]$M_S_nl1/(1-s_obj[[is]]$M_S_nl1)*s_obj[[is]]$M_H_nl1)))%*%s_obj[[is]]$M_X_nlp+
      s_obj[[is]]$M_tX_nip%*%diag(c((s_obj[[is]]$M_S1_ni1*s_obj[[is]]$M_S2_ni1/s_obj[[is]]$M_S1mS2_ni1^2*(s_obj[[is]]$M_H2_ni1-s_obj[[is]]$M_H1_ni1)^2+(s_obj[[is]]$M_S1_ni1*s_obj[[is]]$M_H1_ni1-s_obj[[is]]$M_S2_ni1*s_obj[[is]]$M_H2_ni1)/s_obj[[is]]$M_S1mS2_ni1)))%*%s_obj[[is]]$M_X_nip
    
    H[1:p,(p+sum(m_each[1:is])+1):(p+sum(m_each[(1:is) + 1]))] =
      s_obj[[is]]$M_tX_nop%*%diag(c(s_obj[[is]]$M_mu_no1))%*%s_obj[[is]]$M_Psi_nom+
      s_obj[[is]]$M_tX_nrp%*%diag(c(s_obj[[is]]$M_mu_nr1))%*%s_obj[[is]]$M_Psi_nrm+
      s_obj[[is]]$M_tX_nlp%*%diag(c((s_obj[[is]]$M_S_nl1/(1-s_obj[[is]]$M_S_nl1)^2*s_obj[[is]]$M_H_nl1-s_obj[[is]]$M_S_nl1/(1-s_obj[[is]]$M_S_nl1))*s_obj[[is]]$M_mu_nl1))%*%s_obj[[is]]$M_Psi_nlm+
      s_obj[[is]]$M_tX_nip%*%diag(c(s_obj[[is]]$M_S1_ni1*s_obj[[is]]$M_S2_ni1/s_obj[[is]]$M_S1mS2_ni1^2*(s_obj[[is]]$M_H2_ni1-s_obj[[is]]$M_H1_ni1)*s_obj[[is]]$M_mu_ni1))%*%(s_obj[[is]]$M_Psi2_nim-s_obj[[is]]$M_Psi1_nim)+
      s_obj[[is]]$M_tX_nip%*%diag(c(s_obj[[is]]$M_S1_ni1/s_obj[[is]]$M_S1mS2_ni1*s_obj[[is]]$M_mu_ni1))%*%s_obj[[is]]$M_Psi1_nim-
      s_obj[[is]]$M_tX_nip%*%diag(c(s_obj[[is]]$M_S2_ni1/s_obj[[is]]$M_S1mS2_ni1*s_obj[[is]]$M_mu_ni1))%*%s_obj[[is]]$M_Psi2_nim
    
    H[(p+sum(m_each[1:is])+1):(p+sum(m_each[(1:is) + 1])),1:p] = t(H[1:p,(p+sum(m_each[1:is])+1):(p+sum(m_each[(1:is) + 1]))])
    
    H[(p+sum(m_each[1:is])+1):(p+sum(m_each[(1:is) + 1])),(p+sum(m_each[1:is])+1):(p+sum(m_each[(1:is) + 1]))] =
      s_obj[[is]]$M_tpsi_nom%*%diag(c(1/s_obj[[is]]$M_h0_no1^2))%*%s_obj[[is]]$M_psi_nom+
      s_obj[[is]]$M_tPsi_nlm%*%diag(c(s_obj[[is]]$M_S_nl1/(1-s_obj[[is]]$M_S_nl1)^2*s_obj[[is]]$M_mu_nl1^2))%*%s_obj[[is]]$M_Psi_nlm+
      (s_obj[[is]]$M_tPsi2_nim-s_obj[[is]]$M_tPsi1_nim)%*%diag(c(s_obj[[is]]$M_S1_ni1*s_obj[[is]]$M_S2_ni1/s_obj[[is]]$M_S1mS2_ni1^2*s_obj[[is]]$M_mu_ni1^2))%*%(s_obj[[is]]$M_Psi2_nim-s_obj[[is]]$M_Psi1_nim)
    
  }
  
  s_correction = 1
  
  #M_corr_ll = cbind(rbind(diag(rep(1,p)),s_correction*matrix(rep(M_thetatilde_m1,p)*rep(-M_beta_p1,each=m),ncol=p)),
  #                  rbind(matrix(0,ncol=m,nrow=p),diag(rep(s_correction,m))))
  #M_corr_ll[!pos,] = 0
  

  M_corr_ll = diag(1, p+m_full)
# 
  M_Rstar_ll_full = matrix(0, p+m_full, p+m_full)
  for(is in 1:ns){
    M_Rstar_ll_full[(p + (is-1)*6 + 1:s_obj[[is]]$m), (p + (is-1)*6 + 1:s_obj[[is]]$m)] = s_obj[[is]]$M_Rstar_ll[-c(1:p), -c(1:p)]
    
  }
  
  # 
  M_2    = H+2*s_lambda*M_Rstar_ll_full
  #Q = matrix(0,nrow(X0),p+m_full)
  
  #for(is in 1:ns){
  #  if(ctypeTF[1]){Q[which(ctype[,1] == TRUE & strata == is),1:p] = as.matrix(rep(-s_obj[[is]]$M_H_nr1,p)*s_obj[[is]]$M_X_nrp)}
  #  if(ctypeTF[2]){Q[which(ctype[,2] == TRUE & strata == is),1:p] = rep((1-s_obj[[is]]$M_H_no1),p)*s_obj[[is]]$M_X_nop}
  #  if(ctypeTF[3]){Q[which(ctype[,3] == TRUE & strata == is),1:p] = as.numeric(rep(s_obj[[is]]$M_S_nl1*s_obj[[is]]$M_H_nl1/(1-s_obj[[is]]$M_S_nl1),p))*s_obj[[is]]$M_X_nlp}
  #  if(ctypeTF[4]){Q[which(ctype[,4] == TRUE & strata == is),1:p] = rep((s_obj[[is]]$M_H2_ni1*s_obj[[is]]$M_S2_ni1-s_obj[[is]]$M_H1_ni1*s_obj[[is]]$M_S1_ni1)/s_obj[[is]]$M_S1mS2_ni1,p)*s_obj[[is]]$M_X_nip}
    
  #  if(ctypeTF[1]){Q[which(ctype[,1] == TRUE & strata == is),(p+sum(m_each[1:is])+1):(p+sum(m_each[(1:is) + 1]))] = - as.numeric(s_obj[[is]]$M_mu_nr1)*s_obj[[is]]$M_Psi_nrm }
  #  if(ctypeTF[2]){Q[which(ctype[,2] == TRUE & strata == is),(p+sum(m_each[1:is])+1):(p+sum(m_each[(1:is) + 1]))] = as.numeric(1/s_obj[[is]]$M_h0_no1)*s_obj[[is]]$M_psi_nom -  as.numeric(s_obj[[is]]$M_mu_no1)*s_obj[[is]]$M_Psi_nom}
  #  if(ctypeTF[3]){Q[which(ctype[,3] == TRUE & strata == is),(p+sum(m_each[1:is])+1):(p+sum(m_each[(1:is) + 1]))] = as.numeric(s_obj[[is]]$M_mu_nl1)*s_obj[[is]]$M_Psi_nlm*as.numeric(s_obj[[is]]$M_S_nl1)/as.numeric(1 - s_obj[[is]]$M_S_nl1)}
  #  if(ctypeTF[4]){Q[which(ctype[,4] == TRUE & strata == is),(p+sum(m_each[1:is])+1):(p+sum(m_each[(1:is) + 1]))] = as.numeric(s_obj[[is]]$M_mu_ni1)*(as.numeric(s_obj[[is]]$M_S1_ni1)*s_obj[[is]]$M_Psi1_nim - as.numeric(s_obj[[is]]$M_S2_ni1)*s_obj[[is]]$M_Psi2_nim)/as.numeric(s_obj[[is]]$M_S1_ni1 - s_obj[[is]]$M_S2_ni1) }
    
    
  
  #}
  
  #Sp = Q-matrix(rep(c(rep(0,p),s_obj[[1]]$TwoLRtheta, s_obj[[2]]$TwoLRtheta),nrow(X0)),nrow(X0),byrow=T)/nrow(X0)
  #Q = t(Sp)%*%Sp
  # #pos   = c(if(noX){FALSE}else{rep(TRUE,p)},M_theta_m1>control$min.theta)
  #pos            = c(if(noX){FALSE}else{rep(TRUE,p)},(M_theta_m1>control$min.theta & apply(H[(p+1):(p+m),(p+1):(p+m)],2,sum)>0))
  Minv_1 = Minv_2 = Hinv = matrix(0,p+m_full,p+m_full)                        
  temp = try(solve(H[pos_full,pos_full]))
  #if(class(temp)!="try-error"){
    Minv_2[pos_full,pos_full] = temp        
    #cov_NuNu_M2QM2  = M_corr_ll%*%(Minv_2%*%Q%*%Minv_2)%*%t(M_corr_ll)
    cov_NuNu_M2HM2  = M_corr_ll%*%(Minv_2%*%H%*%Minv_2)%*%t(M_corr_ll)
    #se.Eta_M2QM2    = sqrt(diag(cov_NuNu_M2QM2))
    se.Eta_M2HM2    = sqrt(diag(cov_NuNu_M2HM2))            
  #}else{
  #  cov_NuNu_M2HM2 = matrix(NA,p+m_full,p+m_full)
  #  se.Eta_M2HM2   = rep(NA,p+m_full)
  #}
  
  temp = try(solve(H[pos_full,pos_full]))
  #if(class(temp)!="try-error"){
    Hinv[pos_full,pos_full] = temp
    cov_NuNu_H    = M_corr_ll%*%Hinv%*%t(M_corr_ll) 
    se.Eta_H      = sqrt(diag(cov_NuNu_H))
  #}else{
  #  cov_NuNu_H  = matrix(NA,p+m_full,p+m_full)
  #  se.Eta_H    = rep(NA,p+m_full)
  #}
  #MQ_2    = Q+2*s_lambda*M_Rstar_ll_full
  #temp2 = try(chol2inv(chol(MQ_2[pos_full,pos_full])),silent=T) 
  #if(class(temp2)!="try-error"){
  #  Minv_1[pos_full,pos_full] = temp2        
  #  cov_NuNu_M2QM2  = M_corr_ll%*%(Minv_1%*%Q%*%Minv_1)%*%t(M_corr_ll)
    #cov_NuNu_M2HM2  = M_corr_ll%*%(Minv_2%*%H%*%Minv_2)%*%t(M_corr_ll)
  #  se.Eta_M2QM2    = sqrt(diag(cov_NuNu_M2QM2))
    #se.Eta_M2HM2    = sqrt(diag(cov_NuNu_M2HM2))            
  #}else{
  #  cov_NuNu_M2QM2  = cov_NuNu_M2HM2 = matrix(NA,p+m_full,p+m_full)
  #  se.Eta_M2QM2    = se.Eta_M2HM2   = rep(NA,p+m_full)
  #}
  
  Theta_full = NULL
  
  for(is in 1:ns){
    Theta_full = c(Theta_full, s_obj[[is]]$M_theta_m1)
    
  }
  
  
  #mx.seNu.l5=as.data.frame(cbind(se.Eta_M2QM2,se.Eta_M2HM2,se.Eta_H))
  mx.seNu.l5=as.data.frame(cbind(se.Eta_M2HM2,se.Eta_H))
  #mx.seNu.l5=as.data.frame(cbind(se.Eta_M2HM2,se.Eta_H))
  #colnames(mx.seNu.l5)=c("M2HM2","H")
  #rownames(mx.seNu.l5)[(p+1):(p+m)] = paste("Theta",1:m,sep="")        
  #rownames(mx.seNu.l5)[(1:p)] = paste("Beta",1:p,sep="")
  fit         = list(coef=list(Beta=c(M_beta_p1),Theta=c(Theta_full)), iter = c(iter,full.iter))
  #fit         = list(Beta=c(M_beta_p1), iter = c(iter,full.iter))
  fit$s_obj   = s_obj
  fit$se      = list(Beta=mx.seNu.l5[1:p,],Theta=mx.seNu.l5[-c(1:p),])
  #fit$covar   = list(M2QM2=cov_NuNu_M2QM2,M2HM2=cov_NuNu_M2HM2,H=cov_NuNu_H)
  fit$covar   = list(M2HM2=cov_NuNu_M2HM2,H=cov_NuNu_H)
  fit$knots   = knots
  fit$control = control
  fit$call    = match.call()
  fit$dim     = list(n = n, n.obs = sum(observed), n.ties = sum(ties), p = p, m = knots$m)
  fit$data    = list(time = y, censoring=y[,3L], X = X, name = data.name)# list(name = data.name)#
  fit$df      = s_obj[[1]]$s_df #FIX
  fit$ploglik = s_lik
  # fit$loglik  = s_lik+s_lambda*thetaRtheta
  class(fit)  = "coxph_mpl"
  fit
}







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
coxph_mpl.control <- function(n.obs=NULL, basis = "msplines", smooth = NULL, max.iter=c(1.5e+2,10000,1e+6), tol=1e-6, 
                              n.knots = NULL, n.events_basis = NULL, range.quant = c(0.075,.9),
                              cover.sigma.quant = .25, cover.sigma.fixed=.25, min.theta = 1e-10,
                              penalty = 2L, order = 3L, kappa = 1/.6, epsilon = c(1e-16,1e-10), ties = "epsilon", seed = NULL){
  basis        = basis.name_mpl(basis)
  max.iter     = c(ifelse(is.null(smooth),ifelse(max.iter[1]>0,as.integer(max.iter[1]),1.5e+2),1L),
                   ifelse(max.iter[2]>0,as.integer(max.iter[2]),7.5e+4),
                   ifelse(length(max.iter)==2,1e+6,
                          ifelse(max.iter[3]>ifelse(max.iter[2]>0,as.integer(max.iter[2]),7.5e+4),
                                 as.integer(max.iter[3]),1e+6)))    
  tol          = ifelse(tol>0 & tol<1,tol,1e-7)    
  order        = ifelse(order>0 & order<6,as.integer(order),3L)    
  min.theta    = ifelse(min.theta>0 & min.theta<1e-3,min.theta,1e-10)    
  penalty      = penalty.order_mpl(penalty,basis,order)
  kappa        = ifelse(kappa>1, kappa, 1/.6)
  cover.sigma.quant  = ifelse(cover.sigma.quant>0 & cover.sigma.quant<0.4,cover.sigma.quant,.75)
  cover.sigma.fixed  = ifelse(cover.sigma.fixed>0 & cover.sigma.fixed<0.4,cover.sigma.fixed,.75)    
  if(all(range.quant<=1) & all(range.quant>=0) & length(range.quant)==2){
    range.quant = range.quant[order(range.quant)]
  }else{range.quant = c(0.075,.9)}
  if(is.null(n.knots)|sum(n.knots)<3|length(n.knots)!=2){
    n.knots    = if(basis!='uniform' & basis!='msplines'){c(0,20)}else{c(0,3)}
  }
  if(!is.null(n.events_basis)){
    n.events_basis = ifelse(n.events_basis<1|n.events_basis>floor(n.obs/2),
                            max(round(3.5*log(n.obs)-7.5),1L),round(n.events_basis))
  }else{n.events_basis = max(round(3.5*log(n.obs)-7.5),1L)}    
  if(!is.null(smooth)){
    smooth = ifelse(smooth<0,0,smooth)
  }else{smooth=0}
  out = list(basis = basis, smooth = smooth, max.iter = max.iter, tol = tol,
             order = order, penalty = penalty, n.knots = n.knots, range.quant = range.quant,
             cover.sigma.quant = cover.sigma.quant, cover.sigma.fixed = cover.sigma.fixed,
             n.events_basis = as.integer(n.events_basis), min.theta = min.theta, ties = ties,
             seed = as.integer(seed), kappa = kappa, epsilon = epsilon)
  class(out) = "coxph_mpl.control"
  out
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
basis.name_mpl <- function(k){
  if(k == "discr"| k == "discretized" | k == "discretised" | k == "unif" | k == "uniform"){"uniform"
  }else{if(k == "m" | k == "msplines" | k == "mspline"){"msplines"
  }else{if(k == "gauss" | k == "gaussian"){"gaussian"
  }else{if(k == "epa" | k == "epanechikov"){"epanechikov"
  }else{stop("Unkown basis choice", call. = FALSE)}}}}
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
penalty.order_mpl <- function(p,basis,order){
  p = as.integer(p)
  switch(basis,
         'uniform'  = ifelse(p>0 & p<3,p,2),
         'gaussian' = ifelse(p>0 & p<3,p,2),
         'msplines' = order-1,
         'epa'      = 2)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
knots_mpl=function(control,events){
  n.events  = length(events)
  range = range(events)  
  if(control$n.knots[2]==0){
    Alpha    = quantile(events,seq(0.1,0.9,length.out=(control$n.knots[1]+2)))
  }else{
    Alpha1   = quantile(events,seq(0,control$range.quant[2],length.out=(control$n.knots[1]+1)))    
    Alpha2   = seq(quantile(events,control$range.quant[2]),range(events)[2],length=control$n.knots[2]+2)
    Alpha    = c(Alpha1,Alpha2[-1])
  }
  n.Alpha  = length(Alpha)
  if(control$basis=="gaussian"){
    Sigma = Delta = rep(0,n.Alpha)
    for(aw in 1:n.Alpha){
      if(aw>1 & aw<(n.Alpha-control$n.knots[2])){
        while(sum(events>(Alpha[aw]-2*Sigma[aw])&events<(Alpha[aw]+2*Sigma[aw]))<(n.events*control$cover.sigma.quant)){
          Sigma[aw] = Sigma[aw] + 0.001}
      }else{Sigma[aw] = control$cover.sigma.fixed*(Alpha[n.Alpha]-Alpha[1])/3}
      Delta[aw]= pnorm((range[2]-Alpha[aw])/Sigma[aw])-
        pnorm((range[1]-control$epsilon[1]-Alpha[aw])/Sigma[aw])
    }
    list(m=n.Alpha, Alpha=Alpha, Sigma=Sigma, Delta=Delta)
  }else{if(control$basis=="msplines"|control$basis=="epanechikov"){
    m = n.Alpha+control$order-2
    list(m=m, Alpha=Alpha, Delta=rep(1,m))
  }else{if(control$basis=="uniform"){
    m          = length(Alpha)-1
    Delta      = Alpha[2L:(m+1L)]-Alpha[1L:m]
    list(m=m,Alpha=Alpha,Delta=Delta)
  }else{stop("Unkown basis choice")}}}
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
basis_mpl = function(x,knots,basis,order,which=c(1,2)){
  which.matrix = rep(T,2)
  which.matrix[-which]=FALSE    
  n        = length(x)
  Alpha    = knots$Alpha
  Delta    = knots$Delta
  n.Alpha  = length(Alpha)
  m        = ifelse(basis=="msplines"|basis=="epanechikov",n.Alpha+order-2,knots$m)
  M_Psi_nm = M_psi_nm = matrix(0,n,m)
  ##
  if(basis=="uniform"){
    u_i = sapply(x,function(y,lim=Alpha[-1L])sum(lim<y)+1L)
    for(i in 1:n){
      M_psi_nm[i,u_i[i]]   = 1
      M_Psi_nm[i,1:u_i[i]] = c(if(u_i[i]>1){Delta[1:(u_i[i]-1)]},
                               x[i]-Alpha[u_i[i]])
    }
    ##
  }else{
    if(basis=="gaussian"){
      Sigma = knots$Sigma
      for(u in 1:m){
        M_psi_nm[,u] =  dnorm((x-Alpha[u])/Sigma[u])/(Sigma[u]*Delta[u])
        M_Psi_nm[,u] = (pnorm((x-Alpha[u])/Sigma[u])-
                          pnorm((Alpha[1]-Alpha[u])/Sigma[u]))/Delta[u]
      }
      ##
    }else{
      seq1n = 1:n
      Alpha_star   = as.numeric(c(rep(Alpha[1],order-1L),Alpha,rep(Alpha[n.Alpha],order-1L)))
      M_psi_nm     = M_Psi_nm = cbind(M_psi_nm,0)        
      if(which.matrix[1]){
        Alpha_star_x = sapply(x,function(y,lim=Alpha[-1L])sum(lim<y)+1L)+order-1L
        if(basis=="msplines"){
          M_psi_nm[(Alpha_star_x-1L)*n+seq1n]=1/(Alpha_star[Alpha_star_x+1]-Alpha_star[Alpha_star_x])
          if(order>1){
            for(ow in 2L:order){
              uw_x = Alpha_star_x-ow+1L
              for(pw in 0:(ow-1L)){
                pos_x = (uw_x+pw-1L)*n+seq1n
                M_psi_nm[pos_x]=(ow/((ow-1)*(Alpha_star[1:m+ow]-Alpha_star[1:m])))[uw_x+pw]*
                  ((x-Alpha_star[uw_x+pw])*M_psi_nm[pos_x]+
                     (Alpha_star[uw_x+pw+ow]-x)*M_psi_nm[pos_x+n])    
              }
            }
          }
        }else{
          uw_x = Alpha_star_x-order+1L
          for(pw in 0:(order-1L)){
            pos_x = (uw_x+pw-1L)*n+seq1n
            pos_1 = (uw_x+pw)==1
            pos_m = (uw_x+pw)==m
            pos_other = pos_1==FALSE & pos_m==FALSE    
            M_psi_nm[pos_x[pos_other]]=(6*(x-Alpha_star[uw_x+pw])*(x-Alpha_star[uw_x+pw+order])/ 
                                          ((Alpha_star[uw_x+pw]-Alpha_star[uw_x+pw+order])^3))[pos_other]
            M_psi_nm[pos_x[pos_1]]=(12*(x-Alpha_star[uw_x+pw+order])*(x-2*Alpha_star[uw_x+pw]+Alpha_star[uw_x+pw+order])/ 
                                      ((2*Alpha_star[uw_x+pw]-2*Alpha_star[uw_x+pw+order])^3))[pos_1]
            M_psi_nm[pos_x[pos_m]]=(12*(x-Alpha_star[uw_x+pw])*(x+Alpha_star[uw_x+pw]-2*Alpha_star[uw_x+pw+order])/ 
                                      ((2*Alpha_star[uw_x+pw]-2*Alpha_star[uw_x+pw+order])^3))[pos_m]
          }
        }
        M_psi_nm = M_psi_nm[,1:m,drop=FALSE]
      }    
      if(which.matrix[2]){
        rank.x   = rank(x)
        x        = x[order(x)]    
        Alpha_x  = sapply(x,function(y,lim=Alpha[-1L])sum(lim<y)+1L)    
        up_u     = cumsum(tabulate(Alpha_x,n.Alpha-1))
        for(uw in 1:(m-order+1)){M_Psi_nm[min(n,up_u[uw]+1):n,uw] = 1}    
        if(basis=="msplines"){
          Alpha_star2 = c(rep(Alpha[1],order),Alpha,rep(Alpha[n.Alpha],order))    
          factor_v    = c((Alpha_star2[(order+2):length(Alpha_star2)]-Alpha_star2[1:(length(Alpha_star2)-order-1)])/
                            (order+1),rep(0,order-1))
          M_psi2_nm   = cbind(basis_mpl(x,knots,basis=basis,order=order+1,which=1),matrix(0,n,order-1))
          pos_xo  = rep((Alpha_x-1L)*n,1)+seq1n
          pos_xo1 = rep(pos_xo,order)+rep(1:order,each=n)*n
          for(ow in 0:(order-1)){
            M_Psi_nm[pos_xo+ow*n] = apply(matrix(M_psi2_nm[pos_xo1+ow*n]*
                                                   factor_v[rep(Alpha_x,order)+rep((1:order)+ow,each=n)],ncol=order),1,sum)
          }
        }else{
          Alpha_star_x = sapply(x,function(y,lim=Alpha[-1L])sum(lim<y)+1L)+order-1L
          uw_x = Alpha_star_x-order+1L
          for(pw in 0:(order-1L)){
            pos_x = (uw_x+pw-1L)*n+seq1n
            pos_1 = (uw_x+pw)==1
            pos_m = (uw_x+pw)==m
            pos_other = pos_1==FALSE & pos_m==FALSE    
            M_Psi_nm[pos_x[pos_other]]=((x-Alpha_star[uw_x+pw])^2*(2*x+Alpha_star[uw_x+pw]-3*Alpha_star[uw_x+pw+order])/ 
                                          ((Alpha_star[uw_x+pw]-Alpha_star[uw_x+pw+order])^3))[pos_other]
            M_Psi_nm[pos_x[pos_1]]=((x-Alpha_star[uw_x+pw])*(x^2-2*x*Alpha_star[uw_x+pw]-2*Alpha_star[uw_x+pw]^2+
                                                               6*Alpha_star[uw_x+pw]*Alpha_star[uw_x+pw+order]-3*Alpha_star[uw_x+pw+order]^2)/ 
                                      (2*(Alpha_star[uw_x+pw]-Alpha_star[uw_x+pw+order])^3))[pos_1]
            M_Psi_nm[pos_x[pos_m]]=((x-Alpha_star[uw_x+pw])^2*(x+2*Alpha_star[uw_x+pw]-3*Alpha_star[uw_x+pw+order])/ 
                                      (2*(Alpha_star[uw_x+pw]-Alpha_star[uw_x+pw+order])^3))[pos_m]
          }
        }
        M_Psi_nm = M_Psi_nm[rank.x,1:m,drop=FALSE]
      }}}
  if(all(which.matrix)){list(psi=M_psi_nm,Psi=M_Psi_nm)
  }else{if(which.matrix[1]){M_psi_nm}else{M_Psi_nm}}
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
penalty_mpl=function(control,knots){
  #
  penalty = control$penalty
  m       = knots$m
  if(control$basis=="uniform"){
    D = diag(m)*c(1,-2)[penalty]
    E = diag(m+1)[-(m+1),-1]*c(-1,1)[penalty]
    B = D+E+list(0,t(E))[[penalty]]
    B[c(m+1,(m-1)*m)[1:penalty]] = c(-1,2)[penalty]
    M_R_mm = t(B)%*%B        
  }else{
    M_R_mm = matrix(0,m,m)
    if(control$basis=="gaussian"){
      int_rij_2.fun=function(x,mu_i,mu_j,sig_i,sig_j,t1,tn){
        K     = 4*(pnorm((t1-mu_i)/sig_i)-pnorm((tn-mu_i)/sig_i))*(pnorm((t1-mu_j)/sig_j)-pnorm((tn-mu_j)/sig_j))
        q1q2a = 4*dnorm(x,mu_i,sig_i)*dnorm(x,mu_j,sig_j)*sig_i*sig_j*2*pi*sqrt(sig_i^2+sig_j^2)*
          ((mu_j-x)*sig_i^6*((x-mu_i)^2-sig_i^2)+
             sig_i^4*sig_j^2*((x-4*mu_j+3*mu_i)*sig_i^2-(x-mu_i)^2*(3*x-4*mu_j+mu_i))-
             sig_i^2*sig_j^4*(x-mu_j)^2*(3*x-4*mu_i+mu_j)+
             sig_j^6*((x+3*mu_j-4*mu_i)*sig_i^2-(x-mu_j)^2*(x-mu_i))+
             sig_j^8*(x-mu_i)
          )
        q1q2b = 2*dnorm(mu_j,mu_i,sqrt(sig_i^2+sig_j^2))*pi*sqrt(sig_i^2+sig_j^2)*sig_i^3*sig_j^3*
          (mu_j^4-4*mu_j^3*mu_i+mu_i^4+6*mu_j^2*(mu_i^2-sig_i^2-sig_j^2)-
             6*mu_i^2*(sig_i^2+sig_j^2)+3*(sig_i^2+sig_j^2)^2-4*mu_i*mu_j*(mu_i^2-3*(sig_i^2+sig_j^2))
          )*
          (2*pnorm(((x-mu_j)*sig_i^2+(x-mu_i)*sig_j^2)/(sig_i*sig_j*sqrt(sig_i^2+sig_j^2)))-1)
        q3    = pi*sig_i^3*sig_j^3*(sig_i^2+sig_j^2)^(9/2)*K
        (q1q2a+q1q2b)/q3
      }
      int_rij_1.fun=function(x,mu_i,mu_j,sig_i,sig_j,t1,tn){
        K  = 4*(pnorm((t1-mu_i)/sig_i)-pnorm((tn-mu_i)/sig_i))*(pnorm((t1-mu_j)/sig_j)-pnorm((tn-mu_j)/sig_j))
        q2 =  dnorm((mu_i-mu_j)/sqrt(sig_i^2+sig_j^2))*2*pi*
          sig_i*sig_j*(sig_i^2+sig_j^2-(mu_i-mu_j)^2)*
          (2*pnorm(((x-mu_j)*sig_i^2+(x-mu_i)*sig_j^2)/(sig_i*sig_j*sqrt(sig_i^2+sig_j^2)))-1)
        q1 =   4*pi*dnorm((x-mu_i)/sig_i)*dnorm((x-mu_j)/sig_j)*
          sqrt(sig_i^2+sig_j^2)*((mu_i-x)*sig_i^2+(mu_j-x)*sig_j^2)
        q3 =  pi*sig_i*sig_j*(sig_i^2+sig_j^2)^(5/2)*K
        (q1+q2)/q3    
      }
      for(i in 1:m){
        for(j in i:m){
          if(penalty==2){
            M_R_mm[i,j] = M_R_mm[j,i] = 
              int_rij_2.fun(knots$Alpha[m],knots$Alpha[i],knots$Alpha[j],knots$Sigma[i],knots$Sigma[j],knots$Alpha[1],knots$Alpha[m])-
              int_rij_2.fun(knots$Alpha[1],knots$Alpha[i],knots$Alpha[j],knots$Sigma[i],knots$Sigma[j],knots$Alpha[1],knots$Alpha[m])
          }else{  M_R_mm[i,j] = M_R_mm[j,i] = 
            int_rij_1.fun(knots$Alpha[m],knots$Alpha[i],knots$Alpha[j],knots$Sigma[i],knots$Sigma[j],knots$Alpha[1],knots$Alpha[m])-
            int_rij_1.fun(knots$Alpha[1],knots$Alpha[i],knots$Alpha[j],knots$Sigma[i],knots$Sigma[j],knots$Alpha[1],knots$Alpha[m])}
        }
      }
    }else{  Alpha      = knots$Alpha
    n.Alpha    = length(Alpha)
    order      = control$order
    Alpha_star = c(rep(Alpha[1],order-1L),Alpha,rep(Alpha[n.Alpha],order-1L))
    if(control$basis=="msplines"){
      seq1n        = 1L:(n.Alpha-1)
      n.Alpha_star = length(Alpha_star)
      Alpha_star_x = sapply(Alpha[-1],function(y,lim=Alpha[-1L])sum(lim<y)+1L)+order-1L                
      M_d2f_mm = matrix(0,n.Alpha-1,n.Alpha+order-1L)
      M_d2f_mm[(Alpha_star_x-1L)*(n.Alpha-1)+seq1n]=1/(Alpha_star[Alpha_star_x+1]-Alpha_star[Alpha_star_x])
      for(ow in 2L:order){
        pw   = 1L:ow 
        uw_x = Alpha_star_x-ow+1L
        for(pw in 0:(ow-1L)){
          M_d2f_mm[(uw_x+pw-1L)*(n.Alpha-1)+seq1n]=
            (ow/(Alpha_star[1:(n.Alpha+ow)+ow]-Alpha_star[1:(n.Alpha+ow)]))[uw_x+pw]*
            (M_d2f_mm[(uw_x+pw-1L)*(n.Alpha-1)+seq1n]-M_d2f_mm[(uw_x+pw)*(n.Alpha-1)+seq1n])
        }
      }
      M_d2f_mm = M_d2f_mm[,1:m,drop=FALSE]
      for(uw in 1:m){
        for(vw in uw:m){
          M_R_mm[uw,vw] = M_R_mm[vw,uw] = 
            sum((M_d2f_mm[,uw]*M_d2f_mm[,vw])*(Alpha[-1]-Alpha[-n.Alpha]))
        }
      }
    }else{    for(uw in 1:m){
      f_u = ifelse(uw==1|uw==m,4,1)
      for(vw in uw:m){
        if(Alpha_star[vw]<Alpha_star[uw+order]){
          f_v = ifelse(vw==1|vw==m,4,1)
          M_R_mm[uw,vw] = M_R_mm[vw,uw] = (144*(Alpha_star[uw+order]-Alpha_star[vw]))/
            (f_v*f_u*(Alpha_star[uw+order]-Alpha_star[uw])^3*(Alpha_star[vw+order]-Alpha_star[vw])^3)
        }
      }
    }
    }
    }}
  M_R_mm
}



.nf = function(file){
  temp=file
  vect.colw=seq(1,dim(temp)[2])[sapply(temp,class)=="factor"]
  if(length(vect.colw)>0){
    for(colw in 1:length(vect.colw)){temp[,vect.colw[colw]]=as.character(temp[,vect.colw[colw]])}
  }
  temp
}

.ac = function(...){as.character(...)}

.p = function(...,sep=""){paste(...,sep=sep)}

.an = function(...){as.numeric(...)}



compute_likelihood <- function(x){
  ll=sum(log(x$M_mu_no1)+log(x$M_h0_no1)-x$M_H_no1)-
    sum(x$M_H_nr1)+
    sum(log(1-x$M_S_nl1))+
    sum(log(x$M_S1mS2_ni1))-
    x$s_lambda*x$thetaRtheta
  # cat(c(sum(log(x$M_mu_no1)+log(x$M_h0_no1)-x$M_H_no1), sum(x$M_H_nr1),sum(log(1-x$M_S_nl1)),sum(log(x$M_S1mS2_ni1))), "\n")
  # print(str(list(x$M_mu_no1, x$M_h0_no1, x$M_H_no1, x$M_H_nr1, x$M_S_nl1, x$M_S1mS2_ni1)))
  # cat("#############\n")
  ll
}



