gradBeta <- function(s_obj){
  M_gradbeta_p1_components     <- lapply(s_obj, function(ss){
    (ss$M_tX_nop%*%(1-ss$M_H_no1))-
      ss$M_tX_nrp%*%ss$M_H_nr1+
      ss$M_tX_nlp%*%(ss$M_S_nl1*ss$M_H_nl1/(1-ss$M_S_nl1))+
      ss$M_tX_nip%*%((ss$M_H2_ni1*ss$M_S2_ni1-ss$M_H1_ni1*ss$M_S1_ni1)/ss$M_S1mS2_ni1)
  })
  M_gradbeta_p1 <- Reduce(`+`, M_gradbeta_p1_components)
  M_gradbeta_p1
}


hessianBeta <- function(s_obj){
  # M_hessbeta_p1_components <- lapply(s_obj, function(ss){
  #   ss$M_tX_nop%*%diag(c(ss$M_H_no1),n.ctype[2],n.ctype[2])%*%ss$M_X_nop+
  #     ss$M_tX_nrp%*%diag(c(ss$M_H_nr1),n.ctype[1],n.ctype[1])%*%ss$M_X_nrp+
  #     ss$M_tX_nlp%*%diag(c((ss$M_S_nl1/(1-ss$M_S_nl1)^2*ss$M_H_nl1^2-ss$M_S_nl1/(1-ss$M_S_nl1)*ss$M_H_nl1)),n.ctype[3],n.ctype[3])%*%ss$M_X_nlp+            
  #     ss$M_tX_nip%*%diag(c((ss$M_S1_ni1*ss$M_S2_ni1/ss$M_S1mS2_ni1^2*(ss$M_H2_ni1-ss$M_H1_ni1)^2+
  #                             (ss$M_S1_ni1*ss$M_H1_ni1-ss$M_S2_ni1*ss$M_H2_ni1)/ss$M_S1mS2_ni1
  #     )),n.ctype[4],n.ctype[4])%*%ss$M_X_nip   
  
  M_hessbeta_p1_components <- lapply(s_obj, function(ss){
    ss$M_tX_nop%*%diag(c(ss$M_H_no1))%*%ss$M_X_nop+
      ss$M_tX_nrp%*%diag(c(ss$M_H_nr1))%*%ss$M_X_nrp+
      ss$M_tX_nlp%*%diag(c((ss$M_S_nl1/(1-ss$M_S_nl1)^2*ss$M_H_nl1^2-ss$M_S_nl1/(1-ss$M_S_nl1)*ss$M_H_nl1)))%*%ss$M_X_nlp+
      ss$M_tX_nip%*%diag(c((ss$M_S1_ni1*ss$M_S2_ni1/ss$M_S1mS2_ni1^2*(ss$M_H2_ni1-ss$M_H1_ni1)^2+
                              (ss$M_S1_ni1*ss$M_H1_ni1-ss$M_S2_ni1*ss$M_H2_ni1)/ss$M_S1mS2_ni1)))%*%ss$M_X_nip
  })
  
  # M_hessbeta_p1_components <- lapply(s_obj, function(ss){
  #   ss$M_tX_nop%*%diag(c(ss$M_H_no1))%*%ss$M_X_nop+
  #     ss$M_tX_nrp%*%diag(c(ss$M_H_nr1))%*%ss$M_X_nrp+
  #     ss$M_tX_nlp%*%diag(c((ss$M_S_nl1/(1-ss$M_S_nl1)^2*ss$M_H_nl1^2)))%*%ss$M_X_nlp+            #troubling term removed
  #     ss$M_tX_nip%*%diag(c((ss$M_S1_ni1*ss$M_S2_ni1/ss$M_S1mS2_ni1^2*(ss$M_H2_ni1-ss$M_H1_ni1)^2
  #                             # +(ss$M_S1_ni1*ss$M_H1_ni1-ss$M_S2_ni1*ss$M_H2_ni1)/ss$M_S1mS2_ni1
  #                           )))%*%ss$M_X_nip     
  # })
  M_hessbeta_p1 <- Reduce(`+`, M_hessbeta_p1_components)
  M_hessbeta_p1
}


update_s_obj_beta <- function(s_obj, M_beta_p1, thr, thr2){
  for (is in 1:length(s_obj)){
    s_obj[[is]]$M_mu_no1   <- exp(s_obj[[is]]$M_X_nop%*%M_beta_p1)
    s_obj[[is]]$M_mu_nr1   <- exp(s_obj[[is]]$M_X_nrp%*%M_beta_p1)
    s_obj[[is]]$M_mu_nl1   <- exp(s_obj[[is]]$M_X_nlp%*%M_beta_p1)
    s_obj[[is]]$M_mu_ni1   <- exp(s_obj[[is]]$M_X_nip%*%M_beta_p1)
    s_obj[[is]]$M_H_no1  <- s_obj[[is]]$M_H0_no1 * s_obj[[is]]$M_mu_no1
    s_obj[[is]]$M_H_nr1  <- s_obj[[is]]$M_H0_nr1 * s_obj[[is]]$M_mu_nr1
    s_obj[[is]]$M_H_nl1  <- s_obj[[is]]$M_H0_nl1 * s_obj[[is]]$M_mu_nl1
    s_obj[[is]]$M_H1_ni1 <- s_obj[[is]]$M_H01_ni1* s_obj[[is]]$M_mu_ni1
    s_obj[[is]]$M_H2_ni1 <- s_obj[[is]]$M_H02_ni1* s_obj[[is]]$M_mu_ni1
    s_obj[[is]]$M_S_no1  <- exp(-s_obj[[is]]$M_H_no1)
    s_obj[[is]]$M_S_nl1  <- exp(-s_obj[[is]]$M_H_nl1)
    s_obj[[is]]$M_S1_ni1 <- exp(-s_obj[[is]]$M_H1_ni1)
    s_obj[[is]]$M_S2_ni1 <- exp(-s_obj[[is]]$M_H2_ni1)
    # avoid division by 0
    s_obj[[is]]$M_S_nl1[s_obj[[is]]$M_S_nl1 > thr]   <- thr
    s_obj[[is]]$M_S1_ni1[s_obj[[is]]$M_S1_ni1 > thr] <- thr
    s_obj[[is]]$M_S2_ni1[s_obj[[is]]$M_S2_ni1 > thr] <- thr
    s_obj[[is]]$M_S1mS2_ni1           <- s_obj[[is]]$M_S1_ni1-s_obj[[is]]$M_S2_ni1  
    s_obj[[is]]$M_S1mS2_ni1[s_obj[[is]]$M_S1mS2_ni1<thr2] <- thr2
  }
  s_obj
}

update_omega <- function(s_omega, s_obj){
  if(s_omega>=1e-2){
    s_omega = s_omega/s_obj[[1]]$s_kappa
  }else{if(s_omega<1e-2&s_omega>=1e-5){
    s_omega = s_omega*5e-2
  }else{if(s_omega<1e-5){
    s_omega = s_omega*1e-5    
  }}}
  s_omega
}
  
  
update_s_obj_steptheta <- function(s_obj, thr, thr2, s_nu){
  for (is in 1:length(s_obj)){
    s_obj[[is]]$M_gradthetaA_m1     = 
      s_obj[[is]]$M_tpsi_nom%*%(1/s_obj[[is]]$M_h0_no1)+
      s_obj[[is]]$M_tPsi_nlm%*%(s_obj[[is]]$M_S_nl1*s_obj[[is]]$M_mu_nl1/(1-s_obj[[is]]$M_S_nl1))+
      s_obj[[is]]$M_tPsi2_nim%*%(s_obj[[is]]$M_S2_ni1*s_obj[[is]]$M_mu_ni1/(s_obj[[is]]$M_S1mS2_ni1))-
      s_obj[[is]]$TwoLRtheta*(s_obj[[is]]$TwoLRtheta<0)+0.3        
    s_obj[[is]]$M_gradthetaB_m1     = 
      s_obj[[is]]$M_tPsi_nom%*%s_obj[[is]]$M_mu_no1+
      s_obj[[is]]$M_tPsi_nrm%*%s_obj[[is]]$M_mu_nr1+
      s_obj[[is]]$M_tPsi1_nim%*%(s_obj[[is]]$M_S1_ni1*s_obj[[is]]$M_mu_ni1/(s_obj[[is]]$M_S1mS2_ni1))+
      s_obj[[is]]$TwoLRtheta*(s_obj[[is]]$TwoLRtheta>0)+0.3    
    s_obj[[is]]$M_gradtheta_m1 = s_obj[[is]]$M_gradthetaA_m1-s_obj[[is]]$M_gradthetaB_m1
    #
    s_obj[[is]]$M_s_m1         = s_obj[[is]]$M_theta_m1/s_obj[[is]]$M_gradthetaB_m1
    s_obj[[is]]$M_steptheta_p1 = s_obj[[is]]$M_s_m1 * s_obj[[is]]$M_gradtheta_m1
  }
  s_obj
}


update_s_obj_theta <- function(s_obj, thr, thr2, s_omega){
  for (is in 1:length(s_obj)){
    s_obj[[is]]$M_theta_m1     = s_obj[[is]]$M_theta_m1_OLD+s_omega*s_obj[[is]]$M_steptheta_p1
    s_obj[[is]]$M_theta_m1[s_obj[[is]]$M_theta_m1<thr2] = thr2
    s_obj[[is]]$M_h0_no1    = s_obj[[is]]$M_psi_nom%*%s_obj[[is]]$M_theta_m1
    s_obj[[is]]$M_H0_no1    = s_obj[[is]]$M_Psi_nom%*%s_obj[[is]]$M_theta_m1
    s_obj[[is]]$M_H0_nr1    = s_obj[[is]]$M_Psi_nrm%*%s_obj[[is]]$M_theta_m1
    s_obj[[is]]$M_H0_nl1    = s_obj[[is]]$M_Psi_nlm%*%s_obj[[is]]$M_theta_m1
    s_obj[[is]]$M_H01_ni1   = s_obj[[is]]$M_Psi1_nim%*%s_obj[[is]]$M_theta_m1
    s_obj[[is]]$M_H02_ni1   = s_obj[[is]]$M_Psi2_nim%*%s_obj[[is]]$M_theta_m1
    s_obj[[is]]$Rtheta      = s_obj[[is]]$M_R_mm%*%s_obj[[is]]$M_theta_m1
    s_obj[[is]]$thetaRtheta = t(s_obj[[is]]$M_theta_m1)%*%s_obj[[is]]$Rtheta
    s_obj[[is]]$TwoLRtheta  = s_obj[[is]]$s_lambda*2*s_obj[[is]]$Rtheta            
    s_obj[[is]]$M_H_no1  = s_obj[[is]]$M_H0_no1 * s_obj[[is]]$M_mu_no1
    s_obj[[is]]$M_H_nr1  = s_obj[[is]]$M_H0_nr1 * s_obj[[is]]$M_mu_nr1
    s_obj[[is]]$M_H_nl1  = s_obj[[is]]$M_H0_nl1 * s_obj[[is]]$M_mu_nl1
    s_obj[[is]]$M_H1_ni1 = s_obj[[is]]$M_H01_ni1* s_obj[[is]]$M_mu_ni1
    s_obj[[is]]$M_H2_ni1 = s_obj[[is]]$M_H02_ni1* s_obj[[is]]$M_mu_ni1
    s_obj[[is]]$M_S_no1  = exp(-s_obj[[is]]$M_H_no1)
    s_obj[[is]]$M_S_nl1  = exp(-s_obj[[is]]$M_H_nl1)
    s_obj[[is]]$M_S1_ni1 = exp(-s_obj[[is]]$M_H1_ni1)
    s_obj[[is]]$M_S2_ni1 = exp(-s_obj[[is]]$M_H2_ni1)
    # avoid division by 0
    s_obj[[is]]$M_S_nl1[s_obj[[is]]$M_S_nl1==1]   = thr
    s_obj[[is]]$M_S1_ni1[s_obj[[is]]$M_S1_ni1==1] = thr
    s_obj[[is]]$M_S2_ni1[s_obj[[is]]$M_S2_ni1==1] = thr
    s_obj[[is]]$M_S1mS2_ni1           = s_obj[[is]]$M_S1_ni1-s_obj[[is]]$M_S2_ni1  
    s_obj[[is]]$M_S1mS2_ni1[s_obj[[is]]$M_S1mS2_ni1<thr2] = thr2
  }
  #
  s_obj
}
  
  
  