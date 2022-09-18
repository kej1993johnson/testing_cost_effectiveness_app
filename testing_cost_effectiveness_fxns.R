
#2 means without the isolation abandonment rate



# Functions for cost-effectiveness of testing strategy in a university
solve_transmission_eqns<- function(par_table, vacc_inf_correlation = 'independent', test_type = "antigen"){
  par_names<-colnames(par_table)
  for (i in 1:length(par_names)){
    assign(par_names[i], as.double(par_table[i]))
  }
  
  if (vacc_inf_correlation =='independent'){
  # Assume vaccination is independent of being previously infected
  # P(R|V) = P(R|!V) = P(R)
  init_ru = init_rec- init_rec*init_vacc
  init_rv = init_rec*init_vacc
  init_sv = (1-init_rec)*init_vacc 
  }
  else if (vacc_inf_correlation == 'correlated'){
    # P(R|!V) = x*P(R)
    x = odds_inf_not_vac
    prob_rec_given_vac = init_rec/(init_vacc + x*(1-init_vacc))
    prob_rec_given_not_vac = x*prob_rec_given_vac
    init_ru = prob_rec_given_not_vac*(1-init_vacc) #P(R|!V)*P(!V)
    init_su = (1-prob_rec_given_not_vac)*(1-init_vacc) #P(!R|!V)*P(!V)
    init_rv = prob_rec_given_vac*init_vacc # P(R|V)*P(V) 
    init_sv = (1-prob_rec_given_vac)*init_vacc # P(!R|V)*P(V)
    test<-init_ru+init_rv+init_sv+init_su
    p_rec = init_ru + init_rv
    p_vacc = init_rv + init_sv
  }
  #VE = 1-odds_case_vacc/odds_pop_vacc
  odds_pop_vacc = init_vacc/(1-init_vacc)
  odds_case_vacc = odds_pop_vacc*sigma_v
  prob_cases_vacc = odds_case_vacc/(1+odds_case_vacc)
  
  psym_v = sym_red/sigma_v
  avg_time_asym = psym*t_pre_sym + (1-psym)*t_inf
  avg_time_asym_v = psym_v*t_pre_sym + (1-psym_v)*t_inf
  avg_time_sym_undetected = (1-is_sym)*(t_sym_state) + is_sym*t_iso
  delta_au = 1/avg_time_asym
  delta_av = 1/avg_time_asym_v
  delta_s = 1/t_sym_state
  delta_s = 1/avg_time_sym_undetected
  delta_q = 1/t_quarantine
  beta_a = R0*1/t_inf #transmission rate of asymptomatic individuals
  beta_s = R0*1/t_inf # transmission rate of symptomatic individuals 
  #intro_rate = init_prev/7 #  p = new infections per day per capita in surrounding community. 
  tdays<-seq(from = 0, to = duration_school_year, by = 1)
  dates<-seq(from = as.Date("2021-08-20"), by = "day", length.out = duration_school_year+1)
  # define initial conditions of compartments
  Squ =0
  Equ = 0
  Iaqu = 0
  Isqu = 0
  Rqu = 0
  Su = N*(1-init_ru - init_rv - init_sv - init_prev)# (unvaccinated susceptibles = 1-unvaccinated recovered - vaccinated susceptible - currently infected)
  Eu = 0
  Iau = 0.5*N*init_prev*(1-sigma_v)*(1-prob_cases_vacc) # % pre/asymptomatic * percent infected* percent unvaccinated
  Isu =0.5*N*init_prev*(1-prob_cases_vacc) # % symptomatic * percent infeted * percent unvaccinated
  Ru = N*(init_ru)
  Sqv = 0
  Eqv = 0
  Iaqv = 0
  Isqv = 0
  Rqv = 0
  Sv = N*(init_sv)
  Ev = 0
  Iav = 0.5*N*init_prev*prob_cases_vacc
  Isv = 0.5*N*init_prev*prob_cases_vacc
  Rv = N*(init_rv)
  Ntest = Su+Eu+Iau+Isu+Ru+Sv+Ev+Iav+Isv+Rv # should equal N
  I0 = Iau + Isu + Iaqu + Isqu
  Q0 = Squ + Equ + Iaqu + Isqu + Sqv + Eqv + Iaqv + Isqv 
  Nu = Su+Eu+Iau+Isu +Ru
  Nv = Sv + Ev + Iav + Isv + Rv
  RTs = 0
  PCR = 0
  cumI = 0
  cumsI = 0
  
  # Account for differences in timing of the different test types
  if (test_type == 'antigen'){
    inf_rem = 1
  }
  if (test_type == 'PCR'){
    pct_inf_removed = get_pct_infectiousness_removed(V0=3, Vf=6, delta_add=0, t0topeak=2, Vpeak = 9,
    tpeaktof=7, LOD=5, LOI=6, dt=0.01, t_test = 1, test_sensitivity = 0.95)
    inf_rem<-pct_inf_removed$eff_percent_removed
  }
  
  y <- c(Squ = Squ, Equ=Equ, Iaqu=Iaqu, Isqu=Isqu, Rqu = Rqu, Su=Su, Eu=Eu, Iau=Iau, Isu=Isu, Ru=Ru, Sqv=Sqv, 
         Eqv=Eqv, Iaqv=Iaqv, Isqv=Isqv, Rqv = Rqv, Sv=Sv, Ev=Ev, Iav=Iav, Isv=Isv, Rv=Rv, I=I0,Isymdet = Isqu+Isqv, Idet =Isqu+Isqv+Iaqu+Iaqv,
         Q = Q0, Nu = Nu, Nv = Nv, RTs =0, PCR = 0, cumI =0, cumIvax=0, cumIunvax=0, cumsI=0, FPs=0, TPs=0)
  params<-c(delta_au=delta_au, delta_av = delta_av, delta_s=delta_s, delta_q=delta_q, beta_a=beta_a, beta_s=beta_s, gamma=gamma,
            psym = psym, e_v=e_v, sigma_v=sigma_v, Sp=Sp, Se=Se, 
            w_v =w_v, w_u=w_u, f_v=f_v, f_u=f_u, is=is, k=k, N=N, intro_rate = intros_per_week/7,
            prob_cases_vacc = prob_cases_vacc,
            init_vacc= init_vacc, inf_rem = inf_rem)
  model<-function(tdays, y, params){
    with(as.list(c(y,params)), {
      X = Su+Eu+Iau+Isu+Ru+Sv+Ev+Iav+Isv+Rv # active population
      beta_star = beta_a*(Iau + e_v*Iav) + beta_s*(Isu + e_v*Isv) # force of infection
      # quarantined unvaccinated
      dSqu<- Su*w_u*(1-Sp)*is*f_u - Squ*k
      dEqu<-Eu*w_u*(1-Sp)*is*f_u - Equ*k
      dIaqu<-Iau*w_u*Se*is*f_u*inf_rem - delta_q*Iaqu 
      dIsqu<-Isu*w_u*Se*is*f_u*inf_rem + is_sym*delta_s*Isu  - delta_q*Isqu 
      dRqu<- Ru*w_u*(1-Sp)*is*f_u - Rqu*k
      # unvaccinated
      dSu<- - beta_star*Su/X - Su*w_u*(1-Sp)*is*f_u + Squ*(k) - (1-prob_cases_vacc)*intro_rate
      dEu<- beta_star*Su/X- gamma*Eu - Eu*w_u*(1-Sp)*is*f_u + (1-prob_cases_vacc)*intro_rate
      dIau<- gamma*Eu - delta_au*Iau  - Iau*w_u*Se*is*f_u*inf_rem + Equ*k 
      dIsu<- psym*delta_au*Iau - delta_s*Isu - Isu*w_u*Se*is*f_u*inf_rem 
      dRu<- (1-is_sym)*delta_s*Isu + (1-psym)*delta_au*Iau - Ru*w_u*(1-Sp)*is*f_u + Rqu*k  + delta_q*Isqu +delta_q*Iaqu
      # quarantined vaccinated
      dSqv<-Sv*w_v*(1-Sp)*is*f_v - Sqv*k
      dEqv<- Ev*w_v*(1-Sp)*is*f_v - Eqv*k
      dIaqv<- Iav*w_v*Se*is*f_v*inf_rem - delta_q*Iaqv 
      dIsqv<-Isv*w_v*Se*is*f_v*inf_rem +is_sym*delta_s*Isv  - delta_q*Isqv
      dRqv<-Rv*w_v*(1-Sp)*is*f_v - Rqv*k
      # vaccinated
      dSv<- -sigma_v*beta_star*Sv/X - Sv*w_v*(1-Sp)*is*f_v + Sqv*(k)- prob_cases_vacc*intro_rate
      dEv<- sigma_v*beta_star*Sv/X - gamma*Ev - Ev*w_v*(1-Sp)*is*f_v + prob_cases_vacc*intro_rate
      dIav<- gamma*Ev - delta_av*Iav  - Iav*w_v*Se*is*f_v*inf_rem + Eqv*(k) 
      dIsv<-psym_v*delta_av*Iav - delta_s*Isv - Isv*w_v*Se*is*f_v*inf_rem 
      dRv<- delta_s*(1-is_sym)*Isv + (1-psym_v)*delta_av*Iav - Rv*w_v*(1-Sp)*is*f_v + Rqv*k + delta_q*Isqv +delta_q*Iaqv
      dI<- dIau + dIsu + dIaqu + dIsqu + dIav + dIsv + dIaqv + dIsqv
      dIsymdet<-dIsqu + dIsqv
      dIdet<-dIsqu + dIsqv + dIaqu + dIaqv
      dQ<-dSqu + dEqu + dIaqu + dIsqu + dRqu + dSqv + dEqv + dIaqv + dIsqv + dRqv
      dNu<-dSu+dEu+dIau+dIsu +dRu + dSqu+dEqu+dIaqu+dIsqu +dRqu
      dNv<- dSv + dEv + dIav + dIsv + dRv + dSqv + dEqv + dIaqv + dIsqv + dRqv
      dRTs<-Nu*f_u*w_u + Nv*f_v*w_v
      dPCR<- Su*w_u*(1-Sp)*is*f_u + Eu*w_u*(1-Sp)*is*f_u + Iau*w_u*Se*is*f_u + Isu*w_u*Se*is*f_u + is*Isu+
        Ru*w_u*(1-Sp)*is*f_u + Sv*w_v*(1-Sp)*is*f_v + Ev*w_v*(1-Sp)*is*f_v + Iav*w_v*Se*is*f_v + Isv*w_v*Se*is*f_v+
        is*Isv + Rv*w_v*(1-Sp)*is*f_v# all new quarantines (true positives + false positives)
      dcumI<-beta_star*Su/X + sigma_v*beta_star*Sv/X # all new infections 
      dcumIvax<-sigma_v*beta_star*Sv/X
      dcumIunvax<-beta_star*Su/X
      dcumsI<- Isu*w_u*Se*is*f_u + is_sym*Isu + Isv*w_v*Se*is*f_v +is_sym*Isv# all new symptomatic detected infections
      dFP<-Su*w_u*(1-Sp)*is*f_u + Eu*w_u*(1-Sp)*is*f_u + Ru*w_u*(1-Sp)*is*f_u + Sv*w_v*(1-Sp)*is*f_v +
        Ev*w_v*(1-Sp)*is*f_v + Rv*w_v*(1-Sp)*is*f_v  #-> new false positives?
      dTP<- Iau*w_u*Se*f_u + Isu*w_u*Se*f_u + is_sym*delta_s*Isu + Iav*w_v*Se*f_v + Isv*w_v*Se*f_v +
        is_sym*delta_s*Isv

      
      
      list(c(dSqu, dEqu, dIaqu, dIsqu, dRqu, dSu, dEu, dIau, dIsu, dRu, dSqv, dEqv, dIaqv, dIsqv, dRqv, dSv,
             dEv, dIav, dIsv, dRv, dI,dIsymdet, dIdet, dQ, dNu, dNv, dRTs, dPCR, dcumI, dcumIvax, dcumIunvax, dcumsI, dFP, dTP))
    })
  }
  
  out<-ode(y, tdays,model, params)
  #plot(out)
  Ntot<-rowSums(out[,2:21])
  #plot(tdays, Ntot)
  out_df<-data.frame(out)
  # add some columns
  init_imm<-1-(1-init_vacc)*(1-init_rec)
  pct_vacc<-rep(init_vacc, length(tdays))
  pct_imm<-rep(init_imm, length(tdays))
  testing_freq<-rep(f_u, length(tdays))
  prevalence<-rep(init_prev, length(tdays))
  Nvec<-rep(N, length(tdays))
  
  out_df<-cbind(dates, out_df,pct_vacc, testing_freq, prevalence, Nvec)
  
  return(out_df)
}

testing_model<-function(par_table, cost_table, risk_tolerance, vacc_inf_correlation = 'independent', test_type = 'antigen'){
  out_df<-solve_transmission_eqns(par_table, vacc_inf_correlation, test_type)

  # Gets the outputs of interest for the cost-effectiveness model
  # input cost parameters
  
  cost_names<-colnames(cost_table)
  for (i in 1:length(cost_names)){
    assign(cost_names[i], as.double(cost_table[i]))
  }
  par_names<-colnames(par_table)
  for (i in 1:length(par_names)){
    assign(par_names[i], as.double(par_table[i]))
  }
  
  close_thres = case_when(
    risk_tolerance == "CDC red" ~ 100/100000,
    risk_tolerance == "1.5x CDC red" ~ 150/100000,
    risk_tolerance == "2x CDC red" ~200/100000,
    risk_tolerance == "3x CDC red" ~ 300/100000,
    risk_tolerance == "4x CDC red" ~ 400/100000)
  
  thres_red = 100/100000
  thres_15xred = 150/100000
  thres_2xred = 200/100000
  thres_3xred = 300/100000
  thres_4xred = 400/100000
  
 
  
  
  # outputs
  N<-par_table$N
  n_weeks<-par_table$duration_school_year/7
  n_inf = out_df$cumI[nrow(out_df)] # total number infected at end of the year
  n_inf_vax = out_df$cumIvax[nrow(out_df)]
  n_inf_unvax = out_df$cumIunvax[nrow(out_df)]
  max_symdet<-max(out_df$Isymdet)
  n_sym_det<-out_df$cumsI[nrow(out_df)]
  pct_inf<-n_inf/N
  pct_detected = out_df$TPs[nrow(out_df)]/out_df$cumI[nrow(out_df)]
  testing_freq<-par_table$f_u
  n_pos<-out_df$TPs[nrow(out_df)]
  days_of_online = as.numeric(sum((out_df$Isymdet)/N>close_thres))
  cross_thres <- days_of_online>0
  days_above_red = as.numeric(sum((out_df$Isymdet)/N>thres_red))
  days_above_15xred = as.numeric(sum((out_df$Isymdet)/N>thres_15xred))
  days_above_2xred = as.numeric(sum((out_df$Isymdet)/N>thres_2xred))
  days_above_3xred = as.numeric(sum((out_df$Isymdet)/N>thres_3xred))
  days_above_4xred = as.numeric(sum((out_df$Isymdet)/N>thres_4xred))
  cross_red <- days_above_red>0
  cross_15xred <- days_above_15xred>0
  cross_2xred <-days_above_2xred>0
  cross_3xred<-days_above_3xred>0
  cross_4xred<-days_above_4xred>0
  DLL= sum(out_df$Q) # integral of Q (days of learning lost from quarantine only)
  #DLL_tot = sum(out_df$Q) + sum(N-out_df$Q[out_df$Isymdet/N>close_thres]) # plus days online days*
  # the population not in quarantine 
  cost_of_online = days_of_online*online_price
  cost_DLL = DLL*DLL_price # resources spent by students in response to COVID
  cost_PCR <-out_df$PCR[nrow(out_df)]*PCR_price
  cost_isofac<-out_df$TPs[nrow(out_df)]*pct_pos_isofac*isofac_price*7 # total number of positives * percent that need isolation * 
  # * daily cost of isolation facility usage * 7 days isolation period 
  cost_sequencing<-out_df$TPs[nrow(out_df)]*sequencing_price
  pct_per_week_tested<-par_table$f_u*7*par_table$w_u + par_table$f_v*7*par_table$w_v
  pct_per_week_tested_uv<-par_table$f_u*7*par_table$w_u
  pct_per_week_tested_v<-par_table$f_v*7*par_table$w_v
  n_tests_per_week<- (1-par_table$init_vacc)*par_table$N*pct_per_week_tested_uv + 
    (par_table$init_vacc)*par_table$N*pct_per_week_tested_v
  n_weeks <-par_table$duration_school_year/7
  cost_RT<- n_tests_per_week*RT_price*n_weeks

  
  cost_contact_tracing<-out_df$TPs[nrow(out_df)]*contact_tracing_price

  n_PCR<-out_df$PCR[nrow(out_df)]
  n_RT<-out_df$RTs[nrow(out_df)]
  cost_to_UT<- cost_RT + cost_PCR + cost_isofac + cost_contact_tracing + cost_sequencing + cost_of_online
  #testing, contact-tracing, isolation facilities, cost to maintain a PCR lab/sequencing, 
    #cost for staff/faculty/administrator time
  
  
  # to find things averted, need to run with no testing at all
  par_table$f_u<-0
  par_table$f_v<-0
  out_nt_df<-solve_transmission_eqns(par_table)
  
  DLL_nt<-sum(out_nt_df$Q)
  #DLL_tot_nt = sum(out_nt_df$Q) + sum(N-out_nt_df$Q[out_nt_df$Isymdet/N>close_thres]) # plus days online days*
  n_inf_nt<-out_nt_df$cumI[nrow(out_df)]
  pct_inf_nt<-n_inf_nt/N
  days_of_online_nt<-as.numeric(sum((out_nt_df$Isymdet)/N>close_thres))
  #DLL_averted<-DLL_tot_nt-DLL_tot
  DLL_averted<-DLL_nt- DLL
  inf_averted<-n_inf_nt-n_inf
  pct_inf_averted<- (n_inf_nt-n_inf)/n_inf_nt
  days_of_online_averted<-days_of_online_nt - days_of_online
  pct_detected_nt = out_nt_df$TPs[nrow(out_nt_df)]/out_df$cumI[nrow(out_nt_df)]
  
  cost_nt<-out_nt_df$PCR[nrow(out_nt_df)]*PCR_price
  n_PCR_nt<-out_nt_df$PCR[nrow(out_nt_df)]
  incr_cost_DLL<- DLL_averted*DLL_price
  cost_per_DLL_averted<-cost_RT/DLL_averted
  cost_per_DO_averted<-cost_RT/days_of_online_averted
  cost_per_inf_averted<-cost_RT/inf_averted
  cost_testing_per_student<-cost_RT/par_table$N
  
  
  
  pct_vacc<-par_table$init_vacc
  pct_imm<- 1-(1-par_table$init_rec)*(1-par_table$init_vacc)
  init_prev<-par_table$init_prev
  
  summary_df<-data.frame(pct_vacc, testing_freq,init_prev, n_inf, n_inf_vax, n_inf_unvax, max_symdet, pct_inf, pct_detected, 
                         n_pos, DLL, days_of_online, cost_of_online, cost_RT, cost_PCR,
                         cost_DLL,cost_isofac, cost_contact_tracing,cost_sequencing,
                         cost_to_UT,
                         n_PCR, n_RT, n_inf_nt, pct_inf_nt, DLL_nt, days_of_online_nt, DLL_averted,
                         inf_averted, pct_inf_averted, days_of_online_averted, cost_nt, n_PCR_nt, pct_detected_nt, incr_cost_DLL, 
                         cost_per_DLL_averted, cost_per_DO_averted, cost_per_inf_averted, cost_testing_per_student,
                         n_tests_per_week,close_thres,
                         cross_thres, cross_red, cross_15xred, cross_2xred, cross_3xred, cross_4xred)
  return(summary_df)
}
testing_model_w_uncertainty<-function(par_table, cost_table, risk_tolerance, nsamps, par_bounds, vacc_inf_correlation = 'independent',
                                      test_type = 'antigen'){
  par_names<-colnames(par_table)
  for (i in 1:length(par_names)){
    assign(par_names[i], as.double(par_table[i]))
  }
  
  par_bounds_names<-colnames(par_bounds)
  for (i in 1:length(par_bounds_names)){
    assign(par_bounds_names[i],par_bounds[i])
  }
  

  
  
  
  
  # Make the distributions of parameters to sample from:
  #TRIANGULAR
  R0s<-rtri(n = nsamps, min = min(R0vals), max = max(R0vals), mode = R0) 
  init_recs<-rtri(n= nsamps, min = min(init_rec_vals), max = max(init_rec_vals), mode = init_rec)
  psyms<-rtri(n = nsamps, min = min(psym_vals), max = max(psym_vals), mode = psym) 
  sigma_vs<-rtri(n = nsamps, min = min(sigma_v_vals), max = max(sigma_v_vals), mode = sigma_v)
  sym_reds<-rtri(n = nsamps, min = min(sym_red_vals), max = max(sym_red_vals), mode = sym_red) # #UNIFORM
  iss<-runif(n=nsamps, min = min(is_vals), max = max(is_vals)) 
  is_syms<-runif(n= nsamps, min = min(is_sym_vals), max = max(is_sym_vals)) 
  intros<-runif(n = nsamps, min = min(intro_vals), max = max(intro_vals))
  init_prevs<-rtri(n = nsamps, min = min(init_prev_vals), max = max(init_prev_vals), mode = init_prev)
  
  
  # Run the model for each drawn sample
  par_tablei<-par_table
  for (i in 1:nsamps){
    
    par_tablei$R0<-R0s[i]
    par_tablei$init_rec<-init_recs[i]
    par_tablei$psym<-psyms[i]
    par_tablei$sigma_v<-sigma_vs[i]
    par_tablei$sym_red<-sym_reds[i]
    par_tablei$is<-iss[i]
    par_tablei$is_sym<-is_syms[i]
    par_tablei$intros_per_week<-intros[i]
    par_tablei$init_prev<-init_prevs[i]
    
    if (i==1){
      df_t<-solve_transmission_eqns(par_tablei, vacc_inf_correlation, test_type)
      df<-testing_model(par_tablei, cost_table, risk_tolerance, vacc_inf_correlation, test_type)
      sample<-rep(i, nrow(df_t))
      df_t$sample<-sample
      df$sample<- i
    }
    else{
      df_ti<-solve_transmission_eqns(par_tablei, vacc_inf_correlation, test_type)
      dfi<-testing_model(par_tablei, cost_table, risk_tolerance, vacc_inf_correlation, test_type)
      sample<-rep(i, nrow(df_ti))
      df_ti$sample<-sample
      dfi$sample<- i
      df_t<-rbind(df_t, df_ti)
      df<-rbind(df, dfi)
    }
    
  }
  # df_t is the stacked nsamps*length(t)xncols full simulation data frame
  # df is the stacked nsampsxncols summary otuput data frame
  
  
  
  # MAKE AGGREGATED DATAFRAMES
  # aggr_df_t is now length(t)x 3*length(colnames(df_t) because has header for ub, median, lb
  time<-seq(from = 0, to = duration_school_year, by = 1)
  dates<-seq(from = as.Date("2021-08-20"), by = "day", length.out = duration_school_year+1)
  aggr_df_t<-data.frame(time, dates)
  
  
  col_names<-colnames(df_t)
  
  for (i in 3:(length(col_names)-1)){
    df_long<-df_t%>%dplyr::select(col_names[i], time, sample)
    df_wide<-df_long%>%pivot_wider(names_from = sample, values_from = col_names[i])
    val_mat<-as.matrix(df_wide[,2:ncol(df_wide)])
    med<-apply(val_mat , 1 , quantile , probs = 0.5 , na.rm = TRUE )
    lb<-apply(val_mat , 1 , quantile , probs = 0.025 , na.rm = TRUE )
    ub<-apply(val_mat , 1 , quantile , probs = 0.975 , na.rm = TRUE )
    # put these into data frame 
    med_name<-paste0(col_names[i], '_med')
    lb_name<-paste0(col_names[i], '_lb')
    ub_name<-paste0(col_names[i], '_ub')
    aggr_df_t[[med_name]]<-med
    aggr_df_t[[lb_name]]<-lb
    aggr_df_t[[ub_name]]<-ub
  }
  
  
  
  aggr_df<-data.frame() # all the things in df but with ub, lb, and median + the probability of exceeding each threshold
  run<-1
  aggr_df<-data.frame(run) 
  col_names<-colnames(df)
  for (i in (1:(length(col_names)-7))){
    vec<-df[,i]
    med<-quantile(vec, probs = 0.5, na.rm = TRUE)
    lb<-quantile(vec, probs = 0.05, na.rm = TRUE)
    ub<-quantile(vec, probs = 0.95, na.rm = TRUE)
    # put these into data frame 
    med_name<-paste0(col_names[i], '_med')
    lb_name<-paste0(col_names[i], '_lb')
    ub_name<-paste0(col_names[i], '_ub')
    aggr_df[[med_name]]<-med
    aggr_df[[lb_name]]<-lb
    aggr_df[[ub_name]]<-ub
  }
  
  for (i in ((length(col_names)-6):(length(col_names)-1))){
    vec<-df[,i]
    prob<-sum(vec)/nrow(df)
    # put these into data frame 
    prob_name<-paste0(col_names[i], '_prob')
    aggr_df[[prob_name]]<-prob
  }
  out_list<-list(df, df_t, aggr_df, aggr_df_t)
  return(out_list)
}

get_min_testing_per_vacc<-function(df, threshold_prob, vacc_coverage){

for(j in 1:length(vacc_coverage)){
  df_pv<-df[df$pct_vacc_med == vacc_coverage[j],]
  
  if(j==1){
    
    # find all the testing frequencies that have a probability of exceeding the threshold less than 0.05
    df_prob_red<-df_pv[df_pv$cross_red_prob<threshold_prob,]
    df_prob_15xred<-df_pv[df_pv$cross_15xred_prob<threshold_prob,]
    df_prob_2xred<-df_pv[df_pv$cross_2xred_prob<threshold_prob,]
    df_prob_3xred<-df_pv[df_pv$cross_3xred_prob<threshold_prob,]
    df_prob_4xred<-df_pv[df_pv$cross_4xred_prob<threshold_prob,]
   
    
    
    # find the minim
    df_min_red<-df_prob_red[which.min(df_prob_red$testing_freq_med),]
    df_min_15xred<-df_prob_15xred[which.min(df_prob_15xred$testing_freq_med),]
    df_min_2xred<-df_prob_2xred[which.min(df_prob_2xred$testing_freq_med),]
    df_min_3xred<-df_prob_3xred[which.min(df_prob_3xred$testing_freq_med),]
    df_min_4xred<-df_prob_2xred[which.min(df_prob_4xred$testing_freq_med),]
    df_no_test<-df_pv[df_pv$testing_freq_med == 0,]
  }
  else{
    df_prob_redi<-df_pv[df_pv$cross_red_prob<threshold_prob,]
    df_prob_15xredi<-df_pv[df_pv$cross_15xred_prob<threshold_prob,]
    df_prob_2xredi<-df_pv[df_pv$cross_2xred_prob<threshold_prob,]
    df_prob_3xredi<-df_pv[df_pv$cross_3xred_prob<threshold_prob,]
    df_prob_4xredi<-df_pv[df_pv$cross_4xred_prob<threshold_prob,]
   
    
    
    df_min_redi<-df_prob_redi[which.min(df_prob_redi$testing_freq_med),]
    df_min_15xredi<-df_prob_15xredi[which.min(df_prob_15xredi$testing_freq_med),]
    df_min_2xredi<-df_prob_2xredi[which.min(df_prob_2xredi$testing_freq_med),]
    df_min_3xredi<-df_prob_3xredi[which.min(df_prob_3xredi$testing_freq_med),]
    df_min_4xredi<-df_prob_4xredi[which.min(df_prob_4xredi$testing_freq_med),]
    df_no_testi<-df_pv[df_pv$testing_freq_med == 0,]
    
    df_min_red<-rbind(df_min_red, df_min_redi)
    df_min_15xred<-rbind(df_min_15xred, df_min_15xredi)
    df_min_2xred<-rbind(df_min_2xred, df_min_2xredi)
    df_min_3xred<-rbind(df_min_3xred, df_min_3xredi)
    df_min_4xred<-rbind(df_min_4xred, df_min_4xredi)
    df_no_test<-rbind(df_no_test, df_no_testi)
    
    
  }
}
  df_no_test<-df_no_test%>%mutate(thres = "none")
  df_min_red<-df_min_red%>%mutate(thres = "CDC high")
  df_min_15xred<-df_min_15xred%>%mutate(thres = "1.5x CDC high")
  df_min_2xred<-df_min_2xred%>%mutate(thres = "2x CDC high")
  df_min_3xred<-df_min_3xred%>%mutate(thres = "3x CDC high")
  df_min_4xred<-df_min_4xred%>%mutate(thres = "4x CDC high")
  df_comb<-rbind( df_min_red, df_min_15xred, df_min_2xred, df_min_3xred, df_min_4xred, df_no_test)
  
  return(df_comb)
}

get_pct_infectiousness_removed<-function(V0=3, Vf=6, delta_add=0, t0topeak=2, Vpeak = 9,
                                         tpeaktof=7, LOD=5, LOI=6, dt=0.01, t_test, test_sensitivity){
  
  
  t<-seq(0, 20, dt)
  
  
  # Viral load trajectory
  v<-numeric(length(t))
  growth_rate<-(Vpeak-V0)/t0topeak
  decline_rate<-(Vpeak-Vf)/tpeaktof
  v[t<=t0topeak]<-V0+ growth_rate*(t[t<=t0topeak])
  v[t>t0topeak]<- Vpeak-decline_rate*(t[t>t0topeak]-t0topeak)
  v_i<-v[v>LOI]
  t_rel<-t[v>LOI]
  t_inf<-t_rel- t_rel[1]
  duration_inf<-t_inf[length(t_inf)]
  prop_inf_ind_removed<-(duration_inf-t_test)/duration_inf
  
  v_inf<-v_i-LOI 
  # Find the sum of the AUCs at each stage of infetion (here assuming infectious individuals equally likely to be in any stage)
  AUC_remaining<-rep(0, length(v_inf))
  for (i in 1:length(v_inf)){
    t_stage<-t_inf[i]
    v_rem<-v_inf[t_inf>t_stage]
    AUC_remaining[i]<-sum((v_rem)*dt)
  }
  sum_AUC_all_stages<-sum(AUC_remaining)
  
  # Find the remaining AUC at each stage of infection only including those who would have been caught from a test
  # sum them up, this is the numerator
  v_caught<-v_i[t_inf>t_test] - LOI
  t_caught<-t_inf[t_inf>t_test]
  AUC_remaining_test<-rep(0, length(t_caught))
  for (i in 1:length(v_caught)){
    t_stage<-t_caught[i]
    v_rem_test<-v_caught[t_caught>t_stage]
    AUC_remaining_test[i]<-sum((v_rem_test)*dt)
  }
  sum_AUC_post_test<-sum(AUC_remaining_test)
  pct_removed<-sum_AUC_post_test/sum_AUC_all_stages
  
  eff_percent_removed = test_sensitivity*pct_removed
  df<-data.frame(eff_percent_removed, pct_removed, prop_inf_ind_removed)
  return(df)
  
}

get_daily_p_local<-function(test_data, ARRIVAL_DATE){
  # Convert local infections vector into incidence vector
  days_from_arrival<-as.numeric((test_data[length(test_data[,1]), "ReportedDate"]) - as.Date(ARRIVAL_DATE))
  days_vec<-c(0:days_from_arrival)
  N_SAMPLES <-1000
  time_to_symptoms <- rlnorm(N_SAMPLES, meanlog = 1.54, sdlog = 0.47)
  time_to_test <- rexp(N_SAMPLES, rate = 1/2)
  days_after_infection <- round(time_to_symptoms + time_to_test)
  
  # Store whether sample is local or not
  samples_local<-rep(0, nrow(test_data))
  
  
  ## Run samples ##
  for(i in 1:N_SAMPLES){
    infection_times<-days_vec - sample(days_after_infection, size = length(days_vec), replace = TRUE)
    
    # Any negative numbers are importations, positive local
    for (m in 1:nrow(test_data)) {
      if (infection_times[m] >= 0) {
        samples_local[m] <- samples_local[m] + 1
      }
    }
  }
  
  ## Calculate sample proportions for each day 
  daily_p_local<-samples_local/N_SAMPLES
  return(daily_p_local)
}




Rt_fxn_cases <- function(cases, case_data, arrival, last_reliable_day, MEAN_GI,STD_GI){
  daily_p_local<-get_daily_p_local(case_data, arrival)
  #days <- days_vec
  Q2p5 <- c()
  median <- c()
  Q97p5 <- c()
  BINOM_SAMPLES = 100
  RT_SAMPLES = 100
  local_sample_store <- data.frame()
  import_sample_store <- data.frame()
  ## Binomial samples for each day ##
  for (i in 1:length(daily_p_local)) {
    binom_sample <- rbinom(BINOM_SAMPLES, size = cases[i], prob = daily_p_local[i])
    import_cases <- rep(cases[i], BINOM_SAMPLES) - binom_sample
    
    Q2p5 <- append(Q2p5, quantile(binom_sample, 0.025, na.rm = TRUE))
    median <- append(median, quantile(binom_sample, 0.5, na.rm = TRUE))
    Q97p5 <- append(Q97p5, quantile(binom_sample, 0.975, na.rm = TRUE))
    
    # Store samples for Rt calculation
    colname <- toString(i-1)
    if (i == 1) {
      local_sample_store <- data.frame(binom_sample)
      colnames(local_sample_store) <- colname
      import_sample_store <- data.frame(import_cases)
      colnames(import_sample_store) <- colname
    } else {
      local_sample_store[[colname]] <- binom_sample
      import_sample_store[[colname]] <- import_cases
    }
  }
  
  # Dataframe for local vs. imported detected each day
  df <- data.frame(Q2p5, median, Q97p5, daily_p_local)
  df$daily_p_import <- 1 - daily_p_local
  
  
  SAMPLES_FROM_RT_DIST <- 1000
  
  for (i in 1:RT_SAMPLES) {
    sample_row_index <- round(runif(1, min = 1, max = BINOM_SAMPLES))
    local <- as.numeric(local_sample_store[sample_row_index,])
    import <- as.numeric(import_sample_store[sample_row_index,])
    incid <- data.frame(local, import)
    colnames(incid) <- c("local", "imported")
    
    # Calculate Rt
    config <- make_config(list(mean_si = MEAN_GI, std_si = STD_GI))
    Rt <- estimate_R(incid, method = "parametric_si", config = config)
    
    
    # Rename columns
    colnames(Rt$R)[11] <- "upperbound"
    colnames(Rt$R)[5]  <- "lowerbound"
    colnames(Rt$R)[8]  <- "median"
    colnames(Rt$R)[3]  <- "mean"
    colnames(Rt$R)[4]  <- "sd"
    
    # Make dataframe: each col is sample, rows are days
    if (i == 1) {
      Rt_means <- data.frame(Rt$R$t_end)
      Rt_stds <- data.frame(Rt$R$t_end)
      colnames(Rt_means) <- "time"
      colnames(Rt_stds)   <- "time"
    }
    colname <- toString(i)
    Rt_means[[colname]] <- Rt$R$mean
    Rt_stds[[colname]]   <- Rt$R$sd
    
  }
  
  Rt_times <- Rt_means$time
  Rt_medians <- rep(0, length(Rt_times))
  Rt_upperbounds  <- rep(0, length(Rt_times))
  Rt_lowerbounds <- rep(0, length(Rt_times))
  Rt_mean <- rep(0, length(Rt_times))
  Rt_var <- rep(0, length(Rt_times))
  Rt_a <- rep(0, length(Rt_times))
  Rt_b <- rep(0, length(Rt_times))
  
  # Summarize Rt estimates
  for (i in 1:nrow(Rt_means)) {
    Rt_aggregated_samples <- c()
    # For each Rt sample, use mean and std to 
    #   sample gamma, summarize all gamma draws
    #   for each day
    
    for (j in 2:RT_SAMPLES+1) {
      mean <- Rt_means[i, j]
      std  <- Rt_stds[i, j]
      var  <- std ** 2
      b <- mean / var
      a <- mean * b
      Rt_aggregated_samples <- append(
        Rt_aggregated_samples,
        rgamma(SAMPLES_FROM_RT_DIST, shape = a, rate = b)
      )
    }
    Rt_medians[i] <- quantile(Rt_aggregated_samples, .5)
    Rt_upperbounds[i] <- quantile(Rt_aggregated_samples, .975)
    Rt_lowerbounds[i] <- quantile(Rt_aggregated_samples, .025)
    Rt_mean[i]<-mean
    Rt_a[i]<-a
    Rt_b[i]<-b
    Rt_var[i]<-var
  }
  
  tvec = Rt_times - Rt_times[1]
  
  
  dates<-seq(as.Date(arrival), as.Date(last_reliable_day), "days")
  
  Rt_summary <- data.frame(dates)
  Rt_summary['Rt_times']<-Rt_times[1:length(dates)]
  Rt_summary['Rt_medians']<-Rt_medians[1:length(dates)]
  Rt_summary['Rt_lowerbounds']<-Rt_lowerbounds[1:length(dates)]
  Rt_summary['Rt_upperbounds']<-Rt_upperbounds[1:length(dates)]
  Rt_summary['Rt_mean']<-Rt_mean[1:length(dates)]
  Rt_summary['Rt_var']<-Rt_var[1:length(dates)]
  Rt_summary['Rt_a']<-Rt_a[1:length(dates)]
  Rt_summary['Rt_b']<-Rt_b[1:length(dates)]
  Rt_summary<-Rt_summary%>%filter(dates<=last_reliable_day-7)
  return(Rt_summary)
}

  
