
###  Statistics  ########################################
get_stat<-function(data,covar,mediator,exposure,interactionterm,outcome,M=NULL,O=NULL){
  
  
  # Pick only the covariates
  df_cov <- data[ , covar, drop = FALSE]
  
  # Loop over every column and get the mean
  mcv<-sapply(df_cov, mean, na.rm = T)
  
  form_M <- paste(mediator, "~", exposure, "+", paste(covar, collapse = " + "))
  
  if(M==1){
  
    outm.fit<-glm(as.formula(form_M), family=binomial(link=logit), data = data) 
    ssm<-NULL  
  }
  
  if(M==0){
    outm.fit<-lm(as.formula(form_M), data = data)
    ssm<-(summary(outm.fit)$sigma)**2
  }
  
  bcv<-sapply(covar, function(x) coef(summary(outm.fit))[x,"Estimate"])
  b0<-coef(summary(outm.fit))["(Intercept)","Estimate"]
  b1<-coef(summary(outm.fit))[exposure,"Estimate"]
  bcc<-sum(bcv*mcv)
  
  
  #browser()
  form_Y <- paste( outcome, "~", exposure, "+", mediator, "+", interactionterm, "+", paste(covar, collapse = " + "))
  if(O==1){
    outy.fit<-glm(as.formula(form_Y), family=binomial(link=logit), data = data)
  
  }
  if(O==0){
    outy.fit<-glm(as.formula(form_Y), family=gaussian(link=identity), data = data)
    
  }
  
  
  
  t1<-coef(summary(outy.fit))[exposure,"Estimate"]
  t2<-coef(summary(outy.fit))[mediator, "Estimate"]
  t3<-coef(summary(outy.fit))[paste(exposure,mediator,sep=':'), "Estimate"]
  
  return(c( bcc,b0,t1,t2,t3,b1,ssm ))
  
}


###  Bootstrap  ########################################
boot.bMbO<- function(data, indices)
{
  
  outcome<- Y
  exposure<- A
  mediator<- M
  covar<-COVAR
 
  interactionterm<- paste(mediator,exposure,sep='*')
  
  data<-data[indices, ]     
  
  statistics<-get_stat(data,covar,mediator,exposure,interactionterm,outcome,M=1,O=1)
   
  bcc<-statistics[1]
  b0<-statistics[2]
  t1<-statistics[3]
  t2<-statistics[4]
  t3<-statistics[5]
  b1<-statistics[6]
  
  total<- exp(t1*a)*(1+exp(b0+b1*astar+bcc))*(1+exp(b0+b1*a+bcc+t2+t3*a)) / ( exp(t1*astar)*(1+exp(b0+b1*a+bcc))* (1+ exp(b0+b1*astar+bcc+t2+t3*astar)))
  
  cde_comp1<- exp(t1*(a - astar) + t2*mstar + t3*a*mstar)
  cde_comp2<-(1+exp(b0+b1*astar+bcc))
  cde_comp3<-(1+exp(b0+b1*astar+bcc+t2+t3*astar))
  cde_comp4<- exp(t2*mstar+t3*astar*mstar)
  cde_comp5<-(1+exp(b0+b1*astar+bcc))
  cde_comp6<-(1+exp(b0+b1*astar+bcc+t2+t3*astar))
  cde_comp<- cde_comp1*cde_comp2/cde_comp3-cde_comp4*cde_comp5/cde_comp6
  
  intref_comp1<- exp(t1*(a-astar))
  intref_comp2<-(1+exp(b0+b1*astar+bcc+t2+t3*a))
  intref_comp3<-(1+exp(b0+b1*astar+bcc+t2+t3*astar))
  intref_comp4<- exp(t1*(a-astar)+t2*mstar+t3*a*mstar)
  intref_comp5<-(1+exp(b0+b1*astar+bcc))
  intref_comp6<-(1+exp(b0+b1*astar+bcc+t2+t3*astar))
  intref_comp7<- exp(t2*mstar+t3*astar*mstar)
  intref_comp8<-(1+exp(b0+b1*astar+bcc))
  intref_comp9<-(1+exp(b0+b1*astar+bcc+t2+t3*astar))
  intref_comp<-intref_comp1*intref_comp2/intref_comp3-1-intref_comp4*intref_comp5/intref_comp6 + intref_comp7*intref_comp8/intref_comp9
  
  intmed_comp1<- exp(t1*(a-astar))
  intmed_comp2<-(1+exp(b0+b1*a+bcc+t2+t3*a))
  intmed_comp3<-(1+exp(b0+b1*astar+bcc))
  intmed_comp4<-(1+exp(b0+b1*astar+bcc+t2+t3*astar))
  intmed_comp5<-(1+exp(b0+b1*a+bcc))
  intmed_comp6<-(1+exp(b0+b1*a+bcc+t2+t3*astar))
  intmed_comp7<-(1+exp(b0+b1*astar+bcc)) 
  intmed_comp8<-(1+exp(b0+b1*astar+bcc+t2+t3*astar))
  intmed_comp9<-(1+exp(b0+b1*a+bcc))
  intmed_comp10<- exp(t1*(a-astar))
  intmed_comp11<-(1+exp(b0+b1*astar+bcc+t2+t3*a))
  intmed_comp12<-(1+exp(b0+b1*astar+bcc+t2+t3*astar))
  intmed_comp<- intmed_comp1*intmed_comp2*intmed_comp3/(intmed_comp4*intmed_comp5)-intmed_comp6*intmed_comp7/(intmed_comp8*intmed_comp9)-intmed_comp10*intmed_comp11/intmed_comp12 + 1
  
  pie_comp1<-(1+exp(b0+b1*astar+bcc))
  pie_comp2<-(1+exp(b0+b1*a+bcc+t2+t3*astar))
  pie_comp3<-(1+exp(b0+b1*a+bcc))
  pie_comp4<-(1+exp(b0+b1*astar+bcc+t2+t3*astar))
  pie_comp<-(pie_comp1*pie_comp2)/(pie_comp3*pie_comp4) - 1
  
  terr<- cde_comp + intref_comp + intmed_comp + pie_comp
  
  excess_rr_total<- total - 1
  excess_rr_cde<- cde_comp*(total - 1)/terr
  excess_rr_intref<- intref_comp*(total - 1)/terr
  excess_rr_intmed<- intmed_comp*(total - 1)/terr
  excess_rr_pie<- pie_comp*(total - 1)/terr
  
  prop_cde<- cde_comp/terr
  prop_intmed<- intmed_comp/terr
  prop_intref<- intref_comp/terr
  prop_pie<- pie_comp/terr
  prop_med<- (pie_comp + intmed_comp)/terr
  prop_int<- (intmed_comp + intref_comp)/terr
  prop_elm<- (pie_comp + intmed_comp + intref_comp)/terr
  
  ests<- c(total, excess_rr_total, excess_rr_cde, excess_rr_intref, excess_rr_intmed, excess_rr_pie, 
           prop_cde, prop_intref, prop_intmed, prop_pie, prop_med, prop_int, prop_elm)
  
  return(ests)
}

boot.cMcO<- function(data, indices)
{
  data<-data[indices, ]
  
  outcome<- Y
  exposure<- A
  mediator<- M
  covar<-COVAR
  interactionterm<- paste(mediator,exposure,sep='*')
  
  
  statistics<-get_stat(data,covar,mediator,exposure,interactionterm, outcome,M=0,O=0)
  
  bcc<-statistics[1]
  b0<-statistics[2]
  t1<-statistics[3]
  t2<-statistics[4]
  t3<-statistics[5]
  b1<-statistics[6]
  
  cde<- (t1+t3*mstar)*(a-astar)
  intref<- t3*(b0+b1*astar+bcc-mstar)*(a-astar)
  intmed<- t3*b1*(a-astar)*(a-astar)
  pie_comp<-(t2*b1+t3*b1*astar)*(a-astar)
  te<-cde+ intref + intmed + pie_comp
  prop_cde<- cde/te
  prop_intref<- intref/te
  prop_intmed<- intmed/te
  prop_pie<-pie_comp/te
  prop_med<- (pie_comp + intmed)/te
  prop_int<- (intmed + intref)/te
  prop_elm<- (pie_comp + intmed + intref)/te
  
  ests<- c(te, cde, intref, intmed, pie_comp, prop_cde, prop_intref, prop_intmed, prop_pie, prop_med, prop_int, prop_elm)
  
  return(ests)
}

boot.cMbO<- function(data, indices)
{
  data<-data[indices, ]
  
  outcome<- Y
  exposure<- A
  mediator<- M
  covar<-COVAR
  
  interactionterm<- paste(mediator,exposure,sep='*')
  
  
  statistics<-get_stat(data,covar,mediator,exposure,interactionterm, outcome,M=0,O=1)
  
  bcc<-statistics[1]
  b0<-statistics[2]
  t1<-statistics[3]
  t2<-statistics[4]
  t3<-statistics[5]
  b1<-statistics[6]
  ssm<-statistics[7]
   
  
  
  total<-exp(t1 + t2*b1 + (t3*(b0 + b1*astar + b1*a + bcc + t2*ssm)*(a - astar)) + (0.5*(t3*t3)*ssm*((a*a)-(astar*astar))))
  #browser()
  cde_comp1<-t1*(a - astar) + t2*mstar + t3*a*mstar
  cde_comp2<-(t2 + t3*astar)*(b0 + b1*astar + bcc)
  cde_comp3<-0.5*(t2 + t3*astar)*(t2 + t3*astar)*ssm
  cde_comp4<-t2*mstar + t3*astar*mstar
  cde_comp5<-(t2 + t3*astar)*(b0 + b1*astar + bcc)
  cde_comp6<-0.5*(t2 + t3*astar)*(t2 + t3*astar)*ssm
  cde_comp<-exp(cde_comp1 - cde_comp2 - cde_comp3) - exp(cde_comp4 - cde_comp5 - cde_comp6)
  
  
  
  intref_comp1<-(t1 + t3*(b0 + b1*astar + bcc + t2*ssm)*(a - astar)) + (0.5*(t3*t3)*ssm*((a*a)-(astar*astar)))
  intref_comp2<-(t1*(a - astar) + t2*mstar + t3*a*mstar) - ((t2 + t3*astar)*(b0 + b1*astar + bcc)) - (0.5*(t2 + t3*astar)*(t2 + t3*astar)*ssm)
  intref_comp3<-(t2*mstar + t3*astar*mstar) - ((t2 + t3*astar)*(b0 + b1*astar + bcc)) - (0.5*(t2 + t3*astar)*(t2 + t3*astar)*ssm)
  intref_comp<-exp(intref_comp1) - 1  - exp(intref_comp2) + exp(intref_comp3)
  
  intmed_comp1<-(t1 + t2*b1 + t3*(b0 + b1*astar + b1*a + bcc + t2*ssm))*(a - astar) + (0.5*(t3*t3)*ssm*((a*a)-(astar*astar)))
  intmed_comp2<-(t2*b1 + t3*b1*astar)*(a - astar)
  intmed_comp3<-(t1 + t3*(b0 + b1*astar + bcc + t2*ssm))*(a - astar) + (0.5*(t3*t3)*ssm*((a*a)-(astar*astar)))
  intmed_comp<-exp(intmed_comp1) - exp(intmed_comp2) - exp(intmed_comp3) + 1
  
  pie_comp<-exp(((t2*b1 + t3*b1*astar)*(a - astar))) - 1
  
  terr<-cde_comp + intref_comp + intmed_comp + pie_comp
  
  excess_rr_total<-total - 1
  excess_rr_cde<-cde_comp*(total - 1)/terr
  excess_rr_intref<-intref_comp*(total - 1)/terr
  excess_rr_intmed<-intmed_comp*(total - 1)/terr
  excess_rr_pie<-pie_comp*(total - 1)/terr
  
  prop_cde<-cde_comp/terr
  prop_intref<-intref_comp/terr
  prop_intmed<-intmed_comp/terr
  prop_pie<-pie_comp/terr
  prop_med<-(pie_comp + intmed_comp)/terr
  prop_int<-(intmed_comp + intref_comp)/terr
  prop_elm<-(pie_comp + intmed_comp + intref_comp)/terr
  
 
  
  ests<- c(total, excess_rr_total, excess_rr_cde, excess_rr_intref, excess_rr_intmed, excess_rr_pie, 
           prop_cde, prop_intref, prop_intmed, prop_pie, prop_med, prop_int, prop_elm)
   return(ests)
}


boot.bMcO <- function(data, indices)
{
  data<-data[indices, ]
  
  outcome<- Y
  exposure<- A
  mediator<- M
  covar<-COVAR
  interactionterm<- paste(mediator,exposure,sep='*')
  
  
  statistics<-get_stat(data,covar,mediator,exposure,interactionterm, outcome ,M=1,O=0)
  
  bcc<-statistics[1]
  b0<-statistics[2]
  t1<-statistics[3]
  t2<-statistics[4]
  t3<-statistics[5]
  b1<-statistics[6]
  
  cde<-(t1+t3*mstar)*(a-astar)
  
  intref_comp1<- t3*(a-astar)
  intref_comp2<- exp(b0+b1*astar+bcc)
  intref_comp3<- (1+exp(b0+b1*astar+bcc))-mstar
  intref<- intref_comp1*(intref_comp2/intref_comp3)
  
  intmed_comp1<- t3*(a-astar)
  intmed_comp2<- exp(b0+b1*a+bcc)
  intmed_comp3<- (1+exp(b0+b1*a+bcc))
  intmed_comp4<- exp(b0+b1*astar+bcc)
  intmed_comp5<- (1+exp(b0+b1*astar+bcc))
  intmed<- intmed_comp1*(intmed_comp2/intmed_comp3-intmed_comp4/intmed_comp5)
  
  pie_comp1<-(t2+t3*astar)
  pie_comp2<- exp(b0+b1*a+bcc)
  pie_comp3<-(1+exp(b0+b1*a+bcc))
  pie_comp4<- exp(b0+b1*astar+bcc)
  pie_comp5<-(1+exp(b0+b1*astar+bcc))
  pie_comp<- pie_comp1*(pie_comp2/pie_comp3-pie_comp4/pie_comp5)
  
  te<- cde + intref + intmed + pie_comp
  
  prop_cde<- cde/te
  prop_intref<- intref/te
  prop_intmed<- intmed/te
  prop_pie<- pie_comp/te
  prop_med<- (pie_comp + intmed)/te
  prop_int<- (intmed + intref)/te
  prop_elm<- (pie_comp + intmed + intref)/te
  
  ests<- c(te, cde, intref, intmed, pie_comp, prop_cde, prop_intref, prop_intmed, prop_pie, prop_med, prop_int, prop_elm)
  
  return(ests)
}



###  Results  ########################################

save_results<-function(boot_function=NULL,output=NULL, N=NULL){
  
  #set a seed so that can repeat resampling and get the same results
  set.seed(3000)
  system.time(results<- boot(data = data, statistic = boot_function, R = N))
  
  #collecting results so that they can be output into a single file
   
  df_list<-list()
  
  if (outcome==1) {
    
    ESTIMAND<-c(
      "Total Effect Risk Ratio", "Total Excess Relative Risk", "Excess Relative Risk due to CDE", "Excess Relative Risk due to INTref",
      "Excess Relative Risk due to INTmed","Excess Relative Risk due to PIE", "Proportion CDE", "Proportion INTref",
      "Proportion INTmed", "Proportion PIE", "Overall Proportion Mediated",
      "Overall Proportion Attributable to Interaction","Overall Proportion Eliminated"
    )
    
    index=13
    
    }
  else if (outcome==0) {
    
    ESTIMAND<-c(
      "Total Effect", "CDE", "INTref", "INTmed",
      "PIE","Proportion CDE","Proportion INTref",
      "Proportion INTmed", "Proportion PIE", "Overall Proportion Mediated",
      "Overall Proportion Attributable to Interaction","Overall Proportion Eliminated"
    )
    index=12
    
    }
  
  for (i  in 1:index){
    r<-boot.ci(results, type = "basic", index=i) 
    estq<-r[[2]]
    cis<- r[[4]]
    colnames(cis)<-NULL
    lcl<-cis[1,4]
    ucl<-cis[1,5]
    estimand<-ESTIMAND[i]
    rdf<-data.frame(estimand, estq, lcl, ucl)
    names(rdf)<-c("Estimand", "Estimate", "LCL", "UCL")
    df_list[[i]]<-rdf
  }
  
  # Putting all this together
  allresults<-rbindlist(df_list)
  write.csv(allresults, output)
  
}