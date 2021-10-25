calc_mediation<-function(mediator,input_data){
  
 outputdir="/data/gh_gr_agingandobesity_share/life_shared/Results/External_coops/PV240_Dagher_Scholz_headmotion_gwas/Amendment2021/Results/"  
 model=paste0(' # mediator model \n ',
             mediator, ' ~ a*BMI + Age +Sex \n',
             ' # outcome \n ',
             'mFD ~ c*BMI + b*', mediator, ' + Age + Sex \n',
             ' # indirect effect (a*b)\n',
             ' ab := a*b \n',
             ' # total effect \n',
             ' total := c + (a*b) \n',
             ' # direct effect \n',
             ' direct := c \n')
 
 # Run Model
 nrow(input_data)
 #fit <- sem(model, data = input_data) #when running without bootstrap, we have Nobs=913
 fit <- sem(model, data = input_data,se = "bootstrap")

 #Parameter Estimates
 params=parameterestimates(fit, boot.ci.type = "bca.simple", standardized = TRUE)
 write.csv(params, file = paste0(outputdir,"model_results/", mediator, "_coefficients.csv"))
 
 #Fit Measures
 #fitm=fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "rmsea"))
 #write.csv(fitm, file = paste0(outputdir,"model_results/", mediator, "_fitmeasures.csv"))
 
 #Checking normality of residuals of individual models
 lm1=lm(paste0('mFD ~ BMI + Age + Sex +', mediator), data=input_data)
 png(paste0(outputdir, 'diagplots/', mediator,"_outcome.png"))
 diagnostics.plot(lm1)
 dev.off()
 
 vif1=vif(lm1)
 
 lm2=lm(paste0(mediator, ' ~ BMI +  Age +Sex'), data=input_data)
 png(paste0(outputdir, 'diagplots/', mediator, ".png"))
 diagnostics.plot(lm2)
 dev.off()
 

 write.csv(vif1, file = paste0(outputdir,"model_results/", mediator, "_vif.csv"))

 # Check for influential cases (problematic if Cooks distance > 4/n (as referenced in Wikipedia))
 n=nrow(input_data)
 cooksd=cooks.distance(lm1)
 influential_lm1 <- as.numeric(names(cooksd)[(cooksd > 4/n)])
 cooksd=cooks.distance(lm2)
 influential_lm2 <- as.numeric(names(cooksd)[(cooksd > 4/n)])
 
 all_infl=c(influential_lm1, influential_lm2)
 all_infl=all_infl[!duplicated(all_infl)]
 write.csv(length(all_infl), file = paste0(outputdir,"model_results/", mediator, "_#infl_cases.csv"))
 
 data_without_inf=input_data %>% filter(!row_number() %in% all_infl)
 
 fit_w_infl <- sem(model, data = data_without_inf,se = "bootstrap")
 
 #Parameter Estimates
 params=parameterestimates(fit_w_infl, boot.ci.type = "bca.simple", standardized = TRUE)
 write.csv(params, file = paste0(outputdir,"model_results/", mediator, "_coefficients_wo_infl.csv"))
 
 #Fit Measures
 #fitm=fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "rmsea"))
 #write.csv(fitm, file = paste0(outputdir,"model_results/", mediator, "_fitmeasures_wo_infl.csv"))
 
 return(fit)
}