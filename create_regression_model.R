calc_regression<-function(outcome,input_data){
  
 outputdir="/data/gh_gr_agingandobesity_share/life_shared/Results/External_coops/PV240_Dagher_Scholz_headmotion_gwas/Amendment2021/Results/"  
 
 # Run Model
 regress_model=lm(paste0(outcome, '~ BMI + Age + Sex'), data=input_data)

 #Parameter Estimates
 coeff <- as.data.frame(summary(regress_model)$coefficients)
 CI <- as.data.frame(confint(regress_model))
 coeffs <- merge(coeff, CI, all = T, by="row.names")
 write.csv(coeffs, file = paste0(outputdir,"model_results/", outcome, "regression_coefficients.csv"))
 
 #Fit Measures
 null_model=lm(paste0(outcome, '~ Age + Sex'), data=input_data)
 anv=anova(regress_model,null_model)
 
 write.csv(as.data.frame(anv), file = paste0(outputdir,"model_results/", outcome, "_regression_pval.csv"))
 
 #Checking normality of residuals of individual models
 png(paste0(outputdir, 'diagplots/', outcome,"_regression.png"))
 diagnostics.plot(regress_model)
 dev.off()
 
 # Check for influential cases (problematic if Cooks distance > 4/n (as referenced in Wikipedia))
 n=nrow(input_data)
 cooksd=cooks.distance(regress_model)
 influential_regress_model <- as.numeric(names(cooksd)[(cooksd > 4/n)])

 data_without_inf=input_data %>% filter(!row_number() %in% influential_regress_model)
 write.csv(nrow(input_data)-nrow(data_without_inf), file = paste0(outputdir,"model_results/", outcome, "_#infl_cases.csv"))
 
 rm_wo_infl=lm(paste0(outcome, '~ BMI + Age + Sex'), data=data_without_inf)
 
 #Parameter Estimates
 coeff <- as.data.frame(summary(rm_wo_infl)$coefficients)
 CI <- as.data.frame(confint(rm_wo_infl))
 coeffs <- merge(coeff, CI, all = T, by="row.names")
 write.csv(coeffs, file = paste0(outputdir,"model_results/", outcome, "regression_coefficients_wo_infl.csv"))
 
 #Fit Measures
 null_model=lm(paste0(outcome, '~ Age + Sex'), data=data_without_inf)
 anv=anova(rm_wo_infl,null_model)
 
 write.csv(as.data.frame(anv), file = paste0(outputdir,"model_results/", outcome, "_regression_pval_wo_infl.csv"))
 return(regress_model)
}