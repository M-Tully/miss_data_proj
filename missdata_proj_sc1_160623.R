# The very simple model
# 15-06-23  MTully
#-------------------------

library(lubridate)
library(dplyr)
library(mice) # for MI
library(mitml) #for pooling MI
library(ggpubr) # for arranging multiple plots
setwd("~/Documents/Masters/Miss_Data_Proj/Simulation_2023/")

invlogit <- function(x){return(exp(x)/(1+exp(x)))}

simu <- function(n = 306,
                 
                 #  IcE occurrence variables:   TBC !
                 rho_0 = 0.5,
                 
                 #Baseline DCS variables
                 beta_0 = 77, #constant
                 sigmasq_1 = 15, #error dist term
                 
                 #Followup DCS variables
                 beta_6 = 20, #constant
                 beta_7 = 0.29, #effect of baseline DCS
                 beta_8 = -4.7, #intervention effect   !!-----------goal variable------!!
                 sigmasq_2 = 0.7, #error dist term

                 #ice effect:
                 mu_theta = 0.8, #? to vary
                 sigmasq_theta= 0.1, #? to vary
                 
                 change="additive",
                 seed = 123,
                 print=F,
                 seedata=F){
  set.seed(seed)
  
  M = 100*rho_0
  
  #-------------------------------------------#
  #           1. setup dataset                # 
  #-------------------------------------------#
  
  data <- data.frame(id = 1:n,
                     arm=0, #arm assignment
                     y0 = NA, #baseline DCS
                     y1 = NA, #followup DCS
                     p=0 #start with IcE var (p) set to 'not occurred' for all
  )

  #pick half of participant at random to be in arm 1
  stratify_arm <- sample(1:n, floor(n/2), prob=rep(1/n,n))
  data$arm[stratify_arm] <- 1
  
  #simulate first observations
  data$y0 <- beta_0 + rnorm(n, 0, sigmasq_1)
  
  #-------------------------------------------#
  #            2. IcE occurence               # 
  #-------------------------------------------#
  
  #simulate which participants are impacted by IcE (), probabilities are equal for all
  ice_cal <- runif(n) #random value
  ice_cutoff <- rho_0 # cutoff probability of IcE occuring
  data$p <- ifelse(ice_cal<ice_cutoff, 0, 1) 
  
  #simulate outcome scores
  data$y1 <- beta_6 + beta_7*data$y0 + beta_8*data$arm + rnorm(n, 0, sigmasq_2)
  
  #change value can be additive or multiplicative
  change_val <- rnorm(n, mu_theta, sigmasq_theta)
  
  #create blank data$y1_star variable for observed variable
  data$y1_star <- NA
  
  if(change=="additive"){
    data$y1_star <- data$y1 + change_val*data$p }
  
  else if(change=="multipl"){
    data$y1_star <- data$y1 * change_val^(data$p) }
  
  #new dataset for MI analyisis
  data_mi <- data
  
  #set all IcE affected values to missing for the imputation dataset
  data_mi$y1_star[data_mi$p==0]<- NA
  
  #-------------------------------------------#
  #         3. regular missingness            # 
  #-------------------------------------------#

  #-------------------------------------------#
  #         4. complete case                  # 
  #-------------------------------------------# 
  
  data_cc <- data[data$p==0,]
  
  if(print==T){print(summary(data_cc$diff[data_cc$arm==1]))
    print(summary(data_cc$diff[data_cc$arm==0]))}
  if(seedata==T){return(data)}
  else{
    
    #run regression on complete cases
    mod_cc <- lm(y1_star ~ y0 + arm, data_cc)
    # stratification var not in model I assume?
    
    #-------------------------------------------#
    #        5. impute missingness              # 
    #-------------------------------------------#
    
    data_imp_a0 <- mice(data = data_mi[data_mi$arm==0,c("arm", "y0", "y1_star")], m = M, #method = "norm", 
                        maxit = 5, print=F)
    data_imp_a1 <- mice(data = data_mi[data_mi$arm==1,c("arm", "y0", "y1_star")], m = M, #method = "norm", 
                        maxit = 5, print=F)
    #stratified by ARM
    
    
    imp_res <- data.frame(iter = 1:M,
                          est = NA,
                          sd = NA,
                          tval = NA)
    imp_mods <- list()
    for (m in 1:M) {
      
      # Extract the mth imputed dataset
      data_m_a0 <- complete(data_imp_a0, m)
      data_m_a1 <- complete(data_imp_a1, m)
      
      data_m <- rbind(data_m_a0, data_m_a1)
      
      imp_mods[[m]] <- lm(y1_star ~ y0 + arm, data_m)
      mod_m <- summary(lm(y1_star ~ y0 + arm, data_m))$coefficients
      
      if(print==T){ print( mod_m) }
      imp_res$est[m] <-  mod_m[2,1]
      imp_res$sd[m] <-  mod_m[2,2]
      imp_res$tval[m] <-  mod_m[2,3]
      
    }
    
    ## Combine the estimates
    MI_est <- testEstimates(imp_mods, var.comp = TRUE)
    
    if(print==T){print(confint(MI_est))}
    
    results <- list("mod_cc" = mod_cc,
                    "mod_imp" = imp_res,
                    "MI" = MI_est,
                    "data" = data)
    
    return(results)
    
  }}


#test it works:
d <- simu(n = 30)

#function to run the above simulation over and over and collect key results
run_many_simu <- function(n_iter = 100,
                          rho_0 = 0.5,
                          beta_0 = 77, #constant
                          sigmasq_1 = 15, #error dist term
                          beta_6 = 20, #constant
                          beta_7 = 0.29, #effect of baseline DCS
                          beta_8 = -4.7, #intervention effect   !!-----------goal variable------!!
                          sigmasq_2 = 0.7, #error dist term
                          mu_theta = 0.8, #? to vary
                          sigmasq_theta= 0.1, #? to vary
                          change="additive"){
  
  res <- data.frame(iter = 1:n_iter)
  
  print(paste("Simulation is starting at",Sys.time()))
  
  starttime <- Sys.time()
  
  breaks = floor(seq(1,n_iter, length.out=11))
  
  for(i in 1:n_iter){
    if(i %in% breaks){print(paste("Starting sim ", i, " of ", n_iter, " -- ", round(100*max(0,(i-1))/n_iter), "% complete", sep=""))}
    if(i == breaks[2]){tt = difftime(Sys.time(),starttime)
    print(paste("****     Estimated total runtime: ", round(10.9*(tt),3), " ",units(tt) ,"    ****",sep=""))}
    
    output <- simu(seed=i,
                   rho_0 = rho_0,
                   beta_0 = beta_0, #constant
                   sigmasq_1 = sigmasq_1, #error dist term
                   beta_6 = beta_6, #constant
                   beta_7 = beta_7, #effect of baseline DCS
                   beta_8 = beta_8, #intervention effect   !!-----------goal variable------!!
                   sigmasq_2 = sigmasq_2, #error dist term
                   mu_theta = mu_theta, #? to vary
                   sigmasq_theta= sigmasq_theta, #? to vary
                   change=change)
    mi_res <- output[["MI"]]
    cc_res <- output[["mod_cc"]]
    
    print(mi_res$estimates)
    print(confint(mi_res))
    
    res$estim_mi[i] <- mi_res$estimates[3,1]
    res$sd_mi[i] <- mi_res$estimates[3,2]
    res$CI_low_mi[i] <- confint(mi_res)[3,1]
    res$CI_upp_mi[i] <- confint(mi_res)[3,2]
    res$estim_cc[i] <- summary(cc_res)$coefficients[3,1]
    res$sd_cc[i] <- summary(cc_res)$coefficients[3,2]
    res$CI_low_cc[i] <- confint(cc_res)[3,1]
    res$CI_upp_cc[i] <- confint(cc_res)[3,2]
  }
  
  endtime <- Sys.time()
  
  t=difftime(endtime,starttime)
  print(paste("Simulations complete at ", endtime, " total runtime =", round(t,3)," ",units(t) ,sep=""))
  
  return(res)
}

a <- run_many_simu(n_iter = 1000,
                   rho_0 = 0.2)
b <- run_many_simu(n_iter = 1000,
                   rho_0 = 0.4)
c <- run_many_simu(n_iter = 1000,
                   rho_0 = 0.6)

#Plotting ideas:

#Estimates & CI
(combo <- ggarrange(ggplot(a) + 
                      geom_density(aes(x=CI_low_mi), alpha=0.5, fill="lightgreen") + 
                      geom_density(aes(x=CI_upp_mi), alpha=0.5, fill="lightblue") + 
                      geom_density(aes(x=estim_mi), alpha=0.5, fill="grey") + 
                      ggtitle("Pct IcE = 20%") + xlab("") +
                      theme_light() + xlim(-7.2,-2.5) +
                      geom_vline(xintercept=-4.7),
                    ggplot(b) + 
                      geom_density(aes(x=CI_low_mi), alpha=0.5, fill="lightgreen") + 
                      geom_density(aes(x=CI_upp_mi), alpha=0.5, fill="lightblue") + 
                      geom_density(aes(x=estim_mi), alpha=0.5, fill="grey") + 
                      ggtitle("Pct IcE = 40%") + xlab("") +
                      theme_light() + xlim(-7.2,-2.5) +
                      geom_vline(xintercept=-4.7) ,
                    ggplot(c) + 
                      geom_density(aes(x=CI_low_mi), alpha=0.5, fill="lightgreen") + 
                      geom_density(aes(x=CI_upp_mi), alpha=0.5, fill="lightblue") + 
                      geom_density(aes(x=estim_mi), alpha=0.5, fill="grey") + 
                      ggtitle("Pct IcE = 60%") + xlab("") +
                      theme_light() + xlim(-7.2,-2.5) +
                      geom_vline(xintercept=-4.7) #+xlab("Estimate of 'arm' effect")
                    , ncol=1))
ggsave(paste("s1_density_",lubridate::today(),".jpeg"),combo, width=7, height = 9)

c1 <- "#bc4b51"
c2 <- "#5b8e7d"

#coverage?
# cova_cc <- round(100*sum(a$CI_low_mi<-4.7&a$CI_upp_mi>-4.7)/1000)
# cova_mi <- round(100*sum(a$CI_low_cc<-4.7&a$CI_upp_cc>-4.7)/1000)
# 
# covb_cc <- round(100*sum(b$CI_low_mi<-4.7&b$CI_upp_mi>-4.7)/1000)
# covb_mi <- round(100*sum(b$CI_low_cc<-4.7&b$CI_upp_cc>-4.7)/1000)
# 
# covc_cc <- round(100*sum(c$CI_low_mi<-4.7&c$CI_upp_mi>-4.7)/1000)
# covc_mi <- round(100*sum(c$CI_low_cc<-4.7&c$CI_upp_cc>-4.7)/1000)

plota <- ggplot() + geom_boxplot(data=a, aes(x=estim_mi, y="Beta - MI"), fill=c1, alpha=0.5) +
  geom_boxplot(data=a, aes(x=estim_cc, y="Beta - CC"), fill=c2, alpha=0.5) +
  geom_boxplot(data=a, aes(x=CI_low_mi, y="Low Bound - MI"), fill=c1, alpha=0.5) +
  geom_boxplot(data=a, aes(x=CI_low_cc, y="Low Bound - CC"), fill=c2, alpha=0.5) +
  geom_boxplot(data=a, aes(x=CI_upp_mi, y="Upp Bound - MI"), fill=c1, alpha=0.5) +
  geom_boxplot(data=a, aes(x=CI_upp_cc, y="Upp Bound - CC"), fill=c2, alpha=0.5) +
  theme_light() + xlab("") + ylab("") + ggtitle("20% Affected by the IcE") + theme(text = element_text(size = 20)) 
ggsave(paste("s1_20pct-",today(),".jpeg", sep=""), plota, width=8, height=6)
plotb <- ggplot() + geom_boxplot(data=b, aes(x=estim_mi, y="Beta - MI"), fill=c1, alpha=0.5) +
  geom_boxplot(data=b, aes(x=estim_cc, y="Beta - CC"), fill=c2, alpha=0.5) +
  geom_boxplot(data=b, aes(x=CI_low_mi, y="Low Bound - MI"), fill=c1, alpha=0.5) +
  geom_boxplot(data=b, aes(x=CI_low_cc, y="Low Bound - CC"), fill=c2, alpha=0.5) +
  geom_boxplot(data=b, aes(x=CI_upp_mi, y="Upp Bound - MI"), fill=c1, alpha=0.5) +
  geom_boxplot(data=b, aes(x=CI_upp_cc, y="Upp Bound - CC"), fill=c2, alpha=0.5) +
  theme_light() + xlab("") + ylab("") + ggtitle("40% Affected by the IcE") + theme(text = element_text(size = 20)) 
ggsave(paste("s1_40pct-",today(),".jpeg", sep=""), plotb, width=8, height=6)
plotc <- ggplot() + geom_boxplot(data=c, aes(x=estim_mi, y="Beta - MI"), fill=c1, alpha=0.5) +
  geom_boxplot(data=c, aes(x=estim_cc, y="Beta - CC"), fill=c2, alpha=0.5) +
  geom_boxplot(data=c, aes(x=CI_low_mi, y="Low Bound - MI"), fill=c1, alpha=0.5) +
  geom_boxplot(data=c, aes(x=CI_low_cc, y="Low Bound - CC"), fill=c2, alpha=0.5) +
  geom_boxplot(data=c, aes(x=CI_upp_mi, y="Upp Bound - MI"), fill=c1, alpha=0.5) +
  geom_boxplot(data=c, aes(x=CI_upp_cc, y="Upp Bound - CC"), fill=c2, alpha=0.5) +
  theme_light() + xlab("") + ylab("") + ggtitle("60% Affected by the IcE") + theme(text = element_text(size = 20)) 
ggsave(paste("s1_60pct-",today(),".jpeg", sep=""), plotc, width=8, height=6)

#bias:
bias <- ggarrange(ggplot(a) + 
            geom_density(aes(x=c(-4.7-estim_mi)), alpha=0.5, fill="lightgreen") + 
            geom_density(aes(x=-4.7-estim_cc), alpha=0.5, fill="lightblue") + 
            ggtitle("Pct IcE = 20%") +xlab("") +
            theme_light() + xlim(-1,1) +
            geom_vline(xintercept=0),
          ggplot(b) + 
            geom_density(aes(x=-4.7-estim_mi), alpha=0.5, fill="lightgreen") + 
            geom_density(aes(x=-4.7-estim_cc), alpha=0.5, fill="lightblue") + 
            ggtitle("Pct IcE = 40%") +
            theme_light() +xlim(-1,1) + xlab("") +
            geom_vline(xintercept=0) ,
          ggplot(c) + 
            geom_density(aes(x=-4.7-estim_mi), alpha=0.5, fill="lightgreen") + 
            geom_density(aes(x=-4.7-estim_cc), alpha=0.5, fill="lightblue") + 
            ggtitle("Pct IcE = 60%") +
            theme_light() +xlim(-1,1) +#xlab("Bias") +
            geom_vline(xintercept=0) +xlab("Green = MI, Blue = CC"), ncol=1)#, width=5, height = 9)
ggsave(paste("s1_bias-",today(),".jpeg", sep=""), bias, width=4, height=6)

