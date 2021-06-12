#Using Synthetic Data in a Basic Multiverse Analysis -- ANOVA
#Andrew J. Kelly

##--LIBRARIES--##################################################################################
library("tidyverse") #package for data wrangling
library("synthpop") #package for synthesizing the data
library("effectsize") #package for effect sizes

##--READING THE FILE--##################################################################################
#reading the fle
knaster_aov <- read.csv("knaster.csv") #reads file in
knaster_aov #reports the file

#selects relevent columns from the knaster file
knaster_aov <- knaster_aov %>% select(Current.pain, SUMBDI, Physicalfunction, Negativeviewofself) #puts selected columns into data frame
knaster_aov #reports data frame

##--BDI - FILTERING--##################################################################################

#creating severe category
knaster_aov$SUMBDI[knaster_aov$SUMBDI > 28] = "severe"
knaster_aov #reporting to ensure severe category has been established

#creating moderate category
knaster_aov$SUMBDI[knaster_aov$SUMBDI == 20 | 
                     knaster_aov$SUMBDI == 21 |
                     knaster_aov$SUMBDI == 22 | 
                     knaster_aov$SUMBDI == 23 | 
                     knaster_aov$SUMBDI == 24 | 
                     knaster_aov$SUMBDI == 25 | 
                     knaster_aov$SUMBDI == 26 | 
                     knaster_aov$SUMBDI == 27 | 
                     knaster_aov$SUMBDI == 28] = "moderate"
knaster_aov #reporting to ensure moderate category has been established

#creating mild category
knaster_aov$SUMBDI[knaster_aov$SUMBDI == 14 |
                     knaster_aov$SUMBDI == 15 |
                     knaster_aov$SUMBDI == 16 |
                     knaster_aov$SUMBDI == 17 |
                     knaster_aov$SUMBDI == 18 |
                     knaster_aov$SUMBDI == 19] = "mild"
knaster_aov #reporting to ensure mild category has been established

#creating minimal category
knaster_aov$SUMBDI[knaster_aov$SUMBDI == 1 | 
                     knaster_aov$SUMBDI == 2 | 
                     knaster_aov$SUMBDI == 3 | 
                     knaster_aov$SUMBDI == 4 | 
                     knaster_aov$SUMBDI == 5 | 
                     knaster_aov$SUMBDI == 6 | 
                     knaster_aov$SUMBDI == 7 | 
                     knaster_aov$SUMBDI == 8 | 
                     knaster_aov$SUMBDI == 9 | 
                     knaster_aov$SUMBDI == 10 | 
                     knaster_aov$SUMBDI == 11 | 
                     knaster_aov$SUMBDI == 12 | 
                     knaster_aov$SUMBDI == 13] = "minimal"
knaster_aov #reporting to ensure minimal category has been established

##--ANOVA SECTION--##################################################################################

#linear model for CURRENT PAIN and BDI level
BDI_painlevel <- lm(Current.pain~as.factor(SUMBDI), data = knaster_aov)
aov_painlevelxBDI <- anova(BDI_painlevel) #turns the model into an ANOVA model
aov_painlevelxBDI #reports the outcome

aov_painlevelxBDI_es <- eta_squared(aov_painlevelxBDI) #creates the eta^2 for the model
aov_painlevelxBDI_es #reports the effect size

#linear model for PHYSICAL FUNCTION and BDI level
BDI_physicalfunction <- lm(Physicalfunction~as.factor(SUMBDI), data = knaster_aov)
aov_physicalfunctionxBDI <- anova(BDI_physicalfunction) #turns the model into an ANOVA model
aov_physicalfunctionxBDI #reports the outcome

aov_physicalfunctionxBDI_es <- eta_squared(aov_physicalfunctionxBDI) #creates the eta^2 for the model
aov_physicalfunctionxBDI_es #reports the effect size

#linear model for NEGATIVE VIEW OF SELF and BDI level
BDI_negative <- lm(Negativeviewofself~as.factor(SUMBDI), data = knaster_aov)
aov_negativexBDI <- anova(BDI_negative) #turns the model into an ANOVA model
aov_negativexBDI #reports the outcome

aov_negativexBDI_es <- eta_squared(aov_negativexBDI) #creates the eta^2 for the model
aov_negativexBDI_es #reports the effect size

##--CREATING FRAMES--##################################################################################

#correlation frames
pain_fp <- data.frame(pain_fp=c)
pain_eta <- data.frame(pain_eta=c)
phys_fp <- data.frame(phys_fp=c)
phys_eta <- data.frame(phys_eta=c)
neg_fp <- data.frame(neg_fp=c)
neg_eta <- data.frame(neg_eta=c)

##--PRE-LOOP INFO--##################################################################################

synth_loop <- 1
synth_seed <- 10000

##--LOOP--##################################################################################

#while loop - will loop 2500 times
while(synth_loop < 2501){
  
  synth_seed <- 3000 + synth_loop
  syn_data <- (syn(knaster_aov, seed = synth_seed))$syn #creating a synthetic dataset
  syn_data #reports synthetic dataset
  
  ##--SYNTH ANOVA SECTION--##################################################################################
  
  #linear model for CURRENT PAIN and BDI level
  BDI_painlevel_syn <- lm(Current.pain~as.factor(SUMBDI), data = syn_data)
  aov_painlevelxBDI_syn <- anova(BDI_painlevel_syn) #turns the model into an ANOVA model
  aov_painlevelxBDI_syn #reports the outcome
  
  aov_painlevelxBDI_syn_es <- eta_squared(aov_painlevelxBDI_syn) #creates the eta^2 for the model
  aov_painlevelxBDI_syn_es #reports the effect size
  
  #linear model for PHYSICAL FUNCTION and BDI level
  BDI_physicalfunction_syn <- lm(Physicalfunction~as.factor(SUMBDI), data = syn_data)
  aov_physicalfunctionxBDI_syn <- anova(BDI_physicalfunction_syn) #turns the model into an ANOVA model
  aov_physicalfunctionxBDI_syn #reports the outcome
  
  aov_physicalfunctionxBDI_syn_es <- eta_squared(aov_physicalfunctionxBDI_syn) #creates the eta^2 for the model
  aov_physicalfunctionxBDI_syn_es #reports the effect size
  
  #linear model for NEGATIVE VIEW OF SELF and BDI level
  BDI_negative_syn <- lm(Negativeviewofself~as.factor(SUMBDI), data = syn_data)
  aov_negativexBDI_syn <- anova(BDI_negative_syn) #turns the model into an ANOVA model
  aov_negativexBDI_syn #reports the outcome
  
  aov_negativexBDI_syn_es <- eta_squared(aov_negativexBDI) #creates the eta^2 for the model
  aov_negativexBDI_syn_es$Eta2_partial #reports the effect size
  
  ##--LOOP - BINDING--##################################################################################
  pain_fp <- rbind(pain_fp, aov_painlevelxBDI_syn$`Pr(>F)`[1])
  pain_eta <- rbind(pain_eta, aov_painlevelxBDI_syn_es$Eta2_partial)
  phys_fp <- rbind(phys_fp, aov_physicalfunctionxBDI_syn$`Pr(>F)`[1])
  phys_eta <- rbind(phys_eta, aov_physicalfunctionxBDI_syn_es$Eta2_partial)
  neg_fp <- rbind(neg_fp, aov_negativexBDI_syn$`Pr(>F)`[1])
  neg_eta <- rbind(neg_eta, aov_negativexBDI_syn_es$Eta2_partial)
  
  #adds numbers to the loop
  synth_loop <- synth_loop + 1
}

##--WRANGLING AFTER LOOP--##################################################################################

#binds all correlations and p-values into a single dataframe
syn_data_df <- do.call(cbind, list(pain_fp, pain_eta,
                                   phys_fp, phys_eta,
                                   neg_fp, neg_eta))

#renaming the columns
colnames(syn_data_df) <- c("pain_fp", "pain_eta",
                           "phys_fp", "phys_eta",
                           "neg_fp", "neg_eta")

#saving the file
write.csv(syn_data_df, "multi_syn_anova.csv")

##--FINAL WRANGLING--##################################################################################

#first instance
final_df <- read.csv("multi_syn_anova.csv")
filtered_final_df <- final_df %>% filter(pain_fp > 0.05 &
                                           phys_fp < 0.05 & phys_eta >= 0.610 & phys_eta <= 0.750 & #0.61 and 0.75 were chosen as the cutoffs as they are within the 90% CI of the original eta^2
                                           neg_fp < 0.05 & neg_eta >= 0.530 & neg_eta <= 0.700) #0.53 and 0.70 were chosen as the cutoffs as they are within the 90% CI of the original eta^2
filtered_final_df

write.csv(filtered_final_df, "multi_syn_anova_FILTERED.csv")
