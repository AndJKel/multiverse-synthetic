#Using Synthetic Data in a Basic Multiverse Analysis -- T-test
#Andrew J. Kelly

##--LIBRARIES--##################################################################################
library("tidyverse") #package for data wrangling
library("synthpop") #package for synthesizing the data
library("effectsize") #package for effect sizes


##--READING THE FILE--##################################################################################
#reading the fle
knaster_t <- read.csv("knaster.csv") #reads file in
knaster_t #reports the file

#selects relevent columns from the knaster file
knaster_t <- knaster_t %>% select(Current.pain, SUMBDI, Physicalfunction, Negativeviewofself) #puts selected columns into data frame
knaster_t #reports data frame

#duplicates the SUMBDI column
knaster_t <- knaster_t %>% mutate(BDI_clinical = SUMBDI)  
knaster_t #reports the duplication of the data frame

##--CLINICAL/SUBCLINICAL FILTERING--##################################################################################

#filters out the BDI scores which were 20 and above (i.e., clinical)
knaster_t$BDI_clinical[knaster_t$BDI_clinical >= 20] = "clinical"
knaster_t #reports the data frame to ensure filtering is valid

#assigns values below 20 as subclinical -- probably a better way to do this
knaster_t$BDI_clinical[knaster_t$BDI_clinical == 1 | 
                           knaster_t$BDI_clinical == 2 |
                           knaster_t$BDI_clinical == 3 |
                           knaster_t$BDI_clinical == 4 |
                           knaster_t$BDI_clinical == 5 |
                           knaster_t$BDI_clinical == 6 |
                           knaster_t$BDI_clinical == 7 |
                           knaster_t$BDI_clinical == 8 |
                           knaster_t$BDI_clinical == 9 |
                           knaster_t$BDI_clinical == 10 |
                           knaster_t$BDI_clinical == 11 |
                           knaster_t$BDI_clinical == 12 |
                           knaster_t$BDI_clinical == 13 |
                           knaster_t$BDI_clinical == 14 |
                           knaster_t$BDI_clinical == 15 |
                           knaster_t$BDI_clinical == 16 |
                           knaster_t$BDI_clinical == 17 |
                           knaster_t$BDI_clinical == 18 |
                           knaster_t$BDI_clinical == 19] = "subclinical"
knaster_t #reports to the dataframe to ensure all values are valid with their new category

##--T-TEST TESTING ALL PAIN, PHYS, AND NEG--##################################################################################

#filtering out severe BDI
BDI_clinical_cat <- knaster_t %>% filter(BDI_clinical == "clinical")
BDI_clinical_cat #reports the data frame to ensure filtering occurred

BDI_subclinical_cat <- knaster_t %>% filter(BDI_clinical == "subclinical")
BDI_subclinical_cat #reports the data frame to ensure filtering occurred

#t-test for CURRENT PAIN between the clinical and subclinical individuals
pain_clinicaltype <- t.test(BDI_clinical_cat$Current.pain, BDI_subclinical_cat$Current.pain, var.equal = TRUE) #does the t-test
pain_clinicaltype #reports the t-test outcomes
pain_clinicaltype_es <- cohens_d(BDI_clinical_cat$Current.pain, BDI_subclinical_cat$Current.pain) #reports the es
pain_clinicaltype_es #reports the effect size

#t-test for PHYSICAL FUNCTIONING between the clinical and subclinical individuals
phys_clinicaltype <- t.test(BDI_clinical_cat$Physicalfunction, BDI_subclinical_cat$Physicalfunction, var.equal = TRUE) #does the t-test
phys_clinicaltype #reports the t-test outcomes
phys_clinicaltype_es <- cohens_d(BDI_clinical_cat$Physicalfunction, BDI_subclinical_cat$Physicalfunction) #reports the es
phys_clinicaltype_es #reports the effect size

#t-test for current pain between the clinical and subclinical individuals
neg_clinicaltype <- t.test(BDI_clinical_cat$Negativeviewofself, BDI_subclinical_cat$Negativeviewofself, var.equal = TRUE) #does the t-test
neg_clinicaltype #reports the t-test outcomes
neg_clinicaltype_es <- cohens_d(BDI_clinical_cat$Negativeviewofself, BDI_subclinical_cat$Negativeviewofself) #reports the es
neg_clinicaltype_es #reports the effect size

##--HIGH/LOW PAIN, PHYS, AND NEG FILTERING--##################################################################################

#current pain
BDI_clinicalxhighp_cat <- BDI_clinical_cat %>% filter(Current.pain >= 5.00)
BDI_clinicalxhighp_cat #reports the data frame to ensure filtering occurred

BDI_subclinicalxhighp_cat <- BDI_subclinical_cat %>% filter(Current.pain >= 5.00)
BDI_subclinicalxhighp_cat #reports the data frame to ensure filtering occurred


#physical functioning
BDI_clinicalxlowph_cat <- BDI_clinical_cat %>% filter(Physicalfunction <= 10)
BDI_clinicalxlowph_cat #reports the data frame to ensure filtering occurred

BDI_subclinicalxlowph_cat <- BDI_subclinical_cat %>% filter(Physicalfunction <= 10)
BDI_subclinicalxlowph_cat #reports the data frame to ensure filtering occurred

#negative self
BDI_clinicalxlown_cat <- BDI_clinical_cat %>% filter(Negativeviewofself <= 9)
BDI_clinicalxlown_cat #reports the data frame to ensure filtering occurred

BDI_subclinicalxlown_cat <- BDI_subclinical_cat %>% filter(Negativeviewofself <= 9)
BDI_subclinicalxlown_cat #reports the data frame to ensure filtering occurred

##--HIGH/LOW PAIN, PHYS, AND NEG TESTING--##################################################################################


#t-test & effect size testing - LOW pain x sub/clinical type
typexlowpain_t <- t.test(BDI_clinicalxlowp_cat$Current.pain, 
                      BDI_subclinicalxlowp_cat$Current.pain, 
                      var.equal = TRUE) #does the t-test
typexlowpain_t #reports the t-test outcomes

typexlowpain_es <- cohens_d(BDI_clinicalxlowp_cat$Current.pain,
                                BDI_subclinicalxlowp_cat$Current.pain) #reports the effect size
typexlowpain_es #reports the effect size



#t-test & effect size testing - HIGH pain x sub/clinical type
typexhighpain_t <- t.test(BDI_clinicalxhighp_cat$Current.pain, 
                         BDI_subclinicalxhighp_cat$Current.pain, 
                         var.equal = TRUE) #does the t-test
typexhighpain_t #reports the t-test outcomes

typexhighpain_es <- cohens_d(BDI_clinicalxhighp_cat$Current.pain,
                            BDI_subclinicalxhighp_cat$Current.pain) #reports the effect size
typexhighpain_es #reports the effect size



#t-test & effect size testing - LOW phys x sub/clinical type
typexlowphys_t <- t.test(BDI_clinicalxlowph_cat$Physicalfunction, 
                         BDI_subclinicalxlowph_cat$Physicalfunction, 
                         var.equal = TRUE) #does the t-test
typexlowphys_t #reports the t-test outcomes

typexlowphys_es <- cohens_d(BDI_clinicalxlowph_cat$Physicalfunction,
                            BDI_subclinicalxlowph_cat$Physicalfunction) #reports the effect size
typexlowphys_es #reports the effect size

#t-test & effect size testing - LOW pain x sub/clinical type
typexlowneg_t <- t.test(BDI_clinicalxlown_cat$Negativeviewofself, 
                         BDI_subclinicalxlown_cat$Negativeviewofself, 
                         var.equal = TRUE) #does the t-test
typexlowneg_t #reports the t-test outcomes

typexlowneg_es <- cohens_d(BDI_clinicalxlown_cat$Negativeviewofself,
                            BDI_subclinicalxlown_cat$Negativeviewofself) #reports the effect size
typexlowneg_es #reports the effect size

##--CREATING FRAMES--##################################################################################

#correlation frames
pain_p_all <- data.frame(pain_p_all=c)
pain_d_all <- data.frame(pain_d_all=c)
pain_p_high <- data.frame(pain_p_high=c)
pain_d_high <- data.frame(pain_d_high=c)

phys_p_all <- data.frame(phys_p_all=c)
phys_d_all <- data.frame(phys_d_all=c)
phys_p_low <- data.frame(phys_p_low=c)
phys_d_low <- data.frame(phys_d_low=c)

neg_p_all <- data.frame(neg_p_all=c)
neg_d_all <- data.frame(neg_d_all=c)
neg_p_low <- data.frame(neg_p_low=c)
neg_d_low <- data.frame(neg_d_low=c)

##--PRE-LOOP INFO--##################################################################################

synth_loop <- 1
synth_seed <- 2000

##--LOOP--##################################################################################

#while loop - will loop 2500 times
while(synth_loop < 2501){
  
  synth_seed <- 3000 + synth_loop
  syn_data <- (syn(knaster_t, seed = synth_seed))$syn #creating a synthetic dataset
  syn_data #reports synthetic dataset
  
  
  ##--LOOP - T-TEST - FILTERING--##################################################################################
  BDI_clinical_syn <- syn_data %>% filter(BDI_clinical == "clinical")
  BDI_subclinical_syn <- syn_data %>% filter(BDI_clinical == "subclinical")

  #syn_clinicalxlowp_cat <- BDI_clinical_syn %>% filter(Current.pain <= 4.99)
  syn_clinicalxhighp_cat <- BDI_clinical_syn %>% filter(Current.pain >= 5.00)
  syn_subclinicalxlowp_cat <- BDI_subclinical_syn %>% filter(Current.pain <= 4.99)
  syn_subclinicalxhighp_cat <- BDI_subclinical_syn %>% filter(Current.pain >= 5.00)
  syn_clinicalxlowph_cat <- BDI_clinical_syn %>% filter(Physicalfunction <= 10)
  syn_subclinicalxlowph_cat <- BDI_subclinical_syn %>% filter(Physicalfunction <= 10)
  syn_clinicalxlown_cat <- BDI_clinical_syn %>% filter(Negativeviewofself <= 9)
  syn_subclinicalxlown_cat <- BDI_subclinical_syn %>% filter(Negativeviewofself <= 9)
  
  ##--LOOP - T-TESTS--##################################################################################
  pain_clinicaltype_syn <- t.test(BDI_clinical_syn$Current.pain, 
                                  BDI_subclinical_syn$Current.pain, 
                                  var.equal = TRUE) #does the t-test
  phys_clinicaltype_syn <- t.test(BDI_clinical_syn$Physicalfunction, 
                                  BDI_subclinical_syn$Physicalfunction, 
                                  var.equal = TRUE) #does the t-test
  neg_clinicaltype_syn <- t.test(BDI_clinical_syn$Negativeviewofself, 
                                 BDI_subclinical_syn$Negativeviewofself, 
                                 var.equal = TRUE) #does the t-test
  highpain_clinicaltype_syn <- t.test(syn_clinicalxhighp_cat$Current.pain, 
                                syn_subclinicalxhighp_cat$Current.pain, 
                                var.equal = TRUE) #does the t-test
  lowphys_clinicaltype_syn <- t.test(syn_clinicalxlowph_cat$Physicalfunction, 
                               syn_subclinicalxlowph_cat$Physicalfunction, 
                               var.equal = TRUE) #does the t-test
  lowneg_clinicaltype_syn <- t.test(syn_clinicalxlown_cat$Negativeviewofself, 
                              syn_subclinicalxlown_cat$Negativeviewofself, 
                              var.equal = TRUE) #does the t-test
  
  ##--LOOP - T- EFFECT SIZES--##################################################################################
  pain_clinicaltype_es_syn <- cohens_d(BDI_clinical_syn$Current.pain, 
                                       BDI_subclinical_syn$Current.pain) #reports the effect size
  phys_clinicaltype_es_syn <- cohens_d(BDI_clinical_syn$Physicalfunction, 
                                       BDI_subclinical_syn$Physicalfunction) #reports the effect size
  neg_clinicaltype_es_syn <- cohens_d(BDI_clinical_syn$Negativeviewofself, 
                                      BDI_subclinical_syn$Negativeviewofself) #reports the effect size
  highpain_clinicaltype_es_syn <- cohens_d(syn_clinicalxhighp_cat$Current.pain,
                                           syn_subclinicalxhighp_cat$Current.pain) #reports the effect size
  lowphys_clinicaltype_es_syn <- cohens_d(syn_clinicalxlowph_cat$Physicalfunction,
                                          syn_subclinicalxlowph_cat$Physicalfunction) #reports the effect size
  lowneg_clinicaltype_es_syn <- cohens_d(syn_clinicalxlown_cat$Negativeviewofself,
                                         syn_subclinicalxlown_cat$Negativeviewofself) #reports the effect size

  ##--LOOP - BINDING--##################################################################################
  pain_p_all <- rbind(pain_p_all, pain_clinicaltype_syn$p.value)
  pain_d_all <- rbind(pain_d_all, pain_clinicaltype_es_syn$Cohens_d)
  pain_p_high <- rbind(pain_p_high, highpain_clinicaltype_syn$p.value)
  pain_d_high <- rbind(pain_d_high, highpain_clinicaltype_es_syn$Cohens_d)
  
  phys_p_all <- rbind(phys_p_all, phys_clinicaltype_syn$p.value)
  phys_d_all <- rbind(phys_d_all, phys_clinicaltype_es_syn$Cohens_d)
  phys_p_low <- rbind(phys_p_low, lowphys_clinicaltype_syn$p.value)
  phys_d_low <- rbind(phys_d_low, lowphys_clinicaltype_es_syn$Cohens_d)
  
  neg_p_all <- rbind(neg_p_all, neg_clinicaltype_syn$p.value)
  neg_d_all <- rbind(neg_d_all, neg_clinicaltype_es_syn$Cohens_d)
  neg_p_low <- rbind(neg_p_low, lowneg_clinicaltype_syn$p.value)
  neg_d_low <- rbind(neg_d_low, lowneg_clinicaltype_es_syn$Cohens_d)

  #adds numbers to the loop
  synth_loop <- synth_loop + 1
}


##--WRANGLING AFTER LOOP--##################################################################################

#binds all correlations and p-values into a single dataframe
syn_data_df <- do.call(cbind, list(pain_p_all, pain_d_all,
                                   pain_p_high, pain_d_high,
                                   phys_p_all, phys_d_all,
                                   phys_p_low, phys_d_low,
                                   neg_p_all, neg_d_all,
                                   neg_p_low, neg_d_low))

#renaming the columns
colnames(syn_data_df) <- c("pain_p_all", "pain_d_all",
                           "pain_p_high", "pain_d_high",
                           "phys_p_all", "phys_d_all",
                           "phys_p_low", "phys_d_low",
                           "neg_p_all", "neg_d_all",
                           "neg_p_low", "neg_d_low")

#saving the file
write.csv(syn_data_df, "multi_syn_ttest.csv")


##--FINAL WRANGLING--##################################################################################

#first instance
final_df <- read.csv("multi_syn_ttest.csv")
filtered_final_df <- final_df %>% filter(pain_p_all > 0.05 &
                                           pain_p_high < 0.05 & pain_d_high >= 0.500 & pain_d_high <= 0.800 &
                                           phys_p_all < 0.05 & phys_d_all >= 1.800 & phys_d_all <= 2.000 &
                                           phys_p_low < 0.05 & phys_d_low >= 1.200 & phys_d_low <= 1.500 &
                                           neg_p_all < 0.05 & neg_d_all >= 2.200 & neg_d_all <= 2.500 &
                                           neg_p_low < 0.05 & neg_d_low >= 2.200 & neg_d_low <= 2.500)
filtered_final_df

write.csv(filtered_final_df, "multi_syn_ttest_FILTERED.csv")
