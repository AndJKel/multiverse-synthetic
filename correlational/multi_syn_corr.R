#Using Synthetic Data in a Basic Multiverse Analysis - Correlational
#Andrew J. Kelly

##--LIBRARIES--##################################################################################
library("tidyverse") #package for data wrangling
library("synthpop") #package for synthesizing the data
library("effectsize") #package for effect sizes


##--READING THE FILE--##################################################################################
#reading the fle
knaster <- read.csv("knaster.csv") #reads file in
knaster #reports the file

#selects relevent columns from the knaster file
knaster <- knaster %>% select(Current.pain, SUMBDI, Physicalfunction, Negativeviewofself) #puts selected columns into data frame
knaster #reports data frame


##--BDI - FILTERING--##################################################################################

#filtering out severe BDI
BDI_severe <- knaster %>% filter(SUMBDI > 28)
BDI_severe

#filtering out moderate BDI
BDI_moderate <- knaster %>% filter(SUMBDI >= 20 & SUMBDI <= 28)
BDI_moderate

#filtering out clinical BDI
BDI_clinical <- knaster %>% filter(SUMBDI >= 20)
BDI_clinical

#filtering out mild depression BDI
BDI_mild <- knaster %>% filter(SUMBDI >= 14 & SUMBDI <= 19)
BDI_mild

#filtering out minimal depression BDI
BDI_minimal <- knaster %>% filter(SUMBDI < 14)
BDI_minimal

#filtering out subclinical depression BDI
BDI_subclinical <- knaster %>% filter(SUMBDI <= 19)
BDI_subclinical

##--BDI - CORRELATION TESTS--##################################################################################

severe_BDIxPain <- cor.test(BDI_severe$SUMBDI, BDI_severe$Current.pain)
moderate_BDIxPain <- cor.test(BDI_moderate$SUMBDI, BDI_moderate$Current.pain)
clinical_BDIxPain <- cor.test(BDI_clinical$SUMBDI, BDI_clinical$Current.pain)
mild_BDIxPain <- cor.test(BDI_mild$SUMBDI, BDI_mild$Current.pain)
minimal_BDIxPain <- cor.test(BDI_minimal$SUMBDI, BDI_minimal$Current.pain)
subclinical_BDIxPain <- cor.test(BDI_subclinical$SUMBDI, BDI_subclinical$Current.pain)
everyone_BDIxPain <- cor.test(knaster$SUMBDI, knaster$Current.pain)


##--BDI - CREATING FRAMES--##################################################################################

#correlation frames
severe_BDI_cx <- data.frame(severe_BDI_cx=c)
moderate_BDI_cx <- data.frame(moderate_BDI_cx=c)
clinical_BDI_cx <- data.frame(clinical_BDI_cx=c)
mild_BDI_cx <- data.frame(mild_BDI_cx=c)
minimal_BDI_cx <- data.frame(minimal_BDI_cx=c)
subclinical_BDI_cx <- data.frame(subclinical_BDI_cx=c)
everyone_BDI_cx <- data.frame(everyone_BDI_cx=c)

#p-value frames
severe_BDI_p <- data.frame(severe_BDI_p=c)
moderate_BDI_p <- data.frame(moderate_BDI_p=c)
clinical_BDI_p <- data.frame(clinical_BDI_p=c)
mild_BDI_p <- data.frame(mild_BDI_p=c)
minimal_BDI_p <- data.frame(minimal_BDI_p=c)
subclinical_BDI_p <- data.frame(subclinical_BDI_p=c)
everyone_BDI_p <- data.frame(everyone_BDI_p=c)

##--NEGATIVE - FILTERING--##################################################################################

#filtering out severe BDI
na1 <- knaster %>% filter(Negativeviewofself <= 4)
na1

na2 <- knaster %>% filter(Negativeviewofself >= 5 & Negativeviewofself <= 9)
na2

#na3 <- knaster %>% filter(Negativeviewofself >= 10 & Negativeviewofself <= 14)
#na3

#na4 <- knaster %>% filter(Negativeviewofself >= 15 & Negativeviewofself <= 18)
#na4

nb1 <- knaster %>% filter(Negativeviewofself <= 9)
nb1

nb2 <- knaster %>% filter(Negativeviewofself >= 5 & Negativeviewofself <= 14)
nb2

#nb3 <- knaster %>% filter(Negativeviewofself >= 10 & Negativeviewofself <= 18)
#nb3

nc1 <- knaster %>% filter(Negativeviewofself <= 14)
nc1

nc2 <- knaster %>% filter(Negativeviewofself >= 5 & Negativeviewofself <= 18)
nc2

##--NEGATIVE - CORRELATION TESTS--##################################################################################

na1xPain <- cor.test(na1$Negativeviewofself, na1$Current.pain)
na2xPain <- cor.test(na2$Negativeviewofself, na2$Current.pain)
#na3xPain <- cor.test(na3$Negativeviewofself, na3$Current.pain)
#na4xPain <- cor.test(na4$Negativeviewofself, na4$Current.pain)
nb1xPain <- cor.test(nb1$Negativeviewofself, nb1$Current.pain)
nb2xPain <- cor.test(nb2$Negativeviewofself, nb2$Current.pain)
#nb3xPain <- cor.test(nb3$Negativeviewofself, nb3$Current.pain)
nc1xPain <- cor.test(nc1$Negativeviewofself, nc1$Current.pain)
nc2xPain <- cor.test(nc2$Negativeviewofself, nc2$Current.pain)
everyone_NegxPain <- cor.test(knaster$Negativeviewofself, knaster$Current.pain)

##--NEGATIVE - CREATING FRAMES--##################################################################################

#correlation frames
na1_cx <- data.frame(na1_cx=c)
na2_cx <- data.frame(na2_cx=c)
#na3_cx <- data.frame(na3_cx=c)
#na4_cx <- data.frame(na4_cx=c)
nb1_cx <- data.frame(nb1_cx=c)
nb2_cx <- data.frame(nb2_cx=c)
#nb3_cx <- data.frame(nb3_cx=c)
nc1_cx <- data.frame(nc1_cx=c)
nc2_cx <- data.frame(nc2_cx=c)
everyone_Neg_cx <- data.frame(everyone_Neg_cx=c)

#p-value frames
na1_p <- data.frame(na1_p=c)
na2_p <- data.frame(na2_p=c)
#na3_p <- data.frame(na3_p=c)
#na4_p <- data.frame(na4_p=c)
nb1_p <- data.frame(nb1_p=c)
nb2_p <- data.frame(nb2_p=c)
#nb3_p <- data.frame(nb3_p=c)
nc1_p <- data.frame(nc1_p=c)
nc2_p <- data.frame(nc2_p=c)
everyone_Neg_p <- data.frame(everyone_Neg_p=c)

##--SOMATIC - FILTERING--##################################################################################

#filtering out severe BDI
sa1 <- knaster %>% filter(Physicalfunction <= 5)
sa1

sa2 <- knaster %>% filter(Physicalfunction >= 6 & Physicalfunction <= 9)
sa2

sa3 <- knaster %>% filter(Physicalfunction >= 10 & Physicalfunction <= 12)
sa3

#sa4 <- knaster %>% filter(Physicalfunction >= 13 & Physicalfunction <= 16)
#sa4

#sa5 <- knaster %>% filter(Physicalfunction >= 17 & Physicalfunction <= 21)
#sa5

sb1 <- knaster %>% filter(Physicalfunction <= 9)
sb1

sb2 <- knaster %>% filter(Physicalfunction >= 6 & Physicalfunction <= 12)
sb2

sb3 <- knaster %>% filter(Physicalfunction >= 10 & Physicalfunction <= 16)
sb3

#sb4 <- knaster %>% filter(Physicalfunction >= 13 & Physicalfunction <= 21)
#sb4

sc1 <- knaster %>% filter(Physicalfunction <= 12)
sc1

sc2 <- knaster %>% filter(Physicalfunction >= 6 & Physicalfunction <= 16)
sc2

sc3 <- knaster %>% filter(Physicalfunction >= 10 & Physicalfunction <= 21)
sc3

sd1 <- knaster %>% filter(Physicalfunction <= 16)
sd1

sd2 <- knaster %>% filter(Physicalfunction >= 6 & Physicalfunction <= 21)
sd2

##--SOMATIC - CORRELATION TESTS--##################################################################################

sa1xPain <- cor.test(sa1$Physicalfunction, sa1$Current.pain)
sa2xPain <- cor.test(sa2$Physicalfunction, sa2$Current.pain)
sa3xPain <- cor.test(sa3$Physicalfunction, sa3$Current.pain)
#sa4xPain <- cor.test(sa4$Physicalfunction, sa4$Current.pain)
#sa5xPain <- cor.test(sa5$Physicalfunction, sa5$Current.pain)
sb1xPain <- cor.test(sb1$Physicalfunction, sb1$Current.pain)
sb2xPain <- cor.test(sb2$Physicalfunction, sb2$Current.pain)
sb3xPain <- cor.test(sb3$Physicalfunction, sb3$Current.pain)
#sb4xPain <- cor.test(sb4$Physicalfunction, sb4$Current.pain)
sc1xPain <- cor.test(sc1$Physicalfunction, sc1$Current.pain)
sc2xPain <- cor.test(sc2$Physicalfunction, sc2$Current.pain)
sc3xPain <- cor.test(sc3$Physicalfunction, sc3$Current.pain)
sd1xPain <- cor.test(sd1$Physicalfunction, sd1$Current.pain)
sd2xPain <- cor.test(sd2$Physicalfunction, sd2$Current.pain)
everyone_SomxPain <- cor.test(knaster$Physicalfunction, knaster$Current.pain)

##--SOMATIC - CREATING FRAMES--##################################################################################

#correlation frames
sa1_cx <- data.frame(sa1_cx=c)
sa2_cx <- data.frame(sa2_cx=c)
sa3_cx <- data.frame(sa3_cx=c)
#sa4_cx <- data.frame(sa4_cx=c)
#sa5_cx <- data.frame(sa5_cx=c)
sb1_cx <- data.frame(sb1_cx=c)
sb2_cx <- data.frame(sb2_cx=c)
sb3_cx <- data.frame(sb3_cx=c)
#sb4_cx <- data.frame(sb4_cx=c)
sc1_cx <- data.frame(sc1_cx=c)
sc2_cx <- data.frame(sc2_cx=c)
sc3_cx <- data.frame(sc3_cx=c)
sd1_cx <- data.frame(sd1_cx=c)
sd2_cx <- data.frame(sd2_cx=c)
everyone_Som_cx <- data.frame(everyone_Som_cx=c)

#p-value frames
sa1_p <- data.frame(sa1_p=c)
sa2_p <- data.frame(sa2_p=c)
sa3_p <- data.frame(sa3_p=c)
#sa4_p <- data.frame(sa4_p=c)
#sa5_p <- data.frame(sa5_p=c)
sb1_p <- data.frame(sb1_p=c)
sb2_p <- data.frame(sb2_p=c)
sb3_p <- data.frame(sb3_p=c)
#sb4_p <- data.frame(sb4_p=c)
sc1_p <- data.frame(sc1_p=c)
sc2_p <- data.frame(sc2_p=c)
sc3_p <- data.frame(sc3_p=c)
sd1_p <- data.frame(sd1_p=c)
sd2_p <- data.frame(sd2_p=c)
everyone_Som_p <- data.frame(everyone_Som_p=c)

#setting the seeds outside of the loop
synth_loop <- 1
synth_seed <- 7500


##--LOOP--##################################################################################

#while loop - will loop 2500 times
while(synth_loop < 2501){
  
  synth_seed <- 3000 + synth_loop
  syn_data <- (syn(knaster, seed = synth_seed))$syn #creating a synthetic dataset
  syn_data #reports synthetic dataset

  ##--LOOP - BDI - FILTERING--##################################################################################
  syn_BDI_severe <- syn_data %>% filter(SUMBDI > 28)
  syn_BDI_moderate <- syn_data %>% filter(SUMBDI >= 20 & SUMBDI <= 28)
  syn_BDI_clinical <- syn_data %>% filter(SUMBDI >= 20)
  syn_BDI_mild <- syn_data %>% filter(SUMBDI >= 14 & SUMBDI <= 19)
  syn_BDI_minimal <- syn_data %>% filter(SUMBDI < 14)
  syn_BDI_subclinical <- syn_data %>% filter(SUMBDI <= 19)

  ##--LOOP - BDI - CORS--##################################################################################
  syn_severe_BDI_cx <- cor.test(syn_BDI_severe$SUMBDI, syn_BDI_severe$Current.pain)
  syn_moderate_BDI_cx <- cor.test(syn_BDI_moderate$SUMBDI, syn_BDI_moderate$Current.pain)
  syn_clinical_BDI_cx <- cor.test(syn_BDI_clinical$SUMBDI, syn_BDI_clinical$Current.pain)
  syn_mild_BDI_cx <- cor.test(syn_BDI_mild$SUMBDI, syn_BDI_mild$Current.pain)
  syn_minimal_BDI_cx <- cor.test(syn_BDI_minimal$SUMBDI, syn_BDI_minimal$Current.pain)
  syn_subclinical_BDI_cx <- cor.test(syn_BDI_subclinical$SUMBDI, syn_BDI_subclinical$Current.pain)
  syn_everyone_BDI_cx <- cor.test(syn_data$SUMBDI, syn_data$Current.pain)
  
  ##--LOOP - BDI - BINDING--##################################################################################
  severe_BDI_cx <- rbind(severe_BDI_cx, syn_severe_BDI_cx$estimate)
  severe_BDI_p <- rbind(severe_BDI_p, syn_severe_BDI_cx$p.value)
  moderate_BDI_cx <- rbind(moderate_BDI_cx, syn_moderate_BDI_cx$estimate)
  moderate_BDI_p <- rbind(moderate_BDI_p, syn_moderate_BDI_cx$p.value)
  clinical_BDI_cx <- rbind(clinical_BDI_cx, syn_clinical_BDI_cx$estimate)
  clinical_BDI_p <- rbind(clinical_BDI_p, syn_clinical_BDI_cx$p.value)
  mild_BDI_cx <- rbind(mild_BDI_cx, syn_mild_BDI_cx$estimate)
  mild_BDI_p <- rbind(mild_BDI_p, syn_mild_BDI_cx$p.value)
  minimal_BDI_cx <- rbind(minimal_BDI_cx, syn_minimal_BDI_cx$estimate)
  minimal_BDI_p <- rbind(minimal_BDI_p, syn_minimal_BDI_cx$p.value)
  subclinical_BDI_cx <- rbind(subclinical_BDI_cx, syn_subclinical_BDI_cx$estimate)
  subclinical_BDI_p <- rbind(subclinical_BDI_p, syn_subclinical_BDI_cx$p.value)
  everyone_BDI_cx <- rbind(everyone_BDI_cx, syn_everyone_BDI_cx$estimate)
  everyone_BDI_p <- rbind(everyone_BDI_p, syn_everyone_BDI_cx$p.value)
  
  ##--LOOP - NEGATIVE - FILTERING--##################################################################################
  syn_na1 <- syn_data %>% filter(Negativeviewofself <= 4)
  syn_na2 <- syn_data %>% filter(Negativeviewofself >= 5 & Negativeviewofself <= 9)
  #syn_na3 <- syn_data %>% filter(Negativeviewofself >= 10 & Negativeviewofself <= 14)
  #syn_na4 <- syn_data %>% filter(Negativeviewofself >= 15 & Negativeviewofself <= 18)
  syn_nb1 <- syn_data %>% filter(Negativeviewofself <= 9)
  syn_nb2 <- syn_data %>% filter(Negativeviewofself >= 5 & Negativeviewofself <= 14)
  #syn_nb3 <- syn_data %>% filter(Negativeviewofself >= 10 & Negativeviewofself <= 18)
  syn_nc1 <- syn_data %>% filter(Negativeviewofself <= 14)
  syn_nc2 <- syn_data %>% filter(Negativeviewofself >= 5 & Negativeviewofself <= 18)

  ##--LOOP - NEGATIVE - CORS--##################################################################################
  syn_na1xPain <- cor.test(syn_na1$Negativeviewofself, syn_na1$Current.pain)
  syn_na2xPain <- cor.test(syn_na2$Negativeviewofself, syn_na2$Current.pain)
  #syn_na3xPain <- cor.test(syn_na3$Negativeviewofself, syn_na3$Current.pain)
  #syn_na4xPain <- cor.test(syn_na4$Negativeviewofself, syn_na4$Current.pain)
  syn_nb1xPain <- cor.test(syn_nb1$Negativeviewofself, syn_nb1$Current.pain)
  syn_nb2xPain <- cor.test(syn_nb2$Negativeviewofself, syn_nb2$Current.pain)
  #syn_nb3xPain <- cor.test(syn_nb3$Negativeviewofself, syn_nb3$Current.pain)
  syn_nc1xPain <- cor.test(syn_nc1$Negativeviewofself, syn_nc1$Current.pain)
  syn_nc2xPain <- cor.test(syn_nc2$Negativeviewofself, syn_nc2$Current.pain)
  syn_everyone_NegxPain <- cor.test(syn_data$Negativeviewofself, syn_data$Current.pain)
  
  ##--LOOP - NEGATIVE - BINDING--##################################################################################
  na1_cx <- rbind(na1_cx, syn_na1xPain$estimate)
  na1_p <- rbind(na1_p, syn_na1xPain$p.value)
  na2_cx <- rbind(na2_cx, syn_na2xPain$estimate)
  na2_p <- rbind(na2_p, syn_na2xPain$p.value)
  #na3_cx <- rbind(na3_cx, syn_na3xPain$estimate)
  #na3_p <- rbind(na3_p, syn_na3xPain$p.value)
  #na4_cx <- rbind(na4_cx, syn_na4xPain$estimate)
  #na4_p <- rbind(na4_p, syn_na4xPain$p.value)
  nb1_cx <- rbind(nb1_cx, syn_nb1xPain$estimate)
  nb1_p <- rbind(nb1_p, syn_nb1xPain$p.value)
  nb2_cx <- rbind(nb2_cx, syn_nb2xPain$estimate)
  nb2_p <- rbind(nb2_p, syn_nb2xPain$p.value)
  #nb3_cx <- rbind(nb3_cx, syn_nb3xPain$estimate)
  #nb3_p <- rbind(nb3_p, syn_nb3xPain$p.value)
  nc1_cx <- rbind(nc1_cx, syn_nc1xPain$estimate)
  nc1_p <- rbind(nc1_p, syn_nc1xPain$p.value)
  nc2_cx <- rbind(nc2_cx, syn_nc2xPain$estimate)
  nc2_p <- rbind(nc2_p, syn_nc2xPain$p.value)
  everyone_Neg_cx <- rbind(everyone_Neg_cx, everyone_SomxPain$estimate)
  everyone_Neg_p <- rbind(everyone_Neg_p, everyone_SomxPain$p.value)
  
  ##--LOOP - SOMATIC - FILTERING--##################################################################################
  syn_sa1 <- syn_data %>% filter(Physicalfunction <= 5)
  syn_sa2 <- syn_data %>% filter(Physicalfunction >= 6 & Physicalfunction <= 9)
  syn_sa3 <- syn_data %>% filter(Physicalfunction >= 10 & Physicalfunction <= 12)
  #syn_sa4 <- syn_data %>% filter(Physicalfunction >= 13 & Physicalfunction <= 16)
  #syn_sa5 <- syn_data %>% filter(Physicalfunction >= 17 & Physicalfunction <= 21)
  syn_sb1 <- syn_data %>% filter(Physicalfunction <= 9)
  syn_sb2 <- syn_data %>% filter(Physicalfunction >= 6 & Physicalfunction <= 12)
  syn_sb3 <- syn_data %>% filter(Physicalfunction >= 10 & Physicalfunction <= 16)
  #syn_sb4 <- syn_data %>% filter(Physicalfunction >= 13 & Physicalfunction <= 21)
  syn_sc1 <- syn_data %>% filter(Physicalfunction <= 12)
  syn_sc2 <- syn_data %>% filter(Physicalfunction >= 6 & Physicalfunction <= 16)
  syn_sc3 <- syn_data %>% filter(Physicalfunction >= 10 & Physicalfunction <= 21)
  syn_sd1 <- syn_data %>% filter(Physicalfunction <= 16)
  syn_sd2 <- syn_data %>% filter(Physicalfunction >= 6 & Physicalfunction <= 21)
  
  ##--LOOP - SOMATIC - CORS--##################################################################################
  syn_sa1xPain <- cor.test(syn_sa1$Physicalfunction, syn_sa1$Current.pain)
  syn_sa2xPain <- cor.test(syn_sa2$Physicalfunction, syn_sa2$Current.pain)
  syn_sa3xPain <- cor.test(syn_sa3$Physicalfunction, syn_sa3$Current.pain)
  #syn_sa4xPain <- cor.test(syn_sa4$Physicalfunction, syn_sa4$Current.pain)
  #syn_sa5xPain <- cor.test(syn_sa5$Physicalfunction, syn_sa5$Current.pain)
  syn_sb1xPain <- cor.test(syn_sb1$Physicalfunction, syn_sb1$Current.pain)
  syn_sb2xPain <- cor.test(syn_sb2$Physicalfunction, syn_sb2$Current.pain)
  syn_sb3xPain <- cor.test(syn_sb3$Physicalfunction, syn_sb3$Current.pain)
  #syn_sb4xPain <- cor.test(syn_sb4$Physicalfunction, syn_sb4$Current.pain)
  syn_sc1xPain <- cor.test(syn_sc1$Physicalfunction, syn_sc1$Current.pain)
  syn_sc2xPain <- cor.test(syn_sc2$Physicalfunction, syn_sc2$Current.pain)
  syn_sc3xPain <- cor.test(syn_sc3$Physicalfunction, syn_sc3$Current.pain)
  syn_sd1xPain <- cor.test(syn_sd1$Physicalfunction, syn_sd1$Current.pain)
  syn_sd2xPain <- cor.test(syn_sd2$Physicalfunction, syn_sd2$Current.pain)
  syn_everyone_SomxPain <- cor.test(syn_data$Physicalfunction, syn_data$Current.pain)
  
  ##--LOOP - SOMATIC - BINDING--##################################################################################
  sa1_cx <- rbind(sa1_cx, syn_sa1xPain$estimate)
  sa1_p <- rbind(sa1_p, syn_sa1xPain$p.value)
  sa2_cx <- rbind(sa2_cx, syn_sa2xPain$estimate)
  sa2_p <- rbind(sa2_p, syn_sa2xPain$p.value)
  sa3_cx <- rbind(sa3_cx, syn_sa3xPain$estimate)
  sa3_p <- rbind(sa3_p, syn_sa3xPain$p.value)
  #sa4_cx <- rbind(sa4_cx, syn_sa4xPain$estimate)
  #sa4_p <- rbind(sa4_p, syn_sa4xPain$p.value)
  #sa5_cx <- rbind(sa5_cx, syn_sa5xPain$estimate)
  #sa5_p <- rbind(sa5_p, syn_sa5xPain$p.value)
  sb1_cx <- rbind(sb1_cx, syn_sb1xPain$estimate)
  sb1_p <- rbind(sb1_p, syn_sb1xPain$p.value)
  sb2_cx <- rbind(sb2_cx, syn_sb2xPain$estimate)
  sb2_p <- rbind(sb2_p, syn_sb2xPain$p.value)
  sb3_cx <- rbind(sb3_cx, syn_sb3xPain$estimate)
  sb3_p <- rbind(sb3_p, syn_sb3xPain$p.value)
  #sb4_cx <- rbind(sb4_cx, syn_sb4xPain$estimate)
  #sb4_p <- rbind(sb4_p, syn_sb4xPain$p.value)
  sc1_cx <- rbind(sc1_cx, syn_sc1xPain$estimate)
  sc1_p <- rbind(sc1_p, syn_sc1xPain$p.value)
  sc2_cx <- rbind(sc2_cx, syn_sc2xPain$estimate)
  sc2_p <- rbind(sc2_p, syn_sc2xPain$p.value)
  sc3_cx <- rbind(sc3_cx, syn_sc3xPain$estimate)
  sc3_p <- rbind(sc3_p, syn_sc3xPain$p.value)
  sd1_cx <- rbind(sd1_cx, syn_sd1xPain$estimate)
  sd1_p <- rbind(sd1_p, syn_sd1xPain$p.value)
  sd2_cx <- rbind(sd2_cx, syn_sd2xPain$estimate)
  sd2_p <- rbind(sd2_p, syn_sd2xPain$p.value)
  everyone_Som_cx <- rbind(everyone_Som_cx, syn_everyone_SomxPain$estimate)
  everyone_Som_p <- rbind(everyone_Som_p, syn_everyone_SomxPain$p.value)
  
  #adds numbers to the loop
  synth_loop <- synth_loop + 1
  
}

##--WRANGLING AFTER LOOP--##################################################################################

#binds all correlations and p-values into a single dataframe
syn_data_df <- do.call(cbind, list(severe_BDI_cx, severe_BDI_p,
                                 moderate_BDI_cx, moderate_BDI_p,
                                 clinical_BDI_cx, clinical_BDI_p,
                                 mild_BDI_cx, mild_BDI_p,
                                 minimal_BDI_cx, minimal_BDI_p,
                                 subclinical_BDI_cx, subclinical_BDI_p,
                                 everyone_BDI_cx, everyone_BDI_p,
                                 na1_cx, na1_p,
                                 na2_cx, na2_p,
                                 #na3_cx, na3_p,
                                 #na4_cx, na4_p,
                                 nb1_cx, nb1_p,
                                 nb2_cx, nb2_p,
                                 #nb3_cx, nb3_p,
                                 nc1_cx, nc1_p,
                                 nc2_cx, nc2_p,
                                 everyone_Neg_cx, everyone_Neg_p,
                                 sa1_cx, sa1_p,
                                 sa2_cx, sa2_p,
                                 sa3_cx, sa3_p,
                                 #sa4_cx, sa4_p,
                                 #sa5_cx, sa5_p,
                                 sb1_cx, sb1_p,
                                 sb2_cx, sb2_p, 
                                 sb3_cx, sb3_p,
                                 #sb4_cx, sb4_p,
                                 sc1_cx, sc1_p,
                                 sc2_cx, sc2_p,
                                 sc3_cx, sc3_p,
                                 sd1_cx, sd1_p,
                                 sd2_cx, sd2_p,
                                 everyone_Som_cx, everyone_Som_p))


#REMEMBER TO RENAME THE FUCKING COLUMNS
colnames(syn_data_df) <- c("severe_BDI_cx", "severe_BDI_p",
                           "moderate_BDI_cx", "moderate_BDI_p",
                           "clinical_BDI_cx", "clinical_BDI_p",
                           "mild_BDI_cx", "mild_BDI_p",
                           "minimal_BDI_cx", "minimal_BDI_p",
                           "subclinical_BDI_cx", "subclinical_BDI_p",
                           "everyone_BDI_cx", "everyone_BDI_p",
                           "na1_cx", "na1_p",
                           "na2_cx", "na2_p",
                           #"na3_cx, na3_p,
                           #"na4_cx, na4_p,
                           "nb1_cx", "nb1_p",
                           "nb2_cx", "nb2_p",
                           #"nb3_cx, nb3_p,
                           "nc1_cx", "nc1_p",
                           "nc2_cx", "nc2_p",
                           "everyone_Neg_cx", "everyone_Neg_p",
                           "sa1_cx", "sa1_p",
                           "sa2_cx", "sa2_p",
                           "sa3_cx", "sa3_p",
                           #"sa4_cx, sa4_p,
                           #"sa5_cx, sa5_p,
                           "sb1_cx", "sb1_p",
                           "sb2_cx", "sb2_p", 
                           "sb3_cx", "sb3_p",
                           #"sb4_cx, sb4_p,
                           "sc1_cx", "sc1_p",
                           "sc2_cx", "sc2_p",
                           "sc3_cx", "sc3_p",
                           "sd1_cx", "sd1_p",
                           "sd2_cx", "sd2_p",
                           "everyone_Som_cx", "everyone_Som_p")

#saving the file
write.csv(syn_data_df, "multi_syn_corr.csv")


##--FINAL WRANGLING--##################################################################################


#first instance
final_df <- read.csv("multi_syn_corr.csv")
filtered_final_df <- final_df %>% filter(severe_BDI_p >= 0.050 & 
                                           moderate_BDI_cx >= 0.60 & moderate_BDI_cx <= 0.799 & moderate_BDI_p <= 0.050 & 
                                           clinical_BDI_p >= 0.050 & 
                                           mild_BDI_p >= 0.05 &
                                           minimal_BDI_p >= 0.050 & 
                                           subclinical_BDI_p >= 0.050 &
                                           everyone_BDI_p >= 0.050 &
                                           na1_p >= 0.050 &
                                           na2_p >= 0.050 &
                                           #na3_p >= 0.050 &
                                           #na4_p >= 0.00 &
                                           nb1_p >= 0.050 &
                                           nb2_p >= 0.050 &
                                           #nb3_p >= 0.050 &
                                           nc1_p >= 0.050 &
                                           nc2_p >= 0.050 &
                                           everyone_Neg_p >= 0.050 &
                                           sa1_p >= 0.050 &
                                           sa2_cx >= 0.400 & sa2_cx <= 0.599 & sa2_p <= 0.050 &
                                           sa3_p >= 0.050 &
                                           #sa4_p >= 0.050 &
                                           #sa5_p >= 0.000 &
                                           sb1_p >= 0.050 &
                                           sb2_p >= 0.050 &
                                           sb3_p >= 0.050 &
                                           #sb4_p >= 0.050 &
                                           sc1_p >= 0.050 &
                                           sc2_cx >= 0.200 & sc2_cx >= 0.399 & sc2_p <= 0.050 &
                                           sc3_p >= 0.050 &
                                           sd1_cx >= 0.200 & sd1_cx >= 0.399 & sd1_p <= 0.050 &
                                           sd2_cx >= 0.200 & sd2_cx >= 0.399 & sd2_p <= 0.050 &
                                           everyone_Som_cx >= 0.20 & everyone_Som_cx >= 0.399 & everyone_Som_p <= 0.050)

filtered_final_df #report

write.csv(filtered_final_df, "multi_syn_corr_FILTERED.csv") #writes the final set of filtered, eligible datasets (only the correlations and p-values) to a file