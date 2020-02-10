# calculate correlations between model predictions and observations ---------------------------
rm(list=ls()) #remove previous variable assignments

# load case data
load("data/case_data.RData")

# load vector data
load("data/vector_data.RData")

# combine observed case and vector data
obs.data <- merge(cases, vectors, by=c("Site", "Date"), all=T)

# combine modeled data with observed data
rmodels <- read.csv("SEI-SEIR_simulations_THR.csv", head=T, stringsAsFactors=F)
rmodels$Mtot <- rmodels$M1 + rmodels$M2 + rmodels$M3
rmodels$Date <- as.Date(rmodels$Date, "%Y-%m-%d")
merged_data <- merge(rmodels, obs.data, by=c("Site", "Date"))

# correlation site & by rainfall function ----------------------------------------------------
sites <- unique(rmodels$Site)
r_functions <- c("Briere", "Quadratic", "Inverse")
df_rain <- data.frame(matrix(ncol=4, nrow=0))
colnames(df_rain) <- c("Site", "Rain_function", "Aedes_corr", "Arbovirus_corr")

for (j in 1:length(sites)){
  for (k in 1:length(r_functions)){
    aedes <- cor(merged_data$aedes_total[merged_data$Site==sites[j]&merged_data$Rain_function==r_functions[k]], merged_data$Mtot[merged_data$Site==sites[j]&merged_data$Rain_function==r_functions[k]], use = "complete.obs")
    arboviruses <- cor(merged_data$arboviruses_total[merged_data$Site==sites[j]&merged_data$Rain_function==r_functions[k]], merged_data$I[merged_data$Site==sites[j]&merged_data$Rain_function==r_functions[k]], use = "complete.obs")
    df_rain[(nrow(df_rain)+1),] <- c(sites[j], r_functions[k], round(aedes,2), round(arboviruses,2))
  }
}

# identify best fit models based on correlation --------------------------------------------
library(tidyverse)

corr_aedes <- df_rain %>% 
  group_by(Site) %>%
  filter(Aedes_corr == max(as.numeric(Aedes_corr))) %>%
  select(-Arbovirus_corr)

write.csv(corr_aedes, "correlation_best_model_aedes.csv", row.names=F)

corr_arboviruses <- df_rain %>% 
  group_by(Site) %>%
  filter(Arbovirus_corr == max(as.numeric(Arbovirus_corr))) %>%
  select(-Aedes_corr)
corr_arboviruses <- corr_arboviruses[1:8,]

write.csv(corr_arboviruses, "correlation_best_model_arboviruses.csv", row.names=F)

# mean of means
mean(as.numeric(corr_aedes$Aedes_corr))
mean(as.numeric(corr_arboviruses$Arbovirus_corr))

# correlation for other mosquito life stages -------------------------------------------------
corr_vectors <- data.frame(matrix(ncol=6, nrow=0))
colnames(corr_vectors) <- c("Site", "Rain_function", "Pupae", "Late_instars", "Early_instars", "Eggs")
sites <- c("Chulaimbo", "Kisumu", "Msambweni", "Ukunda")

for (j in 1:length(sites)){
  for (k in 1:length(r_functions)){
    pupae <- cor(merged_data$pupae_total[merged_data$Site==sites[j]&merged_data$Rain_function==r_functions[k]], merged_data$Mtot[merged_data$Site==sites[j]&merged_data$Rain_function==r_functions[k]], use = "complete.obs")
    lateInstars <- cor(merged_data$late_instar_total[merged_data$Site==sites[j]&merged_data$Rain_function==r_functions[k]], merged_data$Mtot[merged_data$Site==sites[j]&merged_data$Rain_function==r_functions[k]], use = "complete.obs")
    earlyInstars <- cor(merged_data$early_instar_total[merged_data$Site==sites[j]&merged_data$Rain_function==r_functions[k]], merged_data$Mtot[merged_data$Site==sites[j]&merged_data$Rain_function==r_functions[k]], use = "complete.obs")
    eggs <- cor(merged_data$egg_total[merged_data$Site==sites[j]&merged_data$Rain_function==r_functions[k]], merged_data$Mtot[merged_data$Site==sites[j]&merged_data$Rain_function==r_functions[k]], use = "complete.obs")
    corr_vectors[(nrow(corr_vectors)+1),] <- c(sites[j], r_functions[k], round(pupae,2), round(lateInstars,2), round(earlyInstars,2), round(eggs,2))
  }
}

# subset modeled data based on best rainfall function ---------------------------------------
cases_and_mods <- corr_arboviruses[,c("Site", "Rain_function")] %>% 
  left_join(merged_data[,c("Site", "Date", "time", "I", "arboviruses_total", "Rain_function")], by=c("Site", "Rain_function"))

save(cases_and_mods, file="model_v_data_cases.RData")

# merge vector data and modeled cases
vectors_and_mods <- corr_aedes[,c("Site", "Rain_function")] %>% 
  left_join(merged_data[,c("Site", "Date", "time", "Mtot", "aedes_total", "Rain_function")], by=c("Site", "Rain_function"))

save(vectors_and_mods, file="model_v_data_vectors.RData")

# sign tests --------------------------------------------------------------------------------
# subset vector data
vectors <- vectors_and_mods[,c("Site", "Date", "Mtot", "aedes_total")]
vectors <- vectors[complete.cases(vectors),]

# calculate whether mosquitoes were predicted and observed to increase or decrease
vectors <- mutate(vectors, Mtot_diff = Mtot - lag(Mtot), aedes_diff = aedes_total - lag(aedes_total))
vectors$sign <- ifelse((vectors$Mtot_diff < 0 & vectors$aedes_diff < 0)|(vectors$Mtot_diff > 0 & vectors$aedes_diff > 0), 1, 0)

# determine first date for each site and make sign NA for that time point
earliestDates <- vectors %>% group_by(Site) %>% summarize(minDate = min(Date))
for (i in 1:nrow(earliestDates)){
  rowX <- which(vectors$Site == earliestDates$Site[i] & vectors$Date == earliestDates$minDate[i])
  vectors$sign[rowX] <- NA
}

# sign test
binom.test(sum(vectors$sign==1, na.rm=T), sum(!is.na(vectors$sign)), p = 0.5, alternative = c("two.sided"), conf.level = 0.95)

# subset case data
cases <- cases_and_mods[,c("Site", "Date", "I", "arboviruses_total")]
cases <- cases[complete.cases(cases),]

# calculate whether cases were predicted and observed to increase or decrease
cases <- mutate(cases, I_diff = I - lag(I), cases_diff = arboviruses_total - lag(arboviruses_total))
cases$sign <- ifelse((cases$I_diff < 0 & cases$cases_diff < 0)|(cases$I_diff > 0 & cases$cases_diff > 0), 1, 0)

# determine first date for each site and make sign NA for that time point
earliestDates <- cases %>% group_by(Site) %>% summarize(minDate = min(Date))
for (i in 1:nrow(earliestDates)){
  rowX <- which(cases$Site == earliestDates$Site[i] & cases$Date == earliestDates$minDate[i])
  cases$sign[rowX] <- NA
}

# sign test
binom.test(sum(cases$sign==1, na.rm=T), sum(!is.na(cases$sign)), p = 0.5, alternative = c("two.sided"), conf.level = 0.95)

# ANOVAs ----------------------------------------------------------------------------------
# vectors
vectors$Year <- format(vectors$Date, "%Y")

vectors2 <- vectors %>% 
  group_by(Site, Year) %>% 
  summarize(percAnnualPred = sum(aedes_total, na.rm=T)/sum(Mtot)) #%>%

# check assumption of homogeneity of variance
library(car)
leveneTest(vectors2$percAnnualPred~vectors2$Site)

# anova
anovaVectors <- aov(vectors2$percAnnualPred~vectors2$Site)
summary(anovaVectors)

# check residuals
resVectors <- anovaVectors$residuals
hist(resVectors)

# Tukey's posthoc test to find out which sites differ
TukeyHSD(anovaVectors)

# cases
cases$Year <- format(cases$Date, "%Y")

cases2 <- cases %>% 
  group_by(Site, Year) %>% 
  summarize(percAnnualPred = sum(arboviruses_total, na.rm=T)/sum(I)) %>%
  filter(percAnnualPred < 1)

# check assumption of homogeneity of variance
leveneTest(cases2$percAnnualPred~cases2$Site)

# anova
anovaCases <- aov(cases2$percAnnualPred~cases2$Site)
summary(anovaCases)

# check residuals
resCases <- anovaCases$residuals
hist(resCases)

# Tukey's posthoc test to find out which sites differ
TukeyHSD(anovaCases)
