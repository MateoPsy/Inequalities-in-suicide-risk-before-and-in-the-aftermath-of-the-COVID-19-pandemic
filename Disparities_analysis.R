##########################
#### COVID-19 VG 2025 ####
##### DATA ANALYSIS  #####
###   Matej Rusiňák    ###
##########################

###### Loading pkg ##### 
library(tidyverse)  # data manipulation
library(boot) # bootstrapping - confidence intervals
library(purrr) # functional programming
library(ggplot2) # visualization
library(stringr)
library(marginaleffects)
library("lmtest")
library(sandwich)
library(sjPlot)
library(patchwork) # optional to patch graphs together
library(forestplot) # for forestplots

# just for multiple imputations
library(mice)
library(mitools)
library(lavaan.mi)
library(manymome)

###### LOADING DATASETS ####
merged_data <- read_csv("merged_data.csv") 
merged_data_diag <- read_csv("merged_data_diag.csv")
data_2022 <- read_csv("data_2022.csv")
data_2022$income <- as.factor(data_2022$income)

###### MULTIPLE IMPUTATIONS ####

### Imputing undisclosed income in 2022 ###
# Independent variables
# education
# sex
# age
# work status
# residence

set.seed(12436)
inc_MI <- data_2022
inc_MI$income[inc_MI$income == "no reply"] <- NA
inc_MI$response <- as.numeric(inc_MI$income)



# Estimate an OLS regression

fitols <- lm(response ~ education + sex + age + 
               work_status_orig + region_residence, data=inc_MI)

p_missing <- unlist(lapply(inc_MI, function(x) sum(is.na(x))))/nrow(inc_MI)
sort(p_missing[p_missing > 0], decreasing = TRUE)

# Select out variables that could cause problems in the imputation process
inc_MI <- select(inc_MI, ID, sex:size_place_residence, age:education, work_status:marital_status, response)


imp <- mice(inc_MI, maxit=0)

predM <- imp$predictorMatrix
meth <- imp$method



# With this command, we tell mice to impute the inc_MI data, create 20
# datasets, and use predM as the predictor matrix

imp2 <- mice(inc_MI, m = 20,
             predictorMatrix = predM, 
             method = meth, print =  TRUE, seed = 12436)

plot(imp2)


# preparing for regression (sr ~ under_10k)
filtered_data <- merged_data_diag[merged_data_diag$collect_wave == 2022,]
suicide_risk <- filtered_data$suicide_risk  # as a vector


# Convert to long format, including the original data
long_data <- complete(imp2, action = "long", include = TRUE)
long_data$response_binary <- ifelse(long_data$response == 1, 1, 0)
long_data$suicide_risk <- rep(suicide_risk, times = length(unique(long_data$.imp)))

# Convert the modified long format back into a `mids` object
imp_reconstructed <- as.mids(long_data)
summary(imp_reconstructed)


# Prevalence with CI (Monte Carlo method)

filtered_list <- list()

# Loop through imputed datasets and filter for income == 1
for (i in 1:imp_reconstructed$m) {
  d <- complete(imp_reconstructed, i)        # extract ith completed dataset
  d <- d[d$response == 1, ]      # filter only income == 1
  filtered_list[[i]] <- d
}

# Fit models on each filtered dataset
models <- lapply(filtered_list, function(d) lm(suicide_risk ~ 1, data = d))

# Combine using mitools::MIcombine
pooled <- MIcombine(models)

# Summary with Monte Carlo-based CI
summary(pooled) # robust, missInfo very little impact

est <- summary(pooled, conf.int = TRUE)
prevalence <- est$results
ci_low <- est$`(lower`
ci_high <- est$`upper)`



### Risk Ratio analysis (Delta method)

filtered_list_high <- list()

for (i in 1:imp_reconstructed$m) {
  d <- long_data[long_data$.imp == i, ]  # extract ith imputed dataset
  d <- d[d$response != 1, ]              # higher income group
  filtered_list_high[[i]] <- d
}



# Fit intercept-only model to estimate prevalence of suicide risk
models_high <- lapply(filtered_list_high, function(d) lm(suicide_risk ~ 1, data = d))

# Combine estimates
pooled_high <- MIcombine(models_high)
est_high <- summary(pooled_high, conf.int = TRUE)

prevalence_high <- est_high$results
ci_low_high <- est_high$`(lower`
ci_high_high <- est_high$`upper)`



# Calculate risk ratio
risk_ratio <- prevalence / prevalence_high

# Delta method approximation of CI
se_log_rr <- sqrt((est$se / prevalence)^2 + (est_high$se / prevalence_high)^2)
ci_lower_rr <- exp(log(risk_ratio) - 1.96 * se_log_rr)
ci_upper_rr <- exp(log(risk_ratio) + 1.96 * se_log_rr)


# final estimates
low_inc_statistics <- data.frame(c(est[,1:4], risk_ratio, ci_lower_rr, ci_upper_rr))
names(low_inc_statistics) <- c("prevalence", "se", "prev_low", "prev_high", "risk_ratio", "risk_ratio_low", "risk_ratio_high")

# write_csv(low_inc_statistics, "low_inc_statistics.csv")

###### ANALYSIS OF SR ####

#### SR IN POPULATION ##
model_sr_2017 <- glm(suicide_risk ~ 1, data = merged_data_diag[merged_data_diag$collect_wave == "2017",], weights = weight)
marginaleffects::hypotheses(model_sr_2017)

model_sr_2022 <- glm(suicide_risk ~ 1, data = merged_data_diag[merged_data_diag$collect_wave == "2022",], weights = weight)
marginaleffects::hypotheses(model_sr_2022)


#### SUPPORTING CODE ###
custom_order <- rev(c("major_depressive_episode", "anxiety_disorders",    "alcohol_use_disorders",    "comorbidity",
                      "any_mental_disorder", "Primary_ed",  "Unemployed",    "Economically_inactive",    "Low_income",              
                      "Not_relationship",    "Single_parenthood",  "Social_mobility",  "Young_females",  "Young_males",   "Older_females",    "Older_males"))


display_labels <- c(
  "Young_females" = "Young females",
  "Young_males" = "Young males",
  "Older_females" = "Older females",
  "Older_males" = "Older males",
  "Primary_ed" = "Primary education",
  "Unemployed" = "In unemployment",
  "Economically_inactive" = "Economically inactive",
  "Low_income" = "Low income",
  "Not_relationship" = "Not in a relationship",
  "Single_parenthood" = "Single parenthood",
  "Social_mobility" = "Downward social mobility",
  "major_depressive_episode" = "Major depressive episode",
  "anxiety_disorders" = "Anxiety disorders",
  "alcohol_use_disorders" = "Alcohol use disorders",
  "comorbidity" = "Co-occuring conditions",
  "any_mental_disorder" = "Any mental disorder"
)



#### PREVALENCE ###
prevalence_2017 <- merged_data_diag %>%
  filter(collect_wave == "2017") %>% 
  select(id, suicide_risk, Young_females:major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity) %>% 
  pivot_longer(cols = c(Young_females:major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity),
               names_to = "VG",
               values_to = "value") %>%
  filter(value == 1) %>% 
  group_by(VG) %>% 
  nest() %>% 
  mutate(LOGM = map(.x = data,
                    ~ glm(suicide_risk ~ 1, data = .x, family = binomial(link = "logit"))),
         CIs = map(LOGM,
                   ~ avg_predictions(.x)))  %>% 
  select(CIs)  %>% 
  unnest() %>% 
  mutate(across(c(estimate, conf.low, conf.high), ~ .x * 100),
         stars = case_when(p.value > 0.05 ~ "",
                           p.value <= 0.05 & p.value > 0.01 ~ "*",
                           p.value <= 0.01 & p.value > 0.001 ~ "**",
                           p.value <= 0.001 ~ "***"),
         prevalence = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                             " (", formatC(round(conf.low, 2), format = "f", digits = 2), ", ", 
                             formatC(round(conf.high, 2), format = "f", digits = 2), ")")) %>% 
  mutate(VG = factor(VG, levels = custom_order)) %>%  # apply your order as factor levels
  arrange(VG)  

# adding variables not included in 2017
low_income <- data.frame( VG    =  "Low_income",  term = NA, estimate = NA, std.error = NA,  statistic = NA, p.value = NA, s.value = NA, conf.low = NA, conf.high = NA, prevalence = NA)
Social_mobility <- data.frame( VG    =  "Social_mobility",  term = NA, estimate = NA, std.error = NA,  statistic = NA, p.value = NA, s.value = NA, conf.low = NA, conf.high = NA, prevalence = NA)
prevalence_2017 <- rbind(prevalence_2017[1:7,], low_income, prevalence_2017[8:9,], Social_mobility, prevalence_2017[10:14,])
# ordering dataset
prevalence_2017$VG <- factor(prevalence_2017$VG, levels = rev(custom_order))
prevalence_2017 <- prevalence_2017[order(prevalence_2017$VG), ]


prevalence_2022 <- merged_data_diag %>%
  filter(collect_wave == "2022") %>% 
  select(id, suicide_risk, Young_females:Social_mobility, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity) %>% 
  pivot_longer(cols = c(Young_females:major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity),
               names_to = "VG",
               values_to = "value") %>% 
  filter(value == 1) %>% 
  group_by(VG) %>% 
  nest() %>% 
  mutate(LOGM = map(.x = data,
                    ~ glm(suicide_risk ~ 1, data = .x, family = binomial(link = "logit"))),
         CIs = map(LOGM,
                   ~ avg_predictions(.x))) %>% 
  select(CIs)  %>% 
  unnest() %>% 
  mutate(across(c(estimate, conf.low, conf.high), ~ .x * 100),
         stars = case_when(p.value > 0.05 ~ "",
                           p.value <= 0.05 & p.value > 0.01 ~ "*",
                           p.value <= 0.01 & p.value > 0.001 ~ "**",
                           p.value <= 0.001 ~ "***"),
         prevalence = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                             " (", formatC(round(conf.low, 2), format = "f", digits = 2), ", ", 
                             formatC(round(conf.high, 2), format = "f", digits = 2), ")")) %>%
  mutate(VG = factor(VG, levels = rev(custom_order))) %>%  # apply your order as factor levels
  arrange(VG)  



#### RISK RATIO ###
risk_ratio_2017 <- merged_data_diag %>%
  filter(collect_wave == "2017") %>% 
  select(id, suicide_risk, Young_females:Social_mobility, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity) %>% 
  pivot_longer(cols = c(Young_females:major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity),
               names_to = "VG",
               values_to = "value") %>%
  group_by(VG) %>% 
  nest() %>% 
  mutate(LOGM = map(.x = data,
                    ~ glm(suicide_risk ~ value, data = .x, family = binomial(link = "logit"))),
         CIs = map(.x = LOGM, 
                   ~ marginaleffects::avg_comparisons(.x, comparison = "ratio")))  %>% 
  select(CIs)  %>%  
  unnest() %>% 
  mutate(stars = case_when(p.value > 0.05 ~ "",
                           p.value <= 0.05 & p.value > 0.01 ~ "*",
                           p.value <= 0.01 & p.value > 0.001 ~ "**",
                           p.value <= 0.001 ~ "***"),
         prevalence = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                             " (", formatC(round(conf.low, 2), format = "f", digits = 2), ", ", 
                             formatC(round(conf.high, 2), format = "f", digits = 2), ")")) %>%
  mutate(VG = factor(VG, levels = rev(custom_order))) %>%  # apply your order as factor levels
  arrange(VG) 


risk_ratio_2022 <- merged_data_diag %>%
  filter(collect_wave == "2022") %>% 
  select(id, suicide_risk, Young_females:Social_mobility, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity) %>% 
  pivot_longer(cols = c(Young_females:major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity),
               names_to = "VG",
               values_to = "value") %>%
  group_by(VG) %>% 
  nest() %>% 
  mutate(LOGM = map(.x = data,
                    ~ glm(suicide_risk ~ value, data = .x, family = binomial(link = "logit"))),
         CIs = map(.x = LOGM, 
                   ~ marginaleffects::avg_comparisons(.x, comparison = "ratio")))  %>% 
  select(CIs)  %>%  
  unnest() %>% 
  mutate(
    stars = case_when(p.value > 0.05 ~ "",
                      p.value <= 0.05 & p.value > 0.01 ~ "*",
                      p.value <= 0.01 & p.value > 0.001 ~ "**",
                      p.value <= 0.001 ~ "***"),
    prevalence = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                        " (", formatC(round(conf.low, 2), format = "f", digits = 2), ", ", 
                        formatC(round(conf.high, 2), format = "f", digits = 2), ")")) %>%
  mutate(VG = factor(VG, levels = rev(custom_order))) %>%  # apply your order as factor levels
  arrange(VG) 


#### TIME TRENDS ###
merged_data_diag$collect_wave <- factor(merged_data_diag$collect_wave)         

VG_lm <- function(df){
  lm(suicide_risk ~ collect_wave * value, data = df)
}

VG_glm <- function(df){
  glm(suicide_risk ~ collect_wave * value, family = binomial, data = df)
}

grid = data.frame(value = 0)


VG_ser_lm <- merged_data_diag %>%
  select(id, weight, collect_wave, suicide_risk, Young_females:Economically_inactive, Not_relationship:Single_parenthood, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity) %>% 
  pivot_longer(cols = c(Young_females:Economically_inactive, Not_relationship:Single_parenthood, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity),
               names_to = "VG",
               values_to = "value") %>%
  group_by(VG) %>% 
  nest() %>% 
  mutate(LPM_difndif = map(data, VG_lm),
         glm_mod = map(data, VG_glm),
         vcov = map(.x = LPM_difndif,
                    ~ coeftest(.x, vcov = vcovHC, type = "HC3")),
         confint = map(.x = vcov,
                       ~ broom::tidy(.x, conf.int = TRUE)),
         avg_comparisons = map(.x = glm_mod, 
                               ~ avg_comparisons(.x, variables = c("collect_wave","value"), cross = TRUE)),
         tinyliny = map(.x = glm_mod,
                        ~ avg_comparisons(.x,
                                          variables = "collect_wave",
                                          by = "value")),
         model_tidy = map(LPM_difndif, mice::tidy),
         est = map_dbl(.x = confint,
                       ~ as.data.frame(.x)[4, "estimate"]),
         ci = map(.x = confint,
                  ~ paste0("(",round((.x)[4,"conf.low"], digits = 3), ", ", round((.x)[4,"conf.high"], digits = 3),")")),
         plot = map2(LPM_difndif, confint, function(model, ci){
           
           new_data <- data.frame(
             value = c(0,0,1,1),
             collect_wave = factor(c("2017", "2022", "2017", "2022"), levels = c("2017", "2022"))
           )
           
           new_data$predicted_prob <- predict(model, newdata = new_data, type = "response")
           
           lower <- new_data$predicted_prob - qnorm(0.975)*ci[,"std.error"]
           upper <- new_data$predicted_prob + qnorm(0.975)*ci[,"std.error"]
           
           tibble(
             fitted_sr = new_data$predicted_prob,
             CIlow = unlist(lower),
             CIhigh = unlist(upper),
             value = new_data$value,
             collect_wave = new_data$collect_wave
           )
         })) 

Final_difndif <- VG_ser_lm %>%
  mutate(
    interaction_difndif = map(tinyliny, ~ {
      tryCatch({
        comp <- .x %>%
          select(value, estimate, std.error) %>%
          pivot_wider(names_from = value, values_from = c(estimate, std.error), names_prefix = "value_") %>%
          mutate(
            interaction = estimate_value_1 - estimate_value_0,
            se = sqrt(std.error_value_1^2 + std.error_value_0^2),
            z = interaction / se,
            p = 2 * (1 - pnorm(abs(z))),
            conf.low = interaction - 1.96 * se,
            conf.high = interaction + 1.96 * se
          ) %>%
          select(interaction, se, z, p, conf.low, conf.high)
        comp
      }, error = function(e) {
        tibble(interaction = NA, se = NA, z = NA, p = NA, conf.low = NA, conf.high = NA)
      })
    })
  )  %>% 
  select(interaction_difndif) %>% 
  unnest() %>% 
  mutate(across(c(interaction, conf.low, conf.high), ~ .x * 100),
         stars = case_when(p > 0.05 ~ "",
                           p <= 0.05 & p > 0.01 ~ "*",
                           p <= 0.01 & p > 0.001 ~ "**",
                           p <= 0.001 ~ "***"),
         estimate = paste0(formatC(round(interaction, 2), format = "f", digits = 2), 
                           " (", formatC(round(conf.low, 2), format = "f", digits = 2), ", ", 
                           formatC(round(conf.high, 2), format = "f", digits = 2), ")")) %>%
  mutate(VG = factor(VG, levels = custom_order)) %>%   # apply your order as factor levels
  arrange(VG)  


low_income_dif <- data.frame( VG    =  "Low_income",  interaction = NA, se = NA, z = NA,  p = NA, conf.low = NA, conf.high = NA, stars = NA, estimate = NA)
Social_mobility_dif <- data.frame( VG    =  "Social_mobility",  interaction = NA, se = NA, z = NA,  p = NA, conf.low = NA, conf.high = NA, stars = NA, estimate = NA)
Final_difndif <- rbind(Final_difndif[1:7,], low_income_dif, Final_difndif[8:9,], Social_mobility_dif, Final_difndif[10:14,])
# ordering dataset
Final_difndif$VG <- factor(Final_difndif$VG, levels = rev(custom_order))
Final_difndif <- Final_difndif[order(Final_difndif$VG), ]

###### FINAL TABLE INTO PUBLICATION ####

## checking for order 
factors <- list(prevalence_2017$VG, prevalence_2022$VG, risk_ratio_2017$VG, risk_ratio_2022$VG, Final_difndif$VG)

if (all(sapply(factors, function(f) levels(f) == levels(factors[[1]])))) {
  # All factor levels match in order – continue with your function
  
  P_RD_table <- data.frame(prevalence_2017$VG, prevalence_2017$prevalence, 
                           risk_ratio_2017$prevalence, prevalence_2022$prevalence, risk_ratio_2022$prevalence, Final_difndif$estimate)
  
} else {
  stop("Error: Factor levels do not match in order!")
}
P_RD_table

names(P_RD_table) <- c("VG", "Prevalence rate (95% CI) 2017", "Risk Ratio (95% CI) 2017", "Prevalence rate (95% CI) 2022", "Risk Ratio (95% CI) 2022", "Time Trends Diff-in-diff (95% CI)")

P_RD_table$`Prevalence rate (95% CI) 2022`[P_RD_table$VG == "Low_income"] <- paste0(round(low_inc_statistics$prevalence*100,2), " (", round(low_inc_statistics$prev_low*100, 2),
                                                                                    ", ", round(low_inc_statistics$prev_high*100, 2), ")")

P_RD_table$`Risk Ratio (95% CI) 2022`[P_RD_table$VG == "Low_income"] <- paste0(round(low_inc_statistics$risk_ratio,2), " (", round(low_inc_statistics$risk_ratio_low, 2),
                                                                               ", ", round(low_inc_statistics$risk_ratio_high, 2), ")")

levels(P_RD_table$VG) <- display_labels[levels(P_RD_table$VG)]

# we no longer consider economically inactive as it presents almost half of the population
P_RD_table_final <- P_RD_table[!P_RD_table$VG == "Economically inactive",]

## adding population statistics

summary_table <- read.csv("summary_table.csv")

# good order
summary_table$`Vulnerable_group` <- factor(summary_table$`Vulnerable_group`, levels = rev(custom_order))
summary_table <- summary_table[order(summary_table$`Vulnerable_group`), ]
levels(summary_table$`Vulnerable_group`) <- display_labels[levels(summary_table$`Vulnerable_group`)]

# we don´t need economically inactive
summary_table <- summary_table[!summary_table$`Vulnerable_group` == "Economically inactive", ]
TABLE_2 <- cbind(P_RD_table_final[,1], summary_table[,2], P_RD_table_final[,2:3], summary_table[,3], P_RD_table_final[,4:6])
names(TABLE_2) <- c("Vulnerable group", "N (%) 2017", "Prevalence rate (95% CI) 2017", "Risk ratio (95% CI) 2017",
                    "N (%) 2022", "Prevalence rate (95% CI) 2022", "Risk ratio (95% CI) 2022", "Time trend Diff-in-Diff (95% CI)")


write_csv(TABLE_2, "TABLE_2.csv")

#### SUPPORTING CODE ####
# let´s prepare for table and graph
names_vg <- factor(c("Young females",            "Young males",              "Older females",             
                     "Older males",                "Primary education",        "In unemployment",              
                     "Economically inactive",    "Not in relationship",      "Single parenthood",       
                     "Major depressive episode", "Anxiety disorders",        "Alcohol use disorders", 
                     "Any mental disorder", "Comorbidity"), levels = c("Major depressive episode", "Anxiety disorders", "Alcohol use disorders",  
                                                                       "Comorbidity", "Any mental disorder", "Primary education", "In unemployment",              
                                                                       "Economically inactive",  "Not in relationship",  "Single parenthood", 
                                                                       "Young females",    "Young males",    "Older females", "Older males"), ordered = TRUE)
VG_ser_lm$VG_names <- names_vg

#hey, comorbidity and any mental is just changed in order and I don´t have any mental power to find a way to change it
names_vg_difndif <- factor(c("Young females",            "Young males",              "Older females",             
                             "Older males",                "Primary education",        "In unemployment",              
                             "Economically inactive",    "Not in relationship",      "Single parenthood",       
                             "Major depressive episode", "Anxiety disorders",        "Alcohol use disorders", "Comorbidity", 
                             "Any mental disorder"), levels = rev(c("Major depressive episode", "Anxiety disorders", "Alcohol use disorders",  
                                                                    "Comorbidity", "Any mental disorder", "Primary education", "In unemployment",              
                                                                    "Economically inactive",  "Not in relationship",  "Single parenthood", 
                                                                    "Young females",    "Young males",    "Older females", "Older males")), ordered = TRUE)


###### STRUCTURAL BARRIERS ####

# when looking for help
nrow(data_2022[data_2022$structural_barriers == 1,]) # N how many people experienced it
sum(data_2022$structural_barriers, na.rm = TRUE)/sum(data_2022$help_seeking_12 == 1, na.rm = TRUE)*100 # proportion from those who looked for help (in %)

sum(data_2022$help_seeking_12 == 1, na.rm = TRUE)
# prevalence across data collection
model_SB <- lm(structural_barriers ~ 1, data = merged_data_diag[merged_data_diag$collect_wave == "2022" & merged_data_diag$help_seeking_12 ==1,])
marginaleffects::hypotheses(model_SB)

SB2022 <- merged_data_diag %>%
  filter(collect_wave == "2022" & help_seeking_12 == 1) %>% 
  select(id, structural_barriers, Young_females:Social_mobility, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity) %>% 
  pivot_longer(cols = c(Young_females:Social_mobility, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity),
               names_to = "VG",
               values_to = "value") %>% 
  filter(value == 1) %>% 
  group_by(VG) %>% 
  nest() %>% 
  mutate(LOGM = map(.x = data,
                    ~ glm(structural_barriers ~ 1, data = .x, family = binomial(link = "logit"))),
         CIs = map(LOGM,
                   ~ avg_predictions(.x))) %>% 
  select(CIs)  %>% 
  unnest() %>% 
  mutate(across(c(estimate, std.error, conf.low, conf.high), ~ .x * 100),
         conf.low = pmax(conf.low, 0),
         conf.high = pmin(conf.high, 100),
         stars = case_when(p.value > 0.05 ~ "",
                           p.value <= 0.05 & p.value > 0.01 ~ "*",
                           p.value <= 0.01 & p.value > 0.001 ~ "**",
                           p.value <= 0.001 ~ "***"),
         prevalence = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                             " (", formatC(round(conf.low, 2), format = "f", digits = 2), ", ", 
                             formatC(round(conf.high,2), format = "f", digits = 2), ")")) %>%
  mutate(VG = factor(VG, levels = custom_order, ordered = TRUE)) %>%  # apply your order as factor levels
  arrange(VG)  

SB2022_1 <- SB2022[order(SB2022$VG), ]
SB2022_1$VG <- display_labels[levels(SB2022_1$VG)]
SB2022_1 <- SB2022_1[nrow(SB2022_1):1, ]




SB2022_rr <- merged_data_diag %>%
  filter(collect_wave == "2022" & help_seeking_12 == 1) %>% 
  select(id, structural_barriers, Young_females:Social_mobility, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity) %>% 
  pivot_longer(cols = c(Young_females:Social_mobility, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity),
               names_to = "VG",
               values_to = "value") %>% 
  group_by(VG) %>% 
  nest() %>% 
  mutate(LOGM = map(.x = data,
                    ~ glm(structural_barriers ~ value, data = .x, family = binomial(link = "logit"))),
         CIs = map(LOGM,
                   ~ marginaleffects::avg_comparisons(.x, comparison = "ratio")))  %>% 
  select(CIs)  %>% 
  unnest() %>% 
  mutate(conf.low = pmax(conf.low, 0),
         stars = case_when(p.value > 0.05 ~ "",
                           p.value <= 0.05 & p.value > 0.01 ~ "*",
                           p.value <= 0.01 & p.value > 0.001 ~ "**",
                           p.value <= 0.001 ~ "***"),
         prevalence = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                             " (", formatC(round(conf.low, 2), format = "f", digits = 2), ", ", 
                             formatC(round(conf.high, 2), format = "f", digits = 2), ")")) %>%
  mutate(VG = factor(VG, levels = custom_order)) %>%  # apply your order as factor levels
  arrange(VG)  

SB2022_2 <- SB2022_rr[order(SB2022_rr$VG), ]
SB2022_2$VG <- display_labels[levels(SB2022_2$VG)]
SB2022_2 <- SB2022_2[nrow(SB2022_2):1, ]


## plotting
plot_levels <- c(
  custom_order[1:4],
  "gap1",
  custom_order[5:11],
  "gap2",
  custom_order[12:16]
)




SB2022 <- SB2022 %>%
  mutate(VG = factor(VG, levels = rev(plot_levels)))

SB2022_rrx1 <- SB2022_rr %>%
  mutate(VG = factor(VG, levels = rev(plot_levels)))


# we worked with variable economically inactive but later we dropped it
SB2022 <- SB2022[!SB2022$VG == "Economically_inactive", ]
SB2022_rrx1 <- SB2022_rrx1[!SB2022_rrx1$VG == "Economically_inactive", ]



plot_labels <- setNames(
  c(
    display_labels[4:1], "",
    display_labels[11:8],
    display_labels[6:5], "",
    display_labels[16:12]
  ),
  plot_levels[!plot_levels == "Economically_inactive"]
)



plot_SB <- SB2022 %>% 
  ggplot(aes(x = estimate, y = VG)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.3) +
  scale_y_discrete(limits = plot_levels[!plot_levels == "Economically_inactive"], labels = plot_labels, drop = FALSE) +
  labs(
    title = "Prevalence of Structural barriers (2022)",
    x = "Prevalence rate in % (95% confidence interval)",
    y = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.y = element_blank())


plot_SBrr <- SB2022_rrx1 %>%
  ggplot(aes(x = estimate, y = VG)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.3) +
  scale_y_discrete(limits = plot_levels[!plot_levels == "Economically_inactive"], labels = plot_labels, drop = FALSE) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  labs(
    title = "Risk Ratios of Structural barriers (2022)",
    x = "Risk ratio (95% confidence interval)",
    y = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.y = element_blank())

combined_SB_1 <- plot_SB + plot_SBrr
# ggsave("plot_SB_rr.pdf", combined_SB_1, width = 20, height = 7)

### saving tables
# write.csv(SB2022_1, "SB2022_final.csv")
# write.csv(SB2022_2, "SB2022_rr_final.csv")


###### SENSITIVITY ANALYSIS ####

# prevalences
sr_thoughts_prevalence_2017 <- merged_data_diag %>%
  filter(collect_wave == "2017") %>% 
  select(id, suicide_thoughts, Young_females:Social_mobility, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity) %>% 
  pivot_longer(cols = c(Young_females:major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity),
               names_to = "VG",
               values_to = "value") %>% 
  filter(value == 1) %>% 
  group_by(VG) %>% 
  nest() %>% 
  mutate(LOGM = map(.x = data,
                    ~ glm(suicide_thoughts ~ 1, data = .x, family = binomial(link = "logit"))),
         CIs = map(LOGM,
                   ~ avg_predictions(.x))) %>% 
  select(CIs)  %>% 
  unnest() %>% 
  mutate(across(c(estimate, conf.low, conf.high), ~ .x * 100),
         conf.low = pmax(conf.low, 0),
         conf.high = pmin(conf.high, 100),
         stars = case_when(p.value > 0.05 ~ "",
                           p.value <= 0.05 & p.value > 0.01 ~ "*",
                           p.value <= 0.01 & p.value > 0.001 ~ "**",
                           p.value <= 0.001 ~ "***"),
         prevalence = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                             " (", formatC(round(conf.low, 2), format = "f", digits = 2), ", ", 
                             formatC(round(conf.high, 2), format = "f", digits = 2), ")")) %>%
  mutate(VG = factor(VG, levels = custom_order)) %>%  # apply your order as factor levels
  arrange(VG)  

sr_behavior_prevalence_2017 <- merged_data_diag %>%
  filter(collect_wave == "2017") %>% 
  select(id, suicide_behavior, Young_females:Social_mobility, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity) %>% 
  pivot_longer(cols = c(Young_females:major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity),
               names_to = "VG",
               values_to = "value") %>% 
  filter(value == 1) %>% 
  group_by(VG) %>% 
  nest() %>% 
  mutate(LOGM = map(.x = data,
                    ~ glm(suicide_behavior ~ 1, data = .x, family = binomial(link = "logit"))),
         CIs = map(LOGM,
                   ~ avg_predictions(.x))) %>% 
  select(CIs)  %>% 
  unnest() %>% 
  mutate(across(c(estimate, conf.low, conf.high), ~ .x * 100),
         conf.low = pmax(conf.low, 0),
         conf.high = pmin(conf.high, 100),
         stars = case_when(p.value > 0.05 ~ "",
                           p.value <= 0.05 & p.value > 0.01 ~ "*",
                           p.value <= 0.01 & p.value > 0.001 ~ "**",
                           p.value <= 0.001 ~ "***"),
         prevalence = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                             " (", formatC(round(conf.low, 2), format = "f", digits = 2), ", ", 
                             formatC(round(conf.high, 2), format = "f", digits = 2), ")")) %>%
  mutate(VG = factor(VG, levels = custom_order)) %>%  # apply your order as factor levels
  arrange(VG)  

sr_lifelong_excl_prevalence_2017 <- merged_data_diag %>%
  filter(collect_wave == "2017") %>% 
  select(id, suicide_risk_without_lifelong, Young_females:Social_mobility, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity) %>% 
  pivot_longer(cols = c(Young_females:major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity),
               names_to = "VG",
               values_to = "value") %>% 
  filter(value == 1) %>% 
  group_by(VG) %>% 
  nest() %>% 
  mutate(LOGM = map(.x = data,
                    ~ glm(suicide_risk_without_lifelong ~ 1, data = .x, family = binomial(link = "logit"))),
         CIs = map(LOGM,
                   ~ avg_predictions(.x))) %>% 
  select(CIs)  %>% 
  unnest() %>% 
  mutate(across(c(estimate, conf.low, conf.high), ~ .x * 100),
         conf.low = pmax(conf.low, 0),
         conf.high = pmin(conf.high, 100),
         stars = case_when(p.value > 0.05 ~ "",
                           p.value <= 0.05 & p.value > 0.01 ~ "*",
                           p.value <= 0.01 & p.value > 0.001 ~ "**",
                           p.value <= 0.001 ~ "***"),
         prevalence = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                             " (", formatC(round(conf.low, 2), format = "f", digits = 2), ", ", 
                             formatC(round(conf.high, 2), format = "f", digits = 2), ")")) %>%
  mutate(VG = factor(VG, levels = custom_order)) %>%  # apply your order as factor levels
  arrange(VG) 

# yr2022
sr_thoughts_prevalence_2022 <- merged_data_diag %>%
  filter(collect_wave == "2022") %>% 
  select(id, suicide_thoughts, Young_females:Social_mobility, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity) %>% 
  pivot_longer(cols = c(Young_females:major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity),
               names_to = "VG",
               values_to = "value") %>% 
  filter(value == 1) %>% 
  group_by(VG) %>% 
  nest() %>% 
  mutate(LOGM = map(.x = data,
                    ~ glm(suicide_thoughts ~ 1, data = .x, family = binomial(link = "logit"))),
         CIs = map(LOGM,
                   ~ avg_predictions(.x))) %>% 
  select(CIs)  %>% 
  unnest() %>% 
  mutate(across(c(estimate, conf.low, conf.high), ~ .x * 100),
         conf.low = pmax(conf.low, 0),
         conf.high = pmin(conf.high, 100),
         stars = case_when(p.value > 0.05 ~ "",
                           p.value <= 0.05 & p.value > 0.01 ~ "*",
                           p.value <= 0.01 & p.value > 0.001 ~ "**",
                           p.value <= 0.001 ~ "***"),
         prevalence = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                             " (", formatC(round(conf.low, 2), format = "f", digits = 2), ", ", 
                             formatC(round(conf.high, 2), format = "f", digits = 2), ")")) %>%
  mutate(VG = factor(VG, levels = custom_order)) %>%  # apply your order as factor levels
  arrange(VG)  

sr_behavior_prevalence_2022 <- merged_data_diag %>%
  filter(collect_wave == "2022") %>% 
  select(id, suicide_behavior, Young_females:Social_mobility, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity) %>% 
  pivot_longer(cols = c(Young_females:major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity),
               names_to = "VG",
               values_to = "value") %>% 
  filter(value == 1) %>% 
  group_by(VG) %>% 
  nest() %>% 
  mutate(LOGM = map(.x = data,
                    ~ glm(suicide_behavior ~ 1, data = .x, family = binomial(link = "logit"))),
         CIs = map(LOGM,
                   ~ avg_predictions(.x))) %>% 
  select(CIs)  %>% 
  unnest() %>% 
  mutate(across(c(estimate, conf.low, conf.high), ~ .x * 100),
         conf.low = pmax(conf.low, 0),
         conf.high = pmin(conf.high, 100),
         stars = case_when(p.value > 0.05 ~ "",
                           p.value <= 0.05 & p.value > 0.01 ~ "*",
                           p.value <= 0.01 & p.value > 0.001 ~ "**",
                           p.value <= 0.001 ~ "***"),
         prevalence = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                             " (", formatC(round(conf.low, 2), format = "f", digits = 2), ", ", 
                             formatC(round(conf.high, 2), format = "f", digits = 2), ")")) %>%
  mutate(VG = factor(VG, levels = custom_order)) %>%  # apply your order as factor levels
  arrange(VG)  

sr_lifelong_excl_prevalence_2022 <- merged_data_diag %>%
  filter(collect_wave == "2022") %>% 
  select(id, suicide_risk_without_lifelong, Young_females:Social_mobility, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity) %>% 
  pivot_longer(cols = c(Young_females:major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity),
               names_to = "VG",
               values_to = "value") %>% 
  filter(value == 1) %>% 
  group_by(VG) %>% 
  nest() %>% 
  mutate(LOGM = map(.x = data,
                    ~ glm(suicide_risk_without_lifelong ~ 1, data = .x, family = binomial(link = "logit"))),
         CIs = map(LOGM,
                   ~ avg_predictions(.x))) %>% 
  select(CIs)  %>% 
  unnest() %>% 
  mutate(across(c(estimate, conf.low, conf.high), ~ .x * 100),
         conf.low = pmax(conf.low, 0),
         conf.high = pmin(conf.high, 100),
         stars = case_when(p.value > 0.05 ~ "",
                           p.value <= 0.05 & p.value > 0.01 ~ "*",
                           p.value <= 0.01 & p.value > 0.001 ~ "**",
                           p.value <= 0.001 ~ "***"),
         prevalence = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                             " (", formatC(round(conf.low, 2), format = "f", digits = 2), ", ", 
                             formatC(round(conf.high, 2), format = "f", digits = 2), ")")) %>%
  mutate(VG = factor(VG, levels = custom_order)) %>%  # apply your order as factor levels
  arrange(VG) 

sr_nonsuic_selfi_prevalence_2022 <- merged_data_diag %>%
  filter(collect_wave == "2022") %>% 
  select(id, Nonsuic_self_i, Young_females:Social_mobility, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity) %>% 
  pivot_longer(cols = c(Young_females:major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity),
               names_to = "VG",
               values_to = "value") %>% 
  filter(value == 1) %>% 
  group_by(VG) %>% 
  nest() %>% 
  mutate(LOGM = map(.x = data,
                    ~ glm(Nonsuic_self_i ~ 1, data = .x, family = binomial(link = "logit"))),
         CIs = map(LOGM,
                   ~ avg_predictions(.x))) %>% 
  select(CIs)  %>% 
  unnest() %>% 
  mutate(across(c(estimate, conf.low, conf.high), ~ .x * 100),
         conf.low = pmax(conf.low, 0),
         conf.high = pmin(conf.high, 100),
         stars = case_when(p.value > 0.05 ~ "",
                           p.value <= 0.05 & p.value > 0.01 ~ "*",
                           p.value <= 0.01 & p.value > 0.001 ~ "**",
                           p.value <= 0.001 ~ "***"),
         prevalence = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                             " (", formatC(round(conf.low, 2), format = "f", digits = 2), ", ", 
                             formatC(round(conf.high, 2), format = "f", digits = 2), ")")) %>%
  mutate(VG = factor(VG, levels = custom_order)) %>%  # apply your order as factor levels
  arrange(VG)  



## Risk Ratios
sr_thoughts_risk_ratio_2017 <- merged_data_diag %>%
  filter(collect_wave == "2017") %>% 
  select(id, suicide_thoughts, Young_females:Social_mobility, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity) %>% 
  pivot_longer(cols = c(Young_females:major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity),
               names_to = "VG",
               values_to = "value") %>%
  group_by(VG) %>% 
  nest() %>% 
  mutate(LOGM = map(.x = data,
                    ~ glm(suicide_thoughts ~ value, data = .x, family = binomial(link = "logit"))),
         CIs = map(.x = LOGM, 
                   ~ marginaleffects::avg_comparisons(.x, comparison = "ratio")))  %>% 
  select(CIs)  %>%  
  unnest() %>% 
  mutate(
    conf.low = pmax(conf.low, 0),
    stars = case_when(p.value > 0.05 ~ "",
                      p.value <= 0.05 & p.value > 0.01 ~ "*",
                      p.value <= 0.01 & p.value > 0.001 ~ "**",
                      p.value <= 0.001 ~ "***"),
    prevalence = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                        " (", formatC(round(conf.low, 2), format = "f", digits = 2), ", ", 
                        formatC(round(conf.high, 2), format = "f", digits = 2), ")")) %>%
  mutate(VG = factor(VG, levels = custom_order)) %>%  # apply your order as factor levels
  arrange(VG) 

sr_behavior_risk_ratio_2017 <- merged_data_diag %>%
  filter(collect_wave == "2017") %>% 
  select(id, suicide_behavior, Young_females:Social_mobility, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity) %>% 
  pivot_longer(cols = c(Young_females:major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity),
               names_to = "VG",
               values_to = "value") %>%
  group_by(VG) %>% 
  nest() %>% 
  mutate(LOGM = map(.x = data,
                    ~ glm(suicide_behavior ~ value, data = .x, family = binomial(link = "logit"))),
         CIs = map(.x = LOGM, 
                   ~ marginaleffects::avg_comparisons(.x, comparison = "ratio")))  %>% 
  select(CIs)  %>%  
  unnest() %>% 
  mutate(
    conf.low = pmax(conf.low, 0),
    stars = case_when(p.value > 0.05 ~ "",
                      p.value <= 0.05 & p.value > 0.01 ~ "*",
                      p.value <= 0.01 & p.value > 0.001 ~ "**",
                      p.value <= 0.001 ~ "***"),
    prevalence = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                        " (", formatC(round(conf.low, 2), format = "f", digits = 2), ", ", 
                        formatC(round(conf.high, 2), format = "f", digits = 2), ")")) %>%
  mutate(VG = factor(VG, levels = custom_order)) %>%  # apply your order as factor levels
  arrange(VG) 

sr_lifelong_excl_risk_ratio_2017 <- merged_data_diag %>%
  filter(collect_wave == "2017") %>% 
  select(id, suicide_risk_without_lifelong, Young_females:Social_mobility, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity) %>% 
  pivot_longer(cols = c(Young_females:major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity),
               names_to = "VG",
               values_to = "value") %>%
  group_by(VG) %>% 
  nest() %>% 
  mutate(LOGM = map(.x = data,
                    ~ glm(suicide_risk_without_lifelong ~ value, data = .x, family = binomial(link = "logit"))),
         CIs = map(.x = LOGM, 
                   ~ marginaleffects::avg_comparisons(.x, comparison = "ratio")))  %>% 
  select(CIs)  %>%  
  unnest() %>% 
  mutate(
    conf.low = pmax(conf.low, 0),
    stars = case_when(p.value > 0.05 ~ "",
                      p.value <= 0.05 & p.value > 0.01 ~ "*",
                      p.value <= 0.01 & p.value > 0.001 ~ "**",
                      p.value <= 0.001 ~ "***"),
    prevalence = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                        " (", formatC(round(conf.low, 2), format = "f", digits = 2), ", ", 
                        formatC(round(conf.high, 2), format = "f", digits = 2), ")")) %>%
  mutate(VG = factor(VG, levels = custom_order)) %>%  # apply your order as factor levels
  arrange(VG) 


# year2022
sr_thoughts_risk_ratio_2022 <- merged_data_diag %>%
  filter(collect_wave == "2022") %>% 
  select(id, suicide_thoughts, Young_females:Social_mobility, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity) %>% 
  pivot_longer(cols = c(Young_females:major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity),
               names_to = "VG",
               values_to = "value") %>%
  group_by(VG) %>% 
  nest() %>% 
  mutate(LOGM = map(.x = data,
                    ~ glm(suicide_thoughts ~ value, data = .x, family = binomial(link = "logit"))),
         CIs = map(.x = LOGM, 
                   ~ marginaleffects::avg_comparisons(.x, comparison = "ratio")))  %>% 
  select(CIs)  %>%  
  unnest() %>% 
  mutate(conf.low = pmax(conf.low, 0),
         stars = case_when(p.value > 0.05 ~ "",
                           p.value <= 0.05 & p.value > 0.01 ~ "*",
                           p.value <= 0.01 & p.value > 0.001 ~ "**",
                           p.value <= 0.001 ~ "***"),
         prevalence = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                             " (", formatC(round(conf.low, 2), format = "f", digits = 2), ", ", 
                             formatC(round(conf.high, 2), format = "f", digits = 2), ")")) %>%
  mutate(VG = factor(VG, levels = custom_order)) %>%  # apply your order as factor levels
  arrange(VG) 

sr_behavior_risk_ratio_2022 <- merged_data_diag %>%
  filter(collect_wave == "2022") %>% 
  select(id, suicide_behavior, Young_females:Social_mobility, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity) %>% 
  pivot_longer(cols = c(Young_females:major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity),
               names_to = "VG",
               values_to = "value") %>%
  group_by(VG) %>% 
  nest() %>% 
  mutate(LOGM = map(.x = data,
                    ~ glm(suicide_behavior ~ value, data = .x, family = binomial(link = "logit"))),
         CIs = map(.x = LOGM, 
                   ~ marginaleffects::avg_comparisons(.x, comparison = "ratio")))  %>% 
  select(CIs)  %>%  
  unnest() %>% 
  mutate(
    conf.low = pmax(conf.low, 0),
    stars = case_when(p.value > 0.05 ~ "",
                      p.value <= 0.05 & p.value > 0.01 ~ "*",
                      p.value <= 0.01 & p.value > 0.001 ~ "**",
                      p.value <= 0.001 ~ "***"),
    prevalence = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                        " (", formatC(round(conf.low, 2), format = "f", digits = 2), ", ", 
                        formatC(round(conf.high, 2), format = "f", digits = 2), ")")) %>%
  mutate(VG = factor(VG, levels = custom_order)) %>%  # apply your order as factor levels
  arrange(VG) 

sr_lifelong_excl_risk_ratio_2022 <- merged_data_diag %>%
  filter(collect_wave == "2022") %>% 
  select(id, suicide_risk_without_lifelong, Young_females:Social_mobility, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity) %>% 
  pivot_longer(cols = c(Young_females:major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity),
               names_to = "VG",
               values_to = "value") %>%
  group_by(VG) %>% 
  nest() %>% 
  mutate(LOGM = map(.x = data,
                    ~ glm(suicide_risk_without_lifelong ~ value, data = .x, family = binomial(link = "logit"))),
         CIs = map(.x = LOGM, 
                   ~ marginaleffects::avg_comparisons(.x, comparison = "ratio")))  %>% 
  select(CIs)  %>%  
  unnest() %>% 
  mutate(
    conf.low = pmax(conf.low, 0),
    stars = case_when(p.value > 0.05 ~ "",
                      p.value <= 0.05 & p.value > 0.01 ~ "*",
                      p.value <= 0.01 & p.value > 0.001 ~ "**",
                      p.value <= 0.001 ~ "***"),
    prevalence = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                        " (", formatC(round(conf.low, 2), format = "f", digits = 2), ", ", 
                        formatC(round(conf.high, 2), format = "f", digits = 2), ")")) %>%
  mutate(VG = factor(VG, levels = custom_order)) %>%  # apply your order as factor levels
  arrange(VG) 

sr_nonsuic_selfi_risk_ratio_2022 <- merged_data_diag %>%
  filter(collect_wave == "2022") %>% 
  select(id, Nonsuic_self_i, Young_females:Social_mobility, major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity) %>% 
  pivot_longer(cols = c(Young_females:major_depressive_episode, anxiety_disorders:any_mental_disorder, comorbidity),
               names_to = "VG",
               values_to = "value") %>%
  group_by(VG) %>% 
  nest() %>% 
  mutate(LOGM = map(.x = data,
                    ~ glm(Nonsuic_self_i ~ value, data = .x, family = binomial(link = "logit"))),
         CIs = map(.x = LOGM, 
                   ~ marginaleffects::avg_comparisons(.x, comparison = "ratio")))  %>% 
  select(CIs)  %>%  
  unnest() %>% 
  mutate(
    conf.low = pmax(conf.low, 0),
    stars = case_when(p.value > 0.05 ~ "",
                      p.value <= 0.05 & p.value > 0.01 ~ "*",
                      p.value <= 0.01 & p.value > 0.001 ~ "**",
                      p.value <= 0.001 ~ "***"),
    prevalence = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                        " (", formatC(round(conf.low, 2), format = "f", digits = 2), ", ", 
                        formatC(round(conf.high, 2), format = "f", digits = 2), ")")) %>%
  mutate(VG = factor(VG, levels = custom_order)) %>%  # apply your order as factor levels
  arrange(VG) 



library(forcats) # mby

plot_levels <- c(
  custom_order[1:4],
  "gap1",
  custom_order[5:11],
  "gap2",
  custom_order[12:16]
)



plot_labels <- setNames(
  c(
    display_labels[4:1], "",
    display_labels[11:5], "",
    display_labels[16:12]
  ),
  plot_levels
)

cb_colors <- c(
  "Thoughts" = "#0072B2" ,
  "Behavior" = "#E69F00",
  "Lifelong Excl." =  "#009E73",
  "Total - Main Analyses" = "#a53606", 
  "NSSI"          = "#b32db5"
)

# combined prevalence_2017
plot_2017 <- bind_rows(
  sr_thoughts_prevalence_2017       %>% mutate(source = "Thoughts"),
  sr_behavior_prevalence_2017       %>% mutate(source = "Behavior"),
  sr_lifelong_excl_prevalence_2017  %>% mutate(source = "Lifelong Excl."),
  prevalence_2017                   %>% mutate(source = "Total - Main Analyses")
) %>%
  mutate(VG = factor(VG, levels = rev(plot_levels))) %>%
  ggplot(aes(x = estimate, y = fct_rev(VG), color = source)) +
  geom_point(position = position_dodge(width = 0.7), size = 2.6) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 position = position_dodge(width = 0.7), height = 0.3) +
  scale_y_discrete(labels = plot_labels, drop = FALSE,  expand = c(0,0)) +
  scale_x_continuous(n.breaks = 10) +
  scale_color_manual(values = cb_colors) + 
  labs(title = "Suicide Risk Prevalence by Vulnerability Group (2017)",
       x = "Prevalence rate in % (95% confidence interval)",
       y = NULL,
       color = "Risk Type") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "bottom")


plot_levels <- c(
  custom_order[1:4],
  "gap1",
  custom_order[5:11],
  "gap2",
  custom_order[12:16]
)



plot_labels <- setNames(
  c(
    display_labels[4:1], "",
    display_labels[11:5], "",
    display_labels[16:12]
  ),
  plot_levels
)


# combined prevalence_2022
plot_2022 <- bind_rows(
  sr_thoughts_prevalence_2022       %>% mutate(source = "Thoughts"),
  sr_behavior_prevalence_2022       %>% mutate(source = "Behavior"),
  sr_lifelong_excl_prevalence_2022  %>% mutate(source = "Lifelong Excl."),
  sr_nonsuic_selfi_prevalence_2022  %>% mutate(source = "NSSI"),
  prevalence_2022                   %>% mutate(source = "Total - Main Analyses")
) %>%
  mutate(VG = factor(VG, levels = rev(plot_levels))) %>%
  ggplot(aes(x = estimate, y = fct_rev(VG), color = source)) +
  geom_point(position = position_dodge(width = 0.7), size = 2.6) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 position = position_dodge(width = 0.7), height = 0.3) +
  scale_y_discrete(labels = plot_labels, drop = FALSE,  expand = c(0,0)) +
  scale_x_continuous(n.breaks = 10) +
  scale_color_manual(values = cb_colors) + 
  labs(title = "Suicide Risk Prevalence by Vulnerability Group (2022)",
       x = "Prevalence rate in % (95% confidence interval)",
       y = NULL,
       color = "Risk Type") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "bottom")




# risk ratio_2017
plot_rr_2017 <- bind_rows(
  sr_thoughts_risk_ratio_2017       %>% mutate(source = "Thoughts"),
  sr_behavior_risk_ratio_2017       %>% mutate(source = "Behavior"),
  sr_lifelong_excl_risk_ratio_2017  %>% mutate(source = "Lifelong Excl."),
  risk_ratio_2017                   %>% mutate(source = "Total - Main Analyses")
) %>%
  mutate(
    across(estimate,
           ~ ifelse(estimate == 1, NA, .x))
  ) %>%
  mutate(VG = factor(VG, levels = rev(plot_levels))) %>%
  ggplot(aes(x = estimate, y = fct_rev(VG), color = source)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_point(position = position_dodge(width = 0.7), size = 2.6) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 position = position_dodge(width = 0.7), height = 0.3) +
  scale_y_discrete(labels = plot_labels, drop = FALSE,  expand = c(0,0)) +
  scale_x_continuous(n.breaks = 10) +
  scale_color_manual(values = cb_colors) + 
  labs(title = "Risk Ratios of Suicide Risk by Vulnerability Group (2017)",
       x = "Risk ratio (95% confidence interval)",
       y = NULL,
       color = "Risk Type") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "bottom")


# risk ratio_2022
plot_rr_2022 <- bind_rows(
  sr_thoughts_risk_ratio_2022       %>% mutate(source = "Thoughts"),
  sr_behavior_risk_ratio_2022       %>% mutate(source = "Behavior"),
  sr_lifelong_excl_risk_ratio_2022  %>% mutate(source = "Lifelong Excl."),
  sr_nonsuic_selfi_risk_ratio_2022  %>% mutate(source = "NSSI"),
  risk_ratio_2022                   %>% mutate(source = "Total - Main Analyses")
) %>%
  mutate(VG = factor(VG, levels = rev(plot_levels))) %>%
  ggplot(aes(x = estimate, y = fct_rev(VG), color = source)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_point(position = position_dodge(width = 0.7), size = 2.6) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 position = position_dodge(width = 0.7), height = 0.3) +
  scale_y_discrete(labels = plot_labels, drop = FALSE,  expand = c(0,0)) +
  scale_x_continuous(n.breaks = 10) +
  scale_color_manual(values = cb_colors) + 
  labs(title = "Risk Ratios of Suicide Risk by Vulnerability Group (2022)",
       x = "Risk ratio (95% confidence interval)",
       y = NULL,
       color = "Risk Type") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "bottom")




# displaying together

combined_prev_plot <- plot_2017 + plot_2022

# ggsave("combined_prev_plot.pdf", combined_prev_plot, width = 20, height = 7)

combined_rr_plot <- plot_rr_2017 + plot_rr_2022

# ggsave("combined_rr_plot.pdf", combined_rr_plot, width = 20, height = 7)


## Tables for sensitivity analysis


display_labels1 <- display_labels
display_labels1 <- display_labels1[!names(display_labels1) %in% c("Low_income", "Social_mobility")]

# PREVALENCES

# list of all data frame names
dfs <- c(
  "sr_thoughts_prevalence_2017",
  "sr_behavior_prevalence_2017",
  "sr_lifelong_excl_prevalence_2017"
)

# apply your transformation to each
for (df_name in dfs) {
  df <- get(df_name)
  
  df$VG <- factor(df$VG)
  df$VG <- display_labels1[levels(df$VG)]
  df <- df[nrow(df):1, ]
  
  assign(df_name, df, envir = .GlobalEnv)
}

dfs <- c(
  "sr_thoughts_prevalence_2022",
  "sr_behavior_prevalence_2022",
  "sr_lifelong_excl_prevalence_2022",
  "sr_nonsuic_selfi_prevalence_2022"
)

# apply your transformation to each
for (df_name in dfs) {
  df <- get(df_name)
  
  df$VG <- factor(df$VG)
  df$VG <- display_labels[levels(df$VG)]
  df <- df[nrow(df):1, ]
  
  assign(df_name, df, envir = .GlobalEnv)
}

sr_thoughts_prevalence_2017
sr_behavior_prevalence_2017
sr_lifelong_excl_prevalence_2017

sensitivity_prevalence_2017 <- data.frame(sr_thoughts_prevalence_2017$VG,
                                          sr_thoughts_prevalence_2017$prevalence,
                                          sr_behavior_prevalence_2017$prevalence,
                                          sr_lifelong_excl_prevalence_2017$prevalence)
names(sensitivity_prevalence_2017) <- c("VG", "SR Thoughts", "SR Behavior", "SR Lifelong Excl.")
# write.csv(sensitivity_prevalence_2017, "sensitivity_prevalence_2017.csv", row.names = FALSE)



sensitivity_prevalence_2022 <- data.frame(sr_thoughts_prevalence_2022$VG,
                                          sr_thoughts_prevalence_2022$prevalence,
                                          sr_behavior_prevalence_2022$prevalence,
                                          sr_lifelong_excl_prevalence_2022$prevalence,
                                          sr_nonsuic_selfi_prevalence_2022$prevalence)
names(sensitivity_prevalence_2022) <- c("VG", "SR Thoughts", "SR Behavior", "SR Lifelong Excl.", "NSSI")
# write.csv(sensitivity_prevalence_2022, "sensitivity_prevalence_2022.csv", row.names = FALSE)

# Risk ratios
sr_lifelong_excl_risk_ratio_2022
# list of all data frame names
dfs <- c(
  "sr_thoughts_risk_ratio_2017",
  "sr_behavior_risk_ratio_2017",
  "sr_lifelong_excl_risk_ratio_2017",
  "sr_thoughts_risk_ratio_2022",
  "sr_behavior_risk_ratio_2022",
  "sr_lifelong_excl_risk_ratio_2022",
  "sr_nonsuic_selfi_risk_ratio_2022"
)

# apply your transformation to each
for (df_name in dfs) {
  df <- get(df_name)
  
  df$VG <- factor(df$VG)
  df$VG <- display_labels[levels(df$VG)]
  df <- df[nrow(df):1, ]
  
  assign(df_name, df, envir = .GlobalEnv)
}


sensitivity_risk_ratios <- data.frame(sr_thoughts_risk_ratio_2017$VG,
                                      sr_thoughts_risk_ratio_2017$prevalence,
                                      sr_behavior_risk_ratio_2017$prevalence,
                                      sr_lifelong_excl_risk_ratio_2017$prevalence,
                                      sr_thoughts_risk_ratio_2022$prevalence,
                                      sr_behavior_risk_ratio_2022$prevalence,
                                      sr_lifelong_excl_risk_ratio_2022$prevalence,
                                      sr_nonsuic_selfi_risk_ratio_2022$prevalence
)
names(sensitivity_risk_ratios) <- c("VG", "SR Thoughts", "SR Behavior", "SR Lifelong Excl.", "SR Thoughts-2022", "SR Behavior-2022", "SR Lifelong Excl.-2022", "NSSI-2022")
# write.csv(sensitivity_risk_ratios, "sensitivity_risk_ratios.csv", row.names = FALSE)


###### FOREST PLOTS COMBINED ####

# empty rows so the graph looks nice
empty_row <- data.frame(VG    =  "",  term = NA, estimate = NA, std.error = NA,  statistic = NA, p.value = NA, s.value = NA, conf.low = NA, conf.high = NA, prevalence = NA, VG_names = NA)
empty_row1 <- data.frame(VG    =  " ",  term = NA, estimate = NA, std.error = NA,  statistic = NA, p.value = NA, s.value = NA, conf.low = NA, conf.high = NA, prevalence = NA, VG_names = NA)

# prevalences
prevalence_2017_plot <- prevalence_2017
prevalence_2017_plot$VG <- display_labels[levels(prevalence_2017_plot$VG)]
prevalence_2017_plot <- rbind(prevalence_2017_plot[1:5,], empty_row, prevalence_2017_plot[6:7,], prevalence_2017_plot[9:12,], empty_row1, prevalence_2017_plot[13:nrow(prevalence_2017_plot),]) # this specific way, so we don´t have economically inactive in there
prevalence_2017_plot$wave <- "2017"

prevalence_2017_plot$VG <- factor(prevalence_2017_plot$VG)


prevalence_2022_plot <- prevalence_2022
prevalence_2022_plot$VG <- display_labels[levels(prevalence_2022_plot$VG)]
prevalence_2022_plot <- rbind(prevalence_2022_plot[1:5,], empty_row, prevalence_2022_plot[6:7,], prevalence_2022_plot[9:12,], empty_row1, prevalence_2022_plot[13:nrow(prevalence_2022_plot),])
prevalence_2022_plot$wave <-  "2022"


# Including Multiple imputation statisitcs
prevalence_2022_plot$estimate[prevalence_2022_plot$VG == "Low income"] <- low_inc_statistics[,1]*100 # due to MI
prevalence_2022_plot$conf.low[prevalence_2022_plot$VG == "Low income"] <- low_inc_statistics[,3]*100 # due to MI
prevalence_2022_plot$conf.high[prevalence_2022_plot$VG == "Low income"] <- low_inc_statistics[,4]*100 # due to MI

prevalence_2022_plot$VG <- factor(prevalence_2022_plot$VG)


# same for risk ratios
risk_ratio_2017_plot <- select(risk_ratio_2017, c(VG, term, estimate, std.error, statistic, p.value, s.value, conf.low, conf.high, prevalence))
risk_ratio_2017_plot$VG <- display_labels[levels(risk_ratio_2017_plot$VG)]
risk_ratio_2017_plot <- rbind(risk_ratio_2017_plot[1:5,], empty_row, risk_ratio_2017_plot[6:7,], risk_ratio_2017_plot[9:12,], empty_row1, risk_ratio_2017_plot[13:nrow(risk_ratio_2017_plot),])
risk_ratio_2017_plot$wave <- "2017"

risk_ratio_2022_plot <- select(risk_ratio_2022, c(VG, term, estimate, std.error, statistic, p.value, s.value, conf.low, conf.high, prevalence))
risk_ratio_2022_plot$VG <- display_labels[levels(risk_ratio_2022_plot$VG)]
risk_ratio_2022_plot <- rbind(risk_ratio_2022_plot[1:5,], empty_row, risk_ratio_2022_plot[6:7,], risk_ratio_2022_plot[9:12,], empty_row1, risk_ratio_2022_plot[13:nrow(risk_ratio_2022_plot),])
risk_ratio_2022_plot$wave <-  "2022"


# Including Multiple imputation statisitcs
risk_ratio_2022_plot$estimate[risk_ratio_2022_plot$VG == "Low income"] <- low_inc_statistics[,5] # due to MI
risk_ratio_2022_plot$conf.low[risk_ratio_2022_plot$VG == "Low income"] <- low_inc_statistics[,6] # due to MI
risk_ratio_2022_plot$conf.high[risk_ratio_2022_plot$VG == "Low income"] <- low_inc_statistics[,7] # due to MI

# and finally dif-in-dif
Final_difndif_plot <- Final_difndif
Final_difndif_plot$VG <- display_labels[levels(Final_difndif_plot$VG)]
Final_difndif_plot <- rbind(Final_difndif_plot[1:5,], empty_row, Final_difndif_plot[6:7,], Final_difndif_plot[9:12,], empty_row1, Final_difndif_plot[13:nrow(Final_difndif_plot),])

### FOREST PLOT PREVALENCE COMBINED ###

combined_prevalence <- bind_rows(prevalence_2017_plot, prevalence_2022_plot)
combined_prevalence$VG <- factor(combined_prevalence$VG, levels = rev(unique(combined_prevalence$VG)))

plot_combined_prevalence <- combined_prevalence %>% 
  ggplot(aes(x = estimate, y = VG, color = wave)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.5) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.3,
    position = position_dodge(width = 0.5)
  ) +
  scale_y_discrete(drop = FALSE) +   # keeps NA rows as gaps
  labs(x = "Prevalence rate in % (95% confidence interval)", y = NULL, title = "Prevalences in 2017 and 2022",
       color = "Year") +
  scale_color_manual(values = c("2017" = "#1f77b4", "2022" = "#ff7f0e")) +
  theme_minimal() +
  scale_x_continuous(n.breaks = 7) +
  theme(
    plot.title = element_text(size = 18, face = "bold", color = "black"), 
    axis.text.y = element_text(hjust = 1, size = 12, color = "black"),
    plot.margin = margin(t = 20, r = 10, b = 20, l = 10, unit = "pt"),
    legend.position = "none"   # <- hides the legend
  )


# ggsave("plot_combined_prevalence.png", width = 9, height = 6, dpi = 300)  

### FOREST PLOT RISK RATIO COMBINED ###


combined_risk_ratio <- bind_rows(risk_ratio_2017_plot, risk_ratio_2022_plot)
combined_risk_ratio$VG <- factor(combined_risk_ratio$VG, levels = rev(unique(combined_risk_ratio$VG)))


plot_combined_risk_ratio <- combined_risk_ratio %>% 
  mutate(
    across(estimate,
           ~ ifelse(estimate == 1, NA, .x))
  ) %>%
  ggplot(aes(x = estimate, y = VG, color = wave)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.5) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.3,
    position = position_dodge(width = 0.5)
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_y_discrete(drop = FALSE) +   # keeps NA rows as gaps
  labs(x = "Risk ratio (95% confidence interval)", y = NULL, title = "Risk ratios in 2017 and 2022",
       color = "Year") +
  scale_color_manual(values = c("2017" = "#1f77b4", "2022" = "#ff7f0e")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold", color = "black"), 
    axis.text.y = element_text(hjust = 1, size = 12, color = "black"),
    plot.margin = margin(t = 20, r = 10, b = 20, l = 10, unit = "pt")
  )

plot_combined_risk_ratio
# ggsave("plot_combined_risk_ratio.png", width = 9, height = 6, dpi = 300)  

cc <- plot_combined_prevalence + plot_combined_risk_ratio
# ggsave("combined_graph_p+rd.png", width = 18, height = 6, dpi = 300)  


###### DIF-IN-DIF ####

Final_difndif_plot$VG <- factor(Final_difndif_plot$VG, levels = rev(unique(Final_difndif_plot$VG)))

forestp_difndif_forrest <- Final_difndif_plot %>% 
  filter(!VG %in% c("Low income", "Downward social mobility")) %>%
  ggplot(aes(y = VG, x = interaction)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Difference in differences (95% confidence interval)", y = NULL, title = "Time trends of SR") +
  theme_minimal() +
  theme(plot.title = element_text(size = 18, face = "bold", color = "black"), 
        axis.text.y = element_text(hjust = 1, size = 12, color = "black"),
        plot.margin = margin(t = 20, r = 40, b = 20, l = 40, unit = "pt")  # adjust these values
  )

forestp_difndif_forrest

# ggsave("forestp_difndif_forrest.png", width = 8, height = 6, dpi = 300)

