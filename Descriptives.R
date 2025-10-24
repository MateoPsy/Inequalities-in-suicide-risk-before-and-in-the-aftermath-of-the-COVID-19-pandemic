##########################
#### COVID-19 VG 2025 ####
####   descriptives   ####
##########################

setwd("Desktop/Data")
###### Loading pkg ##### 
library(tidyverse)  # data manipulation

# loading data set
merged_data <- read_csv("merged_data.csv")
merged_data_diag <- read_csv("merged_data_diag.csv")

###### DESCRIPTIVES #####

## population characteristics for Table 1.

descriptives_by_method_props <- merged_data %>% 
  dplyr::select(collect_wave, collect_method, sex, education, marital_status, work_status, size_place_residence) %>% 
  pivot_longer(-c(collect_wave, collect_method)) %>% 
  count(collect_wave, collect_method, name, value) %>%
  group_by(collect_wave, collect_method, name) %>%
  mutate(prop = round(n/sum(n)*100, 2),
         stat = paste(n, paste0("(", prop, ")"))) %>% 
  dplyr::select(collect_wave, collect_method, name, value, stat) %>% 
  pivot_wider(names_from = c(collect_wave, collect_method), values_from = stat)


descriptives_by_method_age <- merged_data %>% 
  dplyr::select(collect_wave, collect_method, age) %>% 
  group_by(collect_wave, collect_method) %>% 
  summarise(age = paste0(round(mean(age), 2)," (", round(sd(age), 2), ")")) %>% 
  pivot_wider(names_from = c(collect_wave, collect_method), values_from = age) %>% 
  mutate(name = "age", value = "age")

descriptives_by_method <- bind_rows(descriptives_by_method_props, descriptives_by_method_age) %>% 
  mutate(name = factor(as.character(name),
                       levels = c("age", "sex", "education", "marital_status", "work_status", "size_place_residence"),
                       ordered = TRUE)) %>% 
  arrange(name) %>% 
  dplyr::select(name, value, starts_with("2017"), starts_with("2022"))



### Vulnerable groups descriptives for Table 2.

custom_order <- c("Young_females", "Young_males", "Older_females", "Older_males", "Primary_ed", 
                  "Unemployed", "Economically_inactive", "Low_income", "Not_relationship", "Single_parenthood",
                  "Social_mobility", "major_depressive_episode", "anxiety_disorders", "alcohol_use_disorders",
                  "comorbidity", "any_mental_disorder")

# Summarizing data by year
summary_table <- merged_data_diag %>%
  dplyr::select(collect_wave, Young_females, Young_males,
                Older_females, Older_males, Primary_ed, Unemployed, Economically_inactive, 
                Low_income, Not_relationship, Single_parenthood, Social_mobility, major_depressive_episode, 
                anxiety_disorders, alcohol_use_disorders, comorbidity, any_mental_disorder) %>% 
  pivot_longer(-collect_wave, names_to = "Vulnerable_group", values_to = "Is_vulnerable") %>%  # Convert wide to long format
  group_by(collect_wave, Vulnerable_group) %>%
  summarise(
    Cases = sum(Is_vulnerable),                      # Count cases (sum of 1s)
    Percentage = (sum(Is_vulnerable) / n()) * 100    # Calculate percentage
  ) %>%
  ungroup() %>%
  mutate(Percentage = sprintf("%.2f%%", Percentage)) %>%  # Format percentage
  mutate(Year = paste0(Cases, " (", Percentage, ")")) %>%  # Combine cases and percentage
  select(-Percentage) %>%  # Remove separate percentage column
  pivot_wider(names_from = collect_wave, values_from = c(Cases, Year)) %>%  # Convert back to wide format
  mutate(Vulnerable_group = factor(Vulnerable_group, levels = custom_order)) %>%  # Apply custom order
  arrange(Vulnerable_group)

summary_table <- summary_table[c(1,4:5)]

summary_table
### overall - by diagnosis
prevalence_disorders_overall <- merged_data_diag %>% 
  group_by(collect_wave) %>% 
  nest() %>% 
  mutate(survey_design = map(.x = data,
                             ~ survey::svydesign(ids = ~1, # no clusters
                                                 data = .x, 
                                                 weights = .x$weight)), # survey weights - each observation weighted by the inverse of its sampling probability
         diagnoses = map(.x = data,
                         ~ .x %>% dplyr::select(alcohol_dependence, alcohol_abuse, 
                                                major_depressive_episode, 
                                                suicide_risk, 
                                                panic_disorder, agoraphobia, social_phobia, PTSD, GAD,
                                                anxiety_disorders, alcohol_use_disorders, comorbidity, any_mental_disorder)),
         counts = map(.x = diagnoses,
                      ~ .x %>% mutate(across(everything(), ~ sum(.x == 1))) %>% distinct()), # no subgroup with ≤ 5 individuals
         prevalence_mean = map2(.x = diagnoses,
                                .y = survey_design,
                                ~ survey::svymean(~ as.matrix(.x),
                                                  design = .y)),
         prevalence_CI = map(.x = prevalence_mean,
                             ~ confint(.x, 
                                       type = "delta.method", # vahy?
                                       level = 0.95)),
         prevalence = map2(.x = prevalence_mean,
                           .y = prevalence_CI,
                           ~ as.data.frame(cbind(.x, .y)) %>% 
                             rownames_to_column(.) %>% 
                             transmute(diagnosis = str_sub(as.character(rowname), start =  14),
                                       mean =  .x,
                                       lwr_CI = `2.5 %`,
                                       upr_CI = `97.5 %`) %>% 
                             mutate(across(where(is.numeric),
                                           ~ round((.x*100),2)),
                                    estimate = paste0(mean, " (", lwr_CI, ", ", upr_CI, ")")))) 


results_prevalence_disorders_overall <- prevalence_disorders_overall %>% 
  dplyr::select(prevalence) %>% 
  unnest(cols = c(prevalence)) 


table_prevalence_disorders_overall <- results_prevalence_disorders_overall %>% 
  dplyr::select(-c(mean, lwr_CI, upr_CI)) %>% 
  pivot_wider(names_from = "collect_wave",
              values_from = "estimate") %>% 
  mutate(diagnosis = factor(as.character(diagnosis),
                            levels = c("any_mental_disorder", "alcohol_use_disorders", "alcohol_dependence", "alcohol_abuse","major_depressive_episode",
                                       "anxiety_disorders", "social_phobia", "agoraphobia", "panic_disorder", "GAD", "PTSD", "suicide_risk","comorbidity"),
                            ordered = TRUE)) %>% 
  arrange(diagnosis)
