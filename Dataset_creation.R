##########################
#### COVID-19 VG 2025 ####
#### dataset creation ####
##########################


###### Loading pkg ##### 
library(tidyverse)  # data manipulation
library(readxl) # Excel data

###### DATA IMPORT #####
### 2017 ###
data_2017 <- read_delim(file = "Data_2017_raw.csv", delim = ";") %>% 
  mutate(dot = dot.č.)

weights_2017 <- read_delim(file = "Weights_2017.csv", 
                           delim = ";", 
                           locale = locale(decimal_mark = ","))[,c("dot", "weight")]
data_2017_eu <- read_csv("Data_2017_per_european_psychiatry.csv") %>% 
  mutate(dot = ID)

data_2017w <- as_tibble(merge(data_2017, weights_2017, by = "dot"))
data_2017w <- as_tibble(merge(data_2017w, data_2017_eu[,c("dot", "POP")], by = "dot")) %>% 
  rename_with(~ str_replace(.x, "E4([a-z])$", "E4_\\1"), starts_with("E4")) %>%
  rename_with(~ str_replace(.x, "E4_([a-z])$", function(x) {
    letter <- str_sub(x, -1)  # Extract last letter
    number <- match(letter, letters)  # Convert letter to number
    str_replace(x, "[a-z]$", as.character(number))
  }), starts_with("E4")) %>%
  rename_with(~ str_replace(.x, "E1([a-z])", function(x) {
    paste0("E1_", toupper(substr(x, 3, 3))) 
  }), starts_with("E")) %>%
  rename_with(~ str_replace(.x, "A4([a-z])", function(x) {
    paste0("A4_", toupper(substr(x, 3, 3))) 
  }), starts_with("A")) %>%
  rename_with(~ str_replace(.x, "O1([a-z])", function(x) {
    paste0("O1_", toupper(substr(x, 3, 3))) 
  }), starts_with("O")) %>%
  rename_with(~ str_replace(.x, "E2", "E2")) %>%
  mutate(collect_wave = "2017",
         collect_method = "CAPI/PAPI",
         id = dot,
         sex = factor(case_when(Kód_cílové_osoby == 1 ~ Pohlaví_1,
                                Kód_cílové_osoby == 2 ~ Pohlaví_2,
                                Kód_cílové_osoby == 3 ~ Pohlaví_3,
                                Kód_cílové_osoby == 4 ~ Pohlaví_4,
                                Kód_cílové_osoby == 5 ~ Pohlaví_5),
                      levels = c("1","2"), 
                      labels = c("male", "female")),
         age = case_when(Kód_cílové_osoby == 1 ~ Věk_1,
                         Kód_cílové_osoby == 2 ~ Věk_2,
                         Kód_cílové_osoby == 3 ~ Věk_3,
                         Kód_cílové_osoby == 4 ~ Věk_4,
                         Kód_cílové_osoby == 5 ~ Věk_5),
         education = factor(as.character(Ot1),
                            levels = c("1","2","3","4","5"),
                            labels = c("primary", "lower_secondary", "upper_secondary", "universitary", "universitary")),
         work_status = factor(case_when(Ot5_9 == 1 & Ot6 == 1 ~ "self-employed",
                                        Ot5_9 == 1 & Ot6 == 2 ~ "employed",
                                        Ot5_7 == 1 ~ "student",
                                        Ot5_1 == 1 | Ot5_8 == 1 ~ "unemployed", # both intentional and unwanted
                                        Ot5_6 == 1 ~ "retired", 
                                        Ot5_3 == 1 ~ "parental_leave",
                                        TRUE ~ "other")), # in household, sickness leave, disability pension, both employed and self-employed and other
         marital_status = factor(case_when(Ot2 == 1 | Ot2 == 2 | Ot2 == 5 ~ "married|cohabitation", # also married but living apart
                                           Ot2 == 3 ~ "widowed",
                                           Ot2 == 4 ~ "divorced",
                                           Ot2 == 6 ~ "single",
                                           TRUE ~ "other")),
         size_place_residence = factor(case_when(POP %in% 1:2 ~ 1,    
                                                 POP == 3 ~ 2,     
                                                 POP == 4 ~ 3,    
                                                 POP == 5 ~ 4),
                                       levels = c("1","2","3","4"),
                                       labels = c("<5000", "5000-19999", "20000-99999","100000+")),
         psychiatrist = Ot19b_1,
         psychologist = Ot19b_2,
         general_practinioner = Ot19b_3)

## Standardizing names
j2_columns <- grep("^J2", names(data_2017w), value = TRUE)
j3_columns <- grep("^J3", names(data_2017w), value = TRUE)
names(data_2017w)[names(data_2017w) %in% j2_columns] <- paste0("J2_", seq_along(j2_columns))
names(data_2017w)[names(data_2017w) %in% j3_columns] <- paste0("J3_", seq_along(j3_columns))
other_columns <- grep("^A3|^E1|^E4|^O1|^O3|^I4|^I5", names(data_2017w), value = TRUE)
names(data_2017w)[names(data_2017w) %in% other_columns] <- sapply(names(data_2017w)[names(data_2017w) %in% other_columns], function(name) {
  name <- sub("([abcdefghi])", "_\\1", name)
  name <- chartr("abcdefghi", "123456789", name)
  return(name)
}) 


### operationalization of Vulnerable Groups
data_2017w <- data_2017w %>%
  mutate(Young_females = as.numeric(sex == "female" & age <= 29),
         Young_males = as.numeric(sex == "male" & age <= 29),
         Older_females = as.numeric(sex == "female" & age >= 65),
         Older_males = as.numeric(sex == "male" & age >= 65),
         Primary_ed = as.numeric(education == "primary"), # undereducated as primary or sec w/out diploma
         Unemployed = as.numeric(work_status == "unemployed"),
         Economically_inactive = as.numeric(work_status == "unemployed" | work_status == "student" |
                                              work_status == "parental_leave" | work_status == "other"| 
                                              work_status == "retired"),
         Low_income = 0,
         Not_relationship = as.numeric(marital_status == "widowed" | marital_status == "divorced" |marital_status == "single"),
         Single_parenthood = as.numeric(Ot3a == "1" & Not_relationship == "1" & Ot4 >= 1 & age <= 60),
         Social_mobility = 0,
  )
data_2017w$weight <- as.numeric(data_2017w$weight)

### 2022 ###

data_2022 <- read.csv("AZV_2022_CAPI.csv")
data_2022 <- rename_with(data_2022, .fn = ~ str_replace(.x, "Q5_", "SELFI_"),
                         .cols = starts_with("Q5_")) %>% 
  rename_with(~ ifelse(str_detect(.x, "^Q\\d+[A-O]\\d+[A-Z]?(_\\d+)?$"), 
                       str_replace(.x, "^Q\\d+", ""), .x), 
              starts_with("Q")) %>%
  rename_with(~ str_replace(.x, "O1([A-Z])$", "O1_\\1"), starts_with("O1")) %>%
  rename_with(~ str_replace(.x, "E1([A-Z])$", "E1_\\1"), starts_with("E1")) %>%
  rename_with(~ str_replace(.x, "A4([A-Z])$", "A4_\\1"), starts_with("A4")) %>%
  mutate(collect_wave = "2022",
         collect_method = "CAPI/PAPI",
         id = ID,
         weight = VAHA,
         sex = factor(as.character(SEX_B),
                      levels = c("1","2"), 
                      labels = c("male", "female")),
         gender = factor(as.character(SEX),
                         levels = c("1","2", "3", "4"),
                         labels = c("men", "women", "transgender", "non-binary")),
         sexual_orientation = factor(as.character(S_PREF),
                                     levels = c("1","2", "3", "4", "5"),
                                     labels = c("heterosexual", "lesbian", "gay", "bisexual", "other")),
         region_residence = as.factor(REG),
         size_place_residence = factor(as.character(POP),
                                       levels = c("1","2","3","4"),
                                       labels = c("<5000", "5000-19999", "20000-99999","100000+")),
         income = factor(as.character(INC),
                         levels = c("1","2","3","4","5","6","7"),
                         labels = c("0-9k","10-19k","20-29k","30-39k","40-49k","50+k","no reply")),
         age = AGE,
         education = factor(as.character(EDU),
                            levels = c("1","2","3","4"),
                            labels = c("primary", "lower_secondary", "upper_secondary", "universitary")),
         work_status_orig = factor(case_when(EMP == 1 ~ "employed",
                                             EMP == 2 ~ "self-employed",
                                             EMP == 3 ~ "in household",
                                             EMP == 5 ~ "disability pension",
                                             EMP == 6 ~ "sickness leave",
                                             EMP == 7 ~ "student",
                                             EMP == 9 ~ "voluntarily unemployed",  
                                             EMP == 10 ~ "involuntarily unemployed", 
                                             EMP == 8 ~ "retired", 
                                             EMP == 4 ~ "parental_leave")), 
         work_status = factor(case_when(EMP == 1 ~ "employed",
                                        EMP == 2 ~ "self-employed",
                                        EMP == 7 ~ "student",
                                        EMP == 9 | EMP == 10 ~ "unemployed",  # both intentional and unwanted
                                        EMP == 8 ~ "retired", # "not-working retiree"
                                        EMP == 4 ~ "parental_leave",
                                        TRUE ~ "other")), # sickness leave, disability pension, in household
         marital_status = factor(case_when(REL_S == 3 ~ "other", # in relationship, living apart
                                           REL_S == 1 | REL_S == 2 | REL_S == 6 ~ "married|cohabitation", # also married but living apart
                                           REL_S == 4 ~ "widowed",
                                           REL_S == 7 ~ "single",
                                           REL_S == 5 ~ "divorced")),
         marital_status_orig = factor(case_when(REL_S == 3 ~ "in relationship, living apart", # in relationship, living apart
                                                REL_S == 1 ~ "married", 
                                                REL_S == 2 ~ "living with a partner",
                                                REL_S == 6 ~ "married, living apart",
                                                REL_S == 4 ~ "widowed",
                                                REL_S == 7 ~ "single",
                                                REL_S == 5 ~ "divorced")),
         across(.cols = c("SELFI_1", "SELFI_3"),
                ~ 6 - .),
         SELFI_composite = SELFI_1 + SELFI_2 + SELFI_3 + SELFI_4 + SELFI_5, # higher scores indicate higher degree of self-identification with having a mental disorder
         across(.cols = starts_with("Q10_"),
                ~ ifelse(.x %in% c(1, 2, 3, 4, 5, 6), 1, .x)),
         across(.cols = starts_with("Q14_"),
                ~ ifelse(.x %in% c(1, 2, 3, 4, 5, 6), 1, 0)),
         across(.cols = starts_with("Q15_"),
                ~ ifelse(.x %in% c(1, 2, 3, 4), 1, 0))
  ) %>% 
  mutate(structural_barriers = as.numeric(rowSums(dplyr::select(., Q14_1:Q14_5) == 1) >= 1),
         structural_barriers_WO = Q15_2,
         help_seeking_12 = Q9,)



### operationalization of Vulnerable Groups
data_2022 <- data_2022 %>%
  mutate(Young_females = as.numeric(sex == "female" & age <= 29),
         Young_males = as.numeric(sex == "male" & age <= 29),
         Older_females = as.numeric(sex == "female" & age >= 65),
         Older_males = as.numeric(sex == "male" & age >= 65),
         Primary_ed = as.numeric(EDU == 1), # undereducated as primary or sec w/out diploma
         Unemployed = as.numeric(work_status == "unemployed"),
         Economically_inactive = as.numeric(work_status == "unemployed" | work_status == "student" |
                                              work_status == "parental_leave" | work_status == "other"| 
                                              work_status == "retired"),
         Low_income = as.numeric(income == "0-9k"),
         Not_relationship = as.numeric(marital_status == "widowed" | marital_status == "divorced" |marital_status == "single"),
         Single_parenthood = as.numeric(Not_relationship == "1" & CHIL >= 1 & age <= 60),
         Social_mobility = as.numeric(Q78_1 == 1 | Q78_2 == 1 | Q78_3 == 1 | Q78_4 == 1 | INC_LOSS == 2)
  )



names(data_2022$VAHA) <- "weight"


###### MERGING DATASETS ####
merged_data <- bind_rows(data_2017w, data_2022) %>% 
  dplyr::select(id, weight, collect_wave, collect_method, sex, gender, age, education, work_status, work_status_orig, income, marital_status, marital_status_orig, region_residence, size_place_residence, SELFI_composite,
                c(starts_with(c("A1", "A2", "A3", "A4",
                                "E1_", "E2", "E3", "E7",
                                "F1", "F2", 
                                "G1", "G2", "G3", "G4", 
                                "O1_", "O2", "O3", 
                                "I1", "I2", "I3", "I4", "I5", "I6",
                                "J1", "J2", "J3",
                                "C")),
                  E4_1:E4_9, E4_10:E4_13), Q29SP, weight,  # Q30SP nedavam, nakolko to vyplnili len ludia, kt. uviedli 1 pri Q29SP
                # psychiatrist, psychologist, general_practinioner, crisis_intervention, online_therapy, 
                Young_females, Young_males,
                Older_females, Older_males, Primary_ed, Unemployed, Economically_inactive, Low_income, Not_relationship, Single_parenthood, Social_mobility, structural_barriers, structural_barriers_WO, help_seeking_12,
                ends_with(c("_quantity", "_frequency")))


# adding diagnoses based on M.I.N.I. assessment
merged_data_diag <- merged_data %>% 
  mutate(across(.cols = c(starts_with(c("A1", "A2", "A3","E1_", "E2", "E3", "E7","F1", "F2","G1", "G2", "G3", "G4","O1_", "O2", "O3","I1", "I2", "I3", "I4", "I5", "I6","J1", "J2", "J3","C", "Single_parenthood")),
                          E4_1:E4_9, E4_10:E4_13),
                ~ replace(., is.na(.), 0))) %>% 
  mutate(major_depressive_episode = as.numeric((A1 == 2 | A2 == 2) & rowSums(dplyr::select(., A1:A3_7) == 2) >= 5),
         panic_disorder = as.numeric(E1_A == 2 & E1_B == 2 & E2 == 2 & E3 == 2 & rowSums(dplyr::select(., E4_1:E4_13) == 2) >= 4 & E7 == 2),
         agoraphobia = as.numeric(F1 == 2 & F2 == 2),
         social_phobia = as.numeric(rowSums(dplyr::select(., G1:G4) == 2) == 4),
         GAD = as.numeric(O1_A == 2 & O1_B == 2 & O2 == 2 & rowSums(dplyr::select(., O3_1:O3_6) == 2) >= 3),
         PTSD = as.numeric(I1 == 2 & I2 == 2 & I3 == 2 & rowSums(dplyr::select(., starts_with("I4")) == 2) >= 3 & rowSums(dplyr::select(., starts_with("I5")) == 2) >= 2 & I6 == 2),
         Jx = as.numeric(J3_1 == 2 & J3_2 == 2),
         Jy = as.numeric(J2_2 == 2 | J2_3 == 2),
         alcohol_dependence = as.numeric(J1 == 2 & rowSums(dplyr::select(., starts_with("J2")) == 2) >= 3),
         alcohol_abuse = as.numeric(J1 == 2 & rowSums(dplyr::select(., starts_with("J3")) == 2) >= 1),
         #### scoring acc. to MINI - low, medium, or high risk
         C1_rec = (C1 == 2) * 1,
         C2_rec = (C2 == 2) * 2,
         C3_rec = (C3 == 2) * 6,
         C4_rec = (C4 == 2) * 10,
         C5_rec = (C5 == 2) * 10,
         C6_rec = (C6 == 2) * 4,
         suicide_risk = as.numeric(rowSums(dplyr::select(., C1:C6) == 2) >= 1),
         suicide_thoughts = as.numeric(rowSums(dplyr::select(., C1:C4) == 2) >= 1),
         suicide_behavior = as.numeric(rowSums(dplyr::select(., C5:C6) == 2) >= 1),
         suicide_risk_without_lifelong = as.numeric(rowSums(select(., C1:C5) == 2) >= 1),
         Nonsuic_self_i = as.numeric(Q29SP == 1),) %>% 
  mutate(level_suicide_risk = case_when(rowSums(dplyr::select(., C1_rec:C6_rec)) %in% 1:5 ~ 1, # low risk
                                        rowSums(dplyr::select(., C1_rec:C6_rec)) %in% 6:9 ~ 2, # medium risk
                                        rowSums(dplyr::select(., C1_rec:C6_rec)) >= 10 ~ 3, # high risk
                                        TRUE ~ 0),
         anxiety_disorders = as.numeric(rowSums(dplyr::select(.,panic_disorder, agoraphobia, social_phobia, PTSD, GAD) == 1) >= 1), 
         alcohol_use_disorders = as.numeric(rowSums(dplyr::select(., alcohol_dependence, alcohol_abuse) == 1) >= 1),
         any_mental_disorder = as.numeric(major_depressive_episode| anxiety_disorders | alcohol_use_disorders)) %>% 
  mutate(MD_without_alcohol = as.numeric(rowSums(dplyr::select(.,panic_disorder, agoraphobia, social_phobia, PTSD, GAD,
                                                               major_depressive_episode) == 1) >= 1),
         comorbidity = as.numeric(case_when(MD_without_alcohol == 1 & alcohol_use_disorders == 1 ~ 1,
                                            TRUE ~ 0))) 

write_csv(merged_data, "merged_data.csv")
write_csv(merged_data_diag, "merged_data_diag.csv")
write_csv(data_2022, "data_2022.csv")

