
# Load necessary libraries
library(lubridate)
library(dplyr)
library(reticulate)#Allows python calls within R
library(medExtractR)#Use the medExtractR wrapper
library(doParallel)
library(stringr)
library(readr)
library(R.utils)
library(data.table)
library(tidyr)

# Task 3.1 Data Extraction and Patient Identification

# Define TBI-related ICD code prefixes - please check
tbi_icd_prefixes <- c(
  # ICD-9 prefixes
  "800", "801", "803", "804",  # Skull fracture codes
  "850", "851", "852", "853", "854",  # Concussion/intracranial injury
  # ICD-10 prefixes
  "S02",  # Skull fractures
  "S06",  # Intracranial injuries
  "S07",  # Crushing injury of head
  "S09",  # Other head injuries
  "T90"   # Sequelae of head injuries
)
diagnoses_icd<-readRDS("diagnoses_icd.rds")
d_icd_diagnoses<-readRDS("d_icd_diagnoses.rds")
tbi_patients <- diagnoses_icd %>%
  left_join(d_icd_diagnoses, by = c("icd_code", "icd_version")) %>%
  filter(
    (icd_version == 9 & str_starts(icd_code, paste(tbi_icd_prefixes[1:9], collapse = "|"))) |
      (icd_version == 10 & str_starts(icd_code, paste(tbi_icd_prefixes[10:14], collapse = "|")))
  ) %>%
  distinct(subject_id, hadm_id)
rm(diagnoses_icd)

admissions<-readRDS("admissions.rds")

cohort <- tbi_patients %>%
  left_join(admissions, by = c("subject_id", "hadm_id")) %>%
  # Add any additional filters if needed (e.g., exclude newborns)
  filter(admission_type != "NEWBORN")
rm(admissions)

# Summary statistics
cat("Number of unique TBI patients:", n_distinct(cohort$subject_id), "\n")
#Number of unique TBI patients: 9758
cat("Number of TBI admissions:", nrow(cohort), "\n")
#Number of TBI admissions: 10730

# View first few results
glimpse(cohort)

# Save results
write_csv(cohort, "tbi_patients_mimiciv_simplified.csv")

########################################################################################################
# Extract information for all TBI patients from the MIMIC-IV database who received anticoagulation therapy during hospitalisation

meds_list<-c("enoxaparin", "heparin", "lovenox")

# Load and filter prescriptions
prescriptions <- read.csv("prescriptions.csv")
prescriptions <- prescriptions[prescriptions$subject_id %in% cohort$subject_id, ]
prescriptions$drug <- tolower(prescriptions$drug)

# Filter on medication keywords
presc_meds <- prescriptions %>%
  filter(str_detect(drug, regex(paste(meds_list, collapse = "|"), ignore_case = TRUE)))

# View distinct drug names in prescriptions
unique_presc_meds <- presc_meds %>%
  distinct(drug) %>%
  arrange(drug)

print(unique_presc_meds)

# Refine to valid anticoagulant therapies
presc_meds <- presc_meds %>%
  filter(
    str_detect(drug, regex("enoxaparin|lovenox|heparin( sodium)?", ignore_case = TRUE)) &
      !str_detect(drug, regex("flush|lock|topical|desensitization|crrt|impella|priming|lvad|dwell|via|clot directed|intraperitoneal|side port|hemodialysis", ignore_case = TRUE))
  )

# Frequency table of valid drugs
presc_meds %>%
  count(drug, sort = TRUE)

# Clean up
rm(prescriptions)

# presc_meds<-presc_meds %>%
# select("subject_id","hadm_id","starttime","drug")

########
# Get medications from EMAR table
emar <- read.csv("emar.csv")
emar <- emar[emar$subject_id %in% cohort$subject_id,]
emar$medication <- tolower(emar$medication)

emar_meds <- emar %>%
  filter(str_detect(medication, regex(paste(meds_list, collapse = "|"), ignore_case = TRUE)))

# View distinct medication names in EMAR
unique_emar_meds <- emar_meds %>%
  distinct(medication) %>%
  arrange(medication)

print(unique_emar_meds)

emar_meds <- emar_meds %>%
  filter(
    str_detect(medication, regex("enoxaparin|lovenox|heparin( sodium)?", ignore_case = TRUE)) &  # captures 'heparin sodium' too
      !str_detect(medication, regex("flush|lock|topical|desensitization|crrt|impella|priming|lvad|dwell|via|clot directed|intraperitoneal|side port|hemodialysis", ignore_case = TRUE))
  )

emar_meds %>%
  count(medication, sort = TRUE)


rm(emar)

# Standardize EMAR columns to match prescriptions format
if(!is.null(emar_meds) && nrow(emar_meds) > 0){
  emar_meds <- emar_meds %>%
    rename(drug = medication,
           starttime = charttime) %>%
    select(subject_id, hadm_id, starttime, drug)
  emar_meds$prod_strength <- NA
  emar_meds$dose_unit_rx <- NA
  emar_meds$doses_per_24_hrs <- NA
  emar_meds$route <- NA
}

saveRDS(presc_meds,"presc_meds.rds")
saveRDS(emar_meds,"emar_meds.rds")

presc_meds<-presc_meds %>%
  select("subject_id","hadm_id","starttime","drug","prod_strength","dose_unit_rx","doses_per_24_hrs", "route")
# Add EMAR medications if they exist
if(!is.null(emar_meds) && nrow(emar_meds) > 0){
  presc_meds <- rbind(presc_meds, emar_meds)
}
rm(emar_meds)
# presc_meds$starttime<-as.POSIXct(presc_meds$starttime)


#There are lots of synonyms for the same medication in the drugs field so combine the similar ones:
groups <- list()
i <- 1
similarGroups <- function(x, thresh = 0.8){
  grp <- integer(length(x))
  name <- x
  for(i in seq_along(name)){
    if(!is.na(name[i])){
      sim <- agrepl(x[i], x, ignore.case = TRUE, max.distance = 1 - thresh)
      k <- which(sim & !is.na(name))
      grp[k] <- i
      is.na(name) <- k
    }
  }
  grp
}
presc_meds$drug<-gsub("desensitization","",presc_meds$drug)
presc_meds$drug<-gsub("*","",presc_meds$drug)
presc_meds$drug<-str_replace(presc_meds$drug, " \\s*\\([^\\)]+\\)", "")
presc_meds$drug <- stringr::str_replace(presc_meds$drug, '\\*', '') #This is used twice as some have 2 lead asterisks
presc_meds$drug <- stringr::str_replace(presc_meds$drug, '\\*', '')
presc_meds$drug<-gsub("nf ","",presc_meds$drug)
presc_meds$drug <-gsub(".*:","",presc_meds$drug)

write.csv(presc_meds,"tbi_medications.csv")

## List of drugs to create binary columns for
drugs_to_flag <- meds_list

# Convert drug names to lowercase for case-insensitive matching
presc_meds$drug_lower <- tolower(presc_meds$drug)

# Create binary columns for each drug
for (drug in drugs_to_flag) {
  # Clean drug name for column name
  col_name <- gsub("[^a-zA-Z0-9]", "_", tolower(drug))
  col_name <- gsub("_+", "_", col_name)  # Replace multiple underscores with one
  col_name <- gsub("_$|^_", "", col_name)  # Remove leading/trailing underscores

  # Create binary column (1 if drug is found, 0 otherwise)
  presc_meds[[col_name]] <- as.integer(grepl(drug, presc_meds$drug_lower, ignore.case = TRUE))
}

# Remove temporary lowercase column
presc_meds$drug_lower <- NULL


# Aggregate to patient-admission level
patient_meds <- presc_meds %>%
  group_by(subject_id, hadm_id) %>%
  summarise(
    enoxaparin = as.integer(any(enoxaparin == 1, na.rm = TRUE)),
    lovenox = as.integer(any(lovenox == 1, na.rm = TRUE)),
    heparin = as.integer(any(heparin == 1, na.rm = TRUE)),
    .groups = "drop"
  )
patient_meds$Anticoagulant_Therapy<-1
patient_meds$VTEPROPHYLAXISTYPE<-NA
patient_meds$VTEPROPHYLAXISTYPE<-ifelse(patient_meds$enoxaparin==1,"Enoxaparin",NA)
patient_meds$VTEPROPHYLAXISTYPE<-ifelse(patient_meds$lovenox==1,"Lovenox",patient_meds$VTEPROPHYLAXISTYPE)
patient_meds$VTEPROPHYLAXISTYPE<-ifelse(patient_meds$heparin==1,"Heparin",patient_meds$VTEPROPHYLAXISTYPE)

# Merge
presc_meds <- left_join(presc_meds, patient_meds)
presc_meds$Anticoagulant_Therapy<-1
presc_meds$VTEPROPHYLAXISTYPE<-NA
presc_meds$VTEPROPHYLAXISTYPE<-ifelse(presc_meds$enoxaparin==1,"Enoxaparin",NA)
presc_meds$VTEPROPHYLAXISTYPE<-ifelse(presc_meds$lovenox==1,"Lovenox",presc_meds$VTEPROPHYLAXISTYPE)
presc_meds$VTEPROPHYLAXISTYPE<-ifelse(presc_meds$heparin==1,"Heparin",presc_meds$VTEPROPHYLAXISTYPE)


# Save the updated cohort
write.csv(presc_meds, "cohort_medications.csv", row.names = FALSE)
saveRDS(presc_meds, "cohort_medications.rds")

# Load admissions.csv
admissions <- read_csv("admissions.csv")
admissions<-admissions[admissions$hadm_id %in% cohort$hadm_id,]

# Convert admittime and dischtime to datetime (if not already parsed)
# Filter to stays of 1 month (<= 31 days)admissions <- admissions %>%
admissions_filtered <- admissions %>%
  mutate(
    admittime = as_datetime(admittime, tz = "UTC"),
    dischtime = as_datetime(dischtime, tz = "UTC"),
    length_of_stay_days = as.numeric(difftime(dischtime, admittime, units = "days"))
  ) %>%
  filter(length_of_stay_days <= 31)

# View result
print(head(admissions_filtered))

# Convert medication start time to datetime
presc_meds <- presc_meds %>%
  mutate(starttime = as_datetime(starttime, tz = "UTC"))

# Join medication data with admission info
presc_meds_with_admit <- presc_meds %>%
  inner_join(admissions_filtered, by = c("subject_id", "hadm_id"))

# Keep only medication records within first 7 days of admission
presc_meds_firstweek <- presc_meds_with_admit %>%
  filter(starttime >= admittime & starttime <= admittime + days(7))

# Get eligible subject/admission pairs
eligible_subjects <- presc_meds_firstweek %>%
  distinct(subject_id, hadm_id)

# View output
print(presc_meds_firstweek)
cat("Eligible records:", nrow(presc_meds_firstweek), "\n")
cat("Unique patients:", n_distinct(presc_meds_firstweek$subject_id), "\n")

saveRDS(presc_meds_firstweek, "tbi_pats.rds")


# 3.2 Thromboembolic Event Identification

tbi_pats<-readRDS("tbi_pats.rds")

# Load necessary tables (adjust filenames if needed)
diagnoses_icd <- readRDS("diagnoses_icd.rds")
diagnoses_icd_detail<-read.csv("d_icd_diagnoses.csv")
diagnoses_icd<-diagnoses_icd[diagnoses_icd$subject_id %in% tbi_pats$subject_id,]
diagnoses_icd<-dplyr::left_join(diagnoses_icd,diagnoses_icd_detail,by=c("icd_code","icd_version"))

# From ICD-based diagnoses (PE/DVT flags)
# Define regex patterns for PE and DVT ICD codes
# Patterns
pe_icd_pattern <- regex("^4151|^I26", ignore_case = TRUE)
dvt_icd_pattern <- regex("^453|^I82", ignore_case = TRUE)

# Regex-based patterns
pe_regex <- regex("\\bpulmonary embolism\\b", ignore_case = TRUE)
dvt_regex <- regex("deep vein thrombosis|dvt|venous thrombosis|iliac thrombosis|femoral thrombosis|popliteal|peroneal|venous embolism", ignore_case = TRUE)
exclude_phrases <- regex("history of|hx of|personal history", ignore_case = TRUE)


# Flags - ICD
pe_flag <- diagnoses_icd %>% filter(str_detect(icd_code, pe_icd_pattern)) %>% distinct(subject_id, hadm_id) %>% mutate(pe_flag = 1)
dvt_flag <- diagnoses_icd %>% filter(str_detect(icd_code, dvt_icd_pattern)) %>% distinct(subject_id, hadm_id) %>% mutate(dvt_flag = 1)
thrombo_flags <- full_join(pe_flag, dvt_flag, by = c("subject_id", "hadm_id")) %>%
  mutate(
    pe_flag = replace_na(pe_flag, 0),
    dvt_flag = replace_na(dvt_flag, 0),
    thrombo_flag = if_else(pe_flag == 1 | dvt_flag == 1, 1, 0)
  )

# Flags - Regex from long_title
pe_regex_flag <- diagnoses_icd %>% filter(str_detect(long_title, pe_regex) & !str_detect(long_title, exclude_phrases)) %>% distinct(subject_id, hadm_id) %>% mutate(pe_flag_text = 1)
dvt_regex_flag <- diagnoses_icd %>% filter(str_detect(long_title, dvt_regex) & !str_detect(long_title, exclude_phrases)) %>% distinct(subject_id, hadm_id) %>% mutate(dvt_flag_text = 1)
thrombo_regex_flags <- full_join(pe_regex_flag, dvt_regex_flag, by = c("subject_id", "hadm_id")) %>%
  mutate(
    pe_flag_text = replace_na(pe_flag_text, 0),
    dvt_flag_text = replace_na(dvt_flag_text, 0),
    thrombo_flag_text = if_else(pe_flag_text == 1 | dvt_flag_text == 1, 1, 0)
  )

# HCPCS events
hcpcs <- fread("hcpcsevents.csv.gz")
hcpcs_thrombo <- hcpcs %>%
  filter(hcpcs_cd %in% c("G8600", "G8599", "G9060")) %>%
  select(subject_id, hadm_id, thrombo_event_time = chartdate, hcpcs_cd) %>%
  distinct() %>%
  mutate(thrombo_flag_hcpcs = 1)

hcpcs_flags <- hcpcs_thrombo %>%
  distinct(subject_id, hadm_id) %>%
  mutate(thrombo_flag_hcpcs = 1)

#This draws a blank- no indications in hcpcs

# Merge all flags into TBI cohort
tbi_pats <- tbi_pats %>%
  left_join(thrombo_flags, by = c("subject_id", "hadm_id")) %>%
  left_join(thrombo_regex_flags, by = c("subject_id", "hadm_id")) %>%
  left_join(hcpcs_flags, by = c("subject_id", "hadm_id")) %>%
  mutate(across(c(pe_flag, dvt_flag, thrombo_flag, pe_flag_text, dvt_flag_text, thrombo_flag_text, thrombo_flag_hcpcs), ~ replace_na(., 0))) %>%
  mutate(
    thrombo_flag_combined = if_else(thrombo_flag == 1 | thrombo_flag_text == 1 | thrombo_flag_hcpcs == 1, 1, 0),
    thrombo_event_time = if_else(thrombo_flag_combined == 1, admittime, as.POSIXct(NA)),
    thrombo_event_source = case_when(
      thrombo_flag == 1 & thrombo_flag_text == 1 & thrombo_flag_hcpcs == 1 ~ "ICD+Text+HCPCS",
      thrombo_flag == 1 & thrombo_flag_text == 1 ~ "ICD+Text",
      thrombo_flag == 1 ~ "ICD",
      thrombo_flag_text == 1 ~ "Text",
      thrombo_flag_hcpcs == 1 ~ "HCPCS",
      TRUE ~ NA_character_
    )
  )

# Create overlap summary
tbi_pats <- tbi_pats %>% mutate(adm_id = paste(subject_id, hadm_id, sep = "_"))
pe_icd_ids <- tbi_pats %>% filter(pe_flag == 1) %>% pull(adm_id) %>% unique()
pe_text_ids <- tbi_pats %>% filter(pe_flag_text == 1) %>% pull(adm_id) %>% unique()
dvt_icd_ids <- tbi_pats %>% filter(dvt_flag == 1) %>% pull(adm_id) %>% unique()
dvt_text_ids <- tbi_pats %>% filter(dvt_flag_text == 1) %>% pull(adm_id) %>% unique()
thromb_icd_ids <- tbi_pats %>% filter(thrombo_flag == 1) %>% pull(adm_id) %>% unique()
thromb_text_ids <- tbi_pats %>% filter(thrombo_flag_text == 1) %>% pull(adm_id) %>% unique()

# Build summary table
overlap_summary <- tibble(
  category = c("PE", "DVT", "Thrombo"),
  only_icd = c(length(setdiff(pe_icd_ids, pe_text_ids)), length(setdiff(dvt_icd_ids, dvt_text_ids)), length(setdiff(thromb_icd_ids, thromb_text_ids))),
  only_text = c(length(setdiff(pe_text_ids, pe_icd_ids)), length(setdiff(dvt_text_ids, dvt_icd_ids)), length(setdiff(thromb_text_ids, thromb_icd_ids))),
  both = c(length(intersect(pe_icd_ids, pe_text_ids)), length(intersect(dvt_icd_ids, dvt_text_ids)), length(intersect(thromb_icd_ids, thromb_text_ids)))
)

print(overlap_summary)
# A tibble: 3 × 4
#category only_icd only_text  both
#<chr>       <int>     <int> <int>
#  1 PE              1         2    46
#2 DVT            55         7    55
#3 Thrombo        43         8    95

# Determine the Optimal Data Source

# diagnoses_icd	Standardized, high coverage	No event timestamp
# hcpcsevents	Includes event timing	Sparse (only some patients coded)
# noteevents	Best for NLP / regex

#  Use hcpcsevents if available for timing
# Use diagnoses_icd to maximize detection
# Optionally merge both into  tbi_pats

# VTE prophylaxis (anticoagulant therapy) initiation:

vte_prophylaxis_time <- tbi_pats %>%
  group_by(subject_id, hadm_id) %>%
  summarise(vte_start = min(starttime, na.rm = TRUE), .groups = "drop")

tbi_pats <- tbi_pats %>%
  left_join(vte_prophylaxis_time, by = c("subject_id", "hadm_id"))


# 3.3 Covariate Extraction

# Define code patterns for bleeding disorders

# Define regex for ICD-9 bleeding disorders
bleeding_icd9 <- regex("^286|^287", ignore_case = TRUE)

# Define regex for ICD-10 bleeding disorders
bleeding_icd10 <- regex("^D66|^D67|^D68|^D69", ignore_case = TRUE)

# Extract bleeding disorders from both versions
bleeding_disorders <- diagnoses_icd %>%
  filter(
    (icd_version == 9 & str_detect(icd_code, bleeding_icd9)) |
      (icd_version == 10 & str_detect(icd_code, bleeding_icd10))
  ) %>%
  distinct(subject_id, hadm_id) %>%
  mutate(bleeding_disorder_flag = 1)

tbi_pats <- tbi_pats %>%
  left_join(bleeding_disorders, by = c("subject_id", "hadm_id"))


# Pre-hospital Anticoagulation Use
anticoag_list <- c("enoxaparin", "heparin", "lovenox")
pre_hosp_anticoag <- presc_meds %>%
  filter(subject_id %in% tbi_pats$subject_id) %>%
  mutate(starttime = as_datetime(starttime)) %>%
  filter(str_detect(drug, regex(paste(anticoag_list, collapse = "|"), ignore_case = TRUE))) %>%
  inner_join(tbi_pats %>% select(subject_id, hadm_id, admittime) %>% distinct(),
             by = c("subject_id", "hadm_id")) %>%
  filter(starttime < admittime) %>%
  distinct(subject_id, hadm_id) %>%
  mutate(pre_hosp_anticoag_flag = 1)

tbi_pats <- tbi_pats %>%
  left_join(pre_hosp_anticoag, by = c("subject_id", "hadm_id"))

# Count how many patients in tbi_pats had pre-hospital anticoagulation use
tbi_pats %>%
  filter(pre_hosp_anticoag_flag == 1) %>%
  summarise(
    n_records = n(),
    n_unique_patients = n_distinct(subject_id),
    n_unique_admissions = n_distinct(hadm_id)
  )
#  n_records n_unique_patients n_unique_admissions
#1       560                49                  49


# Select just the gender and anchor_age from patients
patients<-read.csv("patients.csv")
patients<-patients[patients$subject_id %in% cohort$subject_id,]  
demographics <- patients %>%
  select(subject_id, gender, anchor_age,dod) %>%
  rename(age = anchor_age)

# Merge into tbi_pats directly
tbi_pats <- tbi_pats %>%
  left_join(demographics, by = "subject_id")

#  Extract Event Resolution Time

# Event resolution was defined as the earliest of hospital discharge (dischtime) or date of death (dod),
# derived from admissions.csv and patients.csv.
# This variable approximates the end of thromboembolic event follow-up during hospitalization.

# Join admissions with patient death data
event_resolution <- admissions %>%
  select(subject_id, hadm_id, dischtime) %>%
  left_join(patients %>% select(subject_id, dod), by = "subject_id") %>%
  left_join(admissions) %>%
  mutate(
    dischtime = as_datetime(dischtime),
    dod = as_datetime(dod),

    # Event resolution time = earliest of discharge or death
    resolution_time = pmin(dischtime, dod, na.rm = TRUE),

    # Label resolution type
    resolution_type = case_when(
      !is.na(dod) & dod <= dischtime ~ "Death",
      TRUE ~ "Discharge"
    )
  ) %>%
  select(subject_id, hadm_id, resolution_time, resolution_type)

tbi_pats <- tbi_pats %>%
  left_join(event_resolution, by = c("subject_id", "hadm_id"))


# 3.4 Data Quality Assessment
# Summarise Missingness and Completeness

# Calculate missingness % per variable
missing_summary <- tbi_pats %>%
  summarise(across(everything(), ~ mean(is.na(.))*100)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "percent_missing") %>%
  arrange(desc(percent_missing))

# View top missing
print(head(missing_summary, 20))

# A tibble: 20 × 2
#variable               percent_missing
#<chr>                            <dbl>
#  1 pre_hosp_anticoag_flag           97.3 
#2 deathtime                        96.4 
#3 thrombo_event_time               93.5 
#4 thrombo_event_source             93.5 
#5 bleeding_disorder_flag           88.0 
#6 doses_per_24_hrs                 75.4 
#7 prod_strength                    72.8 
#8 dose_unit_rx                     72.8 
#9 route                            72.8 
#10 marital_status                   13.3 
#11 edregtime                         5.78
#12 edouttime                         5.78
#13 discharge_location                3.95
#14 subject_id                        0   
#15 hadm_id                           0   
#16 starttime                         0   
#17 drug                              0   
#18 enoxaparin                        0   
#19 heparin                           0   
#20 lovenox                           0 

# Flag Critical Fields for Completeness


# Check if any IDs missing core values
tbi_pats %>%
  filter(is.na(subject_id)| is.na(hadm_id) |is.na(admittime) | is.na(starttime)) %>%
  count()
#0

# Clean Flag Variables
tbi_pats <- tbi_pats %>%
  mutate(
    bleeding_disorder_flag = replace_na(bleeding_disorder_flag, 0),
    pre_hosp_anticoag_flag = replace_na(pre_hosp_anticoag_flag, 0),
    pe_flag = replace_na(pe_flag, 0),
    dvt_flag = replace_na(dvt_flag, 0),
    thrombo_flag = replace_na(thrombo_flag, 0),
    thrombo_flag_icd = replace_na(thrombo_flag, 0)
  )

# Check for Temporal Consistency

# Starttime must be after admittime
tbi_pats %>%
  filter(starttime < admittime) %>%
  count()

# VTE start must be during admission
tbi_pats %>%
  filter(vte_start > dischtime) %>%
  count()

# NOTE: A total of 37 records had anticoagulant therapy start times (vte_start) after the recorded hospital discharge (dischtime).
# These records likely reflect data entry errors, post-discharge prescriptions, or misaligned admission-medication links.
# They may be excluded from the primary analysis to ensure temporal consistency in evaluating in-hospital VTE prophylaxis initiation.

# Flag them for review
tbi_pats <- tbi_pats %>%
  mutate(vte_after_discharge_flag = if_else(vte_start > dischtime, 1, 0))

# Exclude them entirely
  tbi_pats <- tbi_pats %>%
    filter(is.na(vte_start) | vte_start <= dischtime)

saveRDS(tbi_pats, "tbi_pats_cleaned.rds")
write.csv(tbi_pats, "tbi_pats_cleaned.csv", row.names = FALSE)


# EXTRA CODE

# VTE therapy start before admission
tbi_pats %>%
  filter(!is.na(vte_start) & vte_start < admittime) %>%
  count()

# starttime (med administration) before admission
tbi_pats %>%
  filter(!is.na(starttime) & starttime < admittime) %>%
  count()

# starttime after discharge

tbi_pats %>%
  filter(!is.na(starttime) & starttime > dischtime) %>%
  count()

## NOTE: A total of 35 records had medication administration times (starttime) occurring after the recorded hospital discharge (dischtime).
# These may represent discharge prescriptions, outpatient events, or data inconsistencies.
# Consider excluding them from in-hospital analyses to maintain accurate temporal alignment with the admission period.

# Summary statistics of extracted TBI cohort

tbi_pats <- readRDS("tbi_pats_cleaned.rds")

# Basic patient counts
n_patients <- n_distinct(tbi_pats$subject_id)
n_admissions <- n_distinct(tbi_pats$hadm_id)
n_med_admins <- nrow(tbi_pats)

cat("Number of unique patients:", n_patients, "\n")
cat("Number of unique admissions:", n_admissions, "\n")
cat("Total medication administration records:", n_med_admins, "\n\n")

# Gender distribution
tbi_pats %>%
  count(gender) %>%
  mutate(percent = round(n / sum(n) * 100, 1)) %>%
  print()

# Age summary
tbi_pats %>%
  summarise(
    mean_age = round(mean(age, na.rm = TRUE), 1),
    median_age = median(age, na.rm = TRUE),
    min_age = min(age, na.rm = TRUE),
    max_age = max(age, na.rm = TRUE)
  ) %>%
  print()

# Length of stay summary
tbi_pats %>%
  summarise(
    mean_los = round(mean(length_of_stay_days, na.rm = TRUE), 2),
    median_los = median(length_of_stay_days, na.rm = TRUE),
    min_los = min(length_of_stay_days, na.rm = TRUE),
    max_los = max(length_of_stay_days, na.rm = TRUE)
  ) %>%
  print()

# Anticoagulant breakdown
cols_anticoag <- c("enoxaparin", "heparin", "lovenox")
tbi_pats %>%
  summarise(across(all_of(cols_anticoag), ~ sum(. == 1, na.rm = TRUE))) %>%
  print()

# Pre-hospital anticoagulation use
tbi_pats %>%
  count(pre_hosp_anticoag_flag) %>%
  mutate(label = ifelse(pre_hosp_anticoag_flag == 1, "Yes", "No")) %>%
  select(label, n) %>%
  print()

# Thromboembolic event detection by source
tbi_pats %>%
  count(thrombo_event_source) %>%
  mutate(percent = round(n / sum(n) * 100, 1)) %>%
  print()

# Summary complete
cat("\nSummary statistics generated successfully.\n")


