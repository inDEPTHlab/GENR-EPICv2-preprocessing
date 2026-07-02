# Fixing randomization fuck-up 

setwd("/home/s.defina/picklist")

library(foreign)
library(dplyr)

readit <- function(fname, dir = "./files") {
  out <- foreign::read.spss(file.path(dir, paste0(fname, ".sav")), 
                            to.data.frame = TRUE, use.value.labels = TRUE)
  
  # print(summary(out))
  out
}

# list.files("./files")

# Read file prepared by Claudia (27/01/2025): this should have the ID corrected 
# concentration information and sex check (for some samples)
final_selection <- readit("20250127_Serena_SamplestoRandomize") 

# Clean extra samples -----------------------------------------------------------------------

# Exclude low concentration (based on Claudia's variable "exclude")
final_selection <- final_selection[is.na(final_selection$exclude) | final_selection$exclude == 0, ]

# table(final_selection$Period)

# Clean bridge samples ----------------------------------------------------------------------

# check_bridges <- final_selection[final_selection$Period=="Birth",]
# table(table(check_bridges$IDC))

# Remove 1 duplicate bridge sample (processing error on data management side)
final_selection <- final_selection[!duplicated(final_selection$SampleBarcode_FINAL), ]

# -------------------------------------------------------------------------------------------

# Add (final) followup information  
final_selection$IDC <- as.factor(final_selection$IDC)

n_timepoints <- table(final_selection$IDC)
table(n_timepoints)

final_selection$followup <- NA 
for (id in names(n_timepoints)) { 
  final_selection[final_selection$IDC %in% c(id), "followup"] <- n_timepoints[id]
}

addmargins(table(final_selection$followup, final_selection$Period))

# Retrieve information about extra (IMMUNO and Jeanette) samples ---------------

final_selection$biobank <- as.factor(gsub(" ", "", final_selection$biobank))

addmargins(table(final_selection$biobank, final_selection$Period))
addmargins(table(final_selection$biobank, final_selection$Period, final_selection$followup))

# final_selection$immuno_sample <- as.factor(ifelse(grepl("IMMUNO", final_selection$SampleBarcode_FINAL), "yes", "no"))

# final_selection$jeanet_sample <- as.factor(ifelse(final_selection$Period == 5 & 
#                                                     final_selection$immuno_sample == "no" &
#                                                     !grepl("DNA2023", final_selection$BoxID), "yes","no"))
# 
# table(final_selection$immuno_sample)
# table(final_selection$jeanet_sample)
# 
# final_selection$extras_sample <- as.factor(ifelse(final_selection$immuno_sample == "yes" |
#                                                     final_selection$genets_sample == "yes", "yes", "no"))
# 
# addmargins(table(final_selection$Period, final_selection$extras_sample))

# Too many samples were selected so subset to only those needed ================
# Get rid of (3792 - 3776 =) 16 samples

# will randomly select those with 2 timepoints in the extra samples (aka IMMUNO or Jeanette)
set.seed(310896)
to_exclude <- sample(final_selection[final_selection$biobank %in% c("BiobankEpidemiologie", "GENR") & 
                                     final_selection$followup < 3,"new_id"], 16)

final_selection <- final_selection[!final_selection$new_id %in% to_exclude, ]

# Recalculate followup in this set 
n_timepoints <- table(final_selection$IDC)
table(n_timepoints)

final_selection$followup <- NA 
for (id in names(n_timepoints)) { 
  final_selection[final_selection$IDC %in% c(id), "followup"] <- n_timepoints[id]
}

# Clean out variables not needed

final_selection_clean <- final_selection %>% 
  transmute(IDC_period = gsub(" ", "", new_id),
            SampleID = gsub(" ", "", SampleID_FINAL),
            SampleBarcode = gsub(" ", "", SampleBarcode_FINAL),
            IDC = IDC, 
            Period = Period, 
            Followup = followup, 
            Biobank = biobank,
            Concentration = Concentrationngul_FINAL)
# Check
any(duplicated(final_selection_clean$IDC_period))
any(duplicated(final_selection_clean$SampleID))
any(duplicated(final_selection_clean$SampleBarcode))

for (var_name in c('Followup', 'Period', 'Biobank', 'IDC')) {
  
  cat(var_name)
  var <- if (var_name == 'IDC') {
    table(final_selection_clean[, var_name]) } else { final_selection_clean[, var_name] }
  
  addmargins(table(var))
}

# Save 
saveRDS(final_selection_clean, "./files/final_selection_clean_280125.rds")

# OTHER COVARIATES ========================================================================

birth_data <- readit("CHILD-ALLGENERALDATA_15012024") %>% 
  # Combine maternal education time points
  mutate(mom_education = case_when(!is.na(EDUCM5) ~ EDUCM5, # GR1075 Education mother based on 5/6 year questionnaire
                                   !is.na(EDUCM3) ~ EDUCM3, # GR1065 Education mother based on 36 months questionnaire 
                                   TRUE ~ EDUCM), # Education level mother: highest education finished (pregnancy?)
  ) %>% 
  select(IDC, IDM, mom_education, 
         "sex" = "GENDER", 
         "birth_gest_age" = "GESTBIR",
         "birth_weight" = "WEIGHT")

maternal_smoking <- readit("MATERNALSMOKING_22112016") %>% 
  select("IDM" = "idm", 
         "mom_smoking" = "SMOKE_ALL")

covariates <- merge(birth_data, maternal_smoking, by = "IDM", all.x = TRUE)
str(covariates)

# Cleanup
rm(birth_data, maternal_smoking)

# Age at blood draw ---------------------------------------------------------------------

age <- list(readit("CHILDAGE5_08022013"), # 6 years
            readit("CHILDAGE_23032023"), # 10 years
            readit("20240801_Serena_agechild13-17") # # 14 to 18 years
) %>% 
  reduce(full_join, by='IDC') %>% 
  select(IDC, 
         "5" = "agey5child",
         "9" = "agechild9",
         "13" = "AGECHILD13",
         "17" = "ageChild17") %>%
  pivot_longer(cols = "5":"17",
               names_to = "Period", 
               values_to = "age")

# Merge =================================================================================

data <- merge(final_selection_clean, covariates, by = "IDC", all.x = TRUE)

data <- merge(data, age, by = c("IDC","Period"), all.x = TRUE)
data$age[data$Period == "Birth"] <- 0

str(data)
summary(data)

any(duplicated(data$IDC_period)) # FALSE 

# ==============================================================================
# Multivariate, reproducible randomization to proactively counter batch effects
# ==============================================================================

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("Omixer", force = TRUE)

library(Omixer, lib.loc="/home/s.defina/R/x86_64-pc-linux-gnu-library/4.4")

n_wells = 96
n_plates = ceiling(nrow(data) / n_wells) # 39 plates ( = 3744 samples ) 
n_on_last_plate <- nrow(data) - n_wells*(n_plates-1) # 32 samples left on plate 40

mask = c(rep(rep(0, n_wells), n_plates-1), # Fill all wells of first 39 plates
         c(rep(0, n_on_last_plate), 
           rep(1, n_wells-n_on_last_plate))) # Wells to be left empty on last plate

rand_vars = c(# "IDC",
  "Period",
  "Followup",
  "Biobank",
  "sex",
  "age",
  "mom_education",
  "mom_smoking", 
  "birth_gest_age", 
  "birth_weight")

n_iternations = 1000

pdf("./omixer_correlations.pdf", width = 18, height = 4)
mixlayout <- omixerRand(data, 
                        sampleId="IDC_period", # sample ID variable
                        block = "IDC",                # Paired sample identifier
                        iterNum = n_iternations,      # Number of layouts to generate
                        wells = n_wells,              # Number of wells on a plate
                        div = "col",                  # Plate subdivisions
                        positional = TRUE,            # Logical indicator of positional batch effects
                        plateNum = n_plates, 
                        # layout,                     # Custom plate layout as data frame
                        mask = mask,
                        # techVars,                   # Technical covariates (?)
                        randVars = rand_vars          # Randomization variables
)

dev.off()

# Note: SPSS turns IDC factor into factor level ID, clean that
mixlayout$IDC = stringr::str_split_i(mixlayout$sampleId, "_", 1)

# Check layout -----------------------------------------------------------------

problem = mixlayout %>% group_by(IDC) %>% count(plate) %>% count(IDC) %>% filter(n > 1)

# PROBLEM: block setting does not work perfectly, there are 23 instances (0.6%) 
# of time point from the same individual on multiple plates

mixlayout = mixlayout[!is.na(mixlayout$block),]

# To fix it, i manually replace positions with individuals with 1 followup time
for (id in problem$IDC) {
  
  id_positions = mixlayout[mixlayout$IDC == id, ]
  id_plates = table(id_positions$plate)
  # print(id_plates)
  
  small_plate_name = names(which.min(id_plates))
  small_plate_count = unname(id_plates[small_plate_name])
  big_plate_name = names(id_plates[!names(id_plates) %in% small_plate_name])
  
  to_replace = mixlayout[mixlayout$plate == small_plate_name & mixlayout$IDC == id,
                  c("plate","well","row", "column", "chip", "chipPos")]
  
  replacements = mixlayout[mixlayout$plate == big_plate_name & mixlayout$Followup == 1, ]
  replacements = replacements[sample(nrow(replacements), small_plate_count),]
  
  replace_positions = replacements[c("plate","well","row", "column", "chip", "chipPos")]
  
  # print(replace_positions)
  # print(to_replace)
  
  # Swap positions
  mixlayout[mixlayout$plate == small_plate_name & mixlayout$IDC == id,
     c("plate","well","row", "column", "chip", "chipPos")] <- replace_positions
  
  mixlayout[mixlayout$sampleId %in% replacements$sampleId, 
     c("plate","well","row", "column", "chip", "chipPos")] <- to_replace
  
}

# Check swaps
mixlayout %>% group_by(IDC) %>% count(plate) %>% count(IDC) %>% filter(n > 1)
mixlayout %>% group_by(plate) %>% count() %>% filter(n != 96)


# Manually re-randomize within plate -------------------------------------------

set.seed(310896)
rerandom = mixlayout %>%
  mutate(position = paste0(row, "_", column, "_", well, "_", chip,"_", chipPos)) %>%
  group_by(plate) %>%
  mutate(new_position = sample(position)) %>%
  tidyr::separate_wider_delim(new_position, "_",
                       names = c("new_row", "new_column","new_well", "new_chip", "new_chipPos")) %>%
  ungroup()  %>%
  as.data.frame()

tech_vars <- c("plate", "new_chip", "new_chipPos", 'chip', 'chipPos')

omixerCorr2 <- function(x, y) {

  ## Convert variables to numeric
  if(class(x) %in% c("character")) x <- as.numeric(forcats::as_factor(x))
  if(class(y) %in% c("character")) y <- as.numeric(forcats::as_factor(y))
  if(class(x) %in% c("Date", "factor")) x <- as.numeric(x)
  if(class(y) %in% c("Date", "factor")) y <- as.numeric(y)

  ## Save correlation estimates and p values
  if(length(unique(x)) < 5 & length(unique(y)) < 5) {
    ## Two categorical variables
    cvStat <- chisq.test(x, y, correct=FALSE)$statistics
    if(is.null(cvStat)) cvStat <- 0
    cv <- sqrt(cvStat/(length(x)*min(length(unique(x)),
                                     length(unique(y)))-1))
    corVal <- as.numeric(cv)
    corP <- chisq.test(x, y)$p.value
  } else {
    ## Otherwise, use Kendall's correlation coefficient and p-value
    corVal <- cor.test(x, y, method="kendall", exact=FALSE)$estimate
    corP <- cor.test(x, y, method="kendall", exact=FALSE)$p.value
  }

  ## Return as a data frame of estimate and p-value
  corTb <- tibble("corVal"=corVal, "corP"=corP)

  return(corTb)
}

corSelect <- lapply(rand_vars, function(y){
  print(y)
  corTbList <- lapply(tech_vars, function(z){
    print(z)
    cor <- omixerCorr2(rerandom[, y], rerandom[, z])
    corTb <- tibble(layoutNum=1, randVars=y, techVars=z,
                    corVal=cor$corVal, corP=cor$corP)
    return(corTb)
  })
  corTb <-bind_rows(corTbList)
  return(corTb)
})

new_corrs = do.call(rbind.data.frame, corSelect)

rerandom_clean <- rerandom %>%
  transmute(IDC_period = sampleId,
            SampleID, SampleBarcode,
            Plate = plate,
            Well = as.integer(new_well),
            Row = new_row,
            Column = as.integer(new_column),
            Chip = as.integer(new_chip),
            Chip_position = as.integer(new_chipPos),
            IDC = as.character(block),
            Period, Followup, Biobank, Concentration, sex, age) %>%
  arrange(Plate, Well)

mixlayout_clean <- mixlayout %>%
  transmute(IDC_period = sampleId,
            SampleID, SampleBarcode,
            Plate = plate,
            Well = as.integer(well),
            Row = row,
            Column = as.integer(column),
            Chip = as.integer(chip),
            Chip_position = as.integer(chipPos),
            IDC, Period, Followup, Biobank, Concentration, sex, age) %>%
  arrange(Plate, Well)

saveRDS(mixlayout, "randomized_sample_omixer_280125.rds")
saveRDS(rerandom_clean, "randomized_sample_custom_290125.rds")
saveRDS(new_corrs, "custom_correlations_290125.rds")

haven::write_sav(rerandom_clean, "HappyMums_randomized_samples_290125.sav")

###### INSPECT -----------------------------------------------------------------
library(dplyr)

rs = haven::read_sav("~/picklist/HappyMums_randomized_samples_290125.sav")

rs %>% group_by(IDC) %>% count(Plate) %>% count(IDC) %>% filter(n > 1)


rs$IDC = as.factor(rs$IDC)

check1 <- rs %>% filter(Plate == 1)

summary(check1)

table(check1$Followup)

length(unique(check1$IDC))

for (f in 1:5) {
  print(check1 %>% filter(Followup==f) %>% count(IDC) %>% filter(n != f))
}
# ------------------------------------------------------------------------------

