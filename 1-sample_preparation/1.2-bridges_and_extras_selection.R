setwd("/home/s.defina/picklist")

library(foreign)
library(dplyr)

readit <- function(fname, dir = "./files") {
  out <- foreign::read.spss(file.path(dir, paste0(fname, ".sav")), 
                            to.data.frame = TRUE, use.value.labels = TRUE)
  
  print(summary(out))
  out
}

list.files("./files")

# dat0 <- readit("../20240624_Serena_DNAplusWholeblood_Birthuntil17")
# s450k_1 <- readit("Selection GENR_450kmeth_April2014")
# s450k_2 <- readit("Selection_GENR_450kmeth_release2_birth_20161215")

s450k_3 <- readit("Selection_GENR_450kmeth_release3_birth_20230608")
sEPIC_1 <- readit("Selection_GENR_MethylEPIC_release1_birth_20230717")

selectd <- readit("20240709_Serena_DNAisolated_Final")

birthewas <- merge(s450k_3, sEPIC_1, by = "IDC", all = TRUE, suffixes=c("_450k","_EPIC"))
alldata <- merge(selectd, birthewas, by = "IDC", all.x = TRUE)

# NOTE: 2 EPIC samples are inconsistent with the EPIC file 
# alldata$test <- ifelse(is.na(alldata$SampleID), 0, 1)
# table(alldata[,c("sample_EPIC", "test")])

# ==============================================================================

dat <- data.frame(IDC = alldata[,"IDC"])

# dat["sample_450k"] <- as.factor(ifelse(grepl("^\\s*$", dat0$Sample_ID_450Kmeth), "no", "yes"))
# dat["sample_EPIC"] <- as.factor(ifelse(grepl("^\\s*$", dat0$SampleID_MethylEpic_1), "no", "yes"))

dat["sample_450k"] <- as.factor(ifelse(is.na(alldata$Sample_ID), "no", "yes"))
dat["sample_EPIC"] <- as.factor(ifelse(is.na(alldata$SampleID), "no", "yes"))

table(dat[,c("sample_450k","sample_EPIC")], useNA="ifany")

dat["DNAm_birth"] <- as.factor(ifelse(dat[,"sample_450k"]=="yes", "450k", 
                               ifelse(dat[,"sample_EPIC"]=="yes", "EPIC", NA)))

summary(dat["DNAm_birth"])

dat["plate_450k"] <- as.factor(as.numeric(alldata[, "Sample_Plate_450k"]))
dat["plate_EPIC"] <- as.factor(as.numeric(alldata[, "Sample_Plate_EPIC"]))

table(dat[,c("plate_450k","sample_450k")], useNA="ifany")
table(dat[,c("plate_EPIC","sample_EPIC")], useNA="ifany")

# table(dat[,c("plate_450k","plate_EPIC")], useNA="ifany")

# Batch 450K
dat["batch_450k"] <- factor(substring(alldata$Sample_ID, 1, 1),
                            levels = c("7","8","1","2"),
                            labels = c("B1_7","B1_8","B2_1","B2_2"))
# Note: EPIC had only one batch

# ==============================================================================
# 2024 isolated 

dat["sample_06y"] <- ifelse(alldata[, "DNA5_2024_final"] == "DNA => 11.1 ng/ul", 1, 0)
dat["sample_10y"] <- ifelse(alldata[, "DNA9_2024_final"] == "DNA => 11.1 ng/ul", 1, 0)
dat["sample_14y"] <- ifelse(alldata[, "DNA13_2024_final"] == "DNA => 11.1 ng/ul", 1, 0)
dat["sample_18y"] <- ifelse(alldata[, "DNA17_2024_final"] == "DNA => 11.1 ng/ul", 1, 0)

table(dat["sample_06y"], useNA = "ifany")
table(dat["sample_10y"], useNA = "ifany")
table(dat["sample_14y"], useNA = "ifany")
table(dat["sample_18y"], useNA = "ifany")

dat["immuno_extra_06y"] <- ifelse(alldata[, "DNA5_2024_final"] == "DNA => 11.1 ng/ul isolated earlier", 1, 0)
dat["immuno_extra_10y"] <- ifelse(alldata[, "DNA9_2024_final"] == "DNA => 11.1 ng/ul isolated earlier", 1, 0)

addmargins(table(dat[c("immuno_extra_06y", "immuno_extra_10y")], useNA = "ifany"))

dat <- dat %>% 
  mutate(immuno_extras = as.factor(case_when(immuno_extra_06y == 1 & immuno_extra_10y == 1 ~ "both",
                                             immuno_extra_06y == 1 ~ "6 years",
                                             immuno_extra_10y == 1 ~ "10 years")))

table(dat["immuno_extras"], useNA = "ifany")

# -------

dat["followup"] <- rowSums(dat[,c("sample_06y","sample_10y","sample_14y","sample_18y")], na.rm = TRUE)

table(dat["followup"], useNA = "ifany")

# Only with DNAm at birth -------------
samp <- dat[!is.na(dat["DNAm_birth"]), ] # 2508 with either EPIC (1115) or 450k (1393)

table(samp["sample_06y"], useNA = "ifany")
table(samp["sample_10y"], useNA = "ifany")
table(samp["sample_14y"], useNA = "ifany")
table(samp["sample_18y"], useNA = "ifany")
table(samp["followup"], useNA = "ifany")

addmargins(table(samp[c("followup","immuno_extras")]))

samp$pick_06y[samp$followup > 0 & samp$immuno_extra_06y == 1] <- 1
samp$pick_10y[samp$followup > 0 & samp$immuno_extra_10y == 1] <- 1

table(samp$pick_06y) + table(samp$pick_10y) # 206 of 195 needed (note 145 have at least 2 more time points)

samp["tmp_followup_withimmuno"] <- rowSums(samp[c("followup","pick_06y","pick_10y")], na.rm = TRUE)
addmargins(table(samp[c("followup", "tmp_followup_withimmuno")]))

# Save immuno_pick -------------------------------------------------------------
outdata <- merge(selectd, samp, by="IDC")

table(outdata$pick_06y)
table(outdata$pick_10y)

# Rename
outdata$pickDNA_5_July2024 <- outdata$pick_06y
outdata$pickDNA_9_July2024 <- outdata$pick_10y

# haven::write_sav(outdata, "20240710_Serena_DNAisolated_Final_immuno_picklist.sav")

# Bridge selection -------------------------------------------------------------

table(samp[samp["tmp_followup_withimmuno"] > 0, c("batch_450k","plate_450k")])

samp %>%
  group_by(plate_450k) %>%
  filter(tmp_followup_withimmuno > 2) %>%
  count() %>% 
  print(n = 28)

samp %>%
  group_by(plate_EPIC) %>%
  filter(tmp_followup_withimmuno > 3) %>%
  count() %>% 
  print(n = 13)

# Select bridges for EPIC ------------------------------------------------------
samp$bridge_EPIC <- NA

for (p in levels(samp$plate_EPIC)){
  p_sub <- samp %>% filter(plate_EPIC == p & tmp_followup_withimmuno > 3)
  if (nrow(p_sub) < 2) stop("Not enough samples, be more lenient with followup")
  
  # randomly select two per sample
  set.seed(310896)
  random_ids <- sample(p_sub$IDC, size=2, replace = FALSE)
  
  samp$bridge_EPIC[samp$IDC %in% random_ids] <- 1
}

# Check
 table(samp[c("bridge_EPIC","plate_EPIC")], useNA = "ifany")
# Sanity
# table(samp[c("bridge_EPIC","sample_EPIC")], useNA = "ifany")
table(samp[c("bridge_EPIC","tmp_followup_withimmuno")], useNA = "ifany")

# Select bridges for 450k ------------------------------------------------------
samp$bridge_450k <- NA

for (p in levels(samp$plate_450k)){
  
  followup_filter <- 3
  p_sub <- samp %>% filter(plate_450k == p & 
                           tmp_followup_withimmuno > followup_filter)
  
  while (nrow(p_sub) < 2) {
    # cat(nrow(p_sub),"\n")
    followup_filter <- followup_filter - 1
    if (followup_filter < 0) stop("Not enough followup points in plate ", p)
    
    p_sub <- samp %>% filter(plate_450k == p & 
                             tmp_followup_withimmuno > followup_filter)
  }
  
  counts <- p_sub %>% count(tmp_followup_withimmuno, sort = FALSE)
  
  set.seed(310896)
  
  if (nrow(counts) == 1) {
    random_ids <- sample(p_sub$IDC, size = 2, replace = FALSE)
    
    samp$bridge_450k[samp$IDC %in% random_ids] <- 1
    
  } else {
    best_id <- p_sub$IDC[p_sub$tmp_followup_withimmuno == max(p_sub$tmp_followup_withimmuno)]
    if (length(best_id) > 1) stop("Something is up man")
    
    add_id <- sample(p_sub$IDC[p_sub$tmp_followup_withimmuno == min(p_sub$tmp_followup_withimmuno)], 
                     size = 1)
    
    samp$bridge_450k[samp$IDC %in% c(best_id, add_id)] <- 1
  }
}

# Check
table(samp[c("bridge_450k","plate_450k")], useNA = "ifany")
# Sanity
# table(samp[c("bridge_450k","sample_450k")], useNA = "ifany")
table(samp[c("bridge_450k","tmp_followup_withimmuno")], useNA = "ifany")

# Put them together
samp$bridges <- rowSums(samp[c("bridge_450k","bridge_EPIC")], na.rm = TRUE)

addmargins(table(samp[c("bridges","tmp_followup_withimmuno")], useNA = "ifany"))

# Save everything

outdata <- merge(selectd, samp, by="IDC")

table(outdata$pick_06y)
table(outdata$pick_10y)
table(outdata$bridges)

# Rename
outdata$pickDNA_5_July2024 <- outdata$pick_06y
outdata$pickDNA_9_July2024 <- outdata$pick_10y

haven::write_sav(outdata, "20240711_Serena_DNAisolated_Final_immuno_bridge_picklist.sav")
