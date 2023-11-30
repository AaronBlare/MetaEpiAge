rm(list=ls())

###############################################
# Installing packages
###############################################
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ChAMP")
BiocManager::install("methylGSA")
BiocManager::install("preprocessCore")
install.packages("devtools")
install.packages("splitstackshape")
devtools::install_github("danbelsky/DunedinPACE")
devtools::install_github("https://github.com/regRCPqn/regRCPqn")
library("ChAMP")
library("preprocessCore")
library("DunedinPACE")
library("regRCPqn")
library(readxl)
library(splitstackshape)
library("reticulate")
pandas <- import("pandas")

###############################################
# Setting variables
###############################################
dataset <- 'GSE120610'
arraytype <- 'EPIC'

dataset_ref <- 'GSE87571'

###############################################
# Setting path
###############################################
path_data <- "D:/YandexDisk/Work/pydnameth/datasets/GPL23976/GSE120610/raw"
path_horvath <- "D:/YandexDisk/Work/pydnameth/draft/10_MetaEPIClock/MetaEpiAge"
path_harm_ref <- "D:/YandexDisk/Work/pydnameth/draft/10_MetaEPIClock/MetaEpiAge/GPL13534/GSE87571/"
path_pc_clocks <- "D:/YandexDisk/Work/pydnameth/datasets/lists/cpgs/PC_clocks/"
path_work <- "D:/YandexDisk/Work/pydnameth/draft/10_MetaEPIClock/MetaEpiAge/GPL23976/GSE120610"
setwd(path_work)

###############################################
# Load annotations
###############################################
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

###############################################
# Import data
###############################################
pd <- as.data.frame(read_excel(paste(path_data,"/controls.xlsx", sep="")))
pd$description <- gsub(" ","", as.character(pd$description))
rownames(pd) <- pd$description
betas <- as.data.frame(read.csv(paste(path_data,"/GSE120610_Matrix_processed.csv", sep="")))
rownames(betas) <- betas$ID_REF
col_odd <- seq_len(ncol(betas)) %% 2
betas <- betas[ , col_odd == 0]
missed_in_betas <- setdiff(row.names(pd), colnames(betas))
missed_in_pheno <- setdiff(colnames(betas), row.names(pd))
betas <- betas[, rownames(pd)]
cpgs_common <- intersect(rownames(betas), rownames(ann450k))
betas <- betas[cpgs_common, ]

###############################################
# Harmonization
###############################################
mvals <- logit2(betas)
mvals <- data.frame(rownames(mvals), mvals)
colnames(mvals)[1] <- "ID_REF"
mvals <- regRCPqnREF(M_data=mvals, ref_path=path_harm_ref, data_name=dataset_ref)
betas <- ilogit2(mvals)

###############################################
# PC clocks
# You need to setup path to 3 files from original repository (https://github.com/MorganLevineLab/PC-Clocks):
# 1) run_calcPCClocks.R
# 2) run_calcPCClocks_Accel.R
# 3) CalcAllPCClocks.RData (very big file but it is nesessary)
# You also need to apply changes from this issue: https://github.com/MorganLevineLab/PC-Clocks/issues/10
###############################################
source(paste(path_pc_clocks, "run_calcPCClocks.R", sep = ""))
source(paste(path_pc_clocks, "run_calcPCClocks_Accel.R", sep = ""))
pheno <- data.frame(
  'Sex' = pd$Sex,
  'Age' = pd$Age,
  'Tissue' = pd$Tissue
)
pheno['Female'] <- 1
pheno$Age <- as.numeric(pheno$Age)
pheno[pheno$Sex == 'M', 'Female'] <- 0
rownames(pheno) <- rownames(pd)
pc_clocks <- calcPCClocks(
  path_to_PCClocks_directory = path_pc_clocks,
  datMeth = t(betas),
  datPheno = pheno,
  column_check = "skip"
)
pc_clocks <- calcPCClocks_Accel(pc_clocks)
pc_ages <- list("PCHorvath1", "PCHorvath2", "PCHannum", "PCHannum", "PCPhenoAge", "PCGrimAge")
for (pc_age in pc_ages) {
  pd[rownames(pd), pc_age] <- pc_clocks[rownames(pd), pc_age]
}

###############################################
# Create data for Horvath's calculator
###############################################
cpgs_horvath_old <- read.csv(
  paste(path_horvath, "/cpgs_horvath_calculator.csv", sep=""),
  header=TRUE
)$CpG
cpgs_horvath_new <- read.csv(
  paste(path_horvath, "/datMiniAnnotation4_fixed.csv", sep=""),
  header=TRUE
)$Name
cpgs_horvath <- intersect(cpgs_horvath_old, rownames(betas))
cpgs_missed <- setdiff(cpgs_horvath_old, rownames(betas))
betas_missed <- matrix(data='NA', nrow=length(cpgs_missed), dim(betas)[2])
rownames(betas_missed) <- cpgs_missed
colnames(betas_missed) <- colnames(betas)
betas_horvath <- rbind(betas[cpgs_horvath, ], betas_missed)
betas_horvath <- data.frame(row.names(betas_horvath), betas_horvath)
colnames(betas_horvath)[1] <- "ProbeID"
write.csv(
  betas_horvath,
  file="betas_horvath.csv",
  row.names=FALSE,
  quote=FALSE
)

pheno_horvath <- data.frame(
  'Sex' = pd$Sex,
  'Age' = pd$Age,
  'Tissue' = pd$Tissue
)
pheno_horvath['Female'] <- 1
pheno_horvath[pheno_horvath$Sex == 'M', 'Female'] <- 0
pheno_horvath$Age <- as.numeric(pheno_horvath$Age)
rownames(pheno_horvath) <- rownames(pd)
pheno_horvath <- data.frame(row.names(pheno_horvath), pheno_horvath[ ,!(names(pheno_horvath) %in% c("Sex"))])
colnames(pheno_horvath)[1] <- "Sample_Name"
write.csv(
  pheno_horvath,
  file="pheno_horvath.csv",
  row.names=FALSE,
  quote=FALSE
)

###############################################
# DunedinPACE
# This dataset contains many missing values, which will be imputed from reference dataset GSE87571
###############################################
path_ref <- "D:/YandexDisk/Work/pydnameth/datasets/GPL13534/GSE87571/raw/idat"
load_ref <- champ.load(
  directory = path_ref,
  arraytype = '450K',
  method = "minfi",
  methValue = "B",
  autoimpute = TRUE,
  filterDetP = TRUE,
  ProbeCutoff = 0.1,
  SampleCutoff = 0.1,
  detPcut = 0.01,
  filterBeads = FALSE,
  beadCutoff = 0.05,
  filterNoCG = FALSE,
  filterSNPs = FALSE,
  filterMultiHit = FALSE,
  filterXY = FALSE,
  force = TRUE
)
betas_ref <- getBeta(preprocessFunnorm(load_ref$rgSet))
cpgs_to_fill <- setdiff(rownames(betas_ref), rownames(betas))
betas[cpgs_to_fill, ] <- 0
values_to_fill <- apply(betas_ref, 1, median, na.rm=T)
values_to_fill <- data.frame(values_to_fill)

n_filled_samples <- 0
for (sample_id in colnames(betas)) {
  betas[cpgs_to_fill, sample_id] <- values_to_fill[cpgs_to_fill, 'values_to_fill']
  n_filled_samples <- n_filled_samples + 1
  print(n_filled_samples)
}


pace <- PACEProjector(betas)
pd['DunedinPACE'] <- pace$DunedinPACE

###############################################
# Save modified pheno
###############################################
write.csv(pd, file = "pheno.csv")
