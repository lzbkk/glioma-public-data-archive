library(tidyverse)
library(survival)
library(survminer)
library(limma)
library(pheatmap)
getwd()
setwd("~/glioma/Data_Bulk_TCGA/Data_Merged")
gene_name <- "LAP3"
expr_tpm <- readRDS(file.path("data_annotated", "expr_tpm_glioma_anno.rds"))
expr_fpkm <- readRDS(file.path("data_annotated", "expr_fpkm_glioma_anno.rds"))
expr_count <- readRDS(file.path("data_annotated", "expr_count_glioma_anno.rds"))
cli <- readRDS(file.path("data_raw", "clinical_glioma.rds"))
# ============ 1. 样本去重：保留每个患者的一个样本 ============
#batch 为barcode 22-25
cli$barcode[1]
cli$batch2 <- substr(cli$barcode,27,28)
substr("TCGA-HT-7468-01A-11R-2027-07",27,28)
cli$shortLetterCode

# LGG中516个原发，其中14个复发，拥有18个复发样本。其中TCGA-DU-6404、TCGA-DU-6407、TCGA-FG-5965、TCGA-TQ-A7RK均有3样本
#GBM队列中c("TCGA-06-0221","TCGA-14-0736", "TCGA-19-0957", "TCGA-14-1402") 均只有复发样本
#800 行 × 114 列 GBM: 284 LGG: 516
x <- cli %>% filter(cohort == "GBM") %>% 
  group_by(patient) %>% 
  filter(any(shortLetterCode == "TR")) %>% 
  ungroup() %>% 
  select(patient, shortLetterCode) %>% 
  count(patient) %>% 
  arrange(n)
x  

cli %>% filter(patient == "TCGA-06-6698") %>% 
  select(patient, shortLetterCode, barcode)
cli %>% filter(cohort == "GBM", shortLetterCode == "NT") %>%
  select(patient, shortLetterCode) 


x <- colnames(cli)    
length(x)
x[grepl("vital", x)] 
getwd()
dim(cli)



# ============ 2. 表达矩阵与临床数据匹配 
cli_qc <- readRDS("/home/lzb/glioma/Data_Bulk_TCGA/Data_Merged/results/Clinical_Field_QC/clinical_glioma_analysis_fields.rds")
dim(cli_qc)
colnames(cli_qc)
View(cli_qc)
