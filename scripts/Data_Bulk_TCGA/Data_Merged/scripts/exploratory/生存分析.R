library(tidyverse)
getwd()
setwd("~/glioma/Data_Bulk_TCGA/Data_Merged")
gene_name <- "LAP3"
expr_tpm <- readRDS(file.path("data_annotated", "expr_tpm_glioma_anno.rds"))
expr_fpkm <- readRDS(file.path("data_annotated", "expr_fpkm_glioma_anno.rds"))
expr_count <- readRDS(file.path("data_annotated", "expr_count_glioma_anno.rds"))
cli <- readRDS(file.path("data_raw", "clinical_glioma.rds"))

identical(base::setdiff(colnames(expr_count), "gene_type"), cli$barcode)
t(expr_count[gene_name,!colnames(expr_count) %in% "gene_type"])

df1 <- expr_count[gene_name,!colnames(expr_count) %in% "gene_type"] %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("barcode")
df2 <- expr_tpm[gene_name,!colnames(expr_count) %in% "gene_type"] %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("barcode")   

df <- expr_fpkm[gene_name,!colnames(expr_count) %in% "gene_type"] %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("barcode") %>% 
  inner_join(cli, by = "barcode") %>% 
  mutate(group = case_when(LAP3 >= median(LAP3) ~ 1,
                   TRUE ~ 0))

mer_df <-function(gene_name, df1, df2){
  df <- df1[gene_name,!colnames(df1) %in% "gene_type", drop = FALSE] %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("barcode") %>% 
    inner_join(df2, by = "barcode") %>% 
    mutate(group = case_when(
      .data[[gene_name]] >= median(.data[[gene_name]]) ~ 1L,
      TRUE ~ 0L))
  return(df)
} 
x1 <- mer_df("LAP3", expr_tpm, cli)
table(cli$shortLetterCode)
class(length(unique(cli$patient)))
x <- cli %>%
  dplyr::count(patient, name = "asd") %>%
  dplyr::filter(asd > 1) %>% 
  dplyr::arrange(desc(asd))
x
dim(x)
x2 <- cli %>%
  dplyr::group_by(patient) %>% 
  dplyr::filter(n() > 1) %>% 
  dplyr::summarise(asd = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(asd)) %>% 
  dplyr::select(patient,asd)
x2
dim(x2)
x3 <- cli %>% 
  group_by(patient) %>% 
  dplyr::filter(n()>1) %>% 
  dplyr::mutate(t = n()) %>% 
  dplyr::arrange(desc(t)) %>% 
  dplyr::select(patient, t, shortLetterCode)
x3
dim(x3)  
x3
table(x3$shortLetterCode)
x4 <- cli %>% 
  group_by(patient) %>% 
  dplyr::filter(patient %in% x2$patient & all(shortLetterCode != "TP"))
  ungroup()
x4$patient
dim(x4)
# 筛出多行患者，给出重复次数，判断生存时间，生存状态是否一致。
x5 <- cli %>% 
  group_by(patient) %>% 
  dplyr::filter(n() > 1) %>% 
  summarise(t = n(), p_statue = paste(shortLetterCode,collapse = ","),
            os = all(vital_status == dplyr::first(vital_status)),
            follow_up = all(days_to_last_follow_up == dplyr::first(days_to_last_follow_up)),
            `Survival..months` = all(`paper_Survival..months.` == dplyr::first(`paper_Survival..months.`)),
            na = if_else(is.na(`Survival..months`), paste(paper_Survival..months., collapse = ","), NA)
                                   ) %>% 
  ungroup() 
dim(x5)
x5
x6 <- x5 %>% dplyr::filter(is.na(Survival..months))
dim(x6)
x6$patient
cli %>% 
  dplyr::select(patient,vital_status,days_to_last_follow_up,paper_Survival..months.) %>% 
  dplyr::filter(patient %in% x6$patient)
sum(x5$na)
x5 %>% dplyr::filter(follow_up == TRUE)
table(x5$Survival..months)  
any(is.na(cli$paper_Survival..months.))
sum(is.na(cli$days_to_last_follow_up))
cli %>% dplyr::filter(is.na(days_to_last_follow_up)) %>% 
  dplyr::select(patient,vital_status,days_to_last_follow_up)
ver <- colnames(cli)[grepl("paper", colnames(cli))]
x7 <- cli %>% 
  dplyr::select(matches("paper"))
dim(x7)    
colnames(x7)
x8 <- cli %>% 
  dplyr::select(-matches("paper"))
dim(x8)
colnames(x8)
set.seed(123)
x9 <- cli %>% 
  dplyr::slice_sample(n = 100)
print(x9)  
cli$shortLetterCode
#多行患者，仅保留缺失值最少的一行，且只保留“TP”患者
x <- cli %>% 
  dplyr::group_by(patient) %>% 
  dplyr::mutate(na_count = rowSums(is.na(across(everything())))) %>% 
  dplyr::filter(shortLetterCode == "TP") %>% 
  dplyr::slice_min(na_count, n = 1, with_ties = FALSE) %>% 
  ungroup()
length(unique(cli$patient))
dim(x)
length(unique(x$patient))
x1 <- cli %>% 
  dplyr::group_by(patient) %>% 
  dplyr::slice(1) %>% 
  ungroup()
dim(x1)
y <- base::setdiff(unique(cli$patient), unique(x$patient))
y
x2 <- cli %>% 
  dplyr::filter(patient %in% base::setdiff(unique(cli$patient), unique(x$patient))) %>% 
  dplyr::select(patient, shortLetterCode, cohort)
x2
length(unique(x2$patient))
table(cli$shortLetterCode)
cli %>% dplyr::count(shortLetterCode, name = "样本数量")
x3 <- cli %>%
  group_by(patient) %>%
  dplyr::count(cohort, name = "条数") %>%
  ungroup()
x3
dim(x3)
cli_uni <- cli %>% 
  dplyr::group_by(patient) %>% 
  dplyr::mutate(na_count = rowSums(is.na(across(everything())))) %>% 
  dplyr::filter(shortLetterCode == "TP") %>% 
  dplyr::slice_max(na_count, n = 1, with_ties = FALSE) %>% 
  dplyr::ungroup()
cli_uni %>% 
  dplyr::count(cohort)
col <- colnames(cli)
col[grepl("type", col, ignore.case = TRUE)]  
cli %>% dplyr::count(sample_type,shortLetterCode,cohort, .drop = FALSE)  
cli %>% dplyr::filter(cohort == "GBM") %>% 
  dplyr::mutate(type = case_when(substr(barcode,14,15) == "01" ~ "TP",
                                 substr(barcode,14,15) == "02" ~ "TR",
                                 substr(barcode,14,15) == "11" ~ "NT",
                                 TRUE ~ NA_character_)) %>% 
  dplyr::select(barcode, type, shortLetterCode)
x <- "TCGA-06-6390-01A-11R-A96S-41"  
substr(x,16)
cli_uni$barcode  
x <- cli_uni %>% 
  dplyr::filter(substr(barcode, 14, 16) != "01A") %>% 
  dplyr::select(barcode, cohort, sample_type, patient)
x2 <- cli %>% 
  dplyr::filter(patient %in% x$patient) %>% 
  dplyr::select(patient, barcode)
x2 %>% dplyr::count(patient)
cli %>% dplyr::filter(patient == "TCGA-02-0026") %>% 
  dplyr::select(patient, barcode)
table(cli_uni$paper_RNAseq)
se_col <- function(name){x <- colnames(cli_uni)
  x[grepl(name, x, ignore.case = TRUE)]}
se_col("TERT")
fi <- function(i){
  cli %>% dplyr::select(all_of(i))}

fi(se_col("tert"))[1:3, ]

cli_uni %>% dplyr::count(paper_RNAseq, cohort)
cli %>% dplyr::filter(cohort == "GBM") %>%
  dplyr::filter(paper_RNAseq == "Yes") %>% 
  dplyr::count()

cli$paper_RNAseq
cli %>% filter(shortLetterCode == "TP") %>%
  slice(n = 1, .drop = TRUE) %>% 
  count(paper_RNAseq, cohort)

cli %>% 
  filter(shortLetterCode == "TP") %>%
  distinct(patient, .keep_all = TRUE) %>% 
  count(paper_RNAseq, cohort, .drop = TRUE)
usethis::edit_r_profile()
colnames(cli)
colnames(expr_tpm)
table(cli$cohort)
dim(cli)
cli$patient
p_id <- substr(cli$barcode,1,12)
identical(p_id, cli$patient)
