setwd('/pub5/xiaoyun/Jobs/J22/temp/Github')

suppressPackageStartupMessages(library(tidyverse))
library(ggsci)

# functions ---------------------------------------------------------------

# statistical test between varibles

var_relate_als <- function(df_1, df_2 = NULL, var_1, var_2, df1_sam_col, df2_sam_col = NULL){
  
  # df_1：sample characteristics (like mutation burden) 
  # df_2：sample clinical features (like age)
  # var_1：interested varible names of df_1
  # var_2：interested varible names of df_2，if df_2 is NULL，then var_2 is also from df_1
  # df1_sam_col：sample colname of df_1
  # df2_sam_col：sample colname of df_2
  
  if (is.null(df_2)) {
    df <- df_1 %>% select(Sample = df1_sam_col, var_1, var_2)
  } else {
    df_1 <- df_1 %>% select(Sample = df1_sam_col, var_1)
    df_2 <- df_2 %>% select(Sample = df2_sam_col, var_2)
    df <- df_1 %>% left_join(df_2, by = "Sample")
  }  
  
  var_type <- df %>% select(-Sample) %>% map_chr(class)
  var_con <- names(var_type)[var_type %in% c("numeric", "integer")]
  
  var1_con <- intersect(var_1, var_con) %>% set_names(.)
  var2_cat <- setdiff(var_2, var_con) %>% set_names(.)
  
  test_coca <- var1_con %>% map(function(sin_var1){
    sin_var1_p <- map_dbl(var2_cat, possibly(function(sin_var2){
      df_sin_var2_nm <- n_distinct(df[[sin_var2]])
      if (df_sin_var2_nm == 2) {
        p <- wilcox.test(df[[sin_var1]] ~ df[[sin_var2]])$p.value
      } else if (df_sin_var2_nm > 2) {
        p <- kruskal.test(df[[sin_var1]] ~ df[[sin_var2]])$p.value
      }
      p
    }, NA))
  })
  
  test_coca_tib <- test_coca %>% 
    bind_cols() %>% 
    mutate(Variable = var2_cat) %>% 
    select(Variable, everything())
  
}

# pairwise comparison using dunn test after kruskal test

var_dunn_als <- function(df_1, df_2 = NULL, var_1, var_2, df1_sam_col, df2_sam_col = NULL){
  
  # df_1：sample characteristics (like mutation burden) 
  # df_2：sample clinical features (like age)
  # var_1：interested varible names of df_1
  # var_2：interested varible names of df_2，if df_2 is NULL，then var_2 is also from df_1
  # df1_sam_col：sample colname of df_1
  # df2_sam_col：sample colname of df_2
  
  if (is.null(df_2)) {
    df <- df_1 %>% select(Sample = df1_sam_col, var_1, var_2)
  } else {
    df_1 <- df_1 %>% select(Sample = df1_sam_col, var_1)
    df_2 <- df_2 %>% select(Sample = df2_sam_col, var_2)
    df <- df_1 %>% left_join(df_2, by = "Sample")
  }  
  
  var_type <- df %>% select(-Sample) %>% map_chr(class)
  var_con <- names(var_type)[var_type %in% c("numeric", "integer")]
  
  var1_con <- intersect(var_1, var_con) %>% set_names(.)
  var2_cat <- setdiff(var_2, var_con) %>% set_names(.)
  
  test_coca <- var1_con %>% map(function(sin_var1){
    
    sin_var1_result <- map(var2_cat, function(sin_var2){
      dunn_result <- dunn.test::dunn.test(df[[sin_var1]], df[[sin_var2]])
      dunn_tib <- bind_rows(dunn_result[2:length(dunn_result)]) %>% 
        mutate(P.adjusted = p.adjust(dunn_result$P, method = "BH"), var1_var2 = str_c(sin_var1, sin_var2, sep = " - "))
    })
    sin_var1_tib <- sin_var1_result %>% bind_rows()
  })
  test_coca_tib <- test_coca %>% bind_rows()
}

# -------------------------------------------------------------------------

load("data/glioma.mut.ccf.RData") # load mutation clonality data
glioma_mut_clonality_data <- com.data
load("data/gliomaClinicalData.RData")# load clinical data, pan.glio.mer.cli.data


# merge tumor location variables
pan.glio.mer.cli.data <- pan.glio.mer.cli.data %>%
  mutate(age_group = as.character(cut_width(age_at_initial_pathologic_diagnosis, width = 10, boundary = 10)), 
         age_group = ifelse(age_group %in% c("[10,20]", "(20,30]"), "<=30", age_group), 
         age_group = ifelse(age_group %in% c("(70,80]", "(80,90]"), ">70", age_group), 
         grade_temp = ifelse(str_detect(histological_type, "GBM"), "G4", NA), 
         neoplasm_histologic_grade = ifelse(is.na(neoplasm_histologic_grade), grade_temp, neoplasm_histologic_grade))

# number of clonal and subclonal mutations of each sample
sam_clonality_num <- glioma_mut_clonality_data %>%
  mutate(bcr_patient_barcode = str_sub(patient, 1, 12)) %>%
  count(bcr_patient_barcode, CI95.timing) %>%
  pivot_wider(names_from = CI95.timing, values_from = n, values_fill = list(n = 0)) %>%
  left_join(pan.glio.mer.cli.data, by = "bcr_patient_barcode") %>%
  select(bcr_patient_barcode:Subclonal, IDH.status, IDH.codel.subtype, Grade.rough, WHO_subtype)

# test the association between clinical features and clonality --------

# clonality vs IDH codel subtypes in all gliomas

glioma_clonality_kw_test <- var_relate_als(sam_clonality_num, var_1 = c("Clonal", "Subclonal"), var_2 = "IDH.codel.subtype", df1_sam_col = "bcr_patient_barcode") %>% bind_cols()
glioma_clonality_dunn_test <- var_dunn_als(sam_clonality_num, var_1 = c("Clonal", "Subclonal"), var_2 = "IDH.codel.subtype", df1_sam_col = "bcr_patient_barcode") %>% bind_rows()
write_tsv(bind_rows(glioma_clonality_kw_test, glioma_clonality_dunn_test), "result/2.Clinical_molecular_feature_association_of_clonality/clonality_diff_among_subtypes.tsv", na = "")

# clonality vs clinical features by IDH codel subtype

clonality_cli_test_list <- sam_clonality_num %>%
  split(.$IDH.codel.subtype) %>%
  map(~bind_rows(var_relate_als(., pan.glio.mer.cli.data, var_1 = c("Clonal", "Subclonal"), 
                                var_2 = c("age_group", "gender", "neoplasm_histologic_grade", 
                                          "laterality", "family_history_of_cancer", "tumor_location", "first_presenting_symptom", 
                                          "MGMT.promoter.status", "ATRX.status", "TERT.promoter.status", "Chr.7.gain.Chr.10.loss"), 
                                "bcr_patient_barcode", "bcr_patient_barcode")))

clonality_cli_test_tib <- names(clonality_cli_test_list) %>%
  map(~mutate(clonality_cli_test_list[[.]], IDH.codel.subtype = .)) %>%
  bind_rows() %>%
  select(IDH.codel.subtype, everything())

write_tsv(clonality_cli_test_tib, "result/2.Clinical_molecular_feature_association_of_clonality/clonality_cli_test_result_by_subtype.tsv", na = "")


# plots -------------------------------------------------------------------

# combine data
comb_plot_data <- sam_clonality_num %>%
  select(bcr_patient_barcode:Subclonal) %>%
  pivot_longer(Clonal:Subclonal, names_to = "CI95.timing", values_to = "mut_num") %>%
  left_join(pan.glio.mer.cli.data, by = "bcr_patient_barcode") %>%
  mutate(IDH.status = str_c("IDH", IDH.status, sep = " "))

# clonality vs IDH codel subtypes in all gliomas

IDH_codel_clo_all_gliomas <- ggplot(filter(comb_plot_data, !is.na(IDH.codel.subtype)), aes(IDH.codel.subtype, mut_num)) +
  geom_point(aes(col = CI95.timing), position = position_jitterdodge(jitter.width = 0.5, dodge.width = 1), shape = 16, alpha = 7/10, size = 4) + 
  geom_boxplot(fill = NA, position = position_dodge(width = 1), outlier.color = NA) +
  labs(x = NULL, y = "Number of mutations", col = NULL) +
  scale_color_manual(values = c("#D45658", "#3A84BB")) +
  #coord_cartesian(ylim = c(0, 100)) +
  facet_grid(. ~ CI95.timing) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 0.5, vjust = 0.5), axis.ticks.x = element_blank(), 
        axis.title.y = element_text(size = 12), strip.background = element_blank(), strip.text = element_blank(), 
        panel.grid = element_blank(), legend.position = "top")
ggsave("result/2.Clinical_molecular_feature_association_of_clonality/IDH_codel_clo_all_gliomas_limited.pdf", IDH_codel_clo_all_gliomas, width = 8, height = 8)
ggsave("result/2.Clinical_molecular_feature_association_of_clonality/IDH_codel_clo_all_gliomas.pdf", IDH_codel_clo_all_gliomas, width = 8, height = 8)

# clonality vs grade in IDHmut-codel

grade_clonality_IDHmut_codel <- ggplot(filter(comb_plot_data, !is.na(neoplasm_histologic_grade), IDH.codel.subtype %in% "IDHmut-codel"), aes(neoplasm_histologic_grade, mut_num)) +
  geom_point(aes(col = CI95.timing), position = position_jitterdodge(jitter.width = 0.5, dodge.width = 1), shape = 16, alpha = 7/10, size = 4) + 
  geom_boxplot(fill = NA, position = position_dodge(width = 1), outlier.color = NA) +
  labs(x = NULL, y = "Number of mutations", col = NULL) +
  scale_color_manual(values = c("#D45658", "#3A84BB")) +
  scale_x_discrete(labels = c("Grade II", "Grade III")) +
  facet_grid(. ~ CI95.timing) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 0.5, vjust = 0.5), axis.ticks.x = element_blank(), 
        axis.title.y = element_text(size = 12), strip.background = element_blank(), strip.text = element_blank(), 
        panel.grid = element_blank(), legend.position = "top")
ggsave("result/2.Clinical_molecular_feature_association_of_clonality/grade_clonality_IDHmut-codel.pdf", grade_clonality_IDHmut_codel, width = 8, height = 8)

# clonality vs age in IDHmut-codel, IDHmut-non-codel and IDHwt gliomas

age_comb_plot_data <- comb_plot_data %>% 
  filter(!is.na(age_group)) %>%
  mutate(age_group = factor(age_group, levels = c("<=30", "(30,40]", "(40,50]", "(50,60]", "(60,70]", ">70")), 
         age_group = fct_recode(age_group, "30-40" = "(30,40]", "40-50" = "(40,50]", "50-60" = "(50,60]", "60-70" = "(60,70]"))

age_clonality <- ggplot(age_comb_plot_data, aes(age_group, mut_num)) +
  geom_point(aes(col = age_group), position = position_jitterdodge(jitter.width = 2, dodge.width = 1), shape = 16, alpha = 5/10) + 
  geom_boxplot(fill = NA, position = position_dodge(width = 1), outlier.color = NA) +
  labs(x = "Age at diagnosis", y = "Number of mutations", fill = NULL) +
  scale_color_nejm() +
#   coord_cartesian(ylim = c(0, 100)) +
  facet_grid(CI95.timing ~ IDH.codel.subtype, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5), axis.ticks.x = element_blank(), 
        axis.title.y = element_text(size = 12), strip.background = element_blank(), 
        panel.grid = element_blank(), legend.position = "none")
ggsave("result/2.Clinical_molecular_feature_association_of_clonality/age_clonality_limited.pdf", age_clonality)
ggsave("result/2.Clinical_molecular_feature_association_of_clonality/age_clonality.pdf", age_clonality)

# clonality vs age in grade IV of IDHwt gliomas

age_clonality_IDHwt_gradeIV <- ggplot(filter(age_comb_plot_data, !is.na(age_group), WHO_subtype %in% "IDHwt_grade_IV"), aes(age_group, mut_num)) +
  geom_point(aes(col = age_group), position = position_jitterdodge(jitter.width = 2, dodge.width = 1), shape = 16, alpha = 5/10) + 
  geom_boxplot(fill = NA, position = position_dodge(width = 1), outlier.color = NA) +
  labs(x = NULL, y = "Number of mutations", fill = NULL) +
  scale_color_nejm() +
  facet_grid(. ~ CI95.timing) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5), axis.ticks.x = element_blank(), 
        axis.title.y = element_text(size = 12), strip.background = element_blank(),  
        panel.grid = element_blank(), legend.position = "none")
ggsave("result/2.Clinical_molecular_feature_association_of_clonality/age_clonality_IDHwt_gradeIV.pdf", age_clonality_IDHwt_gradeIV)




















