setwd('/pub5/xiaoyun/Jobs/J22/temp/Github')


suppressPackageStartupMessages(library(tidyverse))
library(ggsci)


load("data/gliomaClinicalData.RData") # load clinical data, pan.glio.mer.cli.data
load("data/PureGliomaMutData.RData") # load mutation data
load("data/glioma.mut.ccf.RData") # load mutation clonality data
load("data/pan.glioma.cancer.gene.RData") # load cancer gene data
load("data/panCanAtlas.abs.seg.RData") # load absolute copy number data

glioma_mut_clonality_data <- com.data
pure.glioma.mut.data <- pure.glioma.mut.data %>%
  mutate(bcr_patient_barcode = str_sub(Patient, 1, 12))

# add subtype and grade information
glioma_mut_clonality_data <- glioma_mut_clonality_data %>%
  mutate(bcr_patient_barcode = str_sub(patient, 1, 12), driver_mut = ifelse(Hugo_Symbol %in% pan.glioma.cancer.gene, "Yes", "No")) %>%
  left_join(pan.glio.mer.cli.data, by = "bcr_patient_barcode") %>%
  select(patient:driver_mut, IDH.codel.subtype, Grade.rough, WHO_subtype)

# statistics divided by IDH codel subtypes --------------------------------------

# original mutation number of clonality samples 
clonality_sam_ori_mut_stat <- pure.glioma.mut.data %>%
  semi_join(glioma_mut_clonality_data, by = c("Patient" = "patient")) %>%
  left_join(pan.glio.mer.cli.data, by = "bcr_patient_barcode") %>%
  count(IDH.codel.subtype, name = "ori_sam_mut")

# Number of core samples
glioma_core_sam_stat <- pure.glioma.mut.data %>%
  semi_join(panCanAtlas.abs.seg, by = c("Patient" = "SampleID")) %>%
  left_join(pan.glio.mer.cli.data, by = "bcr_patient_barcode") %>%
  count(IDH.codel.subtype, Patient, name = "mut_num") %>%
  count(IDH.codel.subtype, name = "core_sample_num")
glioma_core_sam <- unique(glioma_core_sam_stat$Patient)
save(glioma_core_sam, file = "data/glioma_core_sam.Rdata")

# Number of clonality samples
sam_stat <- glioma_mut_clonality_data %>% 
  count(IDH.codel.subtype, patient, name = "mut_num") %>%
  count(IDH.codel.subtype, name = "clonality_sample_num")
  
# Number of clonal, subclonal and total mutations
clo_sub_stat <- glioma_mut_clonality_data %>%
  group_by(IDH.codel.subtype, CI95.timing) %>%
  summarise(mut_num = n()) %>%
  mutate(clonality_mut_num = sum(mut_num), frac = mut_num / clonality_mut_num) %>%
  left_join(clonality_sam_ori_mut_stat, by = "IDH.codel.subtype") %>%
  mutate(clonality_over_ori_frac = clonality_mut_num / ori_sam_mut)
  
# Number of clonal and subclonal mutations of each sample
sam_clo_sub_stat <- glioma_mut_clonality_data %>%
  group_by(IDH.codel.subtype, CI95.timing, patient) %>%
  summarise(sam_mut_num = n()) %>%
  summarise(sam_med_num = median(sam_mut_num), 
            sam_min_num = min(sam_mut_num), 
            sam_max_num = max(sam_mut_num))

# Using Chisq test to test the association between clonality and cancer gene mutations
driver_clo_relate_test <- glioma_mut_clonality_data %>%
  filter(!Variant_Classification %in% c("3'Flank", "5'Flank", "3'UTR", "5'UTR", "Intron", "RNA", "Silent")) %>%
  split(.$IDH.codel.subtype) %>%
  map_dbl(~chisq.test(.$driver_mut, .$CI95.timing)$p.value) %>%
  tibble(IDH.codel.subtype = names(.), driver_col_p_value = .)

# Fraction of clonal and subclonal mutations of each cancer gene
cancer_gene_sub_stat <- glioma_mut_clonality_data %>%
  filter(!Variant_Classification %in% c("3'Flank", "5'Flank", "3'UTR", "5'UTR", "Intron", "RNA", "Silent"), driver_mut %in% "Yes") %>%
  group_by(IDH.codel.subtype, Hugo_Symbol, CI95.timing) %>%
  summarise(mut_num = n()) %>%
  mutate(total_mut_num = sum(mut_num), frac = mut_num / total_mut_num) %>%
  pivot_wider(id_cols = c("IDH.codel.subtype", "Hugo_Symbol", "total_mut_num"), names_from = CI95.timing, values_from = frac, values_fill = list(frac = 0)) %>%
  group_by(IDH.codel.subtype) %>%
  arrange(desc(Subclonal), .by_group = TRUE)

# Fraction of clonal and subclonal mutations of each IDH1 protein variant type
IDH1_locus_sub_stat <- glioma_mut_clonality_data %>%
  filter(!Variant_Classification %in% c("3'Flank", "5'Flank", "3'UTR", "5'UTR", "Intron", "RNA", "Silent"), Hugo_Symbol %in% "IDH1") %>%
  group_by(IDH.codel.subtype, HGVSp_Short, CI95.timing) %>%
  summarise(mut_num = n()) %>%
  mutate(total_mut_num = sum(mut_num), frac = mut_num / total_mut_num) %>%
  pivot_wider(id_cols = c("IDH.codel.subtype", "HGVSp_Short", "total_mut_num"), names_from = CI95.timing, values_from = frac, values_fill = list(frac = 0))

merge_stat <- bind_rows(glioma_core_sam_stat, sam_stat, clo_sub_stat, cancer_gene_sub_stat, IDH1_locus_sub_stat, sam_clo_sub_stat, driver_clo_relate_test)
write_tsv(merge_stat, "result/1.Overview_of_mutation_clonality/mut_clonality_stat.tsv", na = "")


# statistical plots -------------------------------------------------------

# sam_mut_num_barplot

sam_mut_order <- glioma_mut_clonality_data %>%
  count(patient, sort = TRUE)

sam_mut_num_barplot_data <- glioma_mut_clonality_data %>%
  mutate(patient = factor(patient, levels = sam_mut_order$patient), 
         IDH.codel.subtype = factor(IDH.codel.subtype, levels = c("IDHwt", "IDHmut-non-codel", "IDHmut-codel")))

barplot_basic <- ggplot(sam_mut_num_barplot_data, aes(patient, fill = CI95.timing)) +
  geom_bar() +
  labs(x = NULL, y = "Number of mutations", fill = NULL) +
  scale_fill_nejm() +
  theme_classic() +
  scale_x_discrete(expand = c(0.015, 0)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank()) +
  facet_grid(. ~ IDH.codel.subtype, space = "free_x", scales = "free_x")

barplot_part1 <- barplot_basic +
  scale_y_continuous(expand = c(0.015, 0)) +
  coord_cartesian(ylim = c(0, 540)) +
  theme(strip.text = element_blank(), legend.position = "bottom")

barplot_part2 <- barplot_basic +
  coord_cartesian(ylim = c(680, 850)) +
  labs(y = "") +
  scale_y_continuous(breaks = c(600, 800, 100)) +
  theme(axis.line.x = element_blank(), legend.position = "none")

comb_plot <- gridExtra::grid.arrange(barplot_part2, barplot_part1, heights = c(1/5, 4/5), ncol = 1, nrow = 2)

ggsave("result/1.Overview_of_mutation_clonality/sam_mut_num_barplot.pdf", comb_plot, width = 10, height = 4)

# driver_clonality_frac_barplot

driver_clonality_frac_barplot <- glioma_mut_clonality_data %>%
  mutate(driver_mut = factor(driver_mut, levels = c("Yes", "No"))) %>%
  ggplot(aes(driver_mut, fill = CI95.timing)) +
  geom_bar(position = "fill") +
  labs(x = NULL, y = "Proportion of mutations", fill = NULL) +
  scale_fill_manual(values = c("#D45658", "#3A84BB")) +
  scale_x_discrete(breaks = c("Yes", "No"), labels = c("Drivers", "Others")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank()) +
  facet_grid(. ~ IDH.codel.subtype)
ggsave("result/1.Overview_of_mutation_clonality/driver_clonality_frac_barplot.pdf", driver_clonality_frac_barplot, width = 8, height = 7)

# cancer_gene_sub_frac_scatterpie

scatterpei_plot <- function(scatterpie_data, y_axis_data){
  
  plot <- ggplot(scatterpie_data) +
    scatterpie::geom_scatterpie(aes(pie_locus, mean_ccf, r = sqrt(sqrt(sqrt(mut_sam_num))) * 0.015), data = scatterpie_data, cols = c("Clonal", "Subclonal"), alpha = 0.9, color = NA) +
    ggrepel::geom_text_repel(aes(text_locus, mean_ccf, label = Hugo_Symbol), size = 1.2) +
    geom_vline(xintercept = 0) +
    geom_errorbar(aes(x = 0, y = mean_ccf, ymin = mean_ccf, ymax = mean_ccf), width = 0.05) +
    geom_errorbar(aes(x = 0, y = y_scale_value, ymin = y_scale_value, ymax = y_scale_value), width = 0.02, data = y_axis_data) +
    geom_text(aes(x = -0.025, y = y_scale_value, label = y_scale_value), data = y_axis_data, size = 1) +
    labs(x = "Mean of CCF", y = NULL, fill = NULL) +
    scale_fill_nejm() +
    scale_color_nejm() +
    coord_fixed() +
    theme_classic() +
    theme(legend.position = "none", axis.ticks = element_blank(), axis.text = element_blank(), axis.line.y = element_blank(), axis.title.x = element_text(size = 6), strip.background = element_blank(), strip.text = element_text(size = 6)) +
    facet_grid(. ~ IDH.codel.subtype)
  return(plot)
}

cancer_gene_scatterpie_data <- glioma_mut_clonality_data %>%
  filter(!Variant_Classification %in% c("3'Flank", "5'Flank", "3'UTR", "5'UTR", "Intron", "RNA", "Silent"), driver_mut %in% "Yes") %>%
  group_by(IDH.codel.subtype, Hugo_Symbol) %>%
  summarise(mut_sam_num = n_distinct(patient), mean_ccf = mean(absolute.ccf)) %>%
  arrange(desc(mean_ccf), .by_group = TRUE) %>%
  left_join(sam_stat, by = "IDH.codel.subtype") %>%
  mutate(mut_sam_frac = mut_sam_num / clonality_sample_num, x_locus = 1, pie_locus = x_locus * c(0.1, -0.1), text_locus = x_locus * c(0.25, -0.25)) %>%
  left_join(cancer_gene_sub_stat, by = c("IDH.codel.subtype", "Hugo_Symbol"))

y_axis_data <- tibble(y_scale_value = seq(0.3, 1, by = 0.1))

cancer_gene_scatterpie_selected_data <- cancer_gene_scatterpie_data %>%
  filter(mut_sam_frac > 0.02, mean_ccf > 0.6)

y_axis_selected_data <- tibble(y_scale_value = seq(0.6, 1, by = 0.1))

cancer_gene_scatterpie <- scatterpei_plot(cancer_gene_scatterpie_data, y_axis_data)
cancer_gene_scatterpie_selected <- scatterpei_plot(cancer_gene_scatterpie_selected_data, y_axis_selected_data)

ggsave("result/1.Overview_of_mutation_clonality/cancer_gene_sub_frac_scatterpie.pdf", cancer_gene_scatterpie)
ggsave("result/1.Overview_of_mutation_clonality/cancer_gene_sub_frac_scatterpie_selected.pdf", cancer_gene_scatterpie_selected)

# cancer gene mutation ccf probability density plot

absolute.cancer.cell.fraction <- function(n.alt, depth, purity, local.copy.number){
  f.function <- function (c,purity,local.copy.number){
    return((purity*c) / (2*(1-purity) + purity*local.copy.number))
  }
  x <- dbinom(n.alt,depth, prob=sapply(seq(0.01,1,length.out=100),f.function,purity,local.copy.number))
  if (min(x) == 0) {
    x[length(x)] <- 1
  }
  x.mat <- data.frame(ccf = seq(0.01,1,length.out=100), ccf_prob = x)
  return(x.mat)
}

interested_mutation <- c("TCGA-CS-4943-01:2:209113112:C", "TCGA-DB-A64U-01:2:209113112:C", "TCGA-DU-8161-01:10:89692904:C", "TCGA-S9-A6TW-01:10:89717661:C")

ccf_density_data_list <- panCanAtlas.abs.seg %>%
  distinct(SampleID, purity = `Aberrant Cell Fraction`) %>%
  right_join(glioma_mut_clonality_data, by = c("SampleID" = "patient")) %>%
  mutate(cn = minor_cn + major_cn, depth = ref_counts + var_counts) %>%
  filter(mutation_id %in% interested_mutation) %>%
  split(.$mutation_id) %>%
  map(~absolute.cancer.cell.fraction(.$var_counts, .$depth, .$purity, .$cn))

ccf_density_plot_data <- names(ccf_density_data_list) %>%
  map(~mutate(ccf_density_data_list[[.]], mutation_id = .)) %>%
  bind_rows() %>%
  left_join(glioma_mut_clonality_data, by = "mutation_id") %>%
  mutate(mutation_id = factor(mutation_id, levels = interested_mutation), sam_gene = str_c(bcr_patient_barcode, Hugo_Symbol, sep = " ")) %>%
  select(mutation_id:ccf, sam_gene, CI95.timing)  

ccf_density_plot <- ggplot(ccf_density_plot_data, aes(x = ccf, y = ccf_prob)) +
  geom_line() +
  geom_area(aes(fill = CI95.timing), show.legend = FALSE) +
  labs(x = "Cancer Cell Fraction", y = "") +
  scale_x_continuous(expand = c(0.0065, 0.0065)) +
  scale_fill_nejm() +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), strip.background = element_blank(), strip.text = element_text(size = 12)) +
  facet_wrap(~ sam_gene, ncol = 1, scales = "free_y")
ggsave("result/1.Overview_of_mutation_clonality/ccf_density_plot.pdf", ccf_density_plot, width = 10, height = 6)





