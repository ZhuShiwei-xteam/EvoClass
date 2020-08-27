setwd("/pub5/xiaoyun/Jobs/J22/temp/Github")

suppressPackageStartupMessages(library(tidyverse))
library(ggsci)


# load mutation clonality data
load("data/glioma.mut.ccf.RData")
glioma_mut_clonality_data <- com.data
# load cancer gene data
load("data/pan.glioma.cancer.gene.RData")


# non-silent driver clonality data
glioma_dri_clonality_data <- glioma_mut_clonality_data %>%
  filter(Hugo_Symbol %in% pan.glioma.cancer.gene, !Variant_Classification %in% c("3'Flank", "5'Flank", "3'UTR", "5'UTR", "Intron", "RNA", "Silent")) %>%
  mutate(ccf.0.8 = ifelse(absolute.ccf > 0.8, "Clonal", "Subclonal"))

test_data <- glioma_dri_clonality_data %>%
  mutate(Hugo_Symbol = "Drivers") %>%
  bind_rows(glioma_dri_clonality_data) %>%
  pivot_longer(cols = c(CI95.timing, prob.clonal.timing, ccf.0.8), names_to = "clonality_method", values_to = "clonality")

# CI95.timing vs prob.clonal.timing
CI95_prob_test_result <- test_data %>%
  filter(clonality_method %in% c("CI95.timing", "prob.clonal.timing"), Hugo_Symbol %in% c("Drivers", "IDH1", "PIK3CA", "TP53")) %>%
  split(.$Hugo_Symbol) %>%
  map_dbl(~fisher.test(.$clonality_method, .$clonality)$p.value)
# Drivers(p = 7.512881e-52), IDH1(p = 1.168346e-19), PIK3CA(p = 6.301870e-06), TP53(p = 5.297441e-12) 

# CI95.timing vs ccf.0.8
CI95_ccf0.8_data <- test_data %>% filter(clonality_method %in% c("CI95.timing", "ccf.0.8"), Hugo_Symbol %in% "Drivers")

CI95_ccf0.8_test <- fisher.test(CI95_ccf0.8_data$clonality_method, CI95_ccf0.8_data$clonality)$p.value # p = 0.2449262
CI95_ccf0.8_diff_mut <- mean(glioma_dri_clonality_data$CI95.timing != glioma_dri_clonality_data$ccf.0.8) # 0.05678233



# plot --------------------------------------------------------------------

# driver clonal fraction barplot

clonality_method_frac_barplot <- ggplot(filter(test_data, !clonality_method %in% "ccf.0.8", Hugo_Symbol %in% c("Drivers", "IDH1", "PIK3CA", "TP53")), aes(clonality_method, fill = clonality)) +
  geom_bar(position = "fill") +
  labs(x = NULL, y = "Proportion of mutations", fill = NULL) +
  scale_fill_manual(values = c("#D45658", "#3A84BB")) +
  scale_x_discrete(breaks = c("CI95.timing", "prob.clonal.timing"), labels = c("CI95 based", "Prob based")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
  facet_grid(. ~ Hugo_Symbol)
ggsave("result/1.Overview_of_mutation_clonality/clonality_method_frac_barplot.pdf", clonality_method_frac_barplot, width = 6, height = 5)

# ccf confidence interval forest plot

clonality_method_forest_plot <- glioma_dri_clonality_data %>%
  mutate(CI95_ccf0.8_compare = ifelse(CI95.timing == ccf.0.8, CI95.timing, "inconsistent")) %>%
  ggplot(aes(fct_reorder(mutation_id, absolute.ccf.0.95), ymin = absolute.ccf.0.05, ymax = absolute.ccf.0.95)) +
  geom_linerange(size = 0.05, aes(col = CI95.timing)) + 
  geom_hline(aes(yintercept = 0.8), colour="#BB0000", linetype="dashed") +
  geom_point(aes(y = absolute.ccf, col = CI95_ccf0.8_compare)) +
  annotate("text", x = 20, y = 0.85, hjust = 0, vjust = 1, label = "CCF = 0.8") +
  scale_color_manual(values = c("#D45658", "#C77CFF", "#3A84BB")) +
  labs(x = NULL, y = "Cancer cell fraction", col = NULL) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(), legend.position = c(1, 0), legend.justification = c(1, 0), legend.background = element_rect(fill = NA), legend.key = element_rect(fill = NA))
ggsave("result/1.Overview_of_mutation_clonality/clonality_method_forest_plot.pdf", clonality_method_forest_plot, width = 6, height = 5)
















