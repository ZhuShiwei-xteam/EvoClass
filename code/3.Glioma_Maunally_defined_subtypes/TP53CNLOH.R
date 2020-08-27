setwd('/pub5/xiaoyun/Jobs/J22/temp/Github')

# The proportion of TP53 mutations had copy number neutral LOH.

# Mutational clonality data
load("data/glioma.mut.ccf.RData")
load("data/gliomaClinicalData.RData")
com.data$patient <- substr(com.data$patient, 1, 12)
glioma.mut.ccf <- com.data
tp53.mut.ccf <- subset(glioma.mut.ccf, Hugo_Symbol=='TP53')
tp53.mut.ccf <- merge(x = tp53.mut.ccf, y = pan.glio.mer.cli.data, by.x = 'patient', by.y = 'bcr_patient_barcode')

# > table(tp53.mut.ccf[, c('minor_cn', 'major_cn')]) 
        # major_cn
# minor_cn   1   2   3   4   5   6
       # 0  24  85  11  41   8   2
       # 1 179   4   4   2   1   0
       # 2   0  90   2   0   0   0
# 85/453 = 18.8%

table(subset(tp53.mut.ccf, IDH.codel.subtype == 'IDHmut-codel')[, c('minor_cn', 'major_cn')])
        # major_cn
# minor_cn 1 2
       # 1 4 0
       # 2 0 3

table(subset(tp53.mut.ccf, IDH.codel.subtype == 'IDHmut-non-codel')[, c('minor_cn', 'major_cn')])
        # major_cn
# minor_cn   1   2   3   4   5   6
       # 0   2  56   7  36   8   2
       # 1 128   3   1   2   1   0
       # 2   0  67   0   0   0   0
# 56/313 = 17.8%

table(subset(tp53.mut.ccf, IDH.codel.subtype == 'IDHwt')[, c('minor_cn', 'major_cn')])
        # major_cn
# minor_cn  1  2  3  4
       # 0 22 29  4  5
       # 1 47  1  3  0
       # 2  0 20  2  0
# 29/133 = 21.8%

