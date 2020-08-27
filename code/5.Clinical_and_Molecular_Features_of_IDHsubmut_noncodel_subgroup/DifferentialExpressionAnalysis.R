setwd('/pub5/xiaoyun/Jobs/J22/temp/Github')
# load library
library("AnnotationDbi")
library("org.Hs.eg.db")
library("DESeq2")

# RNA-SEQ raw counts
setwd(paste(getwd(), '/data/OriginalData/', sep=""))
# utils::untar("data/OriginalData/gdc_download_20200721_025916.394623.tar.gz", exdir = ".")

glioma.samples <- read.csv(file = 'gdc_sample_sheet.2020-07-21.tsv', sep = '\t', header = TRUE, stringsAsFactors = FALSE)

glioma.samples$File.Name <- paste(glioma.samples$File.ID, glioma.samples$File.Name, sep = "/")
glioma.samples <- subset(glioma.samples, Sample.Type == 'Primary Tumor')

glioma.samples <- glioma.samples[!duplicated(glioma.samples$Case.ID), ]
rownames(glioma.samples) <- glioma.samples$Case.ID

# molecular subtype
load('data/GliomaMolecularSubtype.RData')

sample.data <- merge(glioma.samples, glioma.molecular.subtype, by = 'row.names')
sample.data <- sample.data[, c('Case.ID', 'File.Name', 'final.subtype', 'Project.ID')]
colnames(sample.data) <- c('Case.ID', 'File.Name', 'Molecular.Subtype', 'Project.ID')


# Differential expression analysis
sample.data <- subset(sample.data, Molecular.Subtype %in% c('IDHcolmut.non.codel', 'IDHsubmut.non.codel'))
sample.data$Molecular.Subtype <- factor(sample.data$Molecular.Subtype, levels = c('IDHcolmut.non.codel', 'IDHsubmut.non.codel'))
sample.data$Project.ID <- factor(sample.data$Project.ID, levels = c('TCGA-LGG', 'TCGA-GBM'))

dds.htsc <- DESeqDataSetFromHTSeqCount(sampleTable = sample.data, directory = ".", design = ~ Project.ID + Molecular.Subtype)

# Genecode gene ID
setwd('/pub5/xiaoyun/Jobs/J22/temp/Github')
gencode.gtf <- read.csv(file = 'data/OriginalData/gencode.v22.annotation.gtf', header = FALSE, comment.char = "#", sep = '\t', stringsAsFactors = FALSE)
gencode.gtf <- subset(gencode.gtf, V3 == 'gene')

protein.coding.gene <- sapply(strsplit(grep('gene_type protein_coding', gencode.gtf$V9, value = TRUE), split = ';'), function(x) x[1])
protein.coding.gene <- substr(protein.coding.gene, 9, 23)

rownames(dds.htsc) <- substr(rownames(dds.htsc), 1, 15)
dds.htsc <- dds.htsc[intersect(rownames(dds.htsc), protein.coding.gene), ]

save(dds.htsc, file = 'result/5.Clinical_and_Molecular_Features_of_IDHsubmut_noncodel_subgroup/DESeqDataSetFromHTSeqCount.RData')
# Pre-filtering the dataset
keep <- rowSums(counts(dds.htsc) == 0) <= ncol(dds.htsc)*0.5
dds.htsc.filter <- dds.htsc[keep, ]

# Running the differential expression pipeline
dds <- DESeq(dds.htsc.filter, parallel=T)

res <- results(dds, contrast=c("Molecular.Subtype", "IDHsubmut.non.codel", "IDHcolmut.non.codel"), alpha = 0.1)

# save data
save(res, file = 'result/5.Clinical_and_Molecular_Features_of_IDHsubmut_noncodel_subgroup/DESeq2DiffExp.RData')


res <- as.data.frame(res)
res$symbol <- mapIds(org.Hs.eg.db, key=rownames(res), column="SYMBOL", keytype="ENSEMBL", multiVals="first")

#pdf(file = 'result/5.Clinical_and_Molecular_Features_of_IDHsubmut_noncodel_subgroup/DEVolcanoPlot.pdf')
# Make a basic volcano plot
#with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,4)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
#with(subset(res, padj<0.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
#with(subset(res, padj<0.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#with(subset(res, padj<0.05 & abs(log2FoldChange)>1), text(log2FoldChange, y = -log10(pvalue), labels = symbol))
#dev.off()

pdf(file = 'result/5.Clinical_and_Molecular_Features_of_IDHsubmut_noncodel_subgroup/DEVolcanoPlot4.pdf')
# Make a basic volcano plot
with(subset(res, padj<0.05 & log2FoldChange > 1), plot(log2FoldChange, -log10(padj), 
 pch=20, main="Volcano plot", xlim=c(1, 4), ylim = c(1, 12)))
with(subset(res, padj<0.05 & log2FoldChange > 1), 
 text(log2FoldChange, y = -log10(padj), labels = symbol, cex = 0.4))
dev.off()





# differentially expressed genes
diff.exp <- as.data.frame(subset(res, (abs(log2FoldChange) > 1)&(padj < 0.05)))

ens.anno <- select(org.Hs.eg.db, keys = rownames(diff.exp), columns = c("GENENAME", "ENTREZID"), keytype='ENSEMBL')
diff.exp <- cbind(diff.exp,  ens.anno)
diff.exp <- diff.exp[order(diff.exp$log2FoldChange), ]

write.table(diff.exp, file = "result/5.Clinical_and_Molecular_Features_of_IDHsubmut_noncodel_subgroup/DiffExpGene.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# preranked gene list for performing GSEA
res.ranks <- res
res.ranks <- res.ranks[!is.na(res.ranks$pvalue), ]
res.ranks$ranks <- sign(res.ranks$log2FoldChange) * -log10(res.ranks$pvalue)


res.ranks <- res.ranks[, c('symbol', 'ranks')]
res.ranks <- res.ranks[!is.na(res.ranks$symbol), ]

write.table(res.ranks, file = "result/5.Clinical_and_Molecular_Features_of_IDHsubmut_noncodel_subgroup/RankedGeneList.rnk", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# gsea result
gsea.res <- read.table(file = 'data/OriginalData/PrerankedGSEAResults.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)

library('ggplot2')

bar.plot <- ggplot(subset(gsea.res, FDR.q.val < 0.05), aes(reorder(NAME, NES), NES)) + geom_col(aes(fill = FDR.q.val < 0.01)) + 
 coord_flip() + labs(x='Pathway', y='Normalized Enrichment Score', title="Hallmark pathways NES from GSEA") + theme_minimal()

ggsave(filename = 'result/5.Clinical_and_Molecular_Features_of_IDHsubmut_noncodel_subgroup/GSEAResBarplot.pdf',plot = bar.plot)

