library(DESeq2)
library(ggplot2)
library(pheatmap)
library(apeglm)

#Load count data

counts <- read.table("D:/feature_counts/gene_counts.clean.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE)

colnames(counts)


metadata <- read.delim(
  "D:/metadata/Samples.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

#Check dimensions of the data
dim(metadata)
colnames(metadata)

#Check SRR_ID col
table(metadata$SRR_ID)

#Remove empty rows and duplicate Ids
metadata <- metadata[metadata$SRR_ID != "", ]
metadata <- metadata[!is.na(metadata$SRR_ID), ]
metadata <- metadata[!duplicated(metadata$SRR_ID), ]

#check dim again
dim(metadata)

#Trim spaces
metadata$SRR_ID <- trimws(metadata$SRR_ID)

#Set rownames to SRR_IDS
rownames(metadata) <- metadata$SRR_ID

#Set rownames==colnames
all(rownames(metadata) == colnames(counts))

#Create DESEQ2 object

metadata$subtype <- factor(metadata$subtype, levels =c("Normal","TNBC","NonTNBC","HER2"))
metadata$subtype <- relevel(metadata$subtype, ref="Normal")

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ subtype
)

#Filter protein_coding 
library(rtracklayer)

gtf <- import("D:/4.reference/gencode.v44.annotation.gtf")

protein_coding <- unique(
  gtf$gene_id[gtf$gene_type == "protein_coding"]
)

dds_pc <- dds[rownames(dds) %in% protein_coding, ]

nrow(dds_pc)
dds_pc <- DESeq(dds_pc)


#filter low count genes

dds_pc <- dds_pc[rowSums(counts(dds_pc)) >= 20, ]
dds_pc <- DESeq(dds_pc)
vsd <- vst(dds_pc, blind = FALSE)
resultsNames(dds_pc)

#Extract results for each subtype

#TNBC VS Normal

res_TNBC <- lfcShrink(
  dds_pc,
  coef = "subtype_TNBC_vs_Normal",
  type = "apeglm"
)

# NonTNBC vs Normal
res_NonTNBC <- lfcShrink(dds_pc,
                         coef = "subtype_NonTNBC_vs_Normal",
                         type = "apeglm")

# HER2 vs Normal
res_HER2 <- lfcShrink(dds_pc,
                      coef = "subtype_HER2_vs_Normal",
                      type = "apeglm")

#Filter significant genes (some genes have padj=NA , to remove them)

sig_TNBC <- subset(res_TNBC,
                   !is.na(padj) &
                     baseMean > 10 &
                     padj < 0.05 &
                     abs(log2FoldChange) > 1)


sig_NonTNBC <- subset(res_NonTNBC,
                      !is.na(padj) &
                        baseMean > 10 &
                        padj < 0.05 &
                        abs(log2FoldChange) > 1)

sig_HER2 <- subset(res_HER2,
                   !is.na(padj) &
                     baseMean > 10 &
                     padj < 0.05 &
                     abs(log2FoldChange) > 1)

#Save as csv 

write.csv(as.data.frame(sig_TNBC),
          "DEG_TNBC_vs_Normal.csv",
          row.names = TRUE)

write.csv(as.data.frame(sig_NonTNBC),
          "DEG_NonTNBC_vs_Normal.csv",
          row.names = TRUE)

write.csv(as.data.frame(sig_HER2),
          "DEG_HER2_vs_Normal.csv",
          row.names = TRUE)


#PCA Plot

pca <- plotPCA(vsd, intgroup = "subtype", returnData = TRUE)

percentVar <- round(100 * attr(pca, "percentVar"))

ggplot(pca, aes(PC1, PC2, color = subtype)) +
  geom_point(size = 4) +
  stat_ellipse(level = 0.95) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  theme_classic()

table(metadata$subtype)

#Heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

pheatmap(sampleDistMatrix)

#Volcano plots

library(ggplot2)

create_volcano <- function(res, title){
  
  res_df <- as.data.frame(res)
  res_df <- res_df[!is.na(res_df$padj), ]
  
  res_df$Significance <- "Not Significant"
  res_df$Significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
  res_df$Significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"
  
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = Significance), alpha = 0.7, size = 1.8) +
    scale_color_manual(values = c(
      "Downregulated" = "blue",
      "Not Significant" = "grey70",
      "Upregulated" = "red"
    )) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_classic() +
    labs(title = title,
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

create_volcano(res_TNBC, "TNBC vs Normal")
create_volcano(res_NonTNBC, "NonTNBC vs Normal")
create_volcano(res_HER2, "HER2 vs Normal")

sum(res_HER2$padj < 0.05, na.rm=TRUE)

sum(sig_HER2$log2FoldChange > 1)
sum(sig_HER2$log2FoldChange < -1)


sum(res_TNBC$padj < 0.05, na.rm=TRUE)

sum(sig_TNBC$log2FoldChange > 1)
sum(sig_TNBC$log2FoldChange < -1)



sum(res_NonTNBC$padj < 0.05, na.rm=TRUE)

sum(sig_NonTNBC$log2FoldChange > 1)
sum(sig_NonTNBC$log2FoldChange < -1)

#TNBC specific genes
tnbc_genes <- rownames(sig_TNBC)
her2_genes <- rownames(sig_HER2)
nontnbc_genes <- rownames(sig_NonTNBC)

tnbc_specific <- setdiff(tnbc_genes,
                         union(her2_genes, nontnbc_genes))

length(tnbc_specific)





