# to run locally:
# SCRATCH1=/mnt/h/scratch
# docker run -it --rm \
# --workdir $HOME \
# -v /mnt/h/scratch/gse263155:$HOME/project \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/h/scratch" \
# -v /mnt/h/scratch:$HOME/scratch \
# p4rkerw/sctools:R4.3.2e R

# run in rstudio server with NAS mount
# SCRATCH1=/mnt/h/scratch
# workdir=/home/rstudio
# docker run -it --rm \
# -p 8888:8787 \
# -e PASSWORD=password \
# -v /mnt/h/scratch/gse263155:$workdir/project \
# -v /mnt/g/reference:$workdir/reference \
# -v $workdir:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/h/scratch" \
# -v /mnt/h/scratch:$workdir/scratch \
# -e DISABLE_AUTH=true \
# p4rkerw/sctools:R4.3.2e
# # navigate browser to localhost:8888

library(tximport)
library(readr)
library(DESeq2)
library(tidyverse)
library(stringr)
library(biomaRt)
library(apeglm)
library(data.table)
library(dplyr)
library(tibble)


# load files
files <- list.files('project', pattern = "quant.sf", recursive=TRUE, full.names = TRUE)

# Load metadata
samples <- str_split(files, pattern = "/", simplify=TRUE)[,2]
samples.df <- data.frame(samples = samples, path = files)

metadata <- read.csv('project/SraRunTable.csv') %>%
  dplyr::filter(Run %in% samples) %>%
  dplyr::arrange(samples) %>%
  dplyr::filter(treatment %in% c("Vehicle","Empagliflozin"))

samples.df <- samples.df[samples %in% metadata$Run,]
files <- files[samples %in% metadata$Run]

coldata <- data.frame(sample = samples.df$samples, genotype = as.factor(metadata$genotype), treatment = as.factor(metadata$treatment))

# Connect to Ensembl BioMart for mouse (Mus musculus)
# ensembl <- useEnsembl(biomart = "ensembl", dataset = "rnorvegicus_gene_ensembl")

ensembl <- useEnsembl(biomart = "ensembl", dataset = "rnorvegicus_gene_ensembl", mirror = "useast")

# Retrieve transcript IDs and corresponding gene IDs
tx2gene <- getBM(attributes = c("ensembl_transcript_id", "mgi_symbol"),
                 mart = ensembl)

# Check the structure
head(tx2gene)

# Load tx2gene mapping
# tx2gene <- read_csv("tx2gene.csv")

# Prepare named vector of file paths
# files <- setNames(samples$path, samples$sample)

# Import transcript abundance and summarize to gene level
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxV = TRUE)

# Create DESeq2 dataset
dds <- DESeqDataSetFromTximport(txi,
                                colData = coldata,
                                design = ~ genotype + treatment)

# Filter low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run differential expression
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c("genotype", "ZSF1 Obese", "ZSF1 Lean"))

# Shrink log fold changes for plotting
resLFC <- lfcShrink(dds, coef="genotype_ZSF1.Obese_vs_ZSF1.Lean", type="apeglm")

# Order by adjusted p-value
resOrdered <- resLFC[order(resLFC$padj), ]

# Save results
write.csv(as.data.frame(resOrdered), file="project/deseq2_results.csv")

# Get results
res <- results(dds, contrast = c("treatment", "Empagliflozin", "Vehicle"))

# Shrink log fold changes for plotting
resLFC <- lfcShrink(dds, coef="treatment_Empagliflozin_vs_Vehicle", type="apeglm")

# Order by adjusted p-value
resOrdered <- resLFC[order(resLFC$padj), ]

# Save results
write.csv(as.data.frame(resOrdered), file="project/deseq2_results.csv")





# Basic volcano plot
res_df <- as.data.frame(resOrdered) %>%
  rownames_to_column("gene") %>%
  mutate(sig = padj < 0.05 & abs(log2FoldChange) > 1)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10(padj)") +
  scale_color_manual(values = c("grey", "red"))

gene <- "Dclk1"

# Extract normalized counts
norm_counts <- counts(dds, normalized = TRUE)

# Build dataframe for plotting
plot_df <- data.frame(
  expression = norm_counts[gene, ],
  condition = colData(dds)$is_adenine,
  sample = colData(dds)$sample
)

# Optional: log-transform for better visualization
plot_df$log_expression <- log2(plot_df$expression + 1)

# Boxplot + jitter
ggplot(plot_df, aes(x = sample, y = log_expression, fill = is_adenine)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  labs(title = paste("Expression of", gene),
       y = "Log2(Normalized counts + 1)",
       x = "Condition") +
  theme_minimal() +
  theme(legend.position = "none")


# Regularized log transformation (good for small datasets)
rld <- rlog(dds, blind = TRUE)

# OR: Variance Stabilizing Transformation (faster, good for larger datasets)
vsd <- vst(dds, blind = TRUE)

library(ggplot2)
library(dplyr)

# Choose gene ID (e.g., from rownames(dds))
gene <- "Dclk1"  # Replace with your gene

# Use rlog or vsd object
expr_rlog <- assay(rld)[gene, ]
expr_vst  <- assay(vsd)[gene, ]

# Build plotting data frame (choose one)
plot_df <- data.frame(
  expression = expr_rlog,         # or expr_vst
  sample = colData(dds)$sample,
  condition = colData(dds)$is_adenine
)

# Plot: boxplot + jitter
ggplot(plot_df, aes(x = sample, y = expression, fill = is_adenine)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  labs(title = paste("Expression of", gene, "(rlog transformed)"),
       y = "rlog(expression)",
       x = "Condition") +
  theme_minimal() +
  theme(legend.position = "none")

install.packages('msigdbr')

library(DESeq2)
library(fgsea)
library(msigdbr)
library(dplyr)
library(tibble)
library(ggplot2)

# Clean DESeq2 results
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj))

# Use log2FoldChange as ranking metric (can also use stat or signed -log10(padj))
ranking <- res_df %>%
  arrange(desc(log2FoldChange)) %>%
  dplyr::select(gene, log2FoldChange) %>%
  deframe()  # Named vector: names = genes, values = ranking score

# Get gene sets (e.g., GO, KEGG, Reactome) for mouse
msigdb <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")  # C5 = GO, H = Hallmark

# Prepare gene set list
gene_sets <- msigdb %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)

set.seed(42)
gsea_res <- fgsea(pathways = gene_sets,
                  stats    = ranking,
                  minSize  = 15,
                  maxSize  = 500,
                  nperm    = 10000)  # Increase for stable p-values

gsea_res %>%
  arrange(padj) %>%
  head(10) %>%
  ggplot(aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = NES > 0)) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score (NES)", title = "Top Enriched Pathways") +
  theme_minimal()

library(clusterProfiler)

BiocManager::install('org.Mm.eg.db')
library(org.Mm.eg.db)

# GSEA using GO Biological Process
gsea_cp <- GSEA(ranking,
                TERM2GENE = msigdb %>% dplyr::select(gs_name, gene_symbol),
                minGSSize = 15,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                verbose = TRUE)

dotplot(gsea_cp, showCategory = 10)  # Top pathways


