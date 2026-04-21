suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(EnhancedVolcano)
  library(dplyr)
  library(ggrepel)
})

load_counts <- function(counts_file) {
  raw <- read.table(counts_file, header=TRUE, sep="\t",
                    skip=1, row.names=1, check.names=FALSE)
  counts <- raw[, 6:ncol(raw)]
  colnames(counts) <- gsub("results/bam/", "", colnames(counts))
  colnames(counts) <- gsub("\\.bam", "", colnames(counts))
  counts <- as.matrix(counts)
  mode(counts) <- "integer"
  counts <- counts[rowSums(counts) > 0, ]
  return(counts)
}

build_coldata <- function(counts, ref, treat) {
  condition <- ifelse(grepl("^tumor", colnames(counts)), treat, ref)
  data.frame(
    sample    = colnames(counts),
    condition = factor(condition, levels=c(ref, treat)),
    row.names = colnames(counts)
  )
}

run_deseq2 <- function(counts, coldata) {
  dds <- DESeqDataSetFromMatrix(counts, coldata, ~condition)
  DESeq(dds)
}

get_results <- function(dds, ref, treat, padj_cut, lfc_cut) {
  res <- results(dds, contrast=c("condition", treat, ref), alpha=padj_cut)
  res <- lfcShrink(dds, contrast=c("condition", treat, ref),
                   res=res, type="ashr")
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(padj) %>%
    mutate(
      significant = !is.na(padj) & padj < padj_cut & abs(log2FoldChange) > lfc_cut,
      direction   = case_when(
        significant & log2FoldChange > 0 ~ "UP",
        significant & log2FoldChange < 0 ~ "DOWN",
        TRUE ~ "NS"
      )
    )
  list(res=res, res_df=res_df)
}

plot_pca <- function(dds) {
  vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
  pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color = condition, label = name)) +
    geom_point(size = 4) +
    ggrepel::geom_text_repel(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_color_manual(values = c("normal" = "#4f9eff", "tumor" = "#ff6b6b")) +
    theme_minimal(base_size = 13) +
    labs(title = "PCA — Sample Clustering")
}

plot_volcano <- function(res_df, padj_cut, lfc_cut) {
  EnhancedVolcano(res_df,
    lab=res_df$gene, x="log2FoldChange", y="padj",
    pCutoff=padj_cut, FCcutoff=lfc_cut,
    title="Tumor vs Normal — Volcano Plot",
    col=c("grey70","grey70","#4f9eff","#ff6b6b"),
    legendPosition="bottom")
}

plot_heatmap <- function(dds, res_df, n=50) {
  vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
  top_genes <- res_df %>% filter(significant) %>%
    slice_min(padj, n=n) %>% pull(gene)
  mat <- assay(vsd)[top_genes, ]
  mat <- mat - rowMeans(mat)
  anno <- as.data.frame(colData(dds)[,"condition",drop=FALSE])
  pheatmap(mat, annotation_col=anno,
    color=colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
    main=paste("Top", n, "DEGs"))
}

summarise_de <- function(res_df, padj_cut, lfc_cut) {
  sig <- res_df %>% filter(significant)
  list(
    total_tested = nrow(res_df),
    total_sig    = nrow(sig),
    up           = sum(sig$direction=="UP"),
    down         = sum(sig$direction=="DOWN"),
    top10_up     = sig %>% filter(direction=="UP") %>%
                   slice_min(padj, n=10) %>%
                   select(gene, log2FoldChange, padj),
    top10_down   = sig %>% filter(direction=="DOWN") %>%
                   slice_min(padj, n=10) %>%
                   select(gene, log2FoldChange, padj)
  )
}
