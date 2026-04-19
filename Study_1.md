library(edgeR)

# -------------------------
# 1️⃣ Set working directory and load filtered counts
# -------------------------
setwd("/Users/viraa/TestingCode/Testing_Monkey")
dge_filtered <- read.csv("filtered_miRNA_counts.csv", row.names = 1, check.names = FALSE)

# Original counts column names
colnames(dge_filtered)

meta <- read.csv("sample_metadata.csv", stringsAsFactors = FALSE)

meta$Sample_ID <- trimws(meta$Sample_ID)
meta$Timepoint <- factor(trimws(meta$Timepoint), levels = c("baseline", "post"))
meta$Sex <- factor(trimws(meta$Sex), levels = c("F", "M"))
meta$Category <- factor(trimws(meta$Category), levels = c("LD", "BD", "HD"))

common_samples <- intersect(colnames(dge_filtered), meta$Sample_ID)

counts_filtered <- dge_filtered[, colnames(dge_filtered) %in% meta$Sample_ID]
meta_filtered <- meta[meta$Sample_ID %in% common_samples, ]

# Reorder metadata to match counts
meta_filtered <- meta_filtered[match(colnames(counts_filtered), meta_filtered$Sample_ID), ]

# Check that all sample IDs match
if(!all(colnames(counts_filtered) == meta_filtered$Sample_ID)) {
  stop("Sample IDs in counts and metadata do not match!")
}

baseline_samples <- meta_filtered$Sample_ID[meta_filtered$Timepoint == "baseline"]

if(length(baseline_samples) == 0) stop("No baseline samples found. Check Timepoint column!")

counts_baseline <- counts_filtered[, baseline_samples]
meta_baseline <- meta_filtered[meta_filtered$Sample_ID %in% baseline_samples, ]

dge_all <- DGEList(counts = counts_filtered)
dge_all <- calcNormFactors(dge_all, method = "TMM")

CPM_all <- cpm(dge_all, log = FALSE)
write.csv(CPM_all, "Normalized_CPM_all_samples.csv")

logCPM_all <- cpm(dge_all, log = TRUE)
write.csv(logCPM_all, "Normalized_logCPM_all_samples.csv")

baseline_samples <- meta_filtered$Sample_ID[meta_filtered$Timepoint == "baseline"]
post_samples     <- meta_filtered$Sample_ID[meta_filtered$Timepoint == "post"]

logCPM_baseline <- logCPM_all[, baseline_samples]
logCPM_post     <- logCPM_all[, post_samples]

write.csv(logCPM_baseline, "Normalized_logCPM_baseline.csv")
write.csv(logCPM_post,     "Normalized_logCPM_post.csv")

meta_baseline$Sex <- factor(meta_baseline$Sex, levels = c("F", "M"))
meta_baseline$Category <- factor(
  meta_baseline$Category,
  levels = c("LD", "BD", "HD")
)

design <- model.matrix(~ Sex, data = meta_baseline)

dge_baseline <- dge_all[, baseline_samples]

meta_baseline <- meta_filtered[
  match(colnames(dge_baseline), meta_filtered$Sample_ID),
]

meta_baseline$Sex <- droplevels(factor(meta_baseline$Sex, levels = c("F", "M")))

design <- model.matrix(~ Sex, data = meta_baseline)

stopifnot(nrow(design) == ncol(dge_baseline))

dge_baseline <- estimateDisp(dge_baseline, design)

# Fit GLM
fit <- glmQLFit(dge_baseline, design)

# Test for sex effect (assuming 'SexM vs F' is coef 2)
qlf <- glmQLFTest(fit, coef = 2)

# Look at top DE miRNAs
topTags(qlf)


pca_input <- t(logCPM_baseline)

pca_res <- prcomp(pca_input, scale. = TRUE)
library(ggplot2)

# Combine PCA scores with metadata
pca_df <- data.frame(
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  Sample = rownames(pca_res$x),
  Sex = meta_baseline$Sex,
  Category = meta_baseline$Category
)

p <- ggplot(pca_df, aes(x = PC1, y = PC2,
                        color = Category,
                        shape = Sex,
                        label = NA)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.6, size = 3) +
  xlab(paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100, 1), "%)")) +
  ylab(paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100, 1), "%)")) +
  theme_minimal() +
  ggtitle("PCA of Baseline miRNAs (logCPM)\nColored by Drinking Category, Shaped by Sex")
p
ggsave("PCA_plot.png", p, width = 8, height = 6, dpi = 300)

meta_baseline$Sample_ID <- trimws(as.character(meta_baseline$Sample_ID))
common <- intersect(colnames(dge_baseline), meta_baseline$Sample_ID)

dge_baseline  <- dge_baseline[, common, keep.lib.sizes = FALSE]
meta_baseline <- meta_baseline[meta_baseline$Sample_ID %in% common, ]

## 4. Align ONCE using rownames (robust method)
rownames(meta_baseline) <- meta_baseline$Sample_ID
meta_baseline <- meta_baseline[colnames(dge_baseline), , drop = FALSE]

stopifnot(identical(rownames(meta_baseline), colnames(dge_baseline)))

## 5. Set factors cleanly
meta_baseline$Category <- droplevels(factor(
  meta_baseline$Category,
  levels = c("LD","BD","HD")
))
meta_baseline$Sex <- droplevels(factor(
  meta_baseline$Sex,
  levels = c("F","M")
))

# Remove VHD sample (n=1, excluded from analysis)
meta_baseline <- meta_baseline[!is.na(meta_baseline$Category), ]
common        <- intersect(colnames(dge_baseline), meta_baseline$Sample_ID)
dge_baseline  <- dge_baseline[, common, keep.lib.sizes = FALSE]
meta_baseline <- meta_baseline[meta_baseline$Sample_ID %in% common, ]
design <- model.matrix(~ Category + Sex, data = meta_baseline)

dge_baseline <- estimateDisp(dge_baseline, design)
fit <- glmQLFit(dge_baseline, design)


run_contrast <- function(fit, design, contrast_str, label) {
  contrast_vec <- makeContrasts(contrasts = contrast_str, levels = design)
  qlf <- glmQLFTest(fit, contrast = contrast_vec)
  tt <- topTags(qlf, n = Inf)$table
  write.csv(tt, file = file.path("DE_outputs", paste0("DE_", label, ".csv")))
  return(tt)
}



results_list <- list(
  BD_vs_LD  = run_contrast(fit, design, "CategoryBD",                    "BD_vs_LD"),
  HD_vs_LD  = run_contrast(fit, design, "CategoryHD",                    "HD_vs_LD"),
  HD_vs_BD  = run_contrast(fit, design, "CategoryHD - CategoryBD",       "HD_vs_BD")
)

qlf_global <- glmQLFTest(fit, coef = 2:3)
write.csv(topTags(qlf_global, n = Inf)$table,
          file = file.path("DE_outputs", "DE_GLOBAL_Category_effect.csv"))


global_de_list <- rownames(topTags(qlf_global, n = Inf)$table)[
  topTags(qlf_global, n = Inf)$table$FDR < 0.05
]
cat("Total global DE miRNAs:", length(global_de_list), "\n\n")

# Step 2: Extract significant miRNAs from pairwise tests (FDR < 0.05)
bd_vs_ld_sig <- rownames(results_list$BD_vs_LD[
  results_list$BD_vs_LD$FDR < 0.05, 
])

hd_vs_ld_sig <- rownames(results_list$HD_vs_LD[
  results_list$HD_vs_LD$FDR < 0.05, 
])

cat("BD vs LD significant (pairwise):", length(bd_vs_ld_sig), "\n")
cat("HD vs LD significant (pairwise):", length(hd_vs_ld_sig), "\n\n")


# Step 7: Save categorized miRNA lists
write.csv(data.frame(miRNA = bd_vs_ld_sig, Category = "BD-specific"), 
          file = file.path("DE_outputs", "BD_specific_miRNAs.csv"), 
          row.names = FALSE)

write.csv(data.frame(miRNA = hd_vs_ld_sig, Category = "HD-specific"), 
          file = file.path("DE_outputs", "HD_specific_miRNAs.csv"), 
          row.names = FALSE)

library(edgeR)
library(pheatmap)

# logCPM from the *same* object used for DE
logCPM <- cpm(dge_baseline, log = TRUE, prior.count = 1)

mirnas <- unique(c(bd_vs_ld_sig, hd_vs_ld_sig))
# pick the 29 miRNAs
genes <- intersect(mirnas, rownames(logCPM))

# order samples by Category
ord <- order(meta_baseline$Category)
logCPM_ord <- logCPM[genes, ord, drop = FALSE]
meta_ord   <- meta_baseline[ord, ]

# row z-score (visualization)
mat_z <- t(scale(t(logCPM_ord)))

# annotate Category - MAKE SURE ROWNAMES MATCH COLUMN NAMES OF mat_z
ann_col <- data.frame(Category = meta_ord$Category)
rownames(ann_col) <- colnames(mat_z)  # <-- CHANGE THIS LINE

# Verify they match
print(all(colnames(mat_z) == rownames(ann_col)))  # Should print TRUE

# Define colors
ann_colors <- list(
  Category = c(LD = "lightcoral", BD = "lightgreen", HD = "lightblue")
)

png("/Users/viraa/TestingCode/Testing_Monkey/miRNA_heatmap.png", width = 8, height = 10, units = "in", res = 300)

pheatmap(mat_z,
         annotation_col = ann_col,
         annotation_colors = ann_colors,
         show_colnames = FALSE,
         cluster_cols = FALSE,
         fontsize_row = 9,
         cellwidth = 14,
         cellheight = 20,
         main = "Differential miRNA Expression",
         color = colorRampPalette(c("darkblue", "white", "red3"))(50),
         margins = c(2, 44))

dev.off()

library(pROC)

# Define your target miRNAs explicitly
mirnas <- unique(c(bd_vs_ld_sig, hd_vs_ld_sig))

# Rebuild df_DE3 from those miRNAs
logCPM_t <- as.data.frame(t(logCPM[mirnas, ]))
logCPM_t$Sample_ID <- rownames(logCPM_t)
df_DEall <- merge(meta_baseline, logCPM_t, by = "Sample_ID")
rownames(df_DEall) <- df_DEall$Sample_ID

# Remove VHD if present and set factors
df_DEall <- df_DEall[df_DEall$Category != "VHD", , drop = FALSE]
df_DEall$Category <- droplevels(factor(df_DEall$Category, levels = c("LD","BD","HD")))
df_DEall$Sex      <- factor(df_DEall$Sex, levels = c("F","M"))


roc_one_vs_all <- function(data, category, mirna) {
  binary <- ifelse(data$Category == category, 1, 0)
  roc_obj <- roc(binary, data[[mirna]], quiet = TRUE)
  as.numeric(auc(roc_obj))
}

categories <- c("LD","BD","HD")

auc_tables <- lapply(categories, function(cat) {
  auc_vec <- sapply(mirnas, function(m) roc_one_vs_all(df_DEall, cat, m))
  df_auc <- data.frame(miRNA = names(auc_vec), AUC = as.numeric(auc_vec))
  df_auc[order(-df_auc$AUC), ]
})

names(auc_tables) <- categories

df_auc_LD <- auc_tables$LD
df_auc_BD <- auc_tables$BD
df_auc_HD <- auc_tables$HD
df_auc_LD 
df_auc_BD
df_auc_HD

write.csv(df_auc_LD, "AUC_LD.csv", row.names = FALSE)
write.csv(df_auc_BD, "AUC_BD.csv", row.names = FALSE)
write.csv(df_auc_HD, "AUC_HD.csv", row.names = FALSE)

#install.packages("corrplot")
library(corrplot)
library(RColorBrewer)

# Extract only the miRNA columns (excluding Category and Sex)
mirna_data <- df_DEall[, mirnas]

# Calculate correlation matrix
cor_matrix <- cor(mirna_data, use = "complete.obs")

library(pheatmap)

pheatmap(cor_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         fontsize = 8,
         main = "miRNA Correlation Heatmap")



write.csv(cor_matrix, "miRNA_correlations.csv")



#install.packages("pROC")
library(pROC)
library(nnet)


final_miRNAs <- c("mml-miR-1911-3p", "mml-miR-455-5p", "mml-miR-592-5p", 
                  "mml-miR-7182-3p", "mml-miR-9-5p")

logCPM_t <- as.data.frame(t(logCPM[final_miRNAs, ]))
logCPM_t$Sample_ID <- rownames(logCPM_t)
df_DE <- merge(meta_baseline, logCPM_t, by = "Sample_ID")
rownames(df_DE) <- df_DE$Sample_ID


# ROC: BD vs All
png("ROC_BD_vs_All.png", width = 10, height = 2, units = "in", res = 300)
par(mfrow=c(1,5), mar=c(3,3,3,1))
for (m in final_miRNAs) {
  binary  <- ifelse(df_DE$Category == "BD", 1, 0)
  roc_obj <- roc(binary, df_DE[[m]], quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  p_val   <- wilcox.test(df_DE[[m]][df_DE$Category == "BD"],
                         df_DE[[m]][df_DE$Category != "BD"], exact = FALSE)$p.value
  plot(roc_obj, main = paste0(m, "\nAUC=", round(auc_val,3), " p=", round(p_val,3)),
       col = "darkgreen", lwd = 2, cex.main = 0.8)
}
dev.off()


# ROC: HD vs All
png("ROC_HD_vs_All.png", width = 10, height = 2, units = "in", res = 300)
par(mfrow=c(1,5), mar=c(3,3,3,1))
for (m in final_miRNAs) {
  binary  <- ifelse(df_DE$Category == "HD", 1, 0)
  roc_obj <- roc(binary, df_DE[[m]], quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  p_val   <- wilcox.test(df_DE[[m]][df_DE$Category == "HD"],
                         df_DE[[m]][df_DE$Category != "HD"], exact = FALSE)$p.value
  plot(roc_obj, main = paste0(m, "\nAUC=", round(auc_val,3), " p=", round(p_val,3)),
       col = "darkblue", lwd = 2, cex.main = 0.8)
}
dev.off()
# Summary table
results <- data.frame()
for (m in final_miRNAs) {
  p_bd <- wilcox.test(df_DE[[m]][df_DE$Category == "BD"],
                      df_DE[[m]][df_DE$Category != "BD"])$p.value
  p_hd <- wilcox.test(df_DE[[m]][df_DE$Category == "HD"],
                      df_DE[[m]][df_DE$Category != "HD"])$p.value
  results <- rbind(results, data.frame(miRNA = m,
                                       p_BD_vs_All = round(p_bd, 4),
                                       p_HD_vs_All = round(p_hd, 4)))
}
View(results)


write.csv(results, "pvalues_miRNAs_vs_All.csv", row.names = FALSE)

library(ggplot2)
library(ggforce)

list_A_name <- "BD"
list_B_name <- "HD"

list_A <- c("mml-miR-502-3p",
            "mml-miR-7173-5p",
            "mml-miR-573",
            "mml-miR-942-5p",
            "mml-miR-320c",
            "mml-miR-1911-3p",
            "mml-miR-147a",
            "mml-miR-7182-5p",
            "mml-miR-6129",
            "mml-miR-144",
            "mml-miR-7182-3p",
            "mml-miR-192-3p",
            "mml-miR-892c-3p",
            "mml-miR-411-3p",
            "mml-miR-450b-3p",
            "mml-miR-548g-5p",
            "mml-miR-4677-3p"
            )

list_B <- c("mml-miR-877-3p",
            "mml-miR-127-3p",
            "mml-miR-27a-3p",
            "mml-miR-455-5p",
            "mml-miR-9-5p",
            "mml-miR-151-5p",
            "mml-miR-147a",
            "mml-miR-607",
            "mml-miR-448",
            "mml-miR-592-5p",
            "mml-miR-942-5p",
            "mml-miR-875-5p",
            "mml-miR-134-5p",
            "mml-miR-215-5p",
            "mml-miR-361-3p",
            "mml-miR-7173-5p",
            "mml-miR-603",
            "mml-miR-374a-5p",
            "mml-miR-216b")

only_A  <- setdiff(list_A, list_B)
only_B  <- setdiff(list_B, list_A)
both_AB <- intersect(list_A, list_B)


r    <- 3.5
cx   <- 0
cy_A <- 2.5
cy_B <- -2.5


venn <- ggplot() +
  
  # ---- Circle A (top, blue) ----
ggforce::geom_circle(aes(x0 = cx, y0 = cy_A, r = r), fill = "lightgreen", colour = "limegreen", linewidth = 1.8, alpha = 0.55) +
  
  # ---- Circle B (bottom, green) ----
ggforce::geom_circle(aes(x0 = cx, y0 = cy_B, r = r), fill = "lightblue", colour = "skyblue", linewidth = 1.8, alpha = 0.55) +
  
  # ---- Circle name labels (outside, bold) ----
annotate("text",x = cx, y = cy_A + r + 0.55,label    = list_A_name,fontface = "bold", size = 9, colour = "black") +
  annotate("text",x = cx, y = cy_B - r - 0.55, label    = list_B_name, fontface = "bold", size = 9, colour = "black") +
  
  # ---- n= unique counts (inside top/bottom of each circle) ----
annotate("text", x = cx, y = cy_A + 1.8, label    = paste0("n = ", length(only_A)), size     = 6, colour = "black", fontface = "italic") +
  annotate("text",x = cx, y = cy_B - 1.8, label    = paste0("n = ", length(only_B)), size     = 6, colour = "black", fontface = "italic") +
  
  # ---- Shared count in overlap zone ----
annotate("text", x = cx, y = 0, label    = paste0("n = ", length(both_AB)), size     = 6, colour = "black", fontface = "bold") +
  
  # ---- Clean theme, no title ----
coord_fixed() + theme_void() + theme(plot.margin = margin(50, 50, 50, 50))
print(venn)

#ggsave( filename = "miRNA_venn_diagram.png", plot     = venn, width    = 8, height   = 11, dpi      = 300, bg       = "white")

cat("\nSaved: miRNA_venn_diagram.png\n")
print(venn)




library(ggplot2)
library(ggrepel)

# -------------------------
# Configuration
# -------------------------
input_dir <- "DE_outputs"        # Where your DE CSV files are
output_dir <- "DE_outputs"       # Where to save new plots
fdr_cut <- 0.05
lfc_cut <- 0.3

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -------------------------
# Function to create volcano plot from CSV
# -------------------------
create_volcano <- function(csv_file, label, fdr_cut = 0.05, lfc_cut = 0.3) {
  
  # Read the DE results
  res <- read.csv(csv_file, stringsAsFactors = FALSE, row.names = 1)  # Make 1st column row names
  res$miRNA <- rownames(res)  # Convert row names to a column called 'miRNA'
  
  # Prepare plot data
  res$FDR_plot <- pmax(res$FDR, .Machine$double.xmin)
  res$neg_log10_FDR <- -log10(res$FDR_plot)
  res$Significant <- ifelse(abs(res$logFC) >= lfc_cut & res$FDR < fdr_cut, 
                            "Significant", "Not significant")
  
  # Count significant miRNAs
  n_sig <- sum(res$Significant == "Significant")
  message(paste(label, ":", n_sig, "significant miRNAs"))
  
  # Create plot
  p <- ggplot(res, aes(x = logFC, y = neg_log10_FDR)) +
    geom_point(aes(color = Significant), size = 7, alpha = 1) +
    scale_color_manual(values = c("Significant" = "red", 
                                  "Not significant" = "grey60")) +
    
    # Threshold lines
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), 
               linetype = "dashed", color = "blue", linewidth = 0.5) +
    geom_hline(yintercept = -log10(fdr_cut), 
               linetype = "dashed", color = "blue", linewidth = 0.5) +
    
    # Non-overlapping labels for significant points
    geom_text_repel( data = subset(res, Significant == "Significant"),
      aes(label = miRNA),
      size = 6,
      box.padding = 0.5,
      point.padding = 0.3,
      max.overlaps = 20,
      segment.color = "grey50",
      segment.size = 0.3,
      min.segment.length = 0
    ) +
    
    # Labels and styling
    labs(
      x = "log2 Fold Change",
      y = "-log10(FDR)",
      title = paste0("Volcano Plot: ", label),
      subtitle = paste0("FDR < ", fdr_cut, " & |logFC| ≥ ", lfc_cut, 
                        " | Adjusted for Sex"),
      color = ""
    ) +
    
    # Axis limits
    xlim(c(-4, 8)) +
    ylim(c(0, 3)) +
    
    # Theme
    theme_classic(base_size = 16) +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 16),
      axis.title = element_text(size = 16, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2)
    )
  
  # Save plot
  output_file <- file.path(output_dir, paste0("Volcano_", label, "_FDR.png"))
  ggsave(
    filename = output_file,
    plot = p,
    width = 14, 
    height = 10, 
    dpi = 300
  )
  
  message(paste("Saved:", output_file))
  
  return(p)
}


# Define your comparisons (must match your CSV filenames)
comparisons <- list(
  list(file = "DE_BD_vs_LD.csv", label = "BD_vs_LD"),
  list(file = "DE_HD_vs_LD.csv", label = "HD_vs_LD"),
  list(file = "DE_HD_vs_BD.csv", label = "HD_vs_BD")
)

# Create plots
plots <- list()
for(comp in comparisons) {
  csv_path <- file.path(input_dir, comp$file)
  
  if(file.exists(csv_path)) {
    plots[[comp$label]] <- create_volcano(
      csv_file = csv_path,
      label = comp$label,
      fdr_cut = fdr_cut,
      lfc_cut = lfc_cut
    )
  } else {
    warning(paste("File not found:", csv_path))
  }
}



# Base R approach

tsv_path <- "/Users/viraa/TestingCode/Testing_Monkey/Study1/microTHD.tsv"
csv_path <- "/Users/viraa/TestingCode/Testing_Monkey/Study1/microTHD.Csv"

df <- read.delim(tsv_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
write.csv(df, csv_path, row.names = FALSE)

setwd("/Users/viraa/TestingCode/Testing_Monkey/Study1")

BiocManager::install(c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "enrichplot",
  "ggplot2"        
))


library(clusterProfiler)
library(org.Hs.eg.db)  
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)
library(tidyverse)
de_targets <- read.csv("microTHD.Csv")$gene_ensembl_id   # adjust col name
de_targets <- unique(na.omit(de_targets))


all_targets <- read.csv("AllUniquemiRNATargets.csv")$gene_ensembl_id  # adjust col name
all_targets <- unique(na.omit(all_targets))


ensg_to_entrez <- function(ensg_vec) {
  mapped <- mapIds(
    org.Hs.eg.db,
    keys      = ensg_vec,
    column    = "ENTREZID",
    keytype   = "ENSEMBL",
    multiVals = "first"       # keep one Entrez per ENSG
  )
  mapped <- mapped[!is.na(mapped)]
  return(mapped)
}

de_entrez  <- ensg_to_entrez(de_targets)
all_entrez <- ensg_to_entrez(all_targets)


ego <- enrichGO(
  gene          = de_entrez,
  universe      = all_entrez,   # ← your custom background
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",        # BP, MF, CC, or ALL
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE          # converts IDs → gene symbols in output
)


ekegg <- enrichKEGG(
  gene          = de_entrez,
  universe      = all_entrez,   # ← your custom background
  organism      = "hsa",        # hsa = human, mmu = mouse
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

# Write to file
write.csv(as.data.frame(ego),   "HD_GO_enrichment.csv",   row.names = FALSE)
write.csv(as.data.frame(ekegg), "HD_KEGG_enrichment.csv", row.names = FALSE)

colnames(read_csv("HD_GO_enrichment.csv"))

read_pathway_csv <- function(filepath, group_label) {
  read_csv(filepath, show_col_types = FALSE) %>%
    rename(
      Pathway        = Description,
      Hit_count      = Count,
      FDR            = p.adjust,
      RawPvalue      = pvalue
    ) %>%
    mutate(
      Pathway        = str_remove(Pathway, " \\([^)]+\\)$"),
      FoldEnrichment = as.numeric(FoldEnrichment),
      FDR            = as.numeric(FDR),
      logFDR         = -log10(FDR),
      Comparison     = group_label
    )
}

# Load your files
pathway_BDSlim <- read_pathway_csv("BD_GO_enrichment.csv", "BD")
pathway_HDSlim <- read_pathway_csv("HD_GO_enrichment.csv", "HD")



prep_plot_data <- function(data, min_hits = 50, min_fold = 1.2, top_n = 31) {
  data %>%
    filter(FDR < 0.05) %>%
    filter(Hit_count > min_hits) %>%        # <-- gene count filter
    filter(FoldEnrichment > min_fold) %>%   # <-- fold change filter
    arrange(desc(FoldEnrichment)) %>%       # sort highest to lowest
    slice_head(n = top_n)                   # take top N after sorting
}

plot_data <- bind_rows(
  prep_plot_data(pathway_BDSlim),
  prep_plot_data(pathway_HDSlim)
) %>%
  mutate(Pathway = str_trunc(Pathway, 60),
         Pathway = fct_reorder(Pathway, FoldEnrichment, .fun = max))

# --- Plot ---
ggplot(plot_data, aes(x = Comparison, y = Pathway)) +
  geom_point(aes(size = FoldEnrichment, color = logFDR)) +
  scale_size_continuous(
    range = c(6, 10),
    name = "Fold Enrichment"
  ) +
  scale_color_gradient(
    low  = "darkblue",
    high = "red3",
    name = expression(-log[10](FDR))
  ) +
  labs(
    title = "Enriched GO Biological Process (GO-Slim)",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title         = element_text(size = 14, hjust = 0.5, face = "bold"),
    axis.text.x        = element_text(size = 10, face = "bold"),
    axis.text.y        = element_text(size = 10),
    legend.title       = element_text(size = 10),
    legend.text        = element_text(size = 10),
    legend.position    = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    plot.margin        = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  )


# Save the plot
ggsave("Enriched_GO_Biological_Process_BD_vs_HD.png", 
       width = 9, height = 10, dpi = 300, units = "in")
