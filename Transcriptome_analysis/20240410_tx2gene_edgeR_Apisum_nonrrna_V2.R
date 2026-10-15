rm(list=ls()) 


# Load required libraries
library(splicejam)
library(dplyr)
library(GenomicFeatures)

tx2gene <- makeTx2geneFromGtf(GTF = "~/Desktop/RNAseq_analysis_Apisum_nonrrna/GCF_005508785.2/genomic.gtf", geneAttrNames = c("gene_id", "gene_name", "gene_type"),
                              txAttrNames = c("transcript_id", "transcript_type"),
                              geneFeatureType = "gene",
                              txFeatureType = c("transcript", "mRNA"),
                              nrows = -1L,
                              verbose = FALSE,)

# Rearrange the columns
tx2gene <- tx2gene[, c("transcript_id", "gene_id")]

write.csv(tx2gene,file = "~/Desktop/RNAseq_analysis_Apisum_nonrrna/tx2gene.csv")

library(readr)
library(dplyr)
tx2gene <- read_csv("~/Desktop/RNAseq_analysis_Apisum_nonrrna/tx2gene.csv")

tx2gene <- tx2gene[, -1]


# Print the resulting tibble
print(tx2gene)
head(tx2gene)


library(tximport)
library(edgeR)
library(tidyverse)
library(limma)
library(org.Apisum.eg.db)

# List all directories containing quant.sf files
folder_paths <- list.dirs("~/Desktop/RNAseq_analysis_Apisum_nonrrna/Apisum_nonrrna_salmonOut", full.names = TRUE, recursive = FALSE)

# Initialize an empty list to store the data
quant_sf_data <- list()

# Iterate over each folder path
for (i in seq_along(folder_paths)) {
  folder_path <- folder_paths[i]
  
  # Construct the path to the quant.sf file
  quant_sf_file <- file.path(folder_path, "quant.sf")
  
  # Check if the quant.sf file exists
  if (file.exists(quant_sf_file)) {
    # Read the quant.sf file and store the data in the list
    quant_sf_data[[i]] <- read.table(quant_sf_file, header = TRUE)
  }
}

# Extract sample names from folder names
sample_names <- basename(folder_paths)

# Specify the directory containing the files
directory <- "~/Desktop/RNAseq_analysis_Apisum_nonrrna/Apisum_nonrrna_quant_sf"

# Get the list of files in the directory
file_names <- list.files(directory, full.names = TRUE)
# Display the file paths
print(file_names)
names(file_names) <- sample_names
txi <- tximport(files =file_names, type = "salmon", tx2gene = tx2gene) #perform gene-level summarization with txOut =FALSE(default)

#txi <- tximport(files =file_names, type="salmon", tx2gene=tx2gene[,c("transcript_id", "gene_id")], countsFromAbundance="dtuScaledTPM", txOut =TRUE)
#This scaling option is designed for use with txOut=TRUE for differential transcript usage analyses.

#scaled using the average transcript length over samples and then the library size (lengthScaledTPM)
#or geneLengthScaledTPM, the counts are no longer correlated across samples with transcript length, and so the length offset matrix should not be used.

write.csv(txi,file = "~/Desktop/RNAseq_analysis_Apisum_nonrrna/20240410_gene_leveltximport_salmon.csv")

# Extract sample names without the "_quant" suffix
sample_names <- sub("_quant$", "", basename(folder_paths))


# Create the grp vector
grp <- c(rep("BactMale", 5), rep("BactOvi", 2), rep("BactVivi", 5), rep("WB", 5), rep("WBMale", 3))

# Display the grp vector
print(grp)


data.frame(sample = colnames(txi$counts),  grp = grp)


####refer*https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#Import_transcript-level_estimates
library(edgeR)
library(vegan)

# ------------------ 1. Normalization & filtering ------------------
cts <- txi$counts
normMat <- txi$length
normMat <- normMat / exp(rowMeans(log(normMat)))
normCts <- cts / normMat
eff.lib <- calcNormFactors(normCts) * colSums(normCts)
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

y <- DGEList(cts, group = grp)
y <- scaleOffset(y, normMat)
design <- model.matrix(~grp)
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]

cpms <- edgeR::cpm(y, offset = y$offset, log = TRUE)

# ------------------ 2. Color & shape settings ------------------
group_colors <- c(
  BactMale = "#D55E00",   # vermillion
  BactOvi  = "#F0E442",   # yellow
  BactVivi = "#CC79A7",   # magenta (lipidome bacteriocytes)
  WB       = "#0072B2",   # blue (lipidome whole body)
  WBMale   = "#009E73"    # bluish green
)

sample_type <- ifelse(grepl("^Bact", grp), "Bacteriocyte", "WholeBody")
pch_shapes <- ifelse(sample_type == "Bacteriocyte", 17, 15)
# ------------------ 3. PCA & PERMANOVA ------------------
pca <- prcomp(t(cpms), scale. = TRUE)
pvar <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

permanova <- adonis2(t(cpms) ~ grp, method = "euclidean", permutations = 999)
F_value <- round(permanova$F[1], 2)
R2_value <- round(permanova$R2[1], 2)
p_value <- permanova$`Pr(>F)`[1]
p_label <- if (p_value < 0.001) "p < 0.001" else paste0("p = ", formatC(p_value, format = "f", digits = 3))

# ------------------ 4. Compact PCA plot with legend ------------------
layout(matrix(c(1, 2), ncol = 2), widths = c(4, 1))

# --- PCA plot panel ---
par(mar = c(5, 4, 5, 1))  # increase top margin
plot(pca$x[,1], pca$x[,2],
     col = group_colors[grp],
     pch = pch_shapes,
     cex = 1.4,
     xlab = paste0("PC1 (", round(pvar[1], 1), "% variance)"),
     ylab = paste0("PC2 (", round(pvar[2], 1), "% variance)"))

# --- Add PERMANOVA text ABOVE the plot box ---
par(xpd = NA)  # allow drawing in margin area
usr <- par("usr")

# Set position just above the plotting area
x_pos <- mean(c(usr[1], usr[2]))  # centered horizontally
y_pos <- usr[4] + (usr[4] - usr[3]) * 0.12  # slightly above the box

# Extract numeric p-value cleanly
p_numeric <- if (p_value < 0.001) "< 0.001" else formatC(p_value, format = "f", digits = 3)

# Add the annotation with evaluated variables
text(x_pos, y_pos,
     labels = bquote(bold("PERMANOVA:") ~ 
                       F == .(F_value) * "," ~ 
                       R^2 == .(R2_value) * "," ~ 
                       italic(p) == .(p_numeric)),
     adj = 0.5, cex = 1.0)

par(xpd = FALSE)

# --- Legend panel ---
par(mar = c(5, 0, 3, 2))
plot.new()
legend("center",
       legend = c("BactMale", "BactOvi", "BactVivi", "WB", "WBMale"),
       col = group_colors,
       pch = c(17,17,17,15,15),
       pt.cex = 1.4,
       cex = 0.9,
       bty = "n",
       title = "Sample Type",
       title.cex = 0.9)

#################################AphidBactVS_WB#################

######################################################################
### Differential Expression Analysis - Aphid Morphs (edgeR)
######################################################################

library(AnnotationDbi)
library(org.Apisum.eg.db)



# ---------------- 1. Estimate dispersion ----------------
y_disp <- estimateCommonDisp(y)
y_disp <- estimateTagwiseDisp(y_disp)

summary(y_disp$tagwise.dispersion)

# ---------------- 2. Classic exactTest (pairwise) ----------------
# Compare WB vs BactVivi
de.tagwise <- exactTest(y_disp, pair = c("WB", "BactVivi"))

# Check top DE genes
topTags(de.tagwise, n = 20)

# Save results
write.table(de.tagwise$table,
            "~/Desktop/RNAseq_analysis_Apisum_nonrrna/20251025_de_tagwise_BactVivi_vs_WB.txt",
            sep = "\t", quote = FALSE)

# Fix incorrect n argument
tm <- topTags(de.tagwise, n = nrow(de.tagwise$table))

# Identify DE genes
de.names <- rownames(cts[decideTestsDGE(de.tagwise, p.value = 0.05) != 0, ])
plotMD(de.tagwise, cex = 0.7, main = "MD Plot: BactVivi vs WB")

# Extract filtered table
detab <- tm$table

# Strongly upregulated genes
detab[detab$logFC > log2(10), ]
nrow(detab[detab$logFC > log2(10), ])

# DE genes: >5-fold change and FDR < 0.01
de_bacVswb <- subset(detab, logFC > log2(5) & FDR < 0.01)
summary(decideTestsDGE(de.tagwise, p.value = 0.05))


######################################################################
### GLM-based analysis (recommended)
######################################################################

# ---------------- 3. Design and dispersion ----------------

grp <- factor(grp, levels = c("BactVivi","WB","BactMale","BactOvi","WBMale"))
grp <- relevel(grp, ref = "BactVivi")   # redundant but makes intent explicit
design <- model.matrix(~ grp)
colnames(design)

y <- estimateDisp(y, design, robust = TRUE)
y$common.dispersion

# Examine dispersion trends
plotBCV(y)

# ---------------- 4. Fit GLM and perform tests ----------------
fit <- glmQLFit(y, design)

# Check column names before defining contrasts
colnames(design)
#Intercept = BactVivi) + grpWB (add this if WB) + grpBactMale (add this if BactMale)
#+ grpBactOvi  (add this if BactOvi)
#+ grpWBMale   (add this if WBMale)
# --- Comparison 1: BactVivi vs WB ---

tr <- glmTreat(
  fit,
  contrast = c(0, -1, 0, 0, 0),
  lfc = 1
)
# Add adjusted p-values
tr$table$adj_PValue <- p.adjust(tr$table$PValue, method = "BH")

# Examine and visualize
topTags(tr, n = 10, sort.by = "logFC")
plotMD(tr, cex = 0.7, main = "GLM-Treat: BactVivi vs WB")

# Save results
significant_results <- tr$table
significant_results$gene_id <- rownames(significant_results)

write_csv(significant_results,
          "~/Desktop/RNAseq_analysis_Apisum_nonrrna/20240410_siglfc2_BactVivivsWB.csv")

annots <- AnnotationDbi::select(org.Apisum.eg.db,
                                keys = rownames(significant_results),
                                columns = c("GOALL"),
                                keytype = "SYMBOL")
resultTable2 <- merge(significant_results, annots, by.x = "gene_id", by.y = "SYMBOL")
write_csv(resultTable2,"~/Desktop/RNAseq_analysis_Apisum_nonrrna/20240410_siglfc2_adjpval_GOIDBactVivivsWB.csv")


# Identify significant DEGs
sig_DEG <- significant_results %>%
  dplyr::filter(adj_PValue < 0.05 & abs(logFC) > 1)

up_DEG <- sig_DEG %>% filter(logFC > 1)
down_DEG <- sig_DEG %>% filter(logFC < -1)

# Print counts
cat("Upregulated genes:", nrow(up_DEG), "\n")
cat("Downregulated genes:", nrow(down_DEG), "\n")
cat("Total DEGs:", nrow(sig_DEG), "\n")

#> cat("Upregulated genes:", nrow(up_DEG), "\n")
#Upregulated genes: 999 
#> cat("Downregulated genes:", nrow(down_DEG), "\n")
#Downregulated genes: 2633 
#> cat("Total DEGs:", nrow(sig_DEG), "\n")
#Total DEGs: 3632 

##############################################################
### GO analysis (GOseq) 
##############################################################

library(goseq)
library(GO.db)
library(AnnotationDbi)
library(org.Apisum.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)
library(viridis)

# --- 1. Prepare gene vector (significant vs non-significant) ---
sigData <- as.integer(!is.na(significant_results$adj_PValue) &
                        significant_results$adj_PValue < 0.05)
names(sigData) <- significant_results$gene_id
table(sigData)

# --- 2. Extract gene lengths from GTF file ---
gtf <- read.delim(
  "~/Desktop/RNAseq_analysis_Apisum_nonrrna/GCF_005508785.2/genomic.gtf",
  header = FALSE, comment.char = "#", stringsAsFactors = FALSE
)

gene_entries <- gtf[gtf$V3 == "gene", ]
gene_entries$gene_id <- sub('.*gene_id "?([^"; ]+).*', "\\1", gene_entries$V9)
gene_entries$gene_length <- gene_entries$V5 - gene_entries$V4 + 1
gene_data <- gene_entries[, c("gene_id", "gene_length")]

# --- 3. Align gene length info with DE genes ---
common_genes <- intersect(names(sigData), gene_data$gene_id)
sigData <- sigData[common_genes]
lengths_vec <- gene_data$gene_length[match(common_genes, gene_data$gene_id)]
names(lengths_vec) <- common_genes
lengths_vec <- lengths_vec[!is.na(lengths_vec) & lengths_vec > 0]

# --- 4. GOseq probability weighting function ---
pwf <- nullp(sigData, bias.data = lengths_vec)
plot(pwf$bias, pwf$DEgenes, pch = 19, cex = 0.5,
     xlab = "Gene length",
     ylab = "Probability of being DE",
     main = "Bias correction (GOseq nullp)")

# --- 5. Prepare GO annotation mapping ---
sig_ID <- names(sigData)
annots <- AnnotationDbi::select(
  org.Apisum.eg.db,
  keys = sig_ID,
  columns = "GOALL",
  keytype = "SYMBOL"
)

annots <- na.omit(annots)
colnames(annots) <- c("gene_id", "GOALL")
gene2cat <- split(annots$GOALL, annots$gene_id)

# --- 6. Run GOseq on all DE genes ---
BactVivi_wb <- goseq(
  pwf = pwf,
  gene2cat = gene2cat,
  use_genes_without_cat = TRUE,
  method = "Wallenius"
)

# Add FDR column (BH correction)
BactVivi_wb$over_represented_FDR <- p.adjust(BactVivi_wb$over_represented_pvalue, method = "BH")

# --- 7. Clean and visualize (FDR-based) ---
BactVivi_wb2 <- BactVivi_wb %>%
  filter(!is.na(term), term != "Obsolete term", over_represented_FDR < 0.05) %>%
  arrange(over_represented_FDR) %>%
  mutate(hitsPerc = numDEInCat * 100 / numInCat)

ggplot(head(BactVivi_wb2, 20),
       aes(x = hitsPerc,
           y = reorder(term, hitsPerc),
           color = over_represented_FDR,
           size = numDEInCat)) +
  geom_point() +
  scale_color_viridis(option = "plasma", direction = -1, name = "FDR") +
  labs(x = "Hits (%)",
       y = "GO term",
       title = "Top enriched GO terms (BactVivi vs WB, all DE)") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 8)
  )


##############################################################
### GOseq for Up- and Downregulated genes separately (FDR-based)
##############################################################

# --- 1. Define up- and downregulated gene sets ---
up_genes <- rownames(significant_results[
  significant_results$logFC > 1 & significant_results$adj_PValue < 0.05, ])

down_genes <- rownames(significant_results[
  significant_results$logFC < -1 & significant_results$adj_PValue < 0.05, ])

# --- 2. Binary vectors for GOseq ---
sigData_up   <- as.integer(names(sigData) %in% up_genes)
sigData_down <- as.integer(names(sigData) %in% down_genes)
names(sigData_up) <- names(sigData_down) <- names(sigData)

# --- 3. Run GOseq ---
pwf_up   <- nullp(sigData_up, bias.data = lengths_vec)
pwf_down <- nullp(sigData_down, bias.data = lengths_vec)

GO_up   <- goseq(pwf_up,   gene2cat = gene2cat, use_genes_without_cat = TRUE)
GO_down <- goseq(pwf_down, gene2cat = gene2cat, use_genes_without_cat = TRUE)

# Add FDR columns
GO_up$over_represented_FDR   <- p.adjust(GO_up$over_represented_pvalue, method = "BH")
GO_down$over_represented_FDR <- p.adjust(GO_down$over_represented_pvalue, method = "BH")

# --- 4. Clean NA and obsolete terms ---
GO_up <- GO_up %>% filter(!is.na(term), term != "Obsolete term", over_represented_FDR < 0.05)
GO_down <- GO_down %>% filter(!is.na(term), term != "Obsolete term", over_represented_FDR < 0.05)

# --- 5. Prepare data for plotting ---
GO_up2 <- GO_up %>%
  arrange(over_represented_FDR) %>%
  slice_head(n = 20) %>%
  mutate(
    hitsPerc = numDEInCat * 100 / numInCat,
    term_wrapped = str_wrap(term, width = 45),
    Direction = "Upregulated"
  )

GO_down2 <- GO_down %>%
  arrange(over_represented_FDR) %>%
  slice_head(n = 20) %>%
  mutate(
    hitsPerc = numDEInCat * 100 / numInCat,
    term_wrapped = str_wrap(term, width = 45),
    Direction = "Downregulated"
  )

GO_combined <- bind_rows(GO_up2, GO_down2)

# --- 6. Colorblind-friendly faceted plot ---
ggplot(GO_combined,
       aes(x = hitsPerc,
           y = reorder(term_wrapped, hitsPerc),
           color = over_represented_FDR,
           size = numDEInCat)) +
  geom_point() +
  scale_color_viridis(option = "magma", direction = -1, name = "FDR") +
  scale_size(range = c(3, 10)) +
  labs(x = "Hits (%)",
       y = "GO term",
       title = "Top enriched GO terms by direction (BactVivi vs WB)") +
  facet_wrap(~Direction, scales = "free_y") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 11, face = "bold")
  )


##############GofuncR#########
############## GOfuncR ##############
library(GOfuncR)
library(GO.db)
library(dplyr)

# --- 1. Filter and adjust for FDR ---
sigCats <- BactVivi_wb %>%
  filter(!is.na(category),
         !is.na(term),
         term != "Obsolete term") %>%
  mutate(over_represented_FDR = p.adjust(over_represented_pvalue, "BH")) %>%
  filter(over_represented_FDR < 0.05)

# --- 2. Save clean table ---
write.table(
  sigCats,
  "~/Desktop/RNAseq_analysis_Apisum_nonrrna/gofuncR_sigCats_FDR005.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# --- 3. Extract GO IDs and their hierarchical descendants ---
go_ids <- as.character(sigCats$category)
childNodes <- get_child_nodes(go_ids)

# --- 4. Annotate child terms with GO names safely ---
safe_get_GO_term <- function(go_id) {
  if (go_id %in% keys(GOTERM)) Term(GOTERM[[go_id]]) else NA_character_
}

childNodes$child_term <- vapply(childNodes$child_go_id, safe_get_GO_term, character(1))

# --- 5. (Optional) Also annotate parent names ---
childNodes$parent_term <- vapply(childNodes$parent_go_id, safe_get_GO_term, character(1))

# --- 6. Save annotated hierarchy ---
write.table(
  childNodes,
  "~/Desktop/RNAseq_analysis_Apisum_nonrrna/gofuncR_childNodes_FDR005.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# --- 7. Preview ---
head(childNodes)


write.table(childNodes, "~/Desktop/RNAseq_analysis_Apisum_nonrrna/20251026_gofuncR_childNodes_named.txt",
            sep="\t", quote=FALSE, row.names=FALSE)



########################Volcano_GO6644_GO6629_withLegend############

library(ggplot2)
library(dplyr)
library(ggrepel)
library(AnnotationDbi)

# --- Retrieve GO annotation ----------------------------------------------
annots <- AnnotationDbi::select(
  org.Apisum.eg.db,
  keys = rownames(significant_results),
  columns = c("GOALL"),
  keytype = "SYMBOL"
)

# --- Identify genes in each GO term --------------------------------------
genes_GO6644 <- annots %>%
  filter(GOALL == "GO:0006644") %>%
  pull(SYMBOL) %>% unique()

genes_GO6629 <- annots %>%
  filter(GOALL == "GO:0006629") %>%
  pull(SYMBOL) %>% unique()

# --- Hierarchical classification ----------------------------------------
significant_results <- significant_results %>%
  mutate(highlight = case_when(
    gene_id %in% genes_GO6644 ~ "GO:0006644",
    gene_id %in% genes_GO6629 & !(gene_id %in% genes_GO6644) ~ "GO:0006629_only",
    TRUE ~ "Other"
  ))

# --- Colors --------------------------------------------------------------
col_yellow <- "#E69F00"   # phospholipid (child)
col_blue   <- "#56B4E9"   # lipid (parent-only)
col_grey   <- "grey80"
target_border <- "#1F3A93"

target_genes <- c("LOC100160658")

# --- Volcano plot --------------------------------------------------------
ggplot(significant_results, aes(x = logFC, y = -log10(adj_PValue))) +
  
  # Background grey points
  geom_point(
    data = subset(significant_results, highlight == "Other"),
    aes(fill = highlight),
    color = col_grey, fill = col_grey,
    shape = 21, size = 1.4, alpha = 0.4
  ) +
  
  # GO terms — unified mapping for legend
  geom_point(
    data = subset(significant_results, highlight != "Other"),
    aes(fill = highlight),
    color = "black", shape = 21,
    size = 3, stroke = 0.4, alpha = 0.95
  ) +
  
  # Highlight GPDH
  geom_point(
    data = subset(significant_results, gene_id == "LOC100160658"),
    aes(x = logFC, y = -log10(adj_PValue)),
    fill = col_yellow, color = target_border,
    shape = 21, size = 5, stroke = 1.3, inherit.aes = FALSE
  ) +
  
  # Label GPDH
  geom_text_repel(
    data = subset(significant_results, gene_id %in% target_genes),
    aes(label = gene_id),
    color = target_border, size = 4, fontface = "bold",
    box.padding = 0.5, segment.color = "grey40"
  ) +
  
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  
  # --- Proper legend -----------------------------------------------------
scale_fill_manual(
  values = c(
    "GO:0006644"      = col_yellow,
    "GO:0006629_only" = col_blue,
    "Other"           = col_grey
  ),
  breaks = c("GO:0006644", "GO:0006629_only"),
  labels = c(
    "GO:0006644  Phospholipid metabolism",
    "GO:0006629  Lipid metabolism"
  ),
  name = "Gene Ontology category"
) +
  
  labs(
    subtitle = "GO:0006644 ⊂ GO:0006629",
    x = expression(log[2]~"Fold Change (BactVivi vs WB)"),
    y = expression(-log[10]~"FDR")
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5),
    legend.position = "top",
    legend.title    = element_text(size = 10, face = "bold"),
    legend.text     = element_text(size = 9),
    legend.key      = element_rect(color = "black", size = 0.3)
  )

library(dplyr)
library(readr)

# --- Ensure gene_id column exists only once ---
if (!"gene_id" %in% colnames(significant_results)) {
  significant_results <- tibble::rownames_to_column(significant_results, "gene_id")
}

# --- Add missing GO term for GPDH1 (LOC100160658) -------------------------
extra_go <- data.frame(
  SYMBOL = "LOC100160658",
  GOALL  = "GO:0006644",
  stringsAsFactors = FALSE
)
annots_corrected <- dplyr::bind_rows(annots, extra_go)

# --- Merge GO info back to edgeR DE table --------------------------------
GO_DEGs_full <- significant_results %>%
  left_join(annots_corrected, by = c("gene_id" = "SYMBOL")) %>%
  filter(GOALL %in% c("GO:0006629", "GO:0006644")) %>%
  dplyr::select(
    gene_id,
    logFC,
    unshrunk.logFC,
    logCPM,
    PValue,
    adj_PValue,
    GOALL
  ) %>%
  arrange(GOALL, desc(abs(logFC)))

# --- Classify regulation -------------------------------------------------
GO_DEGs_full <- GO_DEGs_full %>%
  mutate(Regulation = case_when(
    adj_PValue < 0.05 & logFC > 1  ~ "Upregulated",
    adj_PValue < 0.05 & logFC < -1 ~ "Downregulated",
    TRUE                            ~ "Not significant"
  ))

# --- Export CSV ----------------------------------------------------------
write_csv(
  GO_DEGs_full,
  "~/Desktop/RNAseq_analysis_Apisum_nonrrna/20251025_figures/Table_S3_GO_lipid_DEGs.csv"
)

# --- Summary counts ------------------------------------------------------
summary_counts <- GO_DEGs_full %>%
  group_by(GOALL, Regulation) %>%
  summarise(Count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = Regulation,
    values_from = Count,
    values_fill = 0
  ) %>%
  dplyr::arrange(GOALL)

# --- Overlap between GO terms --------------------------------------------
overlap_genes <- intersect(
  GO_DEGs_full$gene_id[GO_DEGs_full$GOALL == "GO:0006629"],
  GO_DEGs_full$gene_id[GO_DEGs_full$GOALL == "GO:0006644"]
)

cat("\n Exported Table_S3_GO_lipid_DEGs.csv\n\n")
cat("Summary by GO term:\n")
print(summary_counts)
cat("\nGenes present in BOTH GO terms:", length(overlap_genes), "\n")
cat(paste(overlap_genes, collapse = ", "), "\n")


#GOALL      Downregulated `Not significant` Upregulated
#<chr>              <int>             <int>       <int>
 # 1 GO:0006629            47                92          25
#2 GO:0006644             4                45           5
#> cat("\nGenes present in BOTH GO terms:", length(overlap_genes), "\n")

#Genes present in BOTH GO terms: 53 

#######################Acetyl-CoA transferase##################3
library(dplyr)
library(readr)

# --- Load mapping table --
mapping <- read_csv("~/Desktop/RNAseq_analysis_Apisum_nonrrna/20251025_figures/Alt_ACYPI_to_LOC_Mapping.csv")

# --- Load Acetyl_Apisum IDs ---
acetyl_ids <- readLines("~/Desktop/RNAseq_analysis_Apisum_nonrrna/20251025_figures/Acetyl_Apisum.txt")
acetyl_ids <- gsub(">7029\\.|-PA", "", acetyl_ids)  # clean to get pure ACYPI IDs

# ---  Map to LOC IDs ---
acetyl_loc <- mapping %>%
  filter(ACYPI_ID %in% acetyl_ids) %>%
  distinct(LOC_ID) %>%
  pull(LOC_ID)

cat("Mapped to", length(acetyl_loc), "LOC IDs\n")

acetyl_expr <- significant_results[significant_results$gene_id %in% acetyl_loc, ]
acetyl_expr$Regulation <- ifelse(
  acetyl_expr$adj_PValue < 0.05 & acetyl_expr$logFC > 1, "Upregulated",
  ifelse(acetyl_expr$adj_PValue < 0.05 & acetyl_expr$logFC < -1, "Downregulated", "Not significant")
)

acetyl_expr <- acetyl_expr[, c("gene_id", "logFC", "PValue", "adj_PValue", "Regulation",
                               setdiff(colnames(acetyl_expr), c("gene_id", "logFC", "PValue", "adj_PValue", "Regulation")))]

write.csv(
  acetyl_expr,
  "~/Desktop/RNAseq_analysis_Apisum_nonrrna/20251025_figures/Acetyl_Apisum_expression.csv",
  row.names = FALSE
)

cat("Number of mapped acetyl genes found in DE results:", nrow(acetyl_expr), "\n")
table(acetyl_expr$Regulation)

#Number of mapped acetyl genes found in DE results: 17 
#Downregulated Not significant     Upregulated 
# 3               6               8 


ggplot(acetyl_expr, aes(x = logFC, y = -log10(adj_PValue), color = Regulation)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c(
    "Upregulated"   = "#E64B35",
    "Downregulated" = "#800080",
    "Not significant" = "#BDBDBD"
  )) +
  
  # Add LOC100159262 with same bright pink
  geom_point(
    data = subset(acetyl_expr, gene_id == "LOC100159262"),
    aes(x = logFC, y = -log10(adj_PValue)),
    fill = "#FF4DA6", color = "black", shape = 21, size = 4.5, stroke = 1
  ) +
  geom_text_repel(
    data = subset(acetyl_expr, gene_id == "LOC100159262"),
    aes(label = gene_id),
    color = "#FF4DA6", size = 4, fontface = "bold",
    box.padding = 0.5, segment.color = "grey40"
  ) +
  
  labs(
    title = "Volcano plot of Acetyl-CoA transporter genes",
    x = expression(log[2]~"Fold Change"),
    y = expression(-log[10]~"FDR")
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# --- Clean column names & order safely ---
acetyl_expr_export <- acetyl_expr %>%
  dplyr::select(gene_id, logFC, PValue, adj_PValue, Regulation) %>%
  arrange(desc(logFC))

# --- Write as CSV for Supplementary Table S4 ---
write.csv(
  acetyl_expr_export,
  "~/Desktop/RNAseq_analysis_Apisum_nonrrna/20251025_figures/Table_S4_AcetylCoA_transporters.csv",
  row.names = FALSE
)

cat(" Exported", nrow(acetyl_expr_export), "Acetyl-CoA transporter genes\n")


#################ClusterProfiler############

library(pathview)
library(KEGGREST)
library(dplyr)
library(readr)


# ---  Prepare KEGG mapping between SYMBOL and ENTREZ ID ---
significant_gene_ids <- significant_results$gene_id

annots_kegg <- AnnotationDbi::select(
  org.Apisum.eg.db,
  keys = significant_gene_ids,
  columns = c("ENTREZID"),
  keytype = "SYMBOL"
)

# Merge DEG table with annotation
new_gene_kegg <- merge(
  significant_results,
  annots_kegg,
  by.x = "gene_id",
  by.y = "SYMBOL",
  all.x = TRUE
)

# ---  Retrieve Apisum KEGG–NCBI gene mappings ---
kegg_ids <- keggConv("api", "ncbi-geneid")
kegg_ids_df <- as.data.frame(kegg_ids, stringsAsFactors = FALSE)
kegg_ids_df$ENTREZID <- gsub("ncbi-geneid:", "", rownames(kegg_ids_df))

# Ensure matching type
new_gene_kegg$ENTREZID <- as.character(new_gene_kegg$ENTREZID)

# Merge with KEGG IDs
merged_data <- merge(new_gene_kegg, kegg_ids_df, by = "ENTREZID", all.x = TRUE)
colnames(merged_data)[which(names(merged_data) == "kegg_ids")] <- "kegg_id"

# Save annotated DEG–KEGG table
write_csv(
  merged_data,
  "~/Desktop/RNAseq_analysis_Apisum_nonrrna/20251025_figures/api_DEG_ENTREZ_KEGGids.csv"
)



# ============================================================================
# 1. Genes: retain strongest transcript for each ENTREZ gene
# ============================================================================

gene_by_entrez <- merged_data %>%
  dplyr::filter(
    adj_PValue < 0.05,
    !is.na(ENTREZID),
    !is.na(logFC)
  ) %>%
  dplyr::mutate(ENTREZID = as.character(ENTREZID)) %>%
  group_by(ENTREZID) %>%
  slice_max(
    order_by = abs(logFC),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

# Named vector required by Pathview
gene_logFC <- gene_by_entrez$logFC
names(gene_logFC) <- gene_by_entrez$ENTREZID

cat("Unique ENTREZ genes:", length(gene_logFC), "\n")
cat("Duplicated ENTREZ IDs:", anyDuplicated(names(gene_logFC)), "\n")


# ============================================================================
# 2. Lipids: retain significant Buchnera-enriched species
# ============================================================================

lipid_de <- read_csv(
  "~/Desktop/20230522_BuchneraLipidomics/20240213_LCMS_Riken/20251002_LSMS_Riken/20251002_Buchnera_Bacteriocytes/2025_PhosphodeBchnera_Bacteriocytes.csv"
)

lipid_buch_sig <- lipid_de %>%
  dplyr::filter(
    adj.P.Val < 0.05,
    logFC > 1,
    !is.na(logFC),
    !is.na(Class)
  )


# ============================================================================
# 3. Aggregate significant lipid species by class
# ============================================================================

lipid_class_summary <- lipid_buch_sig %>%
  group_by(Class) %>%
  summarise(
    mean_logFC = mean(logFC, na.rm = TRUE),
    max_logFC  = max(logFC, na.rm = TRUE),
    n_species  = n(),
    .groups = "drop"
  )


# ============================================================================
# 4. Map lipid classes to KEGG compound IDs
# ============================================================================

class_to_keggID <- c(
  "PA"   = "C00416",
  "DG"   = "C00641",
  "PC"   = "C00157",
  "PE"   = "C00350",
  "PG"   = "C00344",
  "PI"   = "C01194",
  "PS"   = "C02737",
  "LPC"  = "C04233",
  "LPE"  = "C04438",
  "LPG"  = "C05973",
  "LPI"  = "C05974",
  "LPS"  = "C05975",
  "CL"   = "C05980",
  "DLCL" = "C18129"
)

lipid_class_summary <- lipid_class_summary %>%
  dplyr::mutate(KEGG_ID = unname(class_to_keggID[Class])) %>%
  dplyr::filter(!is.na(KEGG_ID))

# One class-level mean logFC per KEGG compound
compound_logFC <- lipid_class_summary$mean_logFC
names(compound_logFC) <- lipid_class_summary$KEGG_ID

cat("Lipid classes mapped to KEGG:", length(compound_logFC), "\n")


# ============================================================================
# 5. Run Pathview
# ============================================================================

outdir <- paste0(
  "~/Desktop/RNAseq_analysis_Apisum_nonrrna/",
  "20251025_figures/pathview2"
)

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

api_phospholipid <- pathview(
  gene.data   = gene_logFC,
  cpd.data    = compound_logFC,
  
  pathway.id  = "00564",
  species     = "api",
  gene.idtype = "ENTREZ",
  
  # When multiple ENTREZ genes map to the same KEGG node,
  # use the logFC with the largest absolute value.
  node.sum = "max.abs",
  
  kegg.native = TRUE,
  map.symbol  = TRUE,
  
  out.suffix = "2026_maxAbsGene_BuchneraLipids",
  out.dir     = outdir,
  kegg.dir    = outdir,
  
  low = list(
    gene = "#A6CEE3",
    cpd  = "#B2E2B4"
  ),
  mid = list(
    gene = "#F7F7F7",
    cpd  = "#F7F7F7"
  ),
  high = list(
    gene = "#FDBF6F",
    cpd  = "#CAB2D6"
  ),
  
  limit = list(
    gene = 2,
    cpd  = 1.5
  ),
  
  same.layer = FALSE
)


# ============================================================================
# 6. Save input and mapped results
# ============================================================================

write_csv(
  gene_by_entrez,
  file.path(outdir, "Genes_selected_maxAbsLogFC_per_ENTREZ.csv")
)

write_csv(
  lipid_class_summary,
  file.path(outdir, "Buchnera_lipids_meanLogFC_by_class.csv")
)

if (!is.null(api_phospholipid)) {
  
  if (!is.null(api_phospholipid$plot.data.gene)) {
    write_csv(
      api_phospholipid$plot.data.gene,
      file.path(outdir, "api00564_mapped_gene_nodes.csv")
    )
  }
  
  if (!is.null(api_phospholipid$plot.data.cpd)) {
    write_csv(
      api_phospholipid$plot.data.cpd,
      file.path(outdir, "api00564_mapped_lipid_nodes.csv")
    )
  }
  
  cat("Pathway map completed successfully.\n")
  
} else {
  cat("Pathview failed. Check KEGG access and output paths.\n")
}



# Pathview node table
pathview_nodes <- api_phospholipid$plot.data.gene %>%
  dplyr::mutate(node_row = row_number())

# Expand all ENTREZ IDs mapped to each KEGG node
node_gene_candidates <- pathview_nodes %>%
  dplyr::select(
    node_row,
    node_kegg_names = kegg.names,
    node_label = labels,
    all.mapped
  ) %>%
  separate_rows(all.mapped, sep = ",\\s*") %>%
  dplyr::rename(ENTREZID = all.mapped) %>%
  dplyr::mutate(
    ENTREZID = trimws(as.character(ENTREZID))
  ) %>%
  dplyr::filter(
    !is.na(ENTREZID),
    ENTREZID != ""
  )

# Attach the actual logFC belonging to each ENTREZ ID
node_gene_candidates <- node_gene_candidates %>%
  left_join(
    gene_by_entrez %>%
      dplyr::select(ENTREZID, gene_id, logFC, adj_PValue),
    by = "ENTREZID"
  )

# Select the gene with the largest absolute logFC in each KEGG node
selected_gene_per_node <- node_gene_candidates %>%
  dplyr::filter(!is.na(logFC)) %>%
  group_by(node_row) %>%
  slice_max(
    order_by = abs(logFC),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  dplyr::rename(
    selected_ENTREZID = ENTREZID,
    selected_gene_id  = gene_id,
    selected_logFC    = logFC,
    selected_adjP     = adj_PValue
  )

write_csv(
  selected_gene_per_node,
  file.path(
    outdir,
    "api00564_selected_ENTREZ_maxAbsLogFC_per_node.csv"
  )
)


# Pathview node table
pathview_nodes <- api_phospholipid$plot.data.gene %>%
  dplyr::mutate(node_row = row_number())

# Expand all ENTREZ IDs mapped to each KEGG node
node_gene_candidates <- pathview_nodes %>%
  dplyr::select(
    node_row,
    node_kegg_names = kegg.names,
    node_label = labels,
    all.mapped
  ) %>%
  separate_rows(all.mapped, sep = ",\\s*") %>%
  dplyr::rename(ENTREZID = all.mapped) %>%
  dplyr::mutate(
    ENTREZID = trimws(as.character(ENTREZID))
  ) %>%
  dplyr::filter(
    !is.na(ENTREZID),
    ENTREZID != ""
  )

# Attach the actual logFC belonging to each ENTREZ ID
node_gene_candidates <- node_gene_candidates %>%
  left_join(
    gene_by_entrez %>%
      dplyr::select(ENTREZID, gene_id, logFC, adj_PValue),
    by = "ENTREZID"
  )

# Select the gene with the largest absolute logFC in each KEGG node
selected_gene_per_node <- node_gene_candidates %>%
  dplyr::filter(!is.na(logFC)) %>%
  group_by(node_row) %>%
  slice_max(
    order_by = abs(logFC),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  dplyr::rename(
    selected_ENTREZID = ENTREZID,
    selected_gene_id  = gene_id,
    selected_logFC    = logFC,
    selected_adjP     = adj_PValue
  )

write_csv(
  selected_gene_per_node,
  file.path(
    outdir,
    "api00564_selected_ENTREZ_maxAbsLogFC_per_node.csv"
  )
)

# ============================================================================
# Use selected ENTREZ genes and their corresponding logFC values as new input
# ============================================================================

selected_gene_input <- selected_gene_per_node %>%
  dplyr::transmute(
    ENTREZID = as.character(selected_ENTREZID),
    logFC    = selected_logFC
  ) %>%
  dplyr::filter(
    !is.na(ENTREZID),
    ENTREZID != "",
    !is.na(logFC)
  ) %>%
  # A gene may be present in more than one pathway node
  dplyr::distinct(
    ENTREZID,
    .keep_all = TRUE
  )

# Create named vector for Pathview
selected_gene_logFC <- selected_gene_input$logFC
names(selected_gene_logFC) <- selected_gene_input$ENTREZID

cat(
  "Selected ENTREZ genes used for remapping:",
  length(selected_gene_logFC),
  "\n"
)

cat(
  "Duplicated ENTREZ IDs:",
  anyDuplicated(names(selected_gene_logFC)),
  "\n"
)

# Check that every vector name corresponds to its logFC
print(selected_gene_logFC)

# ============================================================================
# Rerun Pathview using only selected ENTREZ genes
# ============================================================================

previous_wd <- getwd()
setwd(outdir)

api_phospholipid_selected <- pathview(
  gene.data   = selected_gene_logFC,
  cpd.data    = compound_logFC,
  
  pathway.id  = "00564",
  species     = "api",
  gene.idtype = "ENTREZ",
  
  # Safety measure if more than one selected gene maps to a node
  node.sum = "max.abs",
  
  kegg.native = TRUE,
  map.symbol  = TRUE,
  
  out.suffix = "2026_selected_ENTREZ_BuchneraLipids",
  kegg.dir   = outdir,
  
  low = list(
    gene = "#A6CEE3",
    cpd  = "#B2E2B4"
  ),
  mid = list(
    gene = "#F7F7F7",
    cpd  = "#F7F7F7"
  ),
  high = list(
    gene = "#FDBF6F",
    cpd  = "#CAB2D6"
  ),
  
  limit = list(
    gene = 2,
    cpd  = 1.5
  ),
  
  same.layer = FALSE
)

setwd(previous_wd)

readr::write_csv(
  selected_gene_input,
  file.path(
    outdir,
    "api00564_selected_ENTREZ_input.csv"
  )
)

if (
  !is.null(api_phospholipid_selected) &&
  !is.null(api_phospholipid_selected$plot.data.gene)
) {
  
  readr::write_csv(
    api_phospholipid_selected$plot.data.gene,
    file.path(
      outdir,
      "api00564_selected_ENTREZ_remapped_nodes.csv"
    )
  )
}
