rm(list = ls())

#library(GOSemSim)
library(lipidr)
library(dplyr)
#biopan: https://lipidmaps.org/biopan/pathway_analysis.php?id=fyPWt9HSTeAdkc3s
#https://lipidmaps.org/tool/biopan/doc/step3_build.html
#https://www.youtube.com/watch?v=EfTC_LuR3ak&ab_channel=LIPIDMAPS
#https://lora.metabolomics.fgu.cas.cz/

df <- read.csv('~/Desktop/20230522_BuchneraLipidomics/20240213_LCMS_Riken/20251002_LSMS_Riken/20240521_Riken_data_with_blank_Xbtt.csv',check.names = FALSE)
d = as_lipidomics_experiment(df, logged = FALSE, normalized = FALSE)
print(d)

options(max.print = 100, width = 100)

manual <- annotate_lipids(df, no_match = c("warn", "remove", "ignore"))
write.csv(manual, file="~/Desktop/20230522_BuchneraLipidomics/20240213_LCMS_Riken/20240219_manual_annotation.csv")

l=c("PE 14:0/18:3", "FA 14:0", "FA 16:0", "FA 18:0", "PE 36:2", "PE 18:0/18:2", "PC 18:2/18:2", "FA 18:2", "MG 18:0", "PE 18:2/18:2", "MG 16:0", "PE 18:1/18:2", "FA 18:1", "PC 18:1/18:2", "FA 32:0", "TG 46:1", "DG 18:1/18:1", "PI 18:1/18:2", "PE 16:0/17:1", "PE 20:0/18:2", "TG 46:2", "LPC 18:2", "LPE 18:2", "PC 14:0/18:2", "HexCer 19:2(O2)/16:1(O)", "DG 16:1/16:1", "PI 18:2/18:2", "PE 18:0/18:3", "DG 16:1/18:2", "PE 16:0/16:1", "PC 18:2/18:3", "PE 16:0/18:2", "TG 16:0/18:1/18:1", "PE 16:1/18:1", "PE 18:1/18:1", 
    "FA 28:0", "DG 18:2/18:2", "TG 16:0/18:0/18:1", "DG 18:0/18:2", "PE 16:0/18:1", "Cer 17:1(O2)/22:0", "LPE 18:0", "PE 18:2/18:3", "FA 20:0", "FA 18:3", "FA 16:1", "DG 18:0/18:1", "TG 50:2", "PI 18:1/18:1", "TG 18:0/18:1/18:1", "PG 16:0/17:1", "HexCer 17:1(O2)/20:0(O)", "FA 24:0", "PC 16:0/18:2", "MG 18:1", "TG 54:3", "PC 16:1/18:2", "PG 18:2/18:2", "DG 16:0/18:1", "MG 18:2", "HexCer 17:1(O2)/22:0(O)", "DG 16:0/16:1", "PS 18:2/18:2", "PE 18:0/18:1", "LPE 18:2", "TG 8:0/10:0/10:0", "TG 16:0/18:1/18:2", 
    "TG 14:0/18:1/18:2", "DG 18:1/18:2", "FA 26:0", "TG 50:1", "FA 34:0", "PI 34:2", "TG 14:0/16:0/18:1", "PC 16:1/18:1", "TG 12:0/14:0/18:1", "FA 30:0", "TG 18:1/18:1/18:2", "PE 17:1/18:1", "TG 18:0/18:0/18:1", "PS 18:0/18:2", "PE 32:1", "DG 14:0/14:0", "LPC 18:2", "DG 16:0/18:0", "TG 48:2", "PE 16:1/18:2", "LPE 18:1", "PG 16:0/18:1", "PC 14:0/18:1", "DG 8:0/14:0", "TG 14:0/14:0/14:0", "TG 16:0/18:1/18:3", "DG 20:0/18:2", "PE 20:0/18:3", "PG 16:0/18:2", "PG 34:2", "PG 16:1/18:1", "LPC 18:1", "SM 17:1(O2)/18+A65:0", 
    "PE 14:0/18:2", "PE 32:2", "PE 16:1/16:1", "PG 18:1/18:2", "PC 34:4", "Cer 17:1(O2)/18:0", "LPE 17:1", "PC 18:0/18:2", "TG 12:0/14:0/18:2", "PS 18:1/18:2", "PE 20:0/18:1", "PG 16:0/16:1", "TG 14:0/18:2/18:2", "TG 18:1/18:2/18:2", "LPE 16:0", "TG 8:0/14:0/18:1", "HexCer 17:1(O2)/24:0(O)", "TG 8:0/14:0/18:2", "LPE 16:1", "TG 8:0/8:0/8:0", "DG 16:1/18:3", "DG 18:2/18:3", "SPB 17:0(O2)", "HexCer 17:1(O2)/18:0(O)", "PC 34:1", "FA 22:0", "PC 14:0/18:3", "DG 20:0/18:3", "LPE 20:0", "MG 16:1", "DG 14:0/18:2", 
    "TG 42:1", "LPE 18:1", "SPB 17:1(O2)", "SM 17:1(O2)/20:0", "LPI 18:2", "TG 54:6", "PS 18:2/18:3", "Cer 17:1(O2)/18:0(2OH)", "LPS 18:2", "PC 18:3/18:3", "Cer 17:1(O2)/20:0", "TG 44:0", "PE O-20:0/18:2", "DG 18:0/18:0", "PS 16:0/18:2", "PG 17:1/18:1", "Cer 17:0(O2)/22:0", "Cer 17:1(O2)/20:0(2OH)", "HexCer 18:2(O2)/16:1(O)", "TG 42:2", "Cer 17:0(O2)/18:0", "PG 18:2/18:3", "PE 16:0/19:1", "MG 14:0", "Cer 17:1(O2)/22:0(2OH)", "PE 18:3/18:3", "LPI 16:0", "PE 16:0/16:0", "DG 32:3", "TG 16:0/18:2/18:3", 
    "LPI 18:1", "SM 17:1(O2)/22:0", "PE 14:0/16:0", "Cer 17:1(O2)/24:0(2OH)", "HexCer 17:1(O2)/20:0", "PE 34:4", "PI 16:0/18:1", "FA 41:0", "PE 19:0/18:2", "PI 18:2/18:3", "TG 14:0/16:1/18:2", "FA 24:1", "DG 14:0/18:1", "HexCer 17:0(O2)/22:0(O)", "HBMP 18:2/18:1/18:2", "PE 17:0/18:2", "FA 22:1", "SM 39:1(O3)", "ST 27:1(O);S", "LPC 18:1", "FA 17:0", "HBMP 18:2/18:2/18:2", "PE 14:0/17:1", "FA 28:2", "LPE 16:1", "FA 40:0", "CAR 18:1", "LPC 14:0", "PC 18:1/18:1", "FAHFA 20:0/24:0(O)", "PG 14:0/18:2", "FA 15:0", 
    "LPE 19:1", "TG 10:0/10:0/10:1", "Cer 17:0(O2)/22:0(2OH)", "PI 16:1/18:2", "LPC 16:1", "TG 10:0/10:0/10:0", "PE 32:3", "MG 20:0", "PS 18:0/18:1", "LPE 18:3", "FA 20:1", "PG 14:0/18:1", "PG 18:1/18:1", "TG 14:0/18:2/18:3", "CoQ9", "CAR 20:0", "PE 22:0/18:2", "PC 14:0/16:1", "FA 26:1", "FA 26:2", "SM 37:1(O3)", "TG 18:2/18:2/18:3", "DG 16:0/16:0", "FAHFA 20:0/26:0(O)", "PS 16:1/18:2", "FAHFA 20:0/28:0(O)", "CL 18:1/18:2/18:2/18:2", "PG 16:1/18:2", "CAR 18:2", "PG 16:0/19:1", "PC 20:0/18:2", "LPI 18:1", 
    "LPG 18:1", "LPS 18:1", "Cer 17:1(O2)/24:0", "PS 18:0/16:1", "LPI 18:2", "LPS 16:0", "PI 32:1", "PE O-20:0/18:1", "PC 18:0/18:1", "PE 18:1/20:1", "PC 16:1/16:1", "PE 18:0/17:1", "FA 42:0", "PI 14:0/18:2", "PE 30:1", "PE 33:2", "FAHFA 24:0/26:0(O)", "DG 18:3/18:3", "PG 17:1/19:1", "LPG 18:2", "LPG 17:1", "PE 16:1/26:1", "PE 24:0/18:2", "LPE 18:3", "FA 16:2", "FA 17:1", "LPC 16:0", "CL 18:1/18:2/18:2/18:3", "CL 18:1/18:1/18:2/18:2", "LPE 14:0", "PE O-18:0/18:2", "FAHFA 20:0/22:0(O)", "PS 14:0/18:2", 
    "TG 14:0/16:1/18:3", "FA 39:0", "FAHFA 26:0/26:0(O)", "FA 24:2", "PE O-18:0/18:1", "LPC 16:1", "PE 18:1/19:1", "PE 14:1/16:1", "PC 14:1/16:1", "LPS 18:0", "TG 16:0/18:3/18:3", "FA 38:0", "LPC 18:3", "LPI 18:0", "PS 34:4", "PI 18:0/18:1", "Cer 15:1(O2)/15:0", "FA 22:0(2OH)", "FAHFA 28:0/26:0(O)", "PG 32:2", "LPG 16:0", "PC 16:2/18:3", "Cer 17:0(O2)/20:0", "LPC 18:3", "Cer 17:0(O2)/24:0(2OH)", "FA 36:0", "PC 15:0/18:2", "PE 17:1/18:2", "LPC 18:0", "PG 14:0/16:0", "SL 17:0(O)/15:0", "LPG 18:1", "Cer 17:0(O2)/20:0(2OH)", 
    "Cer 18:0(O)3/24:0(2OH)", "FA 30:2", "PG 14:0/17:1", "PG 18:0/18:1", "SM 41:1(O2)", "TG 45:1", "Cer 16:0(O2)/22:0", "PS 18:3/18:3", "PC 14:0/16:2", "Cer 16:1(O2)/18:0", "PI 15:0/18:2", "PG 16:0/16:0", "PC 16:0/16:1", "LPG 18:2", "CL 16:0/18:1/18:2/18:2", "PE O-22:0/18:2", "LPE 16:2", "PG 34:4", "Cer 16:0(O2)/18:0", "PG 16:1/17:1", "PE 44:2", "FA 33:0", "CL 18:1/18:2/18:3/18:3", "PE O-20:0/18:3", "CL 14:0/18:1/18:2/18:2", "LPG 19:1", "FA 28:1", "FA 19:0", "SM 34:1(O2)", "CL 14:0/18:1/18:1/18:2", "PC 14:1/18:3", 
    "LPG 16:1", "PE O-18:0/18:3", "FA 23:0", "PE 14:1/18:3", "FA 37:0", "DLCL 18:1/18:2", "PI 18:3/18:3", "PC 20:0/18:3", "PC 16:1/26:1", "CL 16:1/18:1/18:2/18:2", "PE O-18:1/18:2", "CL 16:0/18:1/18:1/18:2", "LPS 16:1", "CL 14:0/18:1/18:2/18:3", "PI 34:4", "DLCL 18:1/18:1", "MG 22:0", "PE 16:0/17:2", "PE 26:0/16:1", "Cer 16:0(O2)/20:0", "Cer 17:0(O2)/16:0", "FA 20:2", "PC 16:1/28:1", "Cer 17:1(O2)/16:0(2OH)", "CE 16:1", "PG 18:3/18:3", "Cer 17:1(O2)/16:0", "LPE 16:2", "PE 38:4", "LPG 16:1", "FA 27:0", 
    "CL 14:0/16:0/18:1/18:2", "FA 35:0", "PG 14:0/18:3")


# Remove any remaining (O), (O2), and (O3)
l2 = gsub("\\(O\\)|\\(O2\\)|\\(O3\\)", "", l)
l2 = gsub("\\(O2\\)|\\(O\\)|\\(A65\\)", "", l2)
# replace underscores with spaces
l2 <- gsub("_", " ", l2)

# optional: add dummy chain for single-chain species (so lipidr can parse)
l2 <- gsub("Cer ([0-9]+:[0-9]+)$", "Cer \\1/0:0", l2)
l2 <- gsub("HexCer ([0-9]+:[0-9]+)$", "HexCer \\1/0:0", l2)

# Fix sterols
l2 <- gsub(";", "", l2)    
l2 <- gsub("ST ([0-9]+:[0-9]+)", "ST \\1/0:0", l2) 
l2 <- gsub("ST ([0-9]+:[0-9]+)/0:0S", "ST \\1/0:0", l2) # trailing S fix

# Fix quinones
l2 <- gsub("CoQ9", "CoQ9 0:0/0:0", l2)
# Fix quinones
l2 <- gsub("CoQ9", "CoQ9 0:0/0:0", l2)


# Fix sphingomyelin artifact
l2 <- gsub("\\+A65:0", "", l2)
# Fix sphingomyelin with incomplete fatty acyl
l2 <- gsub("SM ([0-9]+:[0-9]+)/([0-9]+)$", "SM \\1/\\2:0", l2)

l2 <- gsub("PE O-", "PE ", l2)   # ether-linked PE
l2 <- gsub("PC O-", "PC ", l2)
l2 <- gsub("PE P-", "PE ", l2)   # plasmalogen

df$Molecule <- l2
print(l2)
write.csv(l2, file = "~/Desktop/20230522_BuchneraLipidomics/20240213_LCMS_Riken/20251002_LSMS_Riken/20251002_lipidr_l2.csv")

# Annotate lipids using lipidr::annotate_lipids
annotation_result <- lipidr::annotate_lipids(l2)

# View the annotation result
print(annotation_result)


# Annotate the remaining lipid names
annotation_result <- lipidr::annotate_lipids(df$Molecule, no_match = "ignore")

write.csv(annotation_result,file = "~/Desktop/20230522_BuchneraLipidomics/20240213_LCMS_Riken/20251002_LSMS_Riken/20251002_sample_lipid_annotation.csv")


##Sample_annotation

d_annotated <- add_sample_annotation(d, "~/Desktop/20230522_BuchneraLipidomics/20240213_LCMS_Riken/20251002_LSMS_Riken/20240521_sample_annotation.csv")


d_normalized = normalize_pqn(d_annotated, measure = "Area", exclude = "blank", log = TRUE)
d_normalized$SampleType <- as.factor(d_normalized$SampleType)

# Extract the 'Area' data from the 'assays' slot
area_data <- d_normalized@assays@data@listData$Area

# Convert the matrix to a data frame
area_df <- as.data.frame(area_data)

# Add column names to the data frame
colnames(area_df) <- d_normalized@colData@rownames

# Extract the 'Molecule' information from the 'elementMetadata' slot
molecule_info <- d_normalized@elementMetadata@listData$Molecule

# Combine 'Molecule' information with 'Area' data
combined_df <- cbind(Molecule = molecule_info, area_df)


# Write the data frame to a CSV file
write.csv(combined_df, file = "~/Desktop/20230522_BuchneraLipidomics/20240213_LCMS_Riken/20251002_LSMS_Riken/20251002_d_normalized.csv")


##########################################heatmap+clustering_class_molecule##################
##########################################
# Heatmap + clustering by lipid class and molecule
##########################################
library(dplyr)
library(pheatmap)

# 1. Assign classes, preserving every input feature
sample_cols <- setdiff(names(combined_df), "Molecule")

combined_df2 <- combined_df %>%
  mutate(
    feature_id = row_number(),
    class_stub = sub("[[:space:]].*$", "", trimws(Molecule)),
    class_stub = case_when(
      grepl("^CoQ[0-9]*$", class_stub) ~ "CoQ",
      class_stub %in% c("PE-O", "PE-P") ~ "PE",
      TRUE ~ class_stub
    )
  )

# 2. Check class assignments
expected_classes <- c(
  "TG", "DG", "MG", "FA", "FAHFA", "CAR",
  "PC", "PE", "PI", "PS", "PG",
  "LPC", "LPE", "LPI", "LPS", "LPG",
  "CL", "DLCL", "HBMP",
  "Cer", "HexCer", "SM", "SPB", "SL",
  "ST", "CE", "CoQ"
)

unassigned <- combined_df2 %>%
  filter(is.na(class_stub) | !class_stub %in% expected_classes) %>%
  dplyr::select(feature_id, Molecule, class_stub)

print(unassigned)

if (nrow(unassigned) > 0) {
  stop("Review the unrecognized lipid names before class aggregation.")
}

# 3. Exclude CoQ and calculate class means and feature counts
combined_df2 <- combined_df2 %>%
  filter(class_stub != "CoQ")

class_df <- combined_df2 %>%
  group_by(class_stub) %>%
  summarise(
    across(all_of(sample_cols), ~ mean(.x, na.rm = TRUE)),
    n_lipids = n(),
    .groups = "drop"
  )

stopifnot(sum(class_df$n_lipids) == nrow(combined_df2))

print(class_df[, c("class_stub", "n_lipids")])

# 4. Convert sample columns to a matrix
class_mat <- as.matrix(
  class_df %>% dplyr::select(all_of(sample_cols))
)

rownames(class_mat) <- paste0(
  class_df$class_stub, " (", class_df$n_lipids, ")"
)

# Check values before clustering
if (any(!is.finite(class_mat))) {
  stop("Class means contain missing or non-finite values; review before plotting.")
}

if (any(apply(class_mat, 1, sd) == 0)) {
  stop("A lipid class has zero variance and cannot be row-scaled.")
}

# 5. Define sample colors
ann_colors <- list(
  SampleType = c(
    Apisum_WholeBody     = "#0072B2",
    Apisum_Appendages    = "#E69F00",
    Buchnera            = "#009E73",
    Apisum_Bacteriocytes = "#CC79A7"
  )
)

# 6. Match sample annotation to matrix columns
sample_annot <- data.frame(
  SampleType = d_normalized$SampleType,
  row.names = colnames(d_normalized)
)

stopifnot(all(colnames(class_mat) %in% rownames(sample_annot)))

sample_annot <- sample_annot[colnames(class_mat), , drop = FALSE]

# 7. Plot class-level heatmap
pheatmap(
  class_mat,
  scale = "row",
  annotation_col = sample_annot,
  annotation_colors = ann_colors,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  show_rownames = TRUE,
  show_colnames = TRUE,
  na_col = "grey90"
)
##########################################
# Molecule-level heatmap (supplementary)
##########################################

mol_mat <- as.matrix(combined_df[, -1])
rownames(mol_mat) <- combined_df$Molecule

# Optional: replace NA by row mean for molecule-level too
for (i in seq_len(nrow(mol_mat))) {
  nas <- is.na(mol_mat[i, ])
  if (any(nas)) {
    mol_mat[i, nas] <- mean(mol_mat[i, !nas], na.rm = TRUE)
  }
}

pheatmap(mol_mat,
         scale = "row",
         annotation_col = sample_annot,
         annotation_colors = ann_colors,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         show_rownames = FALSE,
         show_colnames = TRUE,
         na_col = "grey90")


##################PCA+permanova###########

#use_interactive_graphics(interactive = TRUE)
library(ggplot2)
library(vegan)

# 1. PCA with lipidr
mvaresults <- mva(
  d_normalized,
  measure   = "Area",
  method    = "PCA",
  group_col = "SampleType"
)

p <- plot_mva(
  mvaresults,
  color_by   = "SampleType",
  components = c(1,2),
  ellipse    = TRUE,
  hotelling  = TRUE
) +
  scale_color_manual(values = ann_colors$SampleType)  # keep your colourblind-safe palette

# 2. Run PERMANOVA
mol_mat <- as.matrix(combined_df[,-1])
rownames(mol_mat) <- combined_df$Molecule
mol_t <- t(mol_mat)

sample_annot <- data.frame(SampleType = d_normalized$SampleType)
rownames(sample_annot) <- colnames(mol_mat)

dist_mat <- dist(mol_t, method="euclidean")
permanova_res <- adonis2(dist_mat ~ SampleType, data=sample_annot, permutations=999)

# 3. Extract stats
pval <- signif(permanova_res$`Pr(>F)`[1], 3)
R2   <- round(permanova_res$R2[1], 3)

pval_label <- ifelse(pval < 0.001, "p < 0.001", paste0("p = ", pval))
p + labs(subtitle = bquote("PERMANOVA: R"^2 * " = " * .(R2) * ", " * .(pval_label)))


####################Differential_analysis_Buchnera_Bacteriocytes##############
#################### BEFORE differential analysis ####################

# 1. Add missing lipid class info into d_normalized metadata
meta <- as.data.frame(d_normalized@elementMetadata)

meta <- meta %>%
  mutate(
    Class = case_when(
      grepl("^Cer", Molecule) ~ "Cer",
      grepl("^SM", Molecule) ~ "SM",
      grepl("^HexCer", Molecule) ~ "HexCer",
      grepl("^PE O-|^PE P-", Molecule) ~ "PE",
      grepl("^CoQ", Molecule) ~ "CoQ",
      grepl("^ST", Molecule) ~ "ST",
      grepl("^SL", Molecule) ~ "SL",
      TRUE ~ Class
    ),
    total_cl = as.numeric(total_cl),
    total_cs = as.numeric(total_cs),
    chain_info = ifelse(is.na(total_cl) | is.na(total_cs), "Non-chain", "Chain"),
    total_cl = ifelse(is.na(total_cl), -1, total_cl),
    total_cs = ifelse(is.na(total_cs), -1, total_cs)
  )

# Convert back to DataFrame and re-assign
library(S4Vectors)
d_normalized@elementMetadata <- as(meta, "DataFrame")

####################run differential analysis ####################
library("ggplot2")
de_results <- de_analysis(
  data = d_normalized,
  Buchnera - Apisum_Bacteriocytes,
  measure = "Area"
)

# --- Step 3. Save clean DE results ---
write.csv(
  de_results,
  file = "~/Desktop/20230522_BuchneraLipidomics/20240213_LCMS_Riken/20251002_LSMS_Riken/20251002_Buchnera_Bacteriocytes/2025_deBchnera_Bacteriocytes.csv",
  row.names = FALSE
)

# --- Step 4. Visualize top 20 enriched lipids in Buchnera ---
top20_buchnera <- de_results %>%
  arrange(desc(logFC)) %>%
  head(21)

ggplot(top20_buchnera, aes(x = reorder(Molecule, logFC), y = logFC)) +
  geom_col(fill = "#2ca02c") +
  labs(
    x = "Lipid species",
    y = "log2 Fold Change"
  ) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(angle= 45, hjust = 1,size = 12),
    axis.title  = element_text(size = 14, face = "bold"),
    plot.title  = element_text(size = 16, face = "bold")
  )

ggplot(top20_buchnera,
       aes(x = logFC, y = reorder(Molecule, logFC))) +
  geom_col(fill = "#2ca02c") +
  labs(
    x = expression(log[2]~fold~change),
    y = "Lipid species"
  ) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold")
  )
buchnera_unique <- de_results %>%
  filter(logFC > 0) %>%
  group_by(Molecule) %>%
  summarise(
    logFC = mean(logFC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(logFC)) %>%
  slice_head(n = 20)

buchnera_unique <- buchnera_unique %>%
  arrange(logFC) %>%
  mutate(Molecule = factor(Molecule, levels = Molecule))

ggplot(buchnera_unique,
       aes(x = Molecule, y = logFC)) +
  geom_col(fill = "#2ca02c") +
  labs(
    x = "Lipid species",
    y = expression(log[2]~fold~change)
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 12
    ),
    axis.title = element_text(
      size = 14,
      face = "bold"
    )
  )

ggplot(buchnera_unique,
       aes(x = logFC, y = reorder(Molecule, logFC))) +
  geom_col(fill = "#2ca02c") +
  labs(
    x = expression(log[2]~fold~change),
    y = "Lipid species"
  ) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold")
  )


# --- Step 5. Trend plots ---
plot_trend(de_results, annotation = c("length"))
plot_trend(de_results, annotation = c("unsat"))

# --- Step 6. Lipid set enrichment analysis (LSEA) ---
enrich_results <- lsea(de_results, rank.by = "logFC")
enrich_results <- na.omit(enrich_results)

# Convert the leadingEdge column for export
enrich_results$leadingEdge <- sapply(enrich_results$leadingEdge, paste, collapse = ",")
write.csv(
  enrich_results,
  file = "~/Desktop/20230522_BuchneraLipidomics/20240213_LCMS_Riken/20251002_LSMS_Riken/20251002_Buchnera_Bacteriocytes/20251002_lsea_enrichBuchnera_Bacteriocytes_results_adjPValue.csv",
  row.names = FALSE
)

# --- Step 7. Volcano plot (no NA warnings anymore) ---
plot_results_volcano(de_results, show.labels = TRUE)

# --- Step 8. Enrichment plots ---
significant.sets <- significant_lipidsets(enrich_results, p.cutoff = 0.05)

plot_enrichment(de_results, significant.sets, annotation = "length")
plot_enrichment(de_results, significant.sets, annotation = "unsat")
plot_enrichment(de_results, significant.sets, annotation = "class")

# --- Step 9. Chain distribution plot ---
plot_chain_distribution(de_results)

#####################sub_phospholipid_analysis_#######################

library(dplyr)
library(tidyr)

# --- 1. Define phospholipid-related classes -------------------------------
phospho_classes <- c("PA", "DG", "PC", "PE", "PG", "PI", "PS",
                     "LPC", "LPE", "LPG", "LPI", "LPS",
                     "CL", "DLCL")

# --- 2. Check what classes exist in your dataset --------------------------
unique(rowData(d_normalized)$Class)

# --- 3. Filter phospholipid rows ------------------------------------------
phospho_idx <- rowData(d_normalized) %>%
  as.data.frame() %>%
  filter(Class %in% phospho_classes) %>%
  rownames()

# --- 4. Create d_phospho subset -------------------------------------------
d_phospho <- d_normalized[phospho_idx, ]

# --- 5. Verify -------------------------------------------------------------
dim(d_phospho)
table(rowData(d_phospho)$Class)

phospho_idx <- rowData(d_normalized) %>%
  as.data.frame() %>%
  filter(Class %in% phospho_classes) %>%
  rownames()

d_phospho <- d_normalized[phospho_idx, ]

unique(rowData(d_phospho)$Class)

# Extract assay data (matrix of Molecules × Samples)
area_mat <- assay(d_phospho, "Area")

# Get lipid class info
class_info <- rowData(d_phospho)$Class

# Convert to long data frame
area_df <- as.data.frame(area_mat)
area_df$Class <- class_info

area_long <- area_df %>%
  pivot_longer(-Class, names_to = "Sample", values_to = "Value")

# Add sample annotation
sample_annot <- as.data.frame(colData(d_phospho))
sample_annot$Sample <- rownames(sample_annot)

# Join with sample metadata
area_long <- left_join(area_long, sample_annot, by = "Sample")

# Summarize total intensity per Class × SampleType
class_summary <- area_long %>%
  group_by(Class, SampleType) %>%
  summarise(Total = sum(Value, na.rm = TRUE), .groups = "drop")


library(readr)

# --- Save area_long (per lipid × sample, tidy format) ---
write_csv(
  area_long,
  "~/Desktop/RNAseq_analysis_Apisum_nonrrna/20251025_figures/20251026_phospho_BuchVsbacter.csv")

########## PCA_PERMANOVA (Phospholipids) ##############
library(lipidr)
library(ggplot2)
library(vegan)
library(dplyr)

# --- 1. Run PCA --------------------------------------------------------
mvaresults_phospho <- mva(
  d_phospho,
  measure   = "Area",
  method    = "PCA",
  group_col = "SampleType"
)

# --- 2. Extract PCA scores ---------------------------------------------
scores <- as.data.frame(mvaresults_phospho$scores)
scores$SampleType <- as.factor(colData(d_phospho)$SampleType)

# --- 3. Identify PCA coordinate columns robustly -----------------------
pc_cols <- grep("^(PC|Comp\\.|Dim)", names(scores), value = TRUE)

if (length(pc_cols) < 2) {
  message("No PC-like columns found. Using first two numeric columns instead.")
  num_cols <- names(scores)[sapply(scores, is.numeric)]
  if (length(num_cols) < 2) stop("No numeric columns found in PCA scores.")
  pc_cols <- num_cols[1:2]
}

pc1_col <- pc_cols[1]
pc2_col <- pc_cols[2]
message("Detected PCA columns: ", pc1_col, " and ", pc2_col)

# --- 4. Compute variance explained manually ----------------------------
mol_mat <- assay(d_phospho, "Area")
mol_t   <- t(mol_mat)
pca_manual <- prcomp(mol_t, scale. = TRUE)
var_explained <- (pca_manual$sdev^2) / sum(pca_manual$sdev^2)
pc1_var <- round(var_explained[1] * 100, 1)
pc2_var <- round(var_explained[2] * 100, 1)

# --- 5. PERMANOVA ------------------------------------------------------
sample_annot <- as.data.frame(colData(d_phospho))
dist_mat <- dist(mol_t, method = "euclidean")
permanova_res <- adonis2(dist_mat ~ SampleType, data = sample_annot, permutations = 999)

pval <- signif(permanova_res$`Pr(>F)`[1], 3)
R2   <- round(permanova_res$R2[1], 3)
pval_label <- ifelse(pval < 0.001, "p < 0.001", paste0("p = ", pval))

# --- 6. Compute convex hulls safely -----------------------------------
scores[[pc1_col]] <- as.numeric(scores[[pc1_col]])
scores[[pc2_col]] <- as.numeric(scores[[pc2_col]])

hulls <- scores %>%
  group_by(SampleType) %>%
  filter(!is.na(!!sym(pc1_col)) & !is.na(!!sym(pc2_col))) %>%
  group_modify(~{
    if (nrow(.x) >= 3) .x[chull(.x[[pc1_col]], .x[[pc2_col]]), ]
    else .x
  }) %>%
  ungroup()

# --- 7. Plot PCA with convex hulls & PERMANOVA results -----------------
ggplot(scores, aes(x = !!sym(pc1_col), y = !!sym(pc2_col), color = SampleType)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_polygon(data = hulls, aes(fill = SampleType), color = NA, alpha = 0.25) +
  theme_bw(base_size = 14) +
  labs(
    x = paste0(pc1_col, " (", pc1_var, "%)"),
    y = paste0(pc2_col, " (", pc2_var, "%)"),
    title = paste0("PERMANOVA R² = ", R2, ", ", pval_label)
  ) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )


##########phospholipid_heatmap+cluestering########
class_heatmap <- class_summary %>%
  pivot_wider(names_from = SampleType, values_from = Total, values_fill = 0)

mat <- as.matrix(class_heatmap[,-1])
rownames(mat) <- class_heatmap$Class

library(pheatmap)
pheatmap(mat, scale = "row", cluster_rows = TRUE, angle_col = 45, cluster_cols = TRUE)


###############sub_phospho_differential_analysis_##############

de_resultsPhospho<- de_analysis(
  data=d_phospho, 
  Buchnera - Apisum_Bacteriocytes,
  measure="Area"
)
head(de_resultsPhospho)

write.csv(de_resultsPhospho, file = "~/Desktop/20230522_BuchneraLipidomics/20240213_LCMS_Riken/20251002_LSMS_Riken/20251002_Buchnera_Bacteriocytes/2025_PhosphodeBchnera_Bacteriocytes.csv", row.names = FALSE)



##############Top20_phospho_Buch_bact################


# Select top 20 enriched in Buchnera (sorted by log2FC)
# --- Ensure numeric and remove any weird attributes --------------------
de_df <- de_resultsPhospho %>%
  mutate(
    logFC = as.numeric(logFC),
    adj.P.Val = as.numeric(adj.P.Val)
  )

# --- Now filter and take top 20 ---------------------------------------
top20_buchnera <- de_df %>%
  filter(!is.na(adj.P.Val), adj.P.Val < 0.05, logFC > 1) %>%   # significant & enriched
  arrange(desc(logFC)) %>%
  slice_head(n = 21)   # safer alternative to slice(1:20)


ggplot(top20_buchnera, aes(x = reorder(Molecule, logFC), y = logFC)) +
  geom_col(fill = "#96cf96") +
  coord_flip() +
  labs(
    x = "Lipid species",
    y = "log2 Fold Change Phospholipids"
  ) +
  theme(
    axis.text.y = element_text(size = 12),   # y-axis tick labels (lipid names)
    axis.text.x = element_text(size = 12),   # x-axis tick labels
    axis.title  = element_text(size = 14, face = "bold"), # axis titles
    plot.title  = element_text(size = 16, face = "bold")  # if you add a title
  )


write.csv(top20_buchnera, file = "~/Desktop/20230522_BuchneraLipidomics/20240213_LCMS_Riken/20251002_LSMS_Riken/20251002_Buchnera_Bacteriocytes/20251006_dePhosphoBchnera_Bacteriocytestop20.csv", row.names = FALSE)

write.csv(top20_buchnera, file = "~/Desktop/20230522_BuchneraLipidomics/20240213_LCMS_Riken/20251002_LSMS_Riken/20251002_Buchnera_Bacteriocytes/20251008_d_phospho_sigAll_buchnera.csv", row.names = FALSE)


######try again######

# Ensure numeric
de_df <- de_df %>%
  mutate(
    logFC = as.numeric(logFC),
    adj.P.Val = as.numeric(adj.P.Val)
  )

# Filter significant Buchnera-enriched lipids,
# then collapse duplicate lipid species BEFORE selecting top 20
top20_buchnera <- de_df %>%
  filter(
    !is.na(adj.P.Val),
    adj.P.Val < 0.05,
    logFC > 1
  ) %>%
  group_by(Molecule) %>%
  summarise(
    logFC = mean(logFC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(logFC)) %>%
  slice_head(n = 20)

# Plot
ggplot(
  top20_buchnera,
  aes(x = logFC, y = reorder(Molecule, logFC))
) +
  geom_col(fill = "#96cf96") +
  labs(
    x = expression(log[2]~fold~change),
    y = "Lipid species"
  ) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title  = element_text(size = 14, face = "bold")
  )





###############################################################
#Module_positive_corr_enrich_Buch_BacteriocytesTop30_plusHull 

# --- 0. Load libraries --------------------------------------------------------
library(Hmisc)       # rcorr for correlation
library(igraph)      # network structure
library(ggraph)      # visualization
library(dplyr)
library(tidyr)
library(viridis)
library(dbscan)
library(ggforce)
library(tibble)
library(dynamicTreeCut)

# --- 1. Prepare lipid abundance matrix ---------------------------------------
mat <- assay(d_phospho, "Area")               # rows = lipids, cols = samples
lipid_names <- rowData(d_phospho)$Molecule

# --- 2. Compute Spearman correlations ----------------------------------------
cor_res <- rcorr(t(mat), type = "spearman")
cor_mat <- cor_res$r
p_mat   <- cor_res$P

# keep a copy for clustering (IMPORTANT)
cor_mat_full <- cor_mat

#Adjust p-values (BH) 

upper_idx <- upper.tri(p_mat, diag = FALSE)

p_vec <- p_mat[upper_idx]
p_adj_vec <- p.adjust(p_vec, method = "BH")

p_mat_adj <- matrix(1, nrow = nrow(p_mat), ncol = ncol(p_mat),
                    dimnames = dimnames(p_mat))
p_mat_adj[upper_idx] <- p_adj_vec
p_mat_adj[lower.tri(p_mat_adj)] <- t(p_mat_adj)[lower.tri(p_mat_adj)]
diag(p_mat_adj) <- 0

# --- 3. Keep only strong, positive, significant correlations -----------------
cor_mat[cor_mat < 0.85 | p_mat_adj > 0.05] <- 0
diag(cor_mat) <- 0
# --- 4. Build correlation network --------------------------------------------
graph <- graph_from_adjacency_matrix(cor_mat, mode = "undirected",
                                     weighted = TRUE, diag = FALSE)
V(graph)$name <- lipid_names
graph <- delete_vertices(graph, degree(graph) == 0)

V(graph)$Class <- rowData(d_phospho)$Class[
  match(V(graph)$name, lipid_names)
]
V(graph)$Class[is.na(V(graph)$Class)] <- "Unknown"
V(graph)$degree <- degree(graph)

# --- 5. Cluster lipids into modules (hierarchical clustering) ----------------
dist_mat <- as.dist(1 - abs(cor_mat_full))
hc <- hclust(dist_mat, method = "average")
module_assignments <- cutreeDynamic(
  dendro = hc,
  distM = as.matrix(dist_mat),
  deepSplit = 2,
  minClusterSize = 5
)    

names(module_assignments) <- lipid_names
module_assignments <- module_assignments[names(module_assignments) %in% V(graph)$name]
V(graph)$Module <- module_assignments[match(V(graph)$name, names(module_assignments))]

# --- 6. Compute module activity (mean abundance per module per sample) -------
module_activity <- lapply(unique(V(graph)$Module), function(m) {
  rows <- names(V(graph))[V(graph)$Module == m]
  tmp <- colMeans(mat[match(rows, lipid_names), , drop = FALSE])
  
  # Convert to data frame explicitly
  df <- as.data.frame(t(tmp))
  df$Module <- m
  df
}) %>% bind_rows()

# Reorder columns to have Module first
module_activity <- module_activity %>%
  relocate(Module)
# --- 7. Combine with sample annotation ---------------------------------------
if ("Sample" %in% colnames(sample_annot)) {
  sample_annot2 <- sample_annot
} else {
  sample_annot2 <- sample_annot %>% rownames_to_column(var = "Sample")
}

module_df_long <- module_activity %>%
  pivot_longer(-Module, names_to = "Sample", values_to = "Value") %>%
  left_join(sample_annot2, by = "Sample")

# --- 8. Subset to Buchnera & bacteriocytes -----------------------------------
module_df_subset <- module_df_long %>%
  filter(SampleType %in% c("Buchnera", "Apisum_Bacteriocytes"))

# --- 9. Welch’s t-test per module --------------------------------------------
test_results_bc_t <- module_df_subset %>%
  group_by(Module) %>%
  summarise(
    ttest = if (length(unique(na.omit(SampleType))) == 2)
      t.test(Value ~ SampleType, var.equal = FALSE)$p.value
    else NA_real_,
    mean_buchnera = mean(Value[SampleType == "Buchnera"], na.rm = TRUE),
    mean_bc       = mean(Value[SampleType == "Apisum_Bacteriocytes"], na.rm = TRUE),
    diff          = mean_buchnera - mean_bc,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(ttest, method = "BH"))

sig_modules <- test_results_bc_t %>%
  filter(p_adj < 0.05) %>%
  pull(Module)

V(graph)$is_sig <- V(graph)$Module %in% sig_modules
V(graph)$highlight <- ifelse(V(graph)$is_sig, "sig", "other")
V(graph)$label <- ""


# --- 1. Start clean from the original object ----------------------------
de_df <- as.data.frame(de_resultsPhospho)

# --- 2. Coerce every column to a simple atomic vector -------------------
de_df[] <- lapply(de_df, function(x) {
  if (inherits(x, "Rle")) as.vector(S4Vectors::runValue(x))
  else if (is.list(x)) unlist(x)
  else as.vector(x)
})

# --- 3. Ensure plain data.frame class ----------------------------------
class(de_df) <- "data.frame"


# --- 10. Integrate differential lipid results --------------------------------
# de_resultsPhospho must contain: Molecule, logFC, adj.P.Val

top30_up <- de_df %>%
  filter(!is.na(adj.P.Val), adj.P.Val < 0.05, logFC > 1) %>%
  arrange(desc(logFC)) %>%
  slice_head(n = 30) %>%
  pull(Molecule)

head(top30_up)


top30_down <- de_df %>%
  filter(!is.na(adj.P.Val), adj.P.Val < 0.05, logFC < -1) %>% ####### -1
  arrange(logFC) %>%          # ascending order: most negative first
  slice_head(n = 30) %>%      # safer than slice(1:30)
  pull(Molecule)


V(graph)$logFC <- de_resultsPhospho$logFC[
  match(V(graph)$name, de_resultsPhospho$Molecule)
]
V(graph)$sigdir <- case_when(
  V(graph)$name %in% top30_up   ~ "Buchnera-enriched",
  V(graph)$name %in% top30_down ~ "Bacteriocyte-enriched",
  TRUE ~ "Other"
)
V(graph)$label <- ifelse(V(graph)$sigdir != "Other", V(graph)$name, "")

# --- 11. Compute layout & hulls ----------------------------------------------
layout_df <- create_layout(graph, layout = "kk")

# Compact layout
center_x <- mean(layout_df$x)
center_y <- mean(layout_df$y)
compression_factor <- 0.6
layout_df$x <- center_x + compression_factor * (layout_df$x - center_x)
layout_df$y <- center_y + compression_factor * (layout_df$y - center_y)
set.seed(42)
layout_df$x <- layout_df$x + rnorm(nrow(layout_df), 0, 0.02)
layout_df$y <- layout_df$y + rnorm(nrow(layout_df), 0, 0.02)



# Select only enriched lipids
enrich_df <- layout_df %>%
  filter(sigdir %in% c("Buchnera-enriched", "Bacteriocyte-enriched")) %>%
  mutate(EnrichmentGroup = ifelse(sigdir == "Buchnera-enriched",
                                  "Buchnera", "Apisum_Bacteriocytes"),
         x_scaled = scale(x),
         y_scaled = scale(y))

# DBSCAN clustering per enrichment group
eps_val <- 0.6
minPts_val <- 3
enrich_df <- enrich_df %>%
  group_by(EnrichmentGroup) %>%
  group_modify(~{
    db <- dbscan(cbind(.x$x_scaled, .x$y_scaled), eps = eps_val, minPts = minPts_val)
    .x$cluster <- db$cluster
    .x
  }) %>%
  ungroup() %>%
  filter(cluster != 0) %>%
  mutate(cluster_id = paste(EnrichmentGroup, cluster, sep = "_"))

# --- 12. Define colors -------------------------------------------------------
modules <- sort(unique(V(graph)$Module))
palette_mod <- viridis(length(modules), option = "turbo")
names(palette_mod) <- modules
enrich_colors <- c("Buchnera" = "#009E73",
                   "Apisum_Bacteriocytes" = "#CC79A7")

# --- 13. Plot final network with hulls ---------------------------------------
p_final <- ggraph(layout_df) +
  # Hull outlines
  geom_mark_hull(
    data = enrich_df,
    aes(x = x, y = y, group = cluster_id, color = EnrichmentGroup),
    fill = NA, size = 0.9, expand = 0.03, alpha = 1
  ) +
  # Edges
  geom_edge_link(aes(width = abs(weight)), alpha = 0.15, colour = "grey70") +
  # Nodes
  geom_node_point(aes(color = factor(Module),
                      size = degree,
                      shape = sigdir,
                      alpha = highlight)) +
  # Labels
  geom_node_text(aes(label = label),
                 repel = TRUE, size = 3, color = "black",
                 fontface = "bold", max.overlaps = Inf) +
  # Color & shape scales
  scale_color_manual(values = c(palette_mod, enrich_colors),
                     breaks = c(names(enrich_colors), names(palette_mod)),
                     name = "Module / Enrichment") +
  scale_shape_manual(values = c("Other" = 16,
                                "Buchnera-enriched" = 17,
                                "Bacteriocyte-enriched" = 15),
                     name = "Differential lipids") +
  scale_alpha_manual(values = c("sig" = 1, "other" = 0.4), guide = "none") +
  scale_size(range = c(1, 4)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.justification = "center",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text  = element_text(size = 8),
    plot.margin  = margin(8, 8, 12,8)
  ) +
  guides(
    color = guide_legend(
      title = "Modules",
      nrow = 3, byrow = TRUE,
      label.position = "bottom",
      label.hjust = 0.1,
      override.aes = list(size = 2)
    ),
    shape = guide_legend(title = "Differential lipids")
  ) 

# --- 14. Print final plot -----------------------------------------------------
# Close graphics devices, then display the plot
graphics.off()
print(p_final)




table(V(graph)$Module)
module_table <- data.frame(
  Molecule = V(graph)$name,
  Module   = V(graph)$Module,
  Class    = V(graph)$Class,
  Degree   = V(graph)$degree,
  Sigdir   = V(graph)$sigdir
)
head(module_table)

# Create output folder
dir.create("~/Desktop/modules_lipid_lists", recursive = TRUE, showWarnings = FALSE)
modules <- sort(unique(V(graph)$Module))

for (m in modules) {
  mod_df <- module_table[module_table$Module == m, ]
  
  if (nrow(mod_df) == 0) next   #skip empty modules
  
  write.table(mod_df,
              file = file.path("~/Desktop/modules_lipid_lists", paste0("module_", m, ".txt")),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
}
sig_modules <- test_results_bc_t %>% filter(p_adj < 0.05) %>% pull(Module)

sig_lipids <- module_table[module_table$Module %in% sig_modules, ]
sig_lipids %>% dplyr::arrange(Module)
module_table %>%
  group_by(Module, Class) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = factor(Module), y = n, fill = Class)) +
  geom_col() +
  theme_minimal() +
  labs(x = "Module", y = "Number of lipids", title = "Lipid class composition per module")


write.csv(sig_lipids, file = "~/Desktop/20230522_BuchneraLipidomics/20240213_LCMS_Riken/20251002_LSMS_Riken/20251002_Buchnera_Bacteriocytes/20251010_sig_modules.csv")
write.csv(module_table, file = "~/Desktop/20230522_BuchneraLipidomics/20240213_LCMS_Riken/20251002_LSMS_Riken/20251002_Buchnera_Bacteriocytes/20251010_module_table.csv")
write.csv(modules, file = "~/Desktop/20230522_BuchneraLipidomics/20240213_LCMS_Riken/20251002_LSMS_Riken/20251002_Buchnera_Bacteriocytes/20251010_modules.csv")




#############DifferentialAnalysis_BuchneraVSwholebody##############

de_results <- de_analysis(
  data = d_normalized,
  Buchnera - Apisum_WholeBody,
  measure = "Area"
)

# --- Step 3. Save clean DE results ---
write.csv(
  de_results,
  file = "~/Desktop/20230522_BuchneraLipidomics/20240213_LCMS_Riken/20251002_LSMS_Riken/20251002_Buchnera_Bacteriocytes/2025_deBchnera_Apisum_WholeBody.csv",
  row.names = FALSE
)



# --- Step 5. Trend plots ---
plot_trend(de_results, annotation = c("length"))
plot_trend(de_results, annotation = c("unsat"))

# --- Step 6. Lipid set enrichment analysis (LSEA) ---
enrich_results <- lsea(de_results, rank.by = "logFC")
enrich_results <- na.omit(enrich_results)

# Convert the leadingEdge column for export
enrich_results$leadingEdge <- sapply(enrich_results$leadingEdge, paste, collapse = ",")
write.csv(
  enrich_results,
  file = "~/Desktop/20230522_BuchneraLipidomics/20240213_LCMS_Riken/20251002_LSMS_Riken/20251002_Buchnera_Bacteriocytes/20251002_lsea_enrichBuchnera_wholebody_results_adjPValue.csv",
  row.names = FALSE
)

# --- Step 7. Volcano plot (no NA warnings anymore) ---
plot_results_volcano(de_results, show.labels = TRUE)

# --- Step 8. Enrichment plots ---
significant.sets <- significant_lipidsets(enrich_results, p.cutoff = 0.05)

plot_enrichment(de_results, significant.sets, annotation = "length")
plot_enrichment(de_results, significant.sets, annotation = "unsat")
plot_enrichment(de_results, significant.sets, annotation = "class")

# --- Step 9. Chain distribution plot ---
plot_chain_distribution(de_results)



#############DifferentialAnalysis_BacteriocyteVSwholebody##############

de_results <- de_analysis(
  data = d_normalized,
  Apisum_Bacteriocytes - Apisum_WholeBody,
  measure = "Area"
)

# --- Step 3. Save clean DE results ---
write.csv(
  de_results,
  file = "~/Desktop/20230522_BuchneraLipidomics/20240213_LCMS_Riken/20251002_LSMS_Riken/20251002_Buchnera_Bacteriocytes/2025_deApisum_Bacteriocytes_Apisum_WholeBody.csv",
  row.names = FALSE
)



# --- Step 5. Trend plots ---
plot_trend(de_results, annotation = c("length"))
plot_trend(de_results, annotation = c("unsat"))

# --- Step 6. Lipid set enrichment analysis (LSEA) ---
enrich_results <- lsea(de_results, rank.by = "logFC")
enrich_results <- na.omit(enrich_results)

# Convert the leadingEdge column for export
enrich_results$leadingEdge <- sapply(enrich_results$leadingEdge, paste, collapse = ",")
write.csv(
  enrich_results,
  file = "~/Desktop/20230522_BuchneraLipidomics/20240213_LCMS_Riken/20251002_LSMS_Riken/20251002_Buchnera_Bacteriocytes/20251002_lsea_enrichApisum_Bacteriocytes_wholebody_results_adjPValue.csv",
  row.names = FALSE
)

# --- Step 7. Volcano plot (no NA warnings anymore) ---
plot_results_volcano(de_results, show.labels = TRUE)

# --- Step 8. Enrichment plots ---
significant.sets <- significant_lipidsets(enrich_results, p.cutoff = 0.05)

plot_enrichment(de_results, significant.sets, annotation = "length")
plot_enrichment(de_results, significant.sets, annotation = "unsat")
plot_enrichment(de_results, significant.sets, annotation = "class")

# --- Step 9. Chain distribution plot ---
plot_chain_distribution(de_results)


#################differentialAnalysis_appendagesVSWB#############


de_results <- de_analysis(
  data = d_normalized,
  Apisum_Appendages - Apisum_WholeBody,
  measure = "Area"
)

# --- Step 3. Save clean DE results ---
write.csv(
  de_results,
  file = "~/Desktop/20230522_BuchneraLipidomics/20240213_LCMS_Riken/20251002_LSMS_Riken/20251002_Buchnera_Bacteriocytes/2025_deApisum_Appendages_Apisum_WholeBody.csv",
  row.names = FALSE
)



# --- Step 5. Trend plots ---
plot_trend(de_results, annotation = c("length"))
plot_trend(de_results, annotation = c("unsat"))

# --- Step 6. Lipid set enrichment analysis (LSEA) ---
enrich_results <- lsea(de_results, rank.by = "logFC")
enrich_results <- na.omit(enrich_results)

# Convert the leadingEdge column for export
enrich_results$leadingEdge <- sapply(enrich_results$leadingEdge, paste, collapse = ",")
write.csv(
  enrich_results,
  file = "~/Desktop/20230522_BuchneraLipidomics/20240213_LCMS_Riken/20251002_LSMS_Riken/20251002_Buchnera_Bacteriocytes/20251002_lsea_enrichApisum_Appendages_wholebody_results_adjPValue.csv",
  row.names = FALSE
)

# --- Step 7. Volcano plot (no NA warnings anymore) ---
plot_results_volcano(de_results, show.labels = TRUE)

# --- Step 8. Enrichment plots ---
significant.sets <- significant_lipidsets(enrich_results, p.cutoff = 0.05)

plot_enrichment(de_results, significant.sets, annotation = "length")
plot_enrichment(de_results, significant.sets, annotation = "unsat")
plot_enrichment(de_results, significant.sets, annotation = "class")

# --- Step 9. Chain distribution plot ---
plot_chain_distribution(de_results)


######count checking ### untargeted 
library(VennDiagram)
library(grid)

# 1. Feature-level values used for the heatmap; CoQ already excluded
venn_mat <- as.matrix(
  combined_df2 %>% dplyr::select(all_of(sample_cols))
)

# Match metadata to sample columns
venn_annot <- sample_annot[colnames(venn_mat), , drop = FALSE]

group_map <- c(
  Apisum_WholeBody     = "Wb",
  Apisum_Appendages    = "App",
  Buchnera            = "B",
  Apisum_Bacteriocytes = "BC"
)

sample_groups <- unname(
  group_map[as.character(venn_annot$SampleType)]
)

print(table(sample_groups, useNA = "ifany"))

# 2. Calculate group means for each individual feature
groups <- c("Wb", "App", "B", "BC")

group_means <- sapply(groups, function(g) {
  idx <- which(sample_groups == g)
  
  if (length(idx) == 0) {
    stop(paste("No samples found for:", g))
  }
  
  rowMeans(venn_mat[, idx, drop = FALSE], na.rm = TRUE)
})

# 3. Apply the original group mean > 0 rule
above_zero <- is.finite(group_means) & group_means > 0

# Preserve duplicate names as separate features, as in the heatmap
feature_ids <- as.character(combined_df2$feature_id)

metabolite_list <- setNames(
  lapply(seq_along(groups), function(j) {
    feature_ids[above_zero[, j]]
  }),
  groups
)

print(lengths(metabolite_list))

shared_ids <- Reduce(intersect, metabolite_list)

cat("Features included:", nrow(venn_mat), "\n")
cat("Features shared across all four groups:",
    length(shared_ids), "\n")

# Inspect group means alongside original names
venn_table <- data.frame(
  feature_id = feature_ids,
  Molecule = combined_df2$Molecule,
  group_means,
  shared_all = feature_ids %in% shared_ids,
  check.names = FALSE
)

# 4. Draw the Venn diagram
venn.plot <- venn.diagram(
  x = metabolite_list,
  category.names = names(metabolite_list),
  filename = NULL,
  col = "transparent",
  fill = c("#0072B2", "#E69F00", "#009E73", "#CC79A7"),
  alpha = 0.50,
  cex = 1.5,
  fontfamily = "sans",
  cat.cex = 1.3,
  cat.fontfamily = "sans",
  cat.col = "black",
  cat.dist = rep(0.05, 4),
  margin = 0.10
)

grid.newpage()
grid.draw(venn.plot)
