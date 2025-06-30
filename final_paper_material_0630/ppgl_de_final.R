library(tidyverse)
library(reticulate)

# Specify the path to the .h5ad file
h5ad_file <- "./ppgl_pbulk.h5ad"

# Read the .h5ad file using reticulate
anndata <- import("anndata")
py_index <- import("pandas")$Index
adata <- anndata$read_h5ad(h5ad_file)

counts <- t(adata$X) # Transpose to convert to column-major order
meta_data <- as.data.frame(adata$obs) # Convert obs to a data frame
vars <- adata$var
vars <-  py_to_r(adata$var_names$to_list())

expression_table = counts
rownames(expression_table) <- vars
expression_table <- expression_table[, meta_data$Healthy == "False"]
meta_data <- meta_data[meta_data$Healthy == "False", ]

meta_data$assignments <- meta_data$broad_cluster

#remove where SC_SN is SC
expression_table <- expression_table[, meta_data$SC_SN != "SC"]
meta_data <- meta_data[meta_data$SC_SN != "SC", ]



table(meta_data$assignments)
#remove those cell_type that have a single entry
singlets = meta_data %>%
  group_by(assignments) %>%
  filter(n() < 2) %>%
  pull(broad_cluster)

expression_table = expression_table[,!meta_data$broad_cluster %in% singlets]
meta_data = meta_data %>%
  filter(!broad_cluster %in% singlets)

table(meta_data$assignments)


###################################################
### Now with Limma-voom below
library(limma)
library(edgeR)
# Ensure valid column names for assignments and exptype
meta_data$assignments <- make.names(meta_data$assignments)

dge = DGEList(counts = expression_table, group = meta_data$assignments)

#remove genes expressed in less than 10% of the samples or with less than 1000 total counts
design <- model.matrix(~ 0 + assignments, data = meta_data)

sum(expression_table["SOX2",])
keep_genes <- filterByExpr(dge, design, min.total.count = 100,min.count=0,min.prop=0)
sum(keep_genes)
keep_genes["SOX2"]

dge <- dge[keep_genes, ]

sox2_values <- dge$counts["SOX2", ] / colSums(dge$counts)
assignments <- meta_data$assignments
#barchart
barplot(sox2_values, names.arg = assignments, las = 2, main = "SOX2 Expression by Assignment", col = "lightblue", ylab = "Expression Level")




#norm
dge <- calcNormFactors(dge, method = "TMM")

fit <- voom(dge, design, plot=TRUE)

fit <- lmFit(fit, design)


# Get all unique assignment levels
assignment_levels <- unique(meta_data$assignments)
all_results <- list()
all_results_sig <- list()

# Iterate through each assignment level as the "given" cluster
for (given_cluster in assignment_levels) {
  # Iterate through all other assignment levels for pairwise comparisons
  for (other_cluster in assignment_levels) {
    if (given_cluster == other_cluster) {
      next # Skip comparison of a cluster with itself
    }
    
    # Create contrast name for the current pairwise comparison
    contrast_name <- paste0(given_cluster, "_vs_", other_cluster)
    
    # Create the contrast string
    contrast_string <- paste0("assignments", given_cluster, " - assignments", other_cluster)
    
    # Create contrast
    contrast <- makeContrasts(contrast_string, levels = design)
    
    # Compute differential expression
    diff_expression <- contrasts.fit(fit, contrast)
    diff_expression <- eBayes(diff_expression, trend = TRUE)
    
    # Get top differentially expressed genes
    top_genes <- topTable(diff_expression, number = Inf, sort.by = "P")
    top_genes$genes <- rownames(top_genes)
    top_genes$comparison <- contrast_name
    all_results[[contrast_name]] <- top_genes
    
    # Filter for significant genes (e.g., adjusted p-value < 0.05)
    top_genes_sig <- top_genes %>%
      filter(adj.P.Val < 0.05)
    all_results_sig[[contrast_name]] <- top_genes_sig
  }
}



# Identify genes that are significantly higher in the given cluster in all but one comparison
consistent_markers <- list()
for (given_cluster in assignment_levels) {
  significant_up_in_most <- c()
  gene_list <- unique(unlist(sapply(all_results_sig, rownames))) # Get all significant genes across all comparisons
  
  for (gene in gene_list) {
    significant_comparisons <- 0
    non_significant_comparisons <- 0
    significant_and_higher <- 0
    
    for (other_cluster in assignment_levels) {
      if (given_cluster == other_cluster) {
        #continue
        next # Skip comparison of a cluster with itself
      }
      comparison_name <- paste0(given_cluster, "_vs_", other_cluster)
      if (comparison_name %in% names(all_results_sig)) {
        if (gene %in% rownames(all_results_sig[[comparison_name]])) {
          significant_comparisons <- significant_comparisons + 1
          logfc <- all_results[[comparison_name]][gene, "logFC"]
          padj <- all_results[[comparison_name]][gene, "adj.P.Val"]
          if  (logfc > 0 & padj < 0.05) {
            significant_and_higher <- significant_and_higher + 1
          }
        } else {
          non_significant_comparisons <- non_significant_comparisons + 1
        }
      }
    }
    
    # Check if the gene is significantly higher in the given cluster in all but one comparison
    if (significant_comparisons == 5 & significant_and_higher == significant_comparisons) {
      significant_up_in_most <- c(significant_up_in_most, gene)
      
    }
  }
  consistent_markers[[given_cluster]] <- significant_up_in_most
}

all_results_df <- do.call(rbind, all_results)
all_results_sig_df <- do.call(rbind, all_results_sig)
# Add a column for the comparison name
consistent_markers_df <- data.frame(
  Gene = unlist(consistent_markers),
  Cluster = rep(names(consistent_markers), sapply(consistent_markers, length)))

consistent_markers_df
#make sure Gene is string
consistent_markers_df$Gene <- as.character(consistent_markers_df$Gene)
#save 
save_fold = "./"
write.csv(consistent_markers_df, paste0(save_fold,"consistent_markers.csv"))

#reset rownames

rownames(all_results_df) <- seq(1, nrow(all_results_df))

write.csv(all_results_df, paste0(save_fold,"all_results.csv"))
write.csv(all_results_sig_df, paste0(save_fold,"all_results_sig.csv"))





#Healthy vs PPGL comparisons

library(tidyverse)
library(reticulate)

# Specify the path to the .h5ad file
h5ad_file <- "./ppgl_pbulk.h5ad"

# Read the .h5ad file using reticulate
anndata <- import("anndata")
py_index <- import("pandas")$Index
adata <- anndata$read_h5ad(h5ad_file)

counts <- t(adata$X) # Transpose to convert to column-major order
meta_data <- as.data.frame(adata$obs) # Convert obs to a data frame
vars <- adata$var
vars <-  py_to_r(adata$var_names$to_list())

expression_table = counts
rownames(expression_table) <- vars

meta_data$assignments <- meta_data$broad_cluster

#remove those cell_type that have a single entry
singlets = meta_data %>%
  group_by(assignments) %>%
  filter(n() < 2) %>%
  pull(broad_cluster)

expression_table = expression_table[,!meta_data$broad_cluster %in% singlets]
meta_data = meta_data %>%
  filter(!broad_cluster %in% singlets)


###################################################
### Now with Limma-voom below
library(limma)
library(edgeR)
# Ensure valid column names for assignments and exptype
meta_data$assignments <- make.names(meta_data$assignments)

dge = DGEList(counts = expression_table, group = meta_data$assignments)


design <- model.matrix(~ 0 + assignments + assignments:Healthy, data = meta_data)

#makenames for terms
colnames(design)<-make.names(colnames(design))

keep_genes <- filterByExpr(dge, design, min.total.count = 100,min.count=0,min.prop=0)

dge <- dge[keep_genes, ]

#norm
dge <- calcNormFactors(dge, method = "TMM")

fit <- voom(dge, design, plot=TRUE)

fit <- lmFit(fit, design)

contrast <- makeContrasts(
  assignmentsSustentacular_cells.HealthyTrue,
  levels = design
)

# Compute differential expression
diff_expression <- contrasts.fit(fit, contrast)
diff_expression <- eBayes(diff_expression, robust=TRUE)

# Get top differentially expressed genes
top_genes <- topTable(diff_expression, number = Inf)
top_genes$genes <- rownames(top_genes)
print(rownames((top_genes)))


contrast <- makeContrasts(
  assignmentsChromaffin_cells.HealthyTrue,
  levels = design
)

# Compute differential expression
diff_expression <- contrasts.fit(fit, contrast)
diff_expression <- eBayes(diff_expression, robust=TRUE)

# Get top differentially expressed genes
top_genes <- topTable(diff_expression, number = Inf)
top_genes$genes <- rownames(top_genes)
print(rownames((top_genes)))


#iterate through
coefs = c("assignmentsChromaffin_cells.HealthyTrue","assignmentsEndothelial_cells.HealthyTrue",  
"assignmentsFibroblasts.HealthyTrue","assignmentsMacrophages.HealthyTrue",        
"assignmentsSustentacular_cells.HealthyTrue","assignmentsT_cells.HealthyTrue")

final_res <- list()
final_res_sig <- list()
for (coef in coefs) {
  contrast <- makeContrasts(contrasts = coef, levels = design)
  
  
  # Compute differential expression
  diff_expression <- contrasts.fit(fit, contrast)
  diff_expression <- eBayes(diff_expression, robust=TRUE)
  
  # Get top differentially expressed genes
  top_genes <- topTable(diff_expression, number = Inf)
  top_genes$genes <- rownames(top_genes)
  top_genes$comparison <- coef
  
  print(rownames((top_genes)))
  final_res[[coef]] <- top_genes
  final_res_sig[[coef]] <- top_genes %>%
    filter(adj.P.Val < 0.05)
}

# Combine all results into a single data frame
final_res_df <- do.call(rbind, final_res)
final_res_sig_df <- do.call(rbind, final_res_sig)
# Save the final results to a CSV file
save_fold = "./"

write.csv(final_res_df, paste0(save_fold,"all_results_tumor_vs_healthy.csv"))
write.csv(final_res_sig_df, paste0(save_fold,"all_results_sig_tumor_vs_healthy.csv"))
























##SOX2+ vs SOX2- comparisons



library(tidyverse)
library(reticulate)

# Specify the path to the .h5ad file
h5ad_file <- "./ppgl_pbulk_chrom_sox2pos_neg.h5ad"

# Read the .h5ad file using reticulate
anndata <- import("anndata")
py_index <- import("pandas")$Index
adata <- anndata$read_h5ad(h5ad_file)

counts <- t(adata$X) # Transpose to convert to column-major order
meta_data <- as.data.frame(adata$obs) # Convert obs to a data frame
vars <- adata$var
vars <-  py_to_r(adata$var_names$to_list())







counts <- t(adata$X) # Transpose to convert to column-major order
meta_data <- as.data.frame(adata$obs) # Convert obs to a data frame
vars <- adata$var
vars <-  py_to_r(adata$var_names$to_list())

expression_table = counts
rownames(expression_table) <- vars

meta_data$assignments <- meta_data$broad_cluster

#remove those cell_type that have a single entry
singlets = meta_data %>%
  group_by(assignments) %>%
  filter(n() < 2) %>%
  pull(broad_cluster)

expression_table = expression_table[,!meta_data$broad_cluster %in% singlets]
meta_data = meta_data %>%
  filter(!broad_cluster %in% singlets)
         
         
         
         
#remove those that are Healthy
expression_table <- expression_table[, meta_data$Healthy == "False"]
meta_data <- meta_data[meta_data$Healthy == "False", ]
dge = DGEList(counts = expression_table, group = meta_data$assignments)

#remove genes expressed in less than 10% of the samples or with less than 1000 total counts
design <- model.matrix(~ 0 + assignments + assignments:Chrom_SOX2_status, data = meta_data)

#makenames for terms
colnames(design)<-make.names(colnames(design))

keep_genes <- filterByExpr(dge, design, min.total.count = 200)

dge <- dge[keep_genes, ]

#norm
dge <- calcNormFactors(dge, method = "TMM")

fit <- voom(dge, design, plot=TRUE)

fit <- lmFit(fit, design)

contrast <- makeContrasts(
  assignmentsChromaffin_cells.Chrom_SOX2_statusLow,
  levels = design
)

# Compute differential expression
diff_expression <- contrasts.fit(fit, contrast)
diff_expression <- eBayes(diff_expression, robust=TRUE)

# Get top differentially expressed genes
top_genes <- topTable(diff_expression, number = Inf)
top_genes$genes <- rownames(top_genes)
print(rownames((top_genes)))



contrast <- makeContrasts(
  assignmentsChromaffin_cells.Chrom_SOX2_statusLow,
  levels = design
)

# Compute differential expression
diff_expression <- contrasts.fit(fit, contrast)
diff_expression <- eBayes(diff_expression, robust=TRUE)

# Get top differentially expressed genes
top_genes <- topTable(diff_expression, number = Inf)
top_genes$genes <- rownames(top_genes)
print(rownames((top_genes)))



coefs = c("assignmentsChromaffin_cells.Chrom_SOX2_statusLow")
         # "assignmentsEndothelial_cells.Chrom_SOX2_statusLow",  
          #"assignmentsFibroblasts.Chrom_SOX2_statusLow","assignmentsMacrophages.Chrom_SOX2_statusLow",        
          #"assignmentsSustentacular_cells.Chrom_SOX2_statusLow","assignmentsT_cells.Chrom_SOX2_statusLow")



final_res <- list()
final_res_sig <- list()
for (coef in coefs) {
  contrast <- makeContrasts(contrasts = coef, levels = design)
  
  
  # Compute differential expression
  diff_expression <- contrasts.fit(fit, contrast)
  diff_expression <- eBayes(diff_expression, robust=TRUE)
  
  # Get top differentially expressed genes
  top_genes <- topTable(diff_expression, number = Inf)
  top_genes$genes <- rownames(top_genes)
  top_genes$comparison <- coef
  
  print(rownames((top_genes)))
  final_res[[coef]] <- top_genes
  final_res_sig[[coef]] <- top_genes %>%
    filter(adj.P.Val < 0.05)
}

# Combine all results into a single data frame
final_res_df <- do.call(rbind, final_res)
final_res_sig_df <- do.call(rbind, final_res_sig)
# Save the final results to a CSV file
save_fold = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/ppgl/"

write.csv(final_res_df, paste0(save_fold,"all_results_chrom_sox2pos_neg.csv"))
write.csv(final_res_sig_df, paste0(save_fold,"all_results_sig_chrom_sox2pos_neg.csv"))








### Markers just in healthy



##### Healthy vs Non-metastatic

library(tidyverse)
library(reticulate)

# Specify the path to the .h5ad file
h5ad_file <- "./ppgl_pbulk_0624.h5ad"

# Read the .h5ad file using reticulate
anndata <- import("anndata")
py_index <- import("pandas")$Index
adata <- anndata$read_h5ad(h5ad_file)

counts <- t(adata$X) # Transpose to convert to column-major order
meta_data <- as.data.frame(adata$obs) # Convert obs to a data frame
vars <- adata$var
vars <-  py_to_r(adata$var_names$to_list())

expression_table = counts
rownames(expression_table) <- vars

meta_data$assignments <- meta_data$broad_cluster

#remove those cell_type that have a single entry
singlets = meta_data %>%
  group_by(assignments) %>%
  filter(n() < 2) %>%
  pull(broad_cluster)

expression_table = expression_table[,!meta_data$broad_cluster %in% singlets]
meta_data = meta_data %>%
  filter(!broad_cluster %in% singlets)


###################################################
### Now with Limma-voom below
library(limma)
library(edgeR)
# Ensure valid column names for assignments and exptype
meta_data$assignments <- make.names(meta_data$assignments)



#if Metastasized_to is not nan change it to 1
meta_data$Metastasis <- ifelse(meta_data$Metastasized_to=="nan", 0, 1)

#keep where its either healthy or metastasis is 1
expression_table <- expression_table[, meta_data$Healthy == "True" ]
meta_data <- meta_data[meta_data$Healthy == "True", ]

dge = DGEList(counts = expression_table, group = meta_data$assignments)


design <- model.matrix(~ 0 + assignments, data = meta_data)

#makenames for terms
colnames(design)<-make.names(colnames(design))

keep_genes <- filterByExpr(dge, design, min.total.count =100,min.count=0,min.prop=0)

dge <- dge[keep_genes, ]

#norm
dge <- calcNormFactors(dge, method = "TMM")

fit <- voom(dge, design, plot=TRUE)

fit <- lmFit(fit, design)

contrast <- makeContrasts(
  assignmentsSustentacular_cells,
  levels = design
)

# Compute differential expression
diff_expression <- contrasts.fit(fit, contrast)
diff_expression <- eBayes(diff_expression, robust=TRUE)

# Get top differentially expressed genes
top_genes <- topTable(diff_expression, number = Inf)
top_genes$genes <- rownames(top_genes)
print(rownames((top_genes)))


contrast <- makeContrasts(
  assignmentsCortex_cells,
  levels = design
)


# Compute differential expression
diff_expression <- contrasts.fit(fit, contrast)
diff_expression <- eBayes(diff_expression, robust=TRUE)

# Get top differentially expressed genes
top_genes <- topTable(diff_expression, number = Inf)
top_genes$genes <- rownames(top_genes)
print(rownames((top_genes)))

cortex_markers <- top_genes %>%
  filter(adj.P.Val < 0.05 & logFC > 0)

#take top 2000
cortex_markers <- cortex_markers %>%
  arrange(desc(logFC)) %>%
  head(2000)



cortex_markers <- cortex_markers$genes

#iterate through
coefs = c("assignmentsChromaffin_cells.HealthyTrue","assignmentsEndothelial_cells.HealthyTrue",  
          "assignmentsFibroblasts.HealthyTrue","assignmentsMacrophages.HealthyTrue",        
          "assignmentsSustentacular_cells.HealthyTrue","assignmentsT_cells.HealthyTrue")

final_res <- list()
final_res_sig <- list()
for (coef in coefs) {
  contrast <- makeContrasts(contrasts = coef, levels = design)
  
  
  # Compute differential expression
  diff_expression <- contrasts.fit(fit, contrast)
  diff_expression <- eBayes(diff_expression, robust=TRUE)
  
  # Get top differentially expressed genes
  top_genes <- topTable(diff_expression, number = Inf)
  top_genes$genes <- rownames(top_genes)
  top_genes$comparison <- coef
  
  print(rownames((top_genes)))
  final_res[[coef]] <- top_genes
  final_res_sig[[coef]] <- top_genes %>%
    filter(adj.P.Val < 0.05)
}

# Combine all results into a single data frame
final_res_df <- do.call(rbind, final_res)
final_res_sig_df <- do.call(rbind, final_res_sig)
# Save the final results to a CSV file
save_fold = "./"

write.csv(final_res_df, paste0(save_fold,"all_results_tumor_vs_healthy.csv"))
write.csv(final_res_sig_df, paste0(save_fold,"all_results_sig_tumor_vs_healthy.csv"))










##### Healthy vs metastatic

library(tidyverse)
library(reticulate)

# Specify the path to the .h5ad file
h5ad_file <- "./ppgl_pbulk_0624.h5ad"

# Read the .h5ad file using reticulate
anndata <- import("anndata")
py_index <- import("pandas")$Index
adata <- anndata$read_h5ad(h5ad_file)

counts <- t(adata$X) # Transpose to convert to column-major order
meta_data <- as.data.frame(adata$obs) # Convert obs to a data frame
vars <- adata$var
vars <-  py_to_r(adata$var_names$to_list())

expression_table = counts
rownames(expression_table) <- vars

meta_data$assignments <- meta_data$broad_cluster

#remove those cell_type that have a single entry
singlets = meta_data %>%
  group_by(assignments) %>%
  filter(n() < 2) %>%
  pull(broad_cluster)

expression_table = expression_table[,!meta_data$broad_cluster %in% singlets]
meta_data = meta_data %>%
  filter(!broad_cluster %in% singlets)


###################################################
### Now with Limma-voom below
library(limma)
library(edgeR)
# Ensure valid column names for assignments and exptype
meta_data$assignments <- make.names(meta_data$assignments)



#if Metastasized_to is not nan change it to 1
meta_data$Metastasis <- ifelse(meta_data$Metastasized_to=="nan", 0, 1)

#keep where its either healthy or metastasis is 1
expression_table <- expression_table[, meta_data$Healthy == "True" | meta_data$Metastasis == 1]
meta_data <- meta_data[meta_data$Healthy == "True" | meta_data$Metastasis == 1, ]

dge = DGEList(counts = expression_table, group = meta_data$assignments)


design <- model.matrix(~ 0 + assignments + assignments:Metastasis, data = meta_data)

#makenames for terms
colnames(design)<-make.names(colnames(design))

keep_genes <- filterByExpr(dge, design, min.total.count = 100,min.count=0,min.prop=0)

dge <- dge[keep_genes, ]

#norm
dge <- calcNormFactors(dge, method = "TMM")

fit <- voom(dge, design, plot=TRUE)

fit <- lmFit(fit, design)

contrast <- makeContrasts(
  assignmentsSustentacular_cells.Metastasis,
  levels = design
)

# Compute differential expression
diff_expression <- contrasts.fit(fit, contrast)
diff_expression <- eBayes(diff_expression, robust=TRUE)

# Get top differentially expressed genes
top_genes <- topTable(diff_expression, number = Inf)
top_genes$genes <- rownames(top_genes)
print(rownames((top_genes)))


contrast <- makeContrasts(
  assignmentsSustentacular_cells.Metastasis,
  levels = design
)

# Compute differential expression
diff_expression <- contrasts.fit(fit, contrast)
diff_expression <- eBayes(diff_expression, robust=TRUE)

# Get top differentially expressed genes
top_genes <- topTable(diff_expression, number = Inf)
top_genes$genes <- rownames(top_genes)
print(rownames((top_genes)))


#iterate through
coefs = c("assignmentsChromaffin_cells.HealthyTrue","assignmentsEndothelial_cells.HealthyTrue",  
          "assignmentsFibroblasts.HealthyTrue","assignmentsMacrophages.HealthyTrue",        
          "assignmentsSustentacular_cells.HealthyTrue","assignmentsT_cells.HealthyTrue")

final_res <- list()
final_res_sig <- list()
for (coef in coefs) {
  contrast <- makeContrasts(contrasts = coef, levels = design)
  
  
  # Compute differential expression
  diff_expression <- contrasts.fit(fit, contrast)
  diff_expression <- eBayes(diff_expression, robust=TRUE)
  
  # Get top differentially expressed genes
  top_genes <- topTable(diff_expression, number = Inf)
  top_genes$genes <- rownames(top_genes)
  top_genes$comparison <- coef
  
  print(rownames((top_genes)))
  final_res[[coef]] <- top_genes
  final_res_sig[[coef]] <- top_genes %>%
    filter(adj.P.Val < 0.05)
}

# Combine all results into a single data frame
final_res_df <- do.call(rbind, final_res)
final_res_sig_df <- do.call(rbind, final_res_sig)
# Save the final results to a CSV file
save_fold = "./"

write.csv(final_res_df, paste0(save_fold,"all_results_tumor_vs_healthy.csv"))
write.csv(final_res_sig_df, paste0(save_fold,"all_results_sig_tumor_vs_healthy.csv"))












##### Healthy vs Non-metastatic

library(tidyverse)
library(reticulate)

# Specify the path to the .h5ad file
h5ad_file <- "./ppgl_pbulk_0624.h5ad"

# Read the .h5ad file using reticulate
anndata <- import("anndata")
py_index <- import("pandas")$Index
adata <- anndata$read_h5ad(h5ad_file)

counts <- t(adata$X) # Transpose to convert to column-major order
meta_data <- as.data.frame(adata$obs) # Convert obs to a data frame
vars <- adata$var
vars <-  py_to_r(adata$var_names$to_list())

expression_table = counts
rownames(expression_table) <- vars

meta_data$assignments <- meta_data$broad_cluster

#remove those cell_type that have a single entry
singlets = meta_data %>%
  group_by(assignments) %>%
  filter(n() < 2) %>%
  pull(broad_cluster)

expression_table = expression_table[,!meta_data$broad_cluster %in% singlets]
meta_data = meta_data %>%
  filter(!broad_cluster %in% singlets)


###################################################
### Now with Limma-voom below
library(limma)
library(edgeR)
# Ensure valid column names for assignments and exptype
meta_data$assignments <- make.names(meta_data$assignments)



#if Metastasized_to is not nan change it to 1
meta_data$Metastasis <- ifelse(meta_data$Metastasized_to=="nan", 0, 1)

#keep where its either healthy or metastasis is 1
expression_table <- expression_table[, meta_data$Healthy == "True" | meta_data$Metastasis == 0]
meta_data <- meta_data[meta_data$Healthy == "True" | meta_data$Metastasis == 0, ]

dge = DGEList(counts = expression_table, group = meta_data$assignments)


design <- model.matrix(~ 0 + assignments + assignments:Healthy, data = meta_data)

#makenames for terms
colnames(design)<-make.names(colnames(design))

keep_genes <- filterByExpr(dge, design, min.total.count = 200,min.count=0,min.prop=0)

dge <- dge[keep_genes, ]

print(dim(dge))
#also do not keep cortex_markers
dge <- dge[!rownames(dge) %in% cortex_markers, ]
print(dim(dge))
#norm
dge <- calcNormFactors(dge, method = "TMM")

fit <- voom(dge, design, plot=TRUE)

fit <- lmFit(fit, design)

contrast <- makeContrasts(
  assignmentsChromaffin_cells.HealthyTrue,
  levels = design
)

# Compute differential expression
diff_expression <- contrasts.fit(fit, contrast)
diff_expression <- eBayes(diff_expression, robust=TRUE)

# Get top differentially expressed genes
top_genes <- topTable(diff_expression, number = Inf)
top_genes$genes <- rownames(top_genes)
print(rownames((top_genes)))



sox2_values <- dge$counts["PORCN", ] / colSums(dge$counts)
assignments <- meta_data$Healthy
#barchart
barplot(sox2_values, names.arg = assignments, las = 2, main = "SOX2 Expression by Assignment", col = "lightblue", ylab = "Expression Level")







contrast <- makeContrasts(
  assignmentsSustentacular_cells.Metastasis,
  levels = design
)


# Compute differential expression
diff_expression <- contrasts.fit(fit, contrast)
diff_expression <- eBayes(diff_expression, robust=TRUE)

# Get top differentially expressed genes
top_genes <- topTable(diff_expression, number = Inf)
top_genes$genes <- rownames(top_genes)
print(rownames((top_genes)))


#iterate through
coefs = c("assignmentsChromaffin_cells.HealthyTrue","assignmentsEndothelial_cells.HealthyTrue",  
          "assignmentsFibroblasts.HealthyTrue","assignmentsMacrophages.HealthyTrue",        
          "assignmentsSustentacular_cells.HealthyTrue","assignmentsT_cells.HealthyTrue")

final_res <- list()
final_res_sig <- list()
for (coef in coefs) {
  contrast <- makeContrasts(contrasts = coef, levels = design)
  
  
  # Compute differential expression
  diff_expression <- contrasts.fit(fit, contrast)
  diff_expression <- eBayes(diff_expression, robust=TRUE)
  
  # Get top differentially expressed genes
  top_genes <- topTable(diff_expression, number = Inf)
  top_genes$genes <- rownames(top_genes)
  top_genes$comparison <- coef
  
  print(rownames((top_genes)))
  final_res[[coef]] <- top_genes
  final_res_sig[[coef]] <- top_genes %>%
    filter(adj.P.Val < 0.05)
}

# Combine all results into a single data frame
final_res_df <- do.call(rbind, final_res)
final_res_sig_df <- do.call(rbind, final_res_sig)
# Save the final results to a CSV file
save_fold = "./"

write.csv(final_res_df, paste0(save_fold,"all_results_tumor_vs_healthy.csv"))
write.csv(final_res_sig_df, paste0(save_fold,"all_results_sig_tumor_vs_healthy.csv"))







#### Metastatic to non-metastatic
library(tidyverse)
library(reticulate)

# Specify the path to the .h5ad file
h5ad_file <- "./ppgl_pbulk_0624.h5ad"

# Read the .h5ad file using reticulate
anndata <- import("anndata")
py_index <- import("pandas")$Index
adata <- anndata$read_h5ad(h5ad_file)

counts <- t(adata$X) # Transpose to convert to column-major order
meta_data <- as.data.frame(adata$obs) # Convert obs to a data frame
vars <- adata$var
vars <-  py_to_r(adata$var_names$to_list())

expression_table = counts
rownames(expression_table) <- vars

meta_data$assignments <- meta_data$broad_cluster

#remove those cell_type that have a single entry
singlets = meta_data %>%
  group_by(assignments) %>%
  filter(n() < 2) %>%
  pull(broad_cluster)

expression_table = expression_table[,!meta_data$broad_cluster %in% singlets]
meta_data = meta_data %>%
  filter(!broad_cluster %in% singlets)


###################################################
### Now with Limma-voom below
library(limma)
library(edgeR)
# Ensure valid column names for assignments and exptype
meta_data$assignments <- make.names(meta_data$assignments)



#if Metastasized_to is not nan change it to 1
meta_data$Metastasis <- ifelse(meta_data$Metastasized_to=="nan", 0, 1)

#keep where its either healthy or metastasis is 1
expression_table <- expression_table[, meta_data$Healthy == "False"]
meta_data <- meta_data[meta_data$Healthy == "False",]

dge = DGEList(counts = expression_table, group = meta_data$assignments)


design <- model.matrix(~ 0 + assignments + assignments:Metastasis, data = meta_data)

#makenames for terms
colnames(design)<-make.names(colnames(design))

keep_genes <- filterByExpr(dge, design, min.total.count = 200,min.count=0,min.prop=0)

dge <- dge[keep_genes, ]

#norm
dge <- calcNormFactors(dge, method = "TMM")

fit <- voom(dge, design, plot=TRUE)

fit <- lmFit(fit, design)

contrast <- makeContrasts(
  assignmentsChromaffin_cells.Metastasis,
  levels = design
)

# Compute differential expression
diff_expression <- contrasts.fit(fit, contrast)
diff_expression <- eBayes(diff_expression, robust=TRUE)

# Get top differentially expressed genes
top_genes <- topTable(diff_expression, number = Inf)
top_genes$genes <- rownames(top_genes)
print(rownames((top_genes)))



write.csv(top_genes , paste0(save_fold,"chromaffin_metastatic.csv"))





