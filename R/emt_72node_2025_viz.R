########## SETUP ############

rm(list=ls())
library(sRACIPE) # from GitHub
library(ggplot2)
library(FNN)
library(ClusterR)
library(nnet)
library(reshape2)
library(dplyr)
library(tidyr)
library(prodlim)
library(circlize)
library(viridis)
library(ComplexHeatmap) # from Bioconductor
library(igraph)
library(keyplayer)
#library(forcats) not needed anymore(?)
library(ggrepel) # for PCA loadings
library(limma) # for DEG statistics
library(cowplot) # for plotting
source("R/utils.R")
source("R/utils_clamping.R")
source("R/scratch.R")

# set up directories
topoName <- "emt_ffctopo_72node_CLAMP_2025"
topoDir <- file.path(getwd(),topoName)
plotDir <- file.path(topoDir,"plots")
dataDir <- file.path(topoDir,"data")

# load topology
topo <- loadTopo(topoName)
nGenes <- length(unique(c(topo$Source, topo$Target)))
genes <- unique(c(topo$Source, topo$Target))

set.seed(1234)
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(3,2,4:8)]

# Racipe setup for later
nModels <- 2000
#nModelsFromPyracipe <- -1
racipe_placeholder <- sracipeSimulate(topo, numModels = nModels, plots = FALSE, genIC = TRUE,
                                      genParams = TRUE, integrate = FALSE, integrateStepSize = 0.01,
                                      simulationTime = 0.1, initialNoise = 0, nNoise = 0, simDet = T, anneal = F
)

racipe <- racipe_placeholder
colOrder <- rownames(sracipeIC(racipe))

# plotting utility
square_plot <- function(p) {
  # Extract the legend
  legend <- cowplot::get_legend(p)
  
  # Remove the legend from the plot
  p_no_legend <- p + theme(legend.position = "none")
  
  # Combine with fixed square plot area
  final_plot <- plot_grid( 
    p_no_legend,
    legend,
    rel_widths = c(1, 0.3),
    nrow = 1)
  
  return(final_plot)
  
}

########## UNTREATED CONDITION (FROM PYRACIPE) ############
# read data from pyracipe
pracipe_dir <- file.path(getwd(),paste0("../pyracipe/emt_ffctopo_72node_5k200ic"))
pyracipeNorm <- read.csv(file.path(pracipe_dir,
                                   paste0("emt_ffctopo_72node_5k200ic",
                                          "_pyracipe_states_2025_combined.csv")), row.names = 1)
rownames(pyracipeNorm) <- 1:nrow(pyracipeNorm)

pyracipe_summary <- read.csv(file.path(pracipe_dir,
                                       paste0("emt_ffctopo_72node_5k200ic",
                                              "_pyracipe_summary_2025_combined.csv")), row.names = 1)
rownames(pyracipe_summary) <- 1:nrow(pyracipe_summary)

pyracipe_metadata <- read.csv(file.path(pracipe_dir,
                                        paste0("emt_ffctopo_72node_5k200ic",
                                               "_pyracipe_metadata_2025_combined.csv")), row.names = 1)
rownames(pyracipe_metadata) <- 1:nrow(pyracipe_metadata)
pyracipe_metadata$Index <- 1:nrow(pyracipe_metadata)

pyracipe_params <- read.csv(file.path(pracipe_dir,
                                      paste0("emt_ffctopo_72node_5k200ic",
                                             "_pyracipe_params_2025_combined.csv")), 
                            check.names = F, row.names = 1)
rownames(pyracipe_params) <- 1:nrow(pyracipe_params)

pyracipe_raw_states <- read.csv(file.path(pracipe_dir,
                                          paste0("emt_ffctopo_72node_5k200ic",
                                                 "_pyracipe_states_unStdized_2025_combined.csv")), row.names = 1)
rownames(pyracipe_raw_states) <- 1:nrow(pyracipe_raw_states)
pyracipe_raw_states <- 2^pyracipe_raw_states
pyracipe_raw_states$Index <- 1:nrow(pyracipe_raw_states)

tmpMeans <- rowMeans(log2(t(pyracipe_raw_states)))[colOrder]
tmpSds <- apply(log2(t(1+pyracipe_raw_states)),1,sd)[colOrder]


# regenerate norm data with pseudocount to make normalization consistent
nCols <- ncol(pyracipeNorm)
pyracipeNorm <- log2(1+pyracipe_raw_states[,colOrder])
pyracipeNorm[,1:nCols] <- sweep(pyracipeNorm[,1:nCols], 2, tmpMeans, FUN = "-") # scale
pyracipeNorm[,1:nCols] <- sweep(pyracipeNorm[,1:nCols], 2, tmpSds, FUN = "/") # scale


########## MULTISTABILITY HISTOGRAM ############
## This will plot a histogram showing the number of states per model in initial simulations
pyracipe_summary_rds_fname <- file.path(dataDir,"pyracipe_summary.Rds")
pyracipe_summary <- readRDS(pyracipe_summary_rds_fname)
summary_hist_df_in <- pyracipe_summary

# Ensure summary has the expected structure
count_df <- summary_hist_df_in %>%
  count(NO_STATES, StateIdentity, name = "Count") %>%  # Count occurrences safely
  arrange(NO_STATES)
count_long <- count_df[which(count_df$NO_STATES > 0),]

# Define colors
count_long$StateIdentity <- factor(count_long$StateIdentity, levels=c("E", "Bistable", "M"))
colors <- c("E" = "blue", "Bistable" = "purple", "M" = "red")

max_bin <- 7
count_long <- count_long %>%
  mutate(NO_STATES = ifelse(NO_STATES >= max_bin, paste0(max_bin, "+"), as.character(NO_STATES)))

# Convert NumStates to factor for proper ordering
count_long$NO_STATES <- factor(count_long$NO_STATES, levels = c(as.character(1:max_bin), paste0(max_bin, "+")))

# Create stacked bar plot
image <- ggplot(count_long, aes(x = factor(NO_STATES), y = Count, fill = StateIdentity)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  labs(x = "States",
       y = "Count",
       fill = "Cluster Identity") +
  theme_sticcc() +
  theme(axis.text = element_text(size=24, inherit.blank = FALSE),
        panel.background = element_rect("white"),
        plot.background = element_rect("white"))
image

stab_hist_fname <- file.path(plotDir,"multistability_hist.pdf")
pdf(stab_hist_fname, height = 7.7, width = 10)
print(square_plot(image))
dev.off()



########## WT PCA (ALL) ############
## PCA scatterplot of unperturbed steady states

# PCA
pca_fname <- file.path(dataDir,"pca_uniques.Rds")
pca <- readRDS(pca_fname)
pca_df <- as.data.frame(pca$x)
pc1_weight <- round(100*summary(pca)$importance[2,1],2)
pc2_weight <- round(100*summary(pca)$importance[2,2],2)
plot_xlab <- paste("PC1 (",pc1_weight,"%)",sep="")
plot_ylab <- paste("PC2 (",pc2_weight,"%)",sep="")

image <- ggplot(pca_df) +
  geom_point(aes(x=PC1, y=PC2), size=3) +
  theme(panel.background = element_rect("white"),
        plot.background = element_rect("white"),
        axis.title = element_text(size=28),
        axis.text = element_text(size=24),
        panel.grid = element_line("gray"),
        axis.line = element_line("black")) +
  labs(x = plot_xlab, y = plot_ylab)
image

pca_plot_fname <- file.path(plotDir, "pca_all_states.pdf")
pdf(pca_plot_fname, width = 10, height = 10)
print(image)
dev.off()

########## PCA LOADINGS ############
## Plot PCA loadings for n top contributing genes

loadings_df <- as.data.frame(pca$rotation)
loadings_df$gene <- rownames(loadings_df)
n <- 30
top_genes <- loadings_df %>%
  dplyr::mutate(mag = abs(PC1) + abs(PC2)) %>%
  dplyr::arrange(desc(mag)) %>%
  dplyr::slice(1:30)

image <- ggplot(loadings_df, aes(x = PC1, y = PC2, label = gene)) +
  geom_point(color = "darkblue", size=4) +
  geom_hline(yintercept=0, color="black", alpha=0.5) +
  geom_vline(xintercept=0, color="black", alpha=0.5) +
  geom_text_repel(data=top_genes, size = 10, hjust = 0, vjust = 1, 
                  box.padding = 0.2,         # padding around points
                  point.padding = 0.1,       # connection between point and label
                  max.overlaps = Inf,
                  min.segment.length = 0) +  
  theme(panel.background = element_rect("white"),
        plot.background = element_rect("white"),
        axis.title = element_text(size=28),
        axis.text = element_text(size=24),
        panel.grid = element_line("gray"),
        axis.line = element_line("black")) +
  labs(x = plot_xlab, y = plot_ylab)
image

loading_plot_fname <- file.path(plotDir, paste0("PCA_loadings_top",n,".pdf"))
pdf(loading_plot_fname, width = 10, height = 10)
print(image)
dev.off()

########## CLUSTER K SELECTION ############
## Line plot of silhouette score vs number of clusters

sil_df_fname <- file.path(dataDir, "silhouette_df.Rds")
sil_df <- readRDS(sil_df_fname)

image <- ggplot(sil_df, aes(x=Clusters, y=Silhouette)) +
  geom_point(size=6) +
  geom_line(linewidth=3) +
  theme(panel.background = element_rect("white"),
        plot.background = element_rect("white"),
        axis.title = element_text(size=22),
        axis.text = element_text(size=18),
        panel.grid = element_line("gray"),
        axis.line = element_line("black")) +
  labs(x="k", y="Silhouette Score")
image  

sil_plot_fname <- file.path(plotDir, "silhouette_k_selection.pdf")
pdf(sil_plot_fname, width = 10, height = 10)
print(image)
dev.off()


########## CLUSTER ASSIGNMENT ############
## PCA scatterplot of all steady states, colored by cluster

clust_all_fname <- file.path(dataDir,"clust_all_2025.Rds") 
clust_full <- readRDS(clust_all_fname)

# PCA
image <- ggplot(pca_df) +
  geom_point(aes(x=PC1, y=PC2, color=factor(clust_full)), alpha=0.5, size=3) +
  theme(panel.background = element_rect("white"),
        plot.background = element_rect("white"),
        axis.title = element_text(size=28),
        axis.text = element_text(size=24),
        panel.grid = element_line("gray"),
        axis.line = element_line("black")) +
  labs(x = plot_xlab, y = plot_ylab, color="Cluster") +
  scale_color_manual(values = cbPalette)
image

cluster_plot_fname <- file.path(plotDir, "pca_all_cluster.pdf")
pdf(cluster_plot_fname, width = 10, height = 10)
print(square_plot(image))
dev.off()


########## BISTABLE PCA ############
## Here we will plot all states of E/M bistable models (colored by cluster),
## as well as all initial conditions for clamping simulations

# First, select E/M bistable models
pyracipe_metadata_rds_fname <- file.path(dataDir,"pyracipe_metadata.Rds")
pyracipe_summary_rds_fname <- file.path(dataDir,"pyracipe_summary.Rds")
pyracipe_metadata <- readRDS(pyracipe_metadata_rds_fname)
pyracipe_summary <- readRDS(pyracipe_summary_rds_fname)

keep_models <- pyracipe_summary[which(pyracipe_summary$NO_STATES == 2 &
                                        pyracipe_summary$StateIdentity == "Bistable"),"MODEL_NO"]
keep_states <- c()
keep_states_e <- c()
for(model in keep_models) {
  model_states_all <- pyracipe_metadata[which(pyracipe_metadata$MODEL_NO == model), "Index"]
  model_state_e <- pyracipe_metadata[which(pyracipe_metadata$MODEL_NO == model &
                                             pyracipe_metadata$Cluster == 1), "Index"]
  keep_states <- c(keep_states, model_states_all)
  keep_states_e <- c(keep_states_e, model_state_e)
}

# first, plot all E/M bistable states
image <- ggplot(pca_df[keep_states,], aes(x=PC1, y=PC2, color=as.factor(clust_full[keep_states]))) +
  geom_point(size=3) +
  guides(color=guide_legend(title = "Cluster")) +
  scale_color_manual(values=cbPalette) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  xlim(-6, 12.5) +
  ylim(-6, 7) + 
  theme_sticcc() +
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color="black"),
        panel.background = element_rect("white"),
        plot.background = element_rect("white"))
image

pca_plot_fname <- file.path(plotDir,"pca_wt_states_bistable.pdf")
pdf(pca_plot_fname, height = 10, width = 10)
print(image)
dev.off()


# plot initial conditions for perturbation sims: E states from bistable models
image <- ggplot(pca_df[keep_states_e,], aes(x=PC1, y=PC2, color=as.factor(clust_full[keep_states_e]))) +
  geom_point(size=3) +
  guides(color=guide_legend(title = "Cluster")) +
  scale_color_manual(values=cbPalette) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  xlim(-6, 12.5) +
  ylim(-6, 7) + 
  theme_sticcc() +
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color="black"),
        panel.background = element_rect("white"),
        plot.background = element_rect("white"))
image

pca_plot_fname <- file.path(plotDir,"pca_wt_states_preSignaling.pdf")
pdf(pca_plot_fname, height = 10, width = 10)
print(image)
dev.off()


## finally, plot all states, and overlay the selected E/M bistable states 
image <- ggplot(pca_df[,], aes(x=PC1, y=PC2, color=as.factor(clust_full))) +
  geom_point(size=3, alpha=0.7) +
  geom_point(data=pca_df[keep_states,], 
             aes(x=PC1, y=PC2, fill=as.factor(clust_full[keep_states])), 
             size=3, color="black", pch=21) +
  guides(color=guide_legend(title = "Cluster"), fill="none") +
  scale_color_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  xlim(-6, 12.5) +
  ylim(-6, 7) + 
  theme_sticcc() +
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color="black"),
        panel.background = element_rect("white"),
        plot.background = element_rect("white")
        )
image

pca_plot_fname <- file.path(plotDir,"pca_wt_states_all_bistableOverlay.pdf")
pdf(pca_plot_fname, height = 10, width = 10)
print(image)
dev.off()



########## CLAMP VALUES ############
# For each gene, plot the distribution of steady state values for E and M clusters respectively
clamp_df_fname <- file.path(dataDir,"clamp_values_2025.Rds")
clamp_df <- readRDS(clamp_df_fname)


# Step 1: Calculate mean expression per gene per cluster
gene_order <- clamp_df %>%
  group_by(Gene, Cluster) %>%
  summarize(mean_expr = mean(Expression), .groups = "drop") %>%
  group_by(Gene) %>%
  summarize(diff = max(mean_expr) - min(mean_expr), .groups = "drop") %>%
  arrange(desc(diff)) %>%
  pull(Gene)

# Step 2: Set Gene as a factor with desired order
clamp_df$Gene <- factor(clamp_df$Gene, levels = gene_order)

# Step 3: Plot
image <- ggplot(data = clamp_df, aes(x = Gene, y = Expression, fill = as.factor(Cluster))) + 
  geom_boxplot() +
  labs(title = "Gene Expression by Cluster", 
       x = "Gene", 
       y = "Expression",
       fill = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
image


clamp_boxplot_fname <- file.path(plotDir,"clamp_val_boxplot.pdf")
pdf(clamp_boxplot_fname, height = 10, width = 10)
print(image)
dev.off()

clamp_boxplot_log_fname <- file.path(plotDir,"clamp_val_boxplot_log.pdf")
pdf(clamp_boxplot_log_fname, height = 10, width = 10)
print(image + scale_y_log10())
dev.off()


## Also plot the same for all states, not just E/M bistable ones
clamp_df_full_fname <- file.path(dataDir,"clamp_values_all_2025.Rds")
clamp_df_full <- readRDS(clamp_df_full_fname)


# Step 1: Calculate mean expression per gene per cluster
gene_order_full <- clamp_df_full %>%
  group_by(Gene, Cluster) %>%
  summarize(mean_expr = mean(Expression), .groups = "drop") %>%
  group_by(Gene) %>%
  summarize(diff = max(mean_expr) - min(mean_expr), .groups = "drop") %>%
  arrange(desc(diff)) %>%
  pull(Gene)

# Step 2: Set Gene as a factor with desired order
clamp_df_full$Gene <- factor(clamp_df_full$Gene, levels = gene_order_full)

# Step 3: Plot
image <- ggplot(data = clamp_df_full, aes(x = Gene, y = Expression, fill = as.factor(Cluster))) + 
  geom_boxplot() +
  labs(title = "Gene Expression by Cluster", 
       x = "Gene", 
       y = "Expression",
       fill = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
image


clamp_boxplot_fname <- file.path(plotDir,"clamp_val_ALL_boxplot.pdf")
pdf(clamp_boxplot_fname, height = 10, width = 10)
print(image)
dev.off()

clamp_boxplot_log_fname <- file.path(plotDir,"clamp_val_ALL_boxplot_log.pdf")
pdf(clamp_boxplot_log_fname, height = 10, width = 10)
print(image + scale_y_log10())
dev.off()


########## COMPARE ALL STATES VS E/M BISTABLE ############
## We notice that DEGs between E and M are different for all states vs just E/M bistable states
## Here we plot the difference in DEG rankings 

## Assign E and M genes
assign_cluster_DEGs <- function(expr, cluster_labels, p_thresh = 0.05) {
  design <- model.matrix(~ cluster_labels)
  fit <- lmFit(expr, design, ref=1)
  fit <- eBayes(fit)
  
  top_table <- topTable(fit, coef = 2, number = Inf)
  top_table$cluster <- "Unassigned"
  top_table[which(top_table$logFC > 0),"cluster"] <- "M"
  top_table[which(top_table$logFC < 0),"cluster"] <- "E"
  return(top_table)
}

cluster_degs_all <- assign_cluster_DEGs(t(pyracipeNorm), factor(clust_full))
cluster_degs_bistable <- assign_cluster_DEGs(t(pyracipeNorm[keep_states,]), factor(clust_full[keep_states]))


# Assign ranks (lower = higher rank)
df <- data.frame(
  gene = genes,
  rank1 = match(genes, rownames(cluster_degs_all)),
  rank2 = match(genes, rownames(cluster_degs_bistable))
)

# Replace NA with lowest possible rank + 1 (i.e., items not in one of the lists)
max_rank <- max(c(df$rank1, df$rank2), na.rm = TRUE)
df$rank1[is.na(df$rank1)] <- max_rank + 1
df$rank2[is.na(df$rank2)] <- max_rank + 1

# Calculate change in rank (positive = moved down, negative = moved up)
df$delta <- df$rank2 - df$rank1

# Sort by magnitude of change for display
df <- df[order(abs(df$delta), decreasing = TRUE), ]

# Optional: keep top N most-changed genes
top_n <- 30
df_plot <- head(df, top_n)

# Plot
image <- ggplot(df_plot, aes(x = reorder(gene, delta), y = delta)) +
  geom_segment(aes(xend = gene, y = 0, yend = delta), color = "gray") +
  geom_point(aes(color = delta > 0), size = 3) +
  scale_y_reverse() +
  scale_color_manual(values = c("forestgreen", "firebrick"),
                     labels = c("Moved Up", "Moved Down"),
                     name = "Direction") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Gene", y = "Change in Rank (Bistable - All)",
       title = "Gene Rank Changes Between 'All' and 'E/M Bistable'") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_rect("white"),
        panel.background = element_rect("white"))
image

em_gene_ranking_fname <- file.path(plotDir, "gene_rank_diffs.pdf")
pdf(em_gene_ranking_fname, height = 10, width = 10)
print(image)
dev.off()

########## GENE EXPRESSION HEATMAP ############
## Here we will plot heatmaps showing gene expression by cluster for 1) all steady states, and 2) E/M bistable model states
ha_df <- data.frame(Cluster=clust_full)

# Create an annotation object for the columns
column_annotation <- HeatmapAnnotation(df = ha_df, 
                                       col=list(Cluster=c("1"=unname(cbPalette[1]),"2"=unname(cbPalette[2]),
                                                          "3"=unname(cbPalette[3]),"4"=unname(cbPalette[4]))))
# Create the heatmap with annotation
top_genes <- rownames(cluster_degs_all)[1:30]
image <- Heatmap(as.matrix(t(pyracipeNorm[,top_genes])), 
                 name = "Expression", 
                 top_annotation = column_annotation,
                 row_names_gp=gpar(fontsize=12),
                 show_column_names = F)
image

wt_hmap_fname <- file.path(plotDir,"wt_hmap_allStates_top30Genes.pdf")
pdf(wt_hmap_fname, height = 10, width = 10)
print(image)
dev.off()

wt_hmap_fname_jpg <- file.path(plotDir,"wt_hmap_allStates_top30Genes.jpg")
jpeg(wt_hmap_fname_jpg, height = 6, width = 6, units="in", res = 300)
print(image)
dev.off()

## Now for just bistable states
ha_df <- data.frame(Cluster=clust_full[keep_states])

# Create an annotation object for the columns
column_annotation <- HeatmapAnnotation(df = ha_df, 
                                       col=list(Cluster=c("1"=unname(cbPalette[1]),"2"=unname(cbPalette[2]),
                                                          "3"=unname(cbPalette[3]),"4"=unname(cbPalette[4]))))
# Create the heatmap with annotation
top_genes <- rownames(cluster_degs_bistable)[1:30]
image <- Heatmap(as.matrix(t(pyracipeNorm[keep_states,top_genes])), 
                 name = "Expression", 
                 top_annotation = column_annotation,
                 row_names_gp=gpar(fontsize=12),
                 show_column_names = F)
image

wt_hmap_fname <- file.path(plotDir,"wt_hmap_bistableStates_top30Genes.pdf")
pdf(wt_hmap_fname, height = 10, width = 10)
print(image)
dev.off()

wt_hmap_fname_jpg <- file.path(plotDir,"wt_hmap_bistableStates_top30Genes.jpg")
jpeg(wt_hmap_fname_jpg, height = 6, width = 6, units="in", res = 300)
print(image)
dev.off()


########## MODEL-WISE ASYMMETRY ############

## Multiple ICs - plot basin boundaries
asym_model_use <- 3 # Using the first E/M bistable model in the list
asym_model_params <- pyracipe_params[asym_model_use,]
asym_model_states <- pyracipe_raw_states[which(pyracipe_metadata$MODEL_NO == asym_model_use),genes]
asym_multIC_fname <- file.path(dataDir, paste0("multIC_asym_model_",asym_model_use,".Rds"))
asym_multIC <- readRDS(asym_multIC_fname)

# transform expected states to PCA space
asym_model_states_norm <- asym_model_states[,genes]
asym_model_states_norm[,genes] <- log2(1+asym_model_states_norm[,genes]) # Log transform
asym_model_states_norm[,genes] <- sweep(asym_model_states_norm[,genes], 2, tmpMeans, FUN = "-") # scale
asym_model_states_norm[,genes] <- sweep(asym_model_states_norm[,genes], 2, tmpSds, FUN = "/") # scale
asym_model_states_pca <- scale(asym_model_states_norm, pca$center, pca$scale) %*% pca$rotation

# identify final steady state clusters
asym_multIC_ss_norm <- as.data.frame(t(assay(asym_multIC)))[,genes]
asym_multIC_ss_norm[,genes] <- log2(1+asym_multIC_ss_norm[,genes]) # Log transform
asym_multIC_ss_norm[,genes] <- sweep(asym_multIC_ss_norm[,genes], 2, tmpMeans, FUN = "-") # scale
asym_multIC_ss_norm[,genes] <- sweep(asym_multIC_ss_norm[,genes], 2, tmpSds, FUN = "/") # scale
asym_multIC_ss_pca <- scale(asym_multIC_ss_norm, pca$center, pca$scale) %*% pca$rotation
# replace w/ GMM clustering?
asym_multIC_ss_clusts <- knn_classifier(asym_multIC_ss_pca[,1:num_pcs_clustering], pca_df[,1:num_pcs_clustering], clust_full, k=25)

# pca-transform ICs
asym_multIC_ic_norm <- as.data.frame(t(sracipeIC(asym_multIC)))[,genes]
asym_multIC_ic_norm[,genes] <- log2(1+asym_multIC_ic_norm[,genes]) # Log transform
asym_multIC_ic_norm[,genes] <- sweep(asym_multIC_ic_norm[,genes], 2, tmpMeans, FUN = "-") # scale
asym_multIC_ic_norm[,genes] <- sweep(asym_multIC_ic_norm[,genes], 2, tmpSds, FUN = "/") # scale
asym_multIC_ic_pca <- scale(asym_multIC_ic_norm, pca$center, pca$scale) %*% pca$rotation

# plot ICs and color by what state they ended up in
image <- ggplot() +
  geom_point(data=asym_multIC_ic_pca, aes(x=PC1, y=PC2, color=asym_multIC_ss_clusts)) +
  geom_point(data=asym_model_states_pca, aes(x=PC1, y=PC2), color="red", size=3) +
  theme_sticcc() +
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color="black"),
        axis.title = element_text(size=22),
        plot.background = element_rect("white"),
        panel.background = element_rect("white")) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  labs(color="Cluster")
image

asym_multIC_hmap_fname <- file.path(plotDir,paste0("asym_multIC_PCA_model=",asym_model_use,".pdf"))
pdf(asym_multIC_hmap_fname, height = 10, width = 10)
print(image)
dev.off()



## Long stochastic simulations
## Plot trajectory as heatmap in PCA space
asym_noise_levels <- c(0.02, 0.04, 0.08, 0.1, 0.2)
for(asym_noise_level in asym_noise_levels) {
  asym_longStoch_fname <- file.path(dataDir, 
                                    paste0("longStoch_asym_model_",asym_model_use,
                                           "_noise=",asym_noise_level,".Rds"))
  asym_longStoch <- readRDS(asym_longStoch_fname)
  
  asym_longStoch_raw <- as.data.frame(t(asym_longStoch@metadata$timeSeries))
  for(gene in genes) {
    asym_longStoch_raw[,gene] <- as.numeric(asym_longStoch_raw[,gene])
  }
  
  
  asym_longStoch_norm <- asym_longStoch_raw
  asym_longStoch_norm[,genes] <- log2(1+asym_longStoch_norm[,genes]) # Log transform
  asym_longStoch_norm[,genes] <- sweep(asym_longStoch_norm[,genes], 2, tmpMeans, FUN = "-") # scale
  asym_longStoch_norm[,genes] <- sweep(asym_longStoch_norm[,genes], 2, tmpSds, FUN = "/") # scale
  asym_longStoch_pca <- scale(asym_longStoch_norm, pca$center, pca$scale) %*% pca$rotation
  
  
  image <- ggplot() +
    geom_point(data=asym_longStoch_pca, aes(x=PC1, y=PC2)) +
    geom_point(data=asym_model_states_pca, aes(x=PC1, y=PC2), color="red", size=3) +
    theme_sticcc() +
    theme(axis.line = element_line(linewidth = 1, color = "black"), 
          axis.ticks = element_line(linewidth = 1, color="black"),
          axis.title = element_text(size=22),
          plot.background = element_rect("white"),
          panel.background = element_rect("white")) +
    xlab(plot_xlab) +
    ylab(plot_ylab) +
    labs(color="Cluster")
  print(image)
  
  asym_longStoch_plot_fname <- file.path(plotDir,paste0("longStoch_asym_model_",asym_model_use,
                                                        "_noise=",asym_noise_level,".pdf"))
  pdf(file = asym_longStoch_plot_fname, width = 10, height = 10)
  print(image)
  dev.off()
  
  
  
  image <- ggplot(data=asym_longStoch_pca, aes(x=PC1, y=PC2)) +
    #geom_point(data=asym_longStoch_pca, aes(x=PC1, y=PC2)) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    scale_fill_viridis_c() +
    geom_point(data=asym_model_states_pca, aes(x=PC1, y=PC2), color="red", size=3) +
    theme_sticcc() +
    theme(axis.line = element_line(linewidth = 1, color = "black"), 
          axis.ticks = element_line(linewidth = 1, color="black"),
          axis.title = element_text(size=22),
          plot.background = element_rect("white"),
          panel.background = element_rect("white")) +
    xlab(plot_xlab) +
    ylab(plot_ylab) +
    labs(color="Cluster")
  print(image)
  
  asym_longStoch_plot_fname <- file.path(plotDir,paste0("longStoch_asym_model_",asym_model_use,
                                                        "_noise=",asym_noise_level,"_density.pdf"))
  pdf(file = asym_longStoch_plot_fname, width = 10, height = 10)
  print(image)
  dev.off()
    
    
    
  
  
}

########## SIGNAL EFFICACY HEATMAP ############
# Plot 2D heatmap at different noise levels
# Note: only works with 2-gene signals included
#resultSet_full <- readRDS()


signal_simTime <- 500
signal_relaxTime <- 50
signal_nGenes <- c(1,2)
signal_noise <- c(0, 0.02, 0.04) 
signal_tcorr <- 10
initClust <- 1
tgtClust <- 2
expName <- paste0("ffctopo_t=",signal_simTime,"_relax_OUnoise=",paste0(signal_noise, collapse = "."),
                  "_tau=",signal_tcorr,"_genes=",paste0(signal_nGenes,collapse = "."),"_CLAMPS_2025")
resultSet_fname <- file.path(topoDir,expName,"result_summary.Rds")
genes_x_transitions_matrix_fname <- file.path(topoDir, expName, "genes_x_transitions_matrix.Rds")
resultSet_full <- readRDS(resultSet_fname)
resultSet <- resultSet_full

selectedNoise <- 0.04
eff_hmap_0.04 <- resultSet_full[which(resultSet_full$Noise == selectedNoise),
                                c("Species 1", "Species 2", "ConversionPct")]

eff_hmap_0.04_matrix <- matrix(nrow = 72, ncol = 72)
rownames(eff_hmap_0.04_matrix) <- genes
colnames(eff_hmap_0.04_matrix) <- genes
for(gene1 in genes) {
  for (gene2 in genes) {
    if(gene1 == gene2) {
      eff_val <- eff_hmap_0.04[which(eff_hmap_0.04$`Species 1` == gene1 & 
                                       is.na(eff_hmap_0.04$`Species 2`)), "ConversionPct"]
    } else if(length(which(eff_hmap_0.04$`Species 1` == gene1 & 
                           eff_hmap_0.04$`Species 2` == gene2)) == 1) {
      eff_val <- eff_hmap_0.04[which(eff_hmap_0.04$`Species 1` == gene1 & 
                                       eff_hmap_0.04$`Species 2` == gene2), "ConversionPct"]
    } else {
      eff_val <- eff_hmap_0.04[which(eff_hmap_0.04$`Species 2` == gene1 & 
                                       eff_hmap_0.04$`Species 1` == gene2), "ConversionPct"]
    }
    if(length(eff_val) == 0) {
      eff_val = NA
    }
    eff_hmap_0.04_matrix[gene1, gene2] <- eff_val
    eff_hmap_0.04_matrix[gene2, gene1] <- eff_val
  }
}


image(eff_hmap_0.04_matrix)
eff_hmap_0.04_matrix[lower.tri(eff_hmap_0.04_matrix)] <- NA


# Define the color mapping using viridis
col_fun <- colorRamp2(c(min(eff_hmap_0.04_matrix, na.rm = TRUE), 
                        mean(eff_hmap_0.04_matrix, na.rm = TRUE), 
                        max(eff_hmap_0.04_matrix, na.rm = TRUE)), 
                      viridis(3))  # Generates 3 colors from viridis

# Generate the heatmap
image <- Heatmap(eff_hmap_0.04_matrix, 
                 name = "% EMT", 
                 col = col_fun, 
                 cluster_rows = FALSE, 
                 cluster_columns = FALSE, 
                 show_row_names = TRUE, 
                 show_column_names = TRUE,
                 row_names_gp = gpar(fontsize=20),
                 column_names_gp = gpar(fontsize=20),
                 column_names_rot = 55,
                 column_names_side = "top",
                 row_title = "Gene 1",
                 row_title_gp = gpar(fontsize=24),
                 column_title = "Gene 2",
                 column_title_gp = gpar(fontsize=24),
                 column_title_side = "bottom")
image

plot_fname <- file.path(plotDir,paste0("figxx_2gene_sigEff_heatmap_noise=",paste0(selectedNoise,collapse = ","),".pdf"))
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()


########## ESTIMATE CONTROL EFFICACY ############
## Since ~15 genes have almost zero difference in expression across E/M bistable states,
## we will use the data from noisy simulations targeting these genes as "control" cases
## this will illuminate the efficacy of noise alone in driving EMT

static_genes <- rownames(cluster_degs_bistable)[which(abs(cluster_degs_bistable$logFC) < 1e-10)]
noise_use <- 0.04
summary(resultSet_full[which(resultSet_full$NumGenes == 1 &
                               resultSet_full$`Species 1` %in% static_genes & 
                               resultSet_full$Noise == noise_use),"ConversionPct"])
summary(resultSet_full[which(resultSet_full$NumGenes == 1 &
                               resultSet_full$`Species 1` %in% static_genes & 
                               resultSet_full$Noise == 0),"ConversionPct"])

########## SIGNAL EFFICACY VS TOPOLOGY ############
## Plot signal efficacy against a network topology metrics of signals: closeness & betweenness centrality, out-degree 
cor(resultSet_full$GroupBetweenCentrality, resultSet_full$ConversionPct)
cor(resultSet_full$GroupBetweenCentrality[which(resultSet_full$Noise == 0)], 
    resultSet_full$ConversionPct[which(resultSet_full$Noise == 0)])
cor(resultSet_full$GroupBetweenCentrality[which(resultSet_full$Noise == 0.04)], 
    resultSet_full$ConversionPct[which(resultSet_full$Noise == 0.04)])

selectedNoise <- 0.04
image <- ggplot(resultSet[which(resultSet$Noise == selectedNoise),], aes(x=GroupBetweenCentrality, y=ConversionPct*100)) +
  geom_point(size=3) +
  theme_sticcc() +
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color="black"),
        axis.title = element_text(size=22),
        plot.background = element_rect("white"),
        panel.background = element_rect("white")) +
  xlab("Signal Group Betweenness Centrality") +
  ylab("Models Undergoing EMT (%)")
image

plot_fname <- file.path(plotDir,paste0("eff_vs_betweenness_noise=",
                                       paste0(selectedNoise, collapse = ","),".pdf"))
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()


image <- ggplot(resultSet[which(resultSet$Noise == selectedNoise),], aes(x=GroupClosenessCentrality, y=ConversionPct)) +
  geom_point(size=3) +
  theme_sticcc() +
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color="black"),
        axis.title = element_text(size=22),
        plot.background = element_rect("white"),
        panel.background = element_rect("white")) +
  xlab("Signal Group Closeness Centrality") +
  ylab("Models Undergoing EMT (%)")
image


plot_fname <- file.path(plotDir,paste0("eff_vs_closeness_noise=",
                                       paste0(selectedNoise, collapse = ","),".pdf"))
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()


image <- ggplot(resultSet[which(resultSet$Noise == selectedNoise),], aes(x=TotalOutDegree, y=ConversionPct)) +
  geom_point(size=3) +
  theme_sticcc() +
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color="black"),
        plot.background = element_rect("white"),
        panel.background = element_rect("white")) +
  xlab("Signal Total Out-Degree") +
  ylab("Models Undergoing EMT (%)")
image

plot_fname <- file.path(plotDir,paste0("eff_vs_outDegree_",
                                       paste0(selectedNoise, collapse = ","),".pdf"))
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()




########## EFFECT OF NOISE BY SIGNAL ############
if(all(selectedNoise == c(0, 0.04, 0.2))) {
  signal_summary_df <- data.frame(Signal=unique(paste0(resultSet_full[,c("Species 1")], "_",
                                                       resultSet_full[,c("Species 2")])),
                                  ConversionPct_D0=NA,
                                  ConvertingModels_D0=NA,
                                  ConversionPct_D0.04=NA,
                                  ConvertingModels_D0.04=NA,
                                  ConversionPct_D0.2=NA,
                                  ConvertingModels_D0.2=NA)
  noises <- c(0, 0.04, 0.2)
  for(sig in signal_summary_df$Signal) {
    signal_summary_df[which(signal_summary_df$Signal==sig),"GroupClosenessCentrality"] <- 
      unique(resultSet_full[which(resultSet_full$SigName_Short == sig),"GroupClosenessCentrality"])
    signal_summary_df[which(signal_summary_df$Signal==sig),"GroupBetweenCentrality"] <- 
      unique(resultSet_full[which(resultSet_full$SigName_Short == sig),"GroupBetweenCentrality"])
    for(noise in noises) {
      signal_summary_df[which(signal_summary_df$Signal==sig),paste0("ConversionPct_D",noise)] <- 
        resultSet_full[which(resultSet_full$Noise == noise & 
                               resultSet_full$SigName_Short == sig),"ConversionPct"]
      signal_summary_df[which(signal_summary_df$Signal==sig),paste0("ConvertingModels_D",noise)] <- 
        resultSet_full[which(resultSet_full$Noise == noise & 
                               resultSet_full$SigName_Short == sig),"ConvertingModels"]
    }
  }
  signal_summary_df$D0.04_Gain_Pct <- signal_summary_df$ConversionPct_D0.04 - signal_summary_df$ConversionPct_D0
  signal_summary_df$D0.04_Gain_Pct_ofRemaining <- (signal_summary_df$ConversionPct_D0.04 - signal_summary_df$ConversionPct_D0) / (1-signal_summary_df$ConversionPct_D0)
  
  ggplot(signal_summary_df) +
    geom_histogram(aes(x=ConvertingModels_D0.04-ConvertingModels_D0))
  
  ggplot(signal_summary_df) +
    geom_point(aes(x=GroupClosenessCentrality, y=D0.04_Gain_Pct))
  
  ggplot(signal_summary_df) +
    geom_point(aes(x=GroupBetweenCentrality, y=D0.04_Gain_Pct))
  
  ggplot(signal_summary_df) +
    geom_point(aes(x=GroupBetweenCentrality, y=D0.04_Gain_Pct_ofRemaining)) +
    theme_sticcc() +
    xlab("Betweenness Centrality") +
    ylab("Signal Efficacy, 0.04 vs 0")
  
  # ggplot(signal_summary_df) +
  #   geom_point(aes(x=OutDegree, y=D0.04_Gain_Pct_ofRemaining)) +
  #   theme_sticcc() +
  #   xlab("Out-Degree") +
  #   ylab("Signal Efficacy, 0.04 vs 0")
  
  
  # Heatmap showing gains from noise by signal
  
  eff_gains_0.04_matrix <- matrix(nrow = 26, ncol = 26)
  rownames(eff_gains_0.04_matrix) <- genes_reordered
  colnames(eff_gains_0.04_matrix) <- genes_reordered
  for(gene1 in genes_reordered) {
    for (gene2 in genes_reordered) {
      comb1 <- paste0(gene1, "_", gene2)
      comb2 <- paste0(gene2, "_", gene1)
      if(gene1 == gene2) {
        eff_val <- signal_summary_df[which(signal_summary_df$Signal == paste0(gene1,"_NA")), "D0.04_Gain_Pct_ofRemaining"]
      } else if(length(which(signal_summary_df$Signal == comb1)) == 1) {
        eff_val <- signal_summary_df[which(signal_summary_df$Signal == comb1), "D0.04_Gain_Pct_ofRemaining"]
      } else {
        eff_val <- signal_summary_df[which(signal_summary_df$Signal == comb2), "D0.04_Gain_Pct_ofRemaining"]
      }
      eff_gains_0.04_matrix[gene1, gene2] <- eff_val
      eff_gains_0.04_matrix[gene2, gene1] <- eff_val
    }
  }
  
  
  
  eff_gains_0.04_matrix[lower.tri(eff_gains_0.04_matrix)] <- NA
  
  # Define the color mapping
  library(circlize)
  library(viridis)
  # Define the color mapping using viridis
  col_fun <- colorRamp2(c(min(eff_gains_0.04_matrix, na.rm = TRUE), 
                          mean(eff_gains_0.04_matrix, na.rm = TRUE), 
                          max(eff_gains_0.04_matrix, na.rm = TRUE)), 
                        viridis(3))  # Generates 3 colors from viridis
  
  
  # Generate the heatmap
  image <- Heatmap(eff_gains_0.04_matrix, 
                   name = "Eff. Increase", 
                   col = col_fun, 
                   cluster_rows = FALSE, 
                   cluster_columns = FALSE, 
                   show_row_names = TRUE, 
                   show_column_names = TRUE,
                   row_names_gp = gpar(fontsize=20),
                   column_names_gp = gpar(fontsize=20),
                   column_names_rot = 55,
                   column_names_side = "top",
                   row_title = "Gene 1",
                   row_title_gp = gpar(fontsize=24),
                   column_title = "Gene 2",
                   column_title_gp = gpar(fontsize=24),
                   column_title_side = "bottom")
  image
  
  plot_fname <- file.path(plotDir,paste0("figS1c_2gene_sigEffFromNoise_heatmap_noise=",paste0(selectedNoise,collapse = ","),".pdf"))
  pdf(plot_fname, height = 10, width = 10)
  print(image)
  dev.off()
  
}



########## CELL-WISE SIGNAL EFFICACY HEATMAP ############
cell_signal_df_noise <- c(0.04) # which noise levels to consider when plotting top signals
numToPlot <- 30 # how many of the top signals to plot
cell_signal_df_fname <- file.path(topoDir,expName,
                                  paste0("cell_signal_df_noise=", 
                                         paste0(signal_noise, collapse = ","),".Rds"))
cell_signal_df <- readRDS(cell_signal_df_fname)
cell_signal_df <- cell_signal_df[,which(resultSet$Noise == cell_signal_df_noise)]


numeric_df <- matrix(as.numeric(factor(as.matrix(cell_signal_df), 
                                       levels = c("Rebellious", "Target->Target", "Init->Init", "Init->Target"))),
                     nrow = nrow(cell_signal_df), ncol = ncol(cell_signal_df))
colnames(numeric_df) <- colnames(cell_signal_df)


# select signals to plot
filterIdx <- order(resultSet$ConversionPct[which(resultSet$Noise == cell_signal_df_noise)], decreasing = T)[1:numToPlot]
filter <- resultSet[which(resultSet$Noise == cell_signal_df_noise)[filterIdx],"SetName"]
# subset for those with available data
filter <- filter[which(filter %in% colnames(numeric_df))]
filterIdx <- filterIdx[which(filter %in% colnames(numeric_df))]
numeric_df <- numeric_df[,filter]


# Convert the categories to factors and then back to numerics to ensure unique coding
numeric_df_factor <- as.numeric(as.factor(numeric_df))

# Create a matrix from the factorized data
matrix_data <- matrix(numeric_df_factor, nrow = nrow(numeric_df), ncol = ncol(numeric_df))
colnames(matrix_data) <- gsub("_noise=0.04", "", colnames(cell_signal_df)[filterIdx])

# Define color palette 
#color_palette <- c("red", "lightblue", "pink", "darkblue") 
color_palette <- c("pink", "darkblue") 





# Perform hierarchical clustering on the columns (Conditions)
model_clusters <- as.factor(cutree(hclust(dist(matrix_data)), k=4))
ha_df <- data.frame(Cluster=model_clusters)

# Create an annotation object for the columns
# Color palette
cbPalette <- palette.colors(palette = "Okabe-Ito")[2:9]
names(cbPalette) <- c(1:8)
column_annotation <- HeatmapAnnotation(df = ha_df, 
                                       col=list(Cluster=c("1"=unname(cbPalette[1]),"2"=unname(cbPalette[2]),
                                                          "3"=unname(cbPalette[3]),"4"=unname(cbPalette[4]))))

# Define a color palette for the heatmap (same as before)
#color_palette <- c("red", "lightblue", "pink", "darkblue") # adjust as needed
color_palette <- c("pink", "darkblue")

# Row annotation
colnames(resultSet)
ra_df <- resultSet[filterIdx,c("TotalOutDegree", "GroupBetweenCentrality")]
colnames(ra_df) <- c("OutDeg","Centrality")
#ra_df$Noise <- factor(ra_df$Noise)
row_annotation <- rowAnnotation(df = ra_df)

# Create the heatmap with annotation
image <- Heatmap(t(matrix_data), 
                 name = "EMT Result", 
                 col = colorRamp2(c(min(matrix_data):max(matrix_data)), color_palette), # using categorical colors
                 #top_annotation = column_annotation,
                 right_annotation = row_annotation,
                 row_names_gp=gpar(fontsize=18),
                 show_row_names = T,
                 show_column_names = F,
                 clustering_method_columns = "ward.D2",
                 clustering_method_rows = "ward.D2")


hmap_fname <- file.path(plotDir,paste0("cellwise_efficacy_heatmap_noise=",
                                       paste0(selectedNoise, collapse=","),"_top30.pdf"))
pdf(hmap_fname, height = 12, width = 12)
print(image)
dev.off()

dev.off()
lgd = Legend(labels = c("MET","M \U2192 M", "E \U2192 E", "EMT"), 
             title = "EMT Result", 
             legend_gp = gpar(fill = color_palette),
             labels_gp = gpar(fontsize = 16),
             title_gp = gpar(fontsize = 16, fontface="bold"),
             size=unit(3,"cm"))
cairo_pdf(file.path(plotDir,"hmap_legend_manual.pdf"), height=6, width=6, family="Arial")
image <- draw(lgd)
dev.off()


# # Plot PCA colored by model clusters
# image <- ggplot(pca_df, aes(x=PC1,y=PC2, color=model_clusters)) +
#   scale_color_manual(values=cbPalette[1:4]) +
#   geom_point(size=3) +
#   guides(color=guide_legend(title = "Model Cluster")) + 
#   theme_sticcc() +
#   theme(axis.line = element_line(linewidth = 1, color = "black"), 
#         axis.ticks = element_line(linewidth = 1, color="black"))
# 
# plot_fname <- file.path(plotDir,
#                         paste0("pca_by_transitionLikelihood_noise=", 
#                                paste0(selectedNoise, collapse = ","),".pdf"))
# pdf(plot_fname, height = 10, width = 10)
# print(image)
# dev.off()




