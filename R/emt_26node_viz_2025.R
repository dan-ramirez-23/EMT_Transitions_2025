########## SETUP ############

rm(list=ls())
library(sRACIPE)        # GRN simulations
library(FNN)            # nearest neighbors
library(ClusterR)
library(nnet)
library(reshape2)       # data management
library(dplyr)          # data management
library(tidyr)          # data management
library(prodlim)        # data management
library(igraph)         # network theory tools
library(limma)          # DEG analysis of simulated expression
library(keyplayer)
library(ggplot2)        # plotting
library(cowplot)        # plotting (extract legend separately)
library(ComplexHeatmap) # plotting
library(circlize)       # plotting
library(viridis)        # plotting
library(tibble)         # plotting

source("R/utils.R")
source("R/utils_clamping.R")
source("R/scratch.R")

# set up directories
topoName <- "emt_bhtopo_26node_CLAMP"
topoDir <- file.path(getwd(),topoName)
plotDir <- file.path(topoDir,"plots_apr2025")
dataDir <- file.path(topoDir,"data")

if(!dir.exists(topoDir)) {
  dir.create(topoDir)
}
if(!dir.exists(dataDir)) {
  dir.create(dataDir)
}
if(!dir.exists(plotDir)) {
  dir.create(plotDir)
}


# load topology
topo <- loadTopo(topoName)
nGenes <- length(unique(c(topo$Source, topo$Target)))
genes_reordered <- c("Foxc2","Zeb1","Klf8","Cdh1","miR101", "Zeb2", "Snai1", "miR141",
                     "Tgfbeta","miR200a","miR200b","miR200c","miR205","miR30c","Snai2",
                     "miR34a","Twist2","miR9","Vim","Twist1","Tcf3","Gsc", "Ovol2", "Grhl2",  "Np63a", "Cldn7")

# seed for reproducibility & color palette for plots
set.seed(1234)
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(3,2,4:8)]


# print topology
#plotNetwork(topo, topoName)


racipe_cont_fname <- file.path(topoDir, "racipe_100IC_continued.Rds")
racipe <- readRDS(racipe_cont_fname)
genes <- rownames(racipe)
unnormData <- t(assay(racipe))
simExp <- assay(racipe, 1)
simExp <- log2(1+simExp)
tmpMeans <- rowMeans(simExp)
tmpSds <- apply(simExp,1,sd)

racipeNorm <- sracipeNormalize(racipe)
racipeData <- as.data.frame(t(assay(racipeNorm)))
exprMat_norm <- racipeData

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


########## MULTISTABILITY HISTOGRAM ############
## This will plot a histogram showing the number of states per model in initial simulations
summary_df_fname <- file.path(dataDir,"state_summary_df.Rds")
summary_df <- readRDS(summary_df_fname)
summary_hist_df_in <- summary_df
# reformatting to match expected col names & values
colnames(summary_hist_df_in) <- c("MODEL_NO", "NO_STATES", "StateIdentity")
summary_hist_df_in[which(summary_hist_df_in$StateIdentity == 1), "StateIdentity"] <- "E"
summary_hist_df_in[which(summary_hist_df_in$StateIdentity == 2), "StateIdentity"] <- "M"
summary_hist_df_in[which(summary_hist_df_in$StateIdentity == "bistable"), "StateIdentity"] <- "Bistable"

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

pca_fname <- file.path(dataDir,"pca_2025.Rds")
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
## Plot PCA loadings for all genes

loadings_df <- as.data.frame(pca$rotation)
loadings_df$gene <- rownames(loadings_df)

image <- ggplot(loadings_df, aes(x = PC1, y = PC2, label = gene)) +
  geom_point(color = "darkblue", size=4) +
  geom_hline(yintercept=0, color="black", alpha=0.5) +
  geom_vline(xintercept=0, color="black", alpha=0.5) +
  geom_text_repel(size = 10, hjust = 0, vjust = 1, 
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

loading_plot_fname <- file.path(plotDir, "PCA_loadings.pdf")
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
num_pcs_clustering <- 15 # First 15 PCs cover ~93% of variance

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

models_selected_fname <- file.path(dataDir,"racipe_bistable_indexMap_2025.Rds")
models_selected <- readRDS(models_selected_fname)

ss_unique_fname <- file.path(dataDir,"ss_unique_df.Rds")
ss_unique <- readRDS(ss_unique_fname)
all_bistable_states <- which(ss_unique$Model %in% models_selected)


# first, plot all E/M bistable states
image <- ggplot(pca_df[all_bistable_states,], aes(x=PC1, y=PC2, color=as.factor(clust_full[all_bistable_states]))) +
  geom_point(size=3) +
  guides(color=guide_legend(title = "Cluster")) +
  scale_color_manual(values=cbPalette) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
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
image <- ggplot(pca_df[all_bistable_states,], aes(x=PC1, y=PC2, color=as.factor(clust_full[all_bistable_states]))) +
  geom_point(size=3) +
  guides(color=guide_legend(title = "Cluster")) +
  scale_color_manual(values=cbPalette) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  xlim(-5, 7) +
  ylim(-4.5, 2.5) + 
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
  geom_point(data=pca_df[which(ss_unique$Model %in% models_selected),], 
             aes(x=PC1, y=PC2, fill=as.factor(clust_full[which(ss_unique$Model %in% models_selected)])), 
             size=3, color="black", pch=21) +
  guides(color=guide_legend(title = "Cluster"), fill="none") +
  scale_color_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
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

cluster_degs_all <- assign_cluster_DEGs(t(ss_unique[,genes]), factor(clust_full))
cluster_degs_bistable <- assign_cluster_DEGs(t(ss_unique[which(ss_unique$Model %in% models_selected),]), 
                                             factor(clust_full[which(ss_unique$Model %in% models_selected)]))


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
image <- Heatmap(as.matrix(t(exprMat_norm[ss_unique$StateIndex,1:nGenes])), 
                 name = "Expression", 
                 top_annotation = column_annotation,
                 row_names_gp=gpar(fontsize=12),
                 show_column_names = F)
image

wt_hmap_fname <- file.path(plotDir,"wt_hmap_allStates.pdf")
pdf(wt_hmap_fname, height = 10, width = 10)
print(image)
dev.off()

wt_hmap_fname_jpg <- file.path(plotDir,"wt_hmap_allStates.jpg")
jpeg(wt_hmap_fname_jpg, height = 6, width = 6, units="in", res = 300)
print(image)
dev.off()

## Now for just bistable states
ha_df <- data.frame(Cluster=clust_full[which(ss_unique$Model %in% models_selected)])

# Create an annotation object for the columns
column_annotation <- HeatmapAnnotation(df = ha_df, 
                                       col=list(Cluster=c("1"=unname(cbPalette[1]),"2"=unname(cbPalette[2]),
                                                          "3"=unname(cbPalette[3]),"4"=unname(cbPalette[4]))))
# Create the heatmap with annotation
image <- Heatmap(as.matrix(t(exprMat_norm[ss_unique[which(ss_unique$Model %in% models_selected),"StateIndex"],genes])), 
                 name = "Expression", 
                 top_annotation = column_annotation,
                 row_names_gp=gpar(fontsize=12),
                 show_column_names = F)
image

wt_hmap_fname <- file.path(plotDir,"wt_hmap_bistableStates.pdf")
pdf(wt_hmap_fname, height = 10, width = 10)
print(image)
dev.off()

wt_hmap_fname_jpg <- file.path(plotDir,"wt_hmap_bistableStates.jpg")
jpeg(wt_hmap_fname_jpg, height = 6, width = 6, units="in", res = 300)
print(image)
dev.off()




########## NOISE-ONLY CONTROL - SYMMETRY ############
ctrl_data_dir <- file.path(dataDir, "noise_only_controls")
ctrl_noise_levels <- c(0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.5, 1, 2)
ctrl_trials_per_noise <- 10
ctrl_tcorr <- 10
ctrl_results_df_fname <- file.path(ctrl_data_dir,
                                   paste0("ctrl_summary_trials=",ctrl_trials_per_noise,
                                          "_noise=",paste0(ctrl_noise_levels, collapse = ","),
                                          "_tcorr=",ctrl_tcorr))
ctrl_results_df <- readRDS(ctrl_results_df_fname)


image <- ggplot(ctrl_results_df, aes(x=as.factor(Noise), y=PctConverting*100, fill=as.factor(Noise))) +
  geom_boxplot() +
  ylab("Models Undergoing EMT (%)") +
  xlab("Noise") +
  ylim(0, 75) +
  guides(fill=guide_legend(title="Noise")) +
  theme_sticcc() +
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color="black"),
        plot.background = element_rect("white"))
image

plot_fname <- file.path(plotDir,"noiseOnlyControl_efficacy_EMT_tau=10.pdf")
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()

image <- ggplot(ctrl_results_df, aes(x=as.factor(Noise), y=PctRebellious*100, fill=as.factor(Noise))) +
  geom_boxplot() +
  ylab("Models Undergoing MET (%)") +
  xlab("Noise") +
  ylim(0, 75) +
  guides(fill=guide_legend(title="Noise")) +
  theme_sticcc() +
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color="black"),
        plot.background = element_rect("white"))
image

plot_fname <- file.path(plotDir,"noiseOnlyControl_efficacy_MET_tau=10.pdf")
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()



ctrl_long <- ctrl_results_df %>%
  select(Noise, PctConverting, PctRebellious) %>%
  pivot_longer(cols = c(PctConverting, PctRebellious),
               names_to = "TransitionType",
               values_to = "Pct") %>%
  mutate(TransitionType = recode(TransitionType,
                                 PctConverting = "EMT",
                                 PctRebellious = "MET"),
         Pct = Pct * 100)

# Create the combined plot
image <- ggplot(ctrl_long, aes(x = as.factor(Noise), y = Pct, fill = TransitionType)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  ylab("Models Undergoing Transition (%)") +
  xlab("Noise") +
  ylim(0, 65) +
  guides(fill = guide_legend(title = "Transition")) +
  theme_sticcc() +
  scale_fill_manual(values = cbPalette[c(6,5)]) +
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color = "black"),
        plot.background = element_rect("white"))
image

# Save the plot
plot_fname <- file.path(plotDir, "noiseOnlyControl_efficacy_combined_tau=10.pdf")
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()


########## MODEL-WISE SYMMETRY ############
ctrl_summary_df_fname <- file.path(ctrl_data_dir,"ctrl_summary_df_modelwise.Rds")
ctrl_summary_df <- readRDS(ctrl_summary_df_fname)

ctrl_summary_df$bias_score <- ctrl_summary_df$EMT_Count - ctrl_summary_df$MET_Count
ctrl_summary_df$bias_ratio <- (ctrl_summary_df$EMT_Count - ctrl_summary_df$MET_Count) / (ctrl_summary_df$EMT_Count + ctrl_summary_df$MET_Count + 1e-6)
bias_by_model <- ctrl_summary_df %>%
  dplyr::group_by(Model) %>%
  dplyr::summarize(num_transitions_total = sum(EMT_Count + MET_Count),
                   bias_score = sum(EMT_Count - MET_Count),
                   bias_ratio = sum(EMT_Count - MET_Count) / sum(EMT_Count + MET_Count))
bias_by_model$direction <- case_when(
  bias_by_model$bias_ratio > 0.3 ~ "EMT-biased",
  bias_by_model$bias_ratio < -0.3 ~ "MET-biased",
  TRUE ~ "Neutral"
)



top_models <- bias_by_model %>%
  dplyr::slice_max(abs(bias_score), n = 20)
top_models


## Heatmap of model-wise bias by noise
mat_df <- ctrl_summary_df %>%
  mutate(bias = EMT_Count - MET_Count) %>%
  select(Model, Noise, bias) %>%
  pivot_wider(names_from = Noise, values_from = bias, values_fill = 0) %>%
  column_to_rownames("Model")
bias_mat <- as.matrix(mat_df)
image <- Heatmap(
  bias_mat,
  name = "EMT - MET",
  col = colorRamp2(c(-max(abs(bias_mat)), 0, max(abs(bias_mat))), c("steelblue", "white", "firebrick")),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  row_title = "Models",
  column_title = "Noise Levels",
  row_title_gp = gpar(fontsize=30),
  column_title_gp = gpar(fontsize=30),
  column_names_gp = gpar(fontsize=28),
  heatmap_legend_param = list(title = "Bias Score")
)
image


modelwise_bias_heatmap_fname <- file.path(plotDir, "modelwise_bias_heatmap.pdf")
pdf(file = modelwise_bias_heatmap_fname, width = 10, height = 10)
print(image)
dev.off()

## Bar chart of top biased models
top_models <- bias_by_model %>%
  dplyr::slice_max(abs(bias_score), n = 20)

ggplot(top_models, aes(x = reorder(Model, bias_score), y = bias_score, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("EMT-biased" = "firebrick", "MET-biased" = "steelblue", "Neutral" = "grey70")) +
  labs(title = "Top Biased Models", x = "Model", y = "EMT – MET Count") +
  theme_minimal()



## Multiple ICs - plot basin boundaries
asym_models_use <- c(144, 589, 601)
asym_model_use <- 601 # Using the first E/M bistable model in the list
asym_model_params <- sracipeParams(racipe)[asym_model_use,]
asym_model_states <- ss_unique[which(ss_unique$Model == asym_model_use),genes]
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
asym_multIC_ss_clusts <- knn_classifier(asym_multIC_ss_pca[,1:num_pcs_clustering], pca_df_full[,1:num_pcs_clustering], clust_full, k=25)

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

########## MODEL BIAS CORRELATIONS ############
# Look at differences in parameters between E-biased and M-biased models
racipe_bistable_fname <- file.path(dataDir,"racipe_bistable_2025.Rds")
racipe_bistable_raw <- readRDS(racipe_bistable_fname)

emt_biased_models <- which(bias_by_model$direction == "EMT-biased")
met_biased_models <- which(bias_by_model$direction == "MET-biased")

emt_params <- sracipeParams(racipe_bistable_raw)[emt_biased_models,]
met_params <- sracipeParams(racipe_bistable_raw)[met_biased_models,]

param_test_results <- data.frame(
  Parameter = character(),
  p_value = numeric(),
  median_EMT = numeric(),
  median_MET = numeric(),
  stringsAsFactors = FALSE
)

# Loop through parameters
for (paramNo in 1:ncol(emt_params)) {
  param_name <- colnames(emt_params)[paramNo]
  emt_vals <- emt_params[, paramNo]
  met_vals <- met_params[, paramNo]
  
  # Perform Wilcoxon test
  test_result <- wilcox.test(emt_vals, met_vals)
  
  # Record results
  param_test_results <- rbind(param_test_results, data.frame(
    Parameter = param_name,
    p_value = test_result$p.value,
    median_EMT = median(emt_vals),
    median_MET = median(met_vals),
    stringsAsFactors = FALSE
  ))
}

# adjust for multiple testing
param_test_results$adj_p_value <- p.adjust(param_test_results$p_value, method = "BH")

# View top hits
param_test_results <- param_test_results %>% arrange(adj_p_value)
print(head(param_test_results))


# Select top k parameters
top_params <- param_test_results %>%
  arrange(adj_p_value) %>%
  slice(1:10) %>%
  pull(Parameter)

# Prepare long-format data for ggplot
emt_df <- emt_params[, top_params, drop=FALSE] %>%
  as.data.frame() %>%
  mutate(Bias = "EMT-biased")

met_df <- met_params[, top_params, drop=FALSE] %>%
  as.data.frame() %>%
  mutate(Bias = "MET-biased")

plot_df <- bind_rows(emt_df, met_df) %>%
  pivot_longer(
    cols = -Bias,
    names_to = "Parameter",
    values_to = "Value"
  )

# Boxplot
ggplot(plot_df, aes(x = Parameter, y = Value, fill = Bias)) +
  geom_boxplot(outlier.size = 0.5, linewidth = 0.3) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  scale_y_log10() +
  scale_fill_manual(values = c("EMT-biased" = "#E64B35", "MET-biased" = "#4DBBD5")) +
  labs(
    x = "Parameter",
    y = "Parameter Value",
    fill = "Model Bias",
    title = paste("Top", 10, "Differential Parameters Between EMT- and MET-biased Models")
  )

########## SIGNAL EFFICACY HEATMAP ############
# Plot 2D heatmap at different noise levels
# Note: only works with 2-gene signals included
signal_simTime <- 500
signal_relaxTime <- 50
signal_nGenes <- c(1,2)
signal_noise <- c(0, 0.04, 0.2) #0.04
signal_tcorr <- 10
initClust <- 1
tgtClust <- 2
expName <- paste0("bhtopo_t=",signal_simTime,"_relax_OUnoise=",paste0(signal_noise, collapse = "."),
                  "_tau=",signal_tcorr,"_genes=",paste0(signal_nGenes,collapse = "."),"_CLAMPS_2025")
resultSet_fname <- file.path(topoDir,expName,"result_summary.Rds")
resultSet_full <- readRDS(resultSet_fname)
selectedNoise <- 0.04
eff_hmap_0.04 <- resultSet_full[which(resultSet_full$Noise == selectedNoise),
                                c("Species 1", "Species 2", "ConversionPct")]

eff_hmap_0.04_matrix <- matrix(nrow = length(genes), ncol = length(genes))
rownames(eff_hmap_0.04_matrix) <- genes_reordered
colnames(eff_hmap_0.04_matrix) <- genes_reordered
for(gene1 in genes_reordered) {
  for (gene2 in genes_reordered) {
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



########## SIGNAL EFFICACY VS TOPOLOGY ############
## Plot signal efficacy against a network topology metrics of signals: closeness & betweenness centrality, out-degree 
resultSet <- resultSet_full
selectedNoise <- 0.04
image <- ggplot(resultSet[which(resultSet$Noise == selectedNoise),], aes(x=GroupBetweenCentrality, y=ConversionPct)) +
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






########## COMPARISON TO BOOLEAN MODEL ############
signal_summary_df_fname <- file.path(dataDir,"signal_summary_df.Rds")
signal_summary_df <- readRDS(signal_summary_df_fname)

## Compare 1-gene signals
sigEffs_1gene_racipe_fname <- file.path(dataDir, "sigEffs_1gene_comparison.Rds")
sigEffs_1gene_racipe <- readRDS(sigEffs_1gene_racipe_fname)


ggplot() +
  geom_point(data=sigEffs_1gene_racipe, aes(x=1:26, y=ConversionPct), color="red") +
  geom_point(data=sigEffs_1gene_racipe, aes(x=1:26, y=ConversionPct_Spin), color="blue")

diff_df <- data.frame(Signal=sigEffs_1gene_racipe$`Species 1`,
                      Diff=sigEffs_1gene_racipe$ConversionPct_Spin - sigEffs_1gene_racipe$ConversionPct)
ggplot(data=diff_df[which(!is.na(diff_df$Diff)),]) +
  geom_bar(stat="identity", aes(x=Signal, fill=Signal, y=Diff)) +
  ylab("Conversion %, Spin - RACIPE") +
  xlab("Signal") +
  theme_sticcc() +
  theme(axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle=90))

ggplot() +
  geom_point(data=sigEffs_1gene_racipe, aes(x=ConversionPct, y=ConversionPct_Spin), size=3) +
  xlab("RACIPE Conversion %") +
  ylab("Spin Conversion %") +
  theme_sticcc() +
  theme(axis.line = element_line(color = "black"))

indices <- which(!is.na(sigEffs_1gene_racipe$ConversionPct_Spin))
cor(sigEffs_1gene_racipe$ConversionPct[indices], sigEffs_1gene_racipe$ConversionPct_Spin[indices], 
    use = "complete.obs", method = "spearman")

## Compare 2-gene signals

sigEffs_2gene_racipe_fname <- file.path(dataDir, "sigEffs_2gene_comparison.Rds")
sigEffs_2gene_racipe <- readRDS(sigEffs_2gene_racipe_fname)


diff_df <- data.frame(Signal=paste0(sigEffs_2gene_racipe$`Species 1`,"_",sigEffs_2gene_racipe$`Species 2`),
                      Rank=order(sigEffs_2gene_racipe$ConversionPct),
                      Diff=sigEffs_2gene_racipe$ConversionPct_Spin - sigEffs_2gene_racipe$ConversionPct)
ggplot(data=diff_df[which(!is.na(diff_df$Diff)),]) +
  geom_point(aes(x=Rank, y=Diff), size=3) +
  ylab("Conversion %, Spin - RACIPE") +
  xlab("Signal Rank") +
  theme_sticcc() +
  theme(axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle=90))


ggplot() +
  geom_point(data=sigEffs_2gene_racipe, aes(x=ConversionPct, y=ConversionPct_Spin), size=3) +
  xlab("RACIPE Conversion %") +
  ylab("Spin Conversion %") +
  theme_sticcc() +
  theme(axis.line = element_line(color = "black"))


## Now combine with info from 5B (betweenness centrality)
sigEffs_2gene_racipe$Signal <- paste0(sigEffs_2gene_racipe$`Species 1`, "_", sigEffs_2gene_racipe$`Species 2`)
signal_summary_compare <- merge(sigEffs_2gene_racipe, signal_summary_df, by="Signal")
image <- ggplot() +
  geom_point(data=signal_summary_compare, aes(x=ConversionPct, y=ConversionPct_Spin, color=GroupBetweenCentrality), size=3) +
  xlab("RACIPE Conversion %") +
  ylab("Spin Conversion %") +
  scale_color_gradient(name="Betweenness") +
  theme_sticcc() +
  theme(axis.line = element_line(color = "black"))
image

plot_fname <- file.path(plotDir,paste0("fig4b_spin_vs_racipe_betweenness_noise=",selectedNoise,".pdf"))
pdf(plot_fname, width = 10, height = 10)
print(image)
dev.off()



cor(signal_summary_compare$ConversionPct, signal_summary_compare$ConversionPct_Spin,
    use = "complete.obs", method = "spearman")
cor(signal_summary_compare$ConversionPct, signal_summary_compare$GroupBetweenCentrality,
    use = "complete.obs", method = "spearman")
cor(signal_summary_compare$ConversionPct_Spin, signal_summary_compare$GroupBetweenCentrality,
    use = "complete.obs", method = "spearman")


indices <- which(!is.na(sigEffs_2gene_racipe$ConversionPct_Spin))
cor(sigEffs_2gene_racipe$ConversionPct[indices], sigEffs_2gene_racipe$ConversionPct_Spin[indices], 
    use = "complete.obs", method = "spearman")



# try removing signals where sigEff_spin is one
# sigEffs_2gene_racipe <- sigEffs_2gene_racipe[which(sigEffs_2gene_racipe$ConversionPct_Spin != 1),]
# cor(sigEffs_2gene_racipe$ConversionPct, sigEffs_2gene_racipe$ConversionPct_Spin, 
#     use = "complete.obs", method = "spearman")

## Bland-altman plot
ranks_racipe <- rank(sigEffs_2gene_racipe$ConversionPct)
ranks_spin <- rank(sigEffs_2gene_racipe$ConversionPct_Spin)

plot(ranks_racipe, ranks_spin)

# Calculate means and differences
means <- (ranks_racipe + ranks_spin) / 2
differences <- ranks_racipe - ranks_spin

spinEffOne <- as.factor(sigEffs_2gene_racipe$ConversionPct_Spin == 1)

# Create a Bland-Altman Plot
bland_altman_plot <- ggplot(data = NULL, aes(x = means, y = differences, color=spinEffOne)) +
  geom_point() +  # Add points
  geom_hline(yintercept = mean(differences), linetype = "dashed", color = "red") +  # Add mean line
  geom_hline(yintercept = mean(differences) + 1.96 * sd(differences), linetype = "dashed", color = "blue") +  # Add upper limit
  geom_hline(yintercept = mean(differences) - 1.96 * sd(differences), linetype = "dashed", color = "blue") +  # Add lower limit
  labs(x = "Average Rank", y = "Difference in Rank", title = "Bland-Altman Plot for Ranking Comparisons") +
  theme_minimal()

# Display the plot
print(bland_altman_plot)




########## SINGLE-MODEL TRAJECTORY ############
## TODO: imports, setup, clean up plots

hd_traj_model <- 10 # 10, 297, 144
hd_traj_simTime <- 100
hd_traj_tcorr <- 10
hd_traj_num_trials <- 10
hd_traj_sig <- "Zeb1_noise=0.04"
attemptNo <- 1
hd_traj_models <- c(10, 297, 144, 589, 601)

for(hd_traj_model in hd_traj_models) {
  modelDataDir <- file.path(dataDir,paste0("trajectories_model",hd_traj_model))
  if(!dir.exists(modelDataDir)) {
    dir.create(modelDataDir)
  }
  for(attemptNo in 1:hd_traj_num_trials) {
    traj_fname <- file.path(modelDataDir,paste0("trajectory_model",hd_traj_model,
                                                "_simTime=",hd_traj_simTime,
                                                "_SIG=",hd_traj_sig,"_tau=",hd_traj_tcorr,"_v",attemptNo,".Rds"))
    traj_pca_fname <- file.path(modelDataDir,paste0("trajectoryPCA_model",hd_traj_model,
                                                    "_simTime=",hd_traj_simTime,
                                                    "_SIG=",hd_traj_sig,"_tau=",hd_traj_tcorr,"_v",attemptNo,".Rds"))
    
    hd_traj_pca <- readRDS(traj_pca_fname)
    hd_traj_df_long <- readRDS(traj_fname)
    hd_traj_df_full <- pivot_wider(
      hd_traj_df_long,
      names_from = "Gene",      # Column with names to become new columns
      values_from = "Expression" # Column with values to fill the new columns
    )
    
    
    modelPlotDir <- file.path(plotDir,paste0("trajectories_model",hd_traj_model))
    if(!dir.exists(modelPlotDir)) {
      dir.create(modelPlotDir)
    }
    
    
    gene_subsets <- list(
      'E Markers' = c("Cdh1","miR141","miR34a","miR101",
                      "miR200a","miR200b","miR200c","Ovol2","Grhl2","Cldn7"),
      'M Markers' = c("Vim","Snai1","Snai2",
                      "Twist1","Twist2","Foxc2","Gsc",
                      "Klf8","Tcf3","Tgfbeta"),
      'Clamp Genes' = c("Zeb1")
      #'Disconnected' = c("Np63a","miR9","miR30c","miR205")
    )
    
    # Assign each gene to a subset
    hd_traj_df_long$Subset <- sapply(hd_traj_df_long$Gene, function(gene) {
      subset_name <- names(gene_subsets)[sapply(gene_subsets, function(subset) gene %in% subset)]
      if (length(subset_name) > 0) subset_name else NA
    })
    
    # Filter out genes not in any subset
    filtered_df <- hd_traj_df_long[!is.na(hd_traj_df_long$Subset),]
    
    # Calculate average and standard error for each subset and time point
    avg_traj_df <- filtered_df %>%
      group_by(Subset, Time) %>%
      summarize(
        Mean_Expression = mean(Expression),
        SE_Expression = sd(Expression) / sqrt(n()),
        .groups = "drop"
      )
    
    # Plot the average trajectory with error bars
    image <- ggplot(avg_traj_df[which(avg_traj_df$Time <= 200),], aes(x = Time, y = Mean_Expression, color = Subset)) +
      geom_line(linewidth=2) +
      theme_sticcc() +
      theme(
        axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color = "black"),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 28),
        plot.background = element_rect("white")
      ) +
      labs(y = "Average Expression", x = "Time", color = "Gene Subset", fill = "Gene Subset")
    
    image
    
    
    
    traj_plot_fname <- file.path(modelPlotDir, paste0("trajectory_model",hd_traj_model,
                                                      "_simTime=",hd_traj_simTime,
                                                      "_SIG=",hd_traj_sig,"_tau=",hd_traj_tcorr,"_v",attemptNo,".pdf"))
    pdf(traj_plot_fname, height = 10, width = 10)
    print(image)
    dev.off()
    
    
    
    
    # Plot trajectory, with model states in red
    hd_traj_model_idx <- as.numeric(gsub("Model","",names(models_selected)[which(models_selected == hd_traj_model)])) # for retrieving relevant state
    model_states <- ss_unique[which(ss_unique$Model == hd_traj_model), ]
    model_states_norm <- model_states
    model_states_norm[,1:length(genes)] <- log2(1+model_states_norm[,1:length(genes)]) # Log transform
    model_states_norm[,1:length(genes)] <- sweep(model_states_norm[,1:length(genes)], 2, tmpMeans, FUN = "-") # scale
    model_states_norm[,1:length(genes)] <- sweep(model_states_norm[,1:length(genes)], 2, tmpSds, FUN = "/") # scale
    model_states_pca <- as.data.frame(predict(pca, model_states_norm[,names(tmpMeans)]))
    model_states_pca
    
    
    # Ensure hd_traj_pca only includes Time <= 100
    hd_traj_pca_filtered <- hd_traj_pca %>%
    filter(Time <= 100)
    # Select points at every 10 time units for markers
    icon_points <- hd_traj_pca_filtered %>%
      filter(Time %% 10 == 0)
    
    # Create the plot
    image <- ggplot(pca_df_full[unique_state_idx_list_all,]) +
      geom_point(aes(x = PC1, y = PC2)) +  # Base points
      # Trajectory points and path
      #geom_point(data = hd_traj_pca_filtered, aes(x = PC1, y = PC2, color = Time)) +
      geom_path(data = hd_traj_pca_filtered, aes(x = PC1, y = PC2, color = Time)) +
      # Highlight specific time points with icons (customize shape if needed)
      geom_point(data = icon_points, aes(x = PC1, y = PC2, color=Time), size = 3) +
      # Model states in red
      geom_point(data = model_states_pca, aes(x = PC1, y = PC2), color = "red", size = 3) +
      theme_sticcc() +
      xlab(plot_xlab) +
      ylab(plot_ylab) +
      theme(axis.line = element_line(linewidth = 1, color = "black"), 
            axis.ticks = element_line(linewidth = 1, color = "black"),
            axis.title = element_text(size = 32),
            axis.text = element_text(size = 28),
            plot.background = element_rect("white")
      ) +
      scale_color_gradient() 
    image
    
    
    traj_plot_pca_fname <- file.path(modelPlotDir, paste0("trajectoryPCA_model",hd_traj_model,
                                                          "_simTime=",hd_traj_simTime,
                                                          "_SIG=",hd_traj_sig,"_tau=",hd_traj_tcorr,"_v",attemptNo,".pdf"))
    pdf(traj_plot_pca_fname, height = 10, width = 10)
    print(image)
    dev.off()
    
  }
}







########## TRANSITION TIME VS NOISE ############
num_trials <- 10
time_trial_noise_levels <- c(0, 0.02, 0.04, 0.08, 0.1, 0.2, 0.5)
time_trial_simTime <- 300
times <- c(seq(2, 30, 2), seq(35, 100, 5), seq(120, time_trial_simTime, 20))

time_trial_resultSet_fname <- file.path(dataDir,"Zeb1_timeTrial_paramSets.Rds")
time_trial_resultSet <- readRDS(time_trial_resultSet_fname)
time_trial_setIDList <- rownames(time_trial_resultSet)
time_trial_sig_names <- time_trial_resultSet$SetName#paste0(time_trial_sig_gene,"_noise=",time_trial_noise_levels)
time_trial_expName_list <- paste0("bhtopo_timeTrial_t=",time_trial_simTime,
                                  "_relax_OUnoise=",time_trial_resultSet[time_trial_setIDList, "Noise"],
                                  "_tau=",signal_tcorr,
                                  "_SIG=",time_trial_resultSet[time_trial_setIDList, "SetName"],
                                  "_runNo=",rep(seq(num_trials),each=length(time_trial_noise_levels)))

## multiSet_times (import)
## times (??)

multiSet_fname <- file.path(dataDir,
                            paste0("Zeb1_transition_time_trials=",num_trials,
                                   "_noise=",paste0(time_trial_noise_levels, collapse = ","),"_",".Rds"))
multiSet_times <- readRDS(multiSet_fname)




## Plot conversions vs time
plot_df_list <- list()
idx <- 1
for(setIDNo in seq_along(time_trial_setIDList)) {
  setID <- time_trial_setIDList[setIDNo]
  setName <- time_trial_resultSet[setID, "SetName"]
  sampleSet_times <- multiSet_times[[setIDNo]]
  
  for(timeID in seq_along(times)) {
    time <- times[timeID]
    newrow <- list(Signal=setName, Time=time, Conversions=sampleSet_times[[timeID]][[1]])
    plot_df_list[[idx]] <- newrow
    idx <- idx+1
  }
  
}

plot_df <- do.call(rbind, lapply(plot_df_list, function(x) as.data.frame(t(x))))
plot_df$Time <- as.numeric(plot_df$Time)
plot_df$Conversions <- as.numeric(plot_df$Conversions)
plot_df$Signal <- as.character(plot_df$Signal)
plot_df$Noise <- as.numeric(gsub("Zeb1_noise=", "", plot_df$Signal))
plot_df$ConversionPct <- plot_df$Conversions / 500 * 100


summary_df <- plot_df %>%
  group_by(Time, Noise, Signal) %>%
  summarise(
    mean_conv = mean(ConversionPct),
    sd_conv = sd(ConversionPct),
    .groups = "drop"
  )

# Step 2: Plot with error bars
image <- ggplot(summary_df[which(summary_df$Noise < 0.5),], aes(x = Time, y = mean_conv, color = as.factor(Noise))) +
  geom_line(linewidth = 1.2) +
  geom_errorbar(aes(ymin = mean_conv - sd_conv, ymax = mean_conv + sd_conv), 
                width = 0.2, linewidth = 0.8) +
  #facet_wrap(~Signal) +  # Optional: only if Signal should be faceted
  theme_sticcc() +
  theme(
    axis.line = element_line(linewidth = 1, color = "black"), 
    axis.ticks = element_line(linewidth = 1, color = "black"),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 24),
    plot.background = element_rect("white")
  ) +
  labs(
    y = "Conversion % (mean ± SD)",
    color = "Noise"
  )

image


transitions_vs_time_vs_noise_fname <- file.path(plotDir,paste0("transitions_vs_noise_sigGene=Zeb1.pdf"))
pdf(file = transitions_vs_time_vs_noise_fname, width = 10, height = 10)
print(image)
dev.off()



