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
all_bistable_states <- which(ss_unique$Model %in% models_selected)


# first, plot all E/M bistable states
image <- ggplot(pca_df_full[all_bistable_states,], aes(x=PC1, y=PC2, color=as.factor(clust_full[all_bistable_states]))) +
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
image <- ggplot(pca_df_full[all_bistable_states,], aes(x=PC1, y=PC2, color=as.factor(clust_full[racipe_bistable_indices]))) +
  geom_point(size=3) +
  guides(color=guide_legend(title = "Cluster")) +
  scale_color_manual(values=cbPalette) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  xlim(-5, 10) +
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
image <- ggplot(pca_df_full[,], aes(x=PC1, y=PC2, color=as.factor(clust_full))) +
  geom_point(size=3, alpha=0.7) +
  geom_point(data=pca_df_full[which(ss_unique$Model %in% models_selected),], 
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



## Compare 1-gene signals
sigEffs_1gene_racipe_fname <- file.path(dataDir, "sigEffs_1gene_comparison.Rds")
sigEffs_1gene_racipe <- readRDS(sigEffs_1gene_racipe_fname)


ggplot() +
  geom_point(data=sigEffs_1gene_racipe, aes(x=1:26, y=ConversionPct), color="red") +
  geom_point(data=sigEffs_1gene_spin, aes(x=1:26, y=ConversionPct), color="blue")

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





########## TRANSITION TIME VS NOISE ############

## setIDList (construct)
## time_trial_resultSet (import)
## multiSet_times (import)
## times (??)

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

image <- ggplot(plot_df[which(plot_df$Signal %in% c("Zeb1_noise=0.04", "Zeb1_noise=0")),], aes(x=Time, y=Conversions, color=Signal)) +
  geom_point() +
  geom_path() +
  ggtitle(paste0(expName,"_20ManualSignals"))

image





