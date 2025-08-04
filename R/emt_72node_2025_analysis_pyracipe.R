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
library(mclust)
library(cluster) # for silhouette index
source("R/utils.R")
source("R/utils_clamping.R")
source("R/scratch.R")

# set up directories
topoName <- "emt_ffctopo_72node_CLAMP_2025"
topoDir <- file.path(getwd(),topoName)
plotDir <- file.path(topoDir,"plots")
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





########## PCA & CLUSTERING UNIQUE STATES ############

# PCA
pca_fname <- file.path(dataDir,"pca_uniques.Rds")
if(!file.exists(pca_fname)) {
  # PCA on full data
  pca <- prcomp(pyracipeNorm[,1:nGenes])
  pca_df <- as.data.frame(pca$x)
  
  # We will flip PC1 to put E states on the left
  pca_df$PC1 <- -1 * pca_df$PC1
  pca$rotation[,1] <- -1 * pca$rotation[,1]
  pca$x[,1] <- -1 * pca$x[,1]
  
  # sanity check
  ggplot(pca_df, aes(x=PC1, y=PC2, color=pyracipeNorm$Ecadherin)) + geom_point()
  
  saveRDS(pca, pca_fname)
} else {
  pca <- readRDS(pca_fname)
  pca_df <- as.data.frame(pca$x)
}


# Cluster all steady states
# Helper function to select k
calc_silhouette <- function(X, k_range) {
  silhouette_means <- c()
  nll_values <- c()
  
  for (k in k_range) {
    gmm_model <- Mclust(X, G = k, verbose = FALSE)
    clusters <- gmm_model$classification
    sil <- silhouette(clusters, dist(X))
    silhouette_means <- c(silhouette_means, mean(sil[, 3]))
  }
  
  df <- data.frame(Clusters=k_range, Silhouette=silhouette_means)
  
  return(df)
  
}


clust_all_fname <- file.path(dataDir,"clust_all_2025.Rds") 
sil_df_fname <- file.path(dataDir, "silhouette_df.Rds")
gmm_fname <- file.path(dataDir, "clust_gmm.Rds")
num_pcs_clustering <- 45 # First 45 PCs cover ~91% of variance
if(!file.exists(clust_all_fname)) { 
  
  # Identify optimal k
  if(!file.exists(sil_df_fname)) {
    k_range <- 2:8
    
    sil_df <- calc_silhouette(pca_df[,1:num_pcs_clustering], k_range)
    saveRDS(sil_df, sil_df_fname)
  }
  
  # Silhouette suggests k = 2, which is used now
  gmm = GMM(pca_df[,1:num_pcs_clustering], 2, 
            dist_mode = "eucl_dist", seed_mode = "random_subset", km_iter = 10,
            em_iter = 10, verbose = F)

  
  # predict centroids, covariance matrix and weights
  clust_full = predict(gmm, newdata = pca_df[,1:num_pcs_clustering])
  
  # Sanity check plots
  ggplot(pca_df, aes(x=PC1, y=PC2, color=as.factor(clust_full))) + geom_point()
  ggplot(pca_df, aes(x=PC1, y=PC2, color=pyracipeNorm$Ecadherin)) + geom_point()
  
  # Im reversing it so the initial state (E) is cluster 1
  revClust <- clust_full
  revClust <- ifelse(clust_full == 1, 2, ifelse(clust_full == 2, 1, clust_full))
  
  perm <- c(2,1)
  gmm$centroids <- gmm$centroids[perm, , drop = FALSE]
  gmm$covariance_matrices <- gmm$covariance_matrices[perm, , drop=FALSE]
  gmm$weights <- gmm$weights[perm]
  gmm$Log_likelihood <- gmm$Log_likelihood[, perm]
  
  saveRDS(gmm, file = gmm_fname)
  saveRDS(revClust, file = clust_all_fname)
} else {
  clust_full <- readRDS(clust_all_fname)
  gmm <- readRDS(gmm_fname)
}


########## SELECT E/M BISTABLE MODELS ############
pyracipe_metadata_rds_fname <- file.path(dataDir,"pyracipe_metadata.Rds")
pyracipe_summary_rds_fname <- file.path(dataDir,"pyracipe_summary.Rds")

if("Cluster" %in% colnames(pyracipe_metadata)) {
  pyracipe_metadata <- readRDS(pyracipe_metadata_rds_fname)
  pyracipe_summary <- readRDS(pyracipe_summary_rds_fname)
} else {
  pyracipe_metadata$Cluster <- clust_full
  pyracipe_metadata$StateIdentity <- NA
  for(modelNo in pyracipe_summary$MODEL_NO) {
    model_clusters <- clust_full[which(pyracipe_metadata == modelNo)]
    
    if(1 %in% model_clusters & 2 %in% model_clusters) {
      pyracipe_metadata[which(pyracipe_metadata$MODEL_NO == modelNo), "StateIdentity"] <- "Bistable"
      pyracipe_summary[which(pyracipe_summary$MODEL_NO == modelNo), "StateIdentity"] <- "Bistable"  
    } else if(1 %in% model_clusters) {
      pyracipe_metadata[which(pyracipe_metadata$MODEL_NO == modelNo), "StateIdentity"] <- "E"  
      pyracipe_summary[which(pyracipe_summary$MODEL_NO == modelNo), "StateIdentity"] <- "E"  
    } else if(2 %in% model_clusters) {
      pyracipe_metadata[which(pyracipe_metadata$MODEL_NO == modelNo), "StateIdentity"] <- "M"  
      pyracipe_summary[which(pyracipe_summary$MODEL_NO == modelNo), "StateIdentity"] <- "M"  
    }
  }
  
  saveRDS(pyracipe_metadata, pyracipe_metadata_rds_fname)
  saveRDS(pyracipe_summary, pyracipe_summary_rds_fname)
}


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

# sanity check: plot PCA with only kept states, and only bistable ones
ggplot(pca_df[keep_states,]) + 
  geom_point(aes(x=PC1, y=PC2, color=factor(clust_full[keep_states]))) +
  xlim(-5, 10) +
  ylim(-6, 7)
ggplot(pca_df[keep_states_e,]) + geom_point(aes(x=PC1, y=PC2)) +
  xlim(-5, 10) +
  ylim(-6, 7)


########## CONVERT PYRACIPE PARAMETERS ############

replacement_map <- c(MPR = "G", DNR = "K", TSH = "TH", HCO = "N", FCH = "FC")
all_replaced_strings <- replace_first_substring(colnames(pyracipe_params), replacement_map)
all_replaced_strings <- gsub("-","_",all_replaced_strings)
colnames(pyracipe_params) <- unname(all_replaced_strings)
paramColOrder <- colnames(sracipeParams(racipe))
pyracipe_params <- pyracipe_params[,paramColOrder]

# Need to reciprocate columns with FC values (429:end)
for(colID in 429:570) {
  g1 <- strsplit(colnames(pyracipe_params)[colID], "_")[[1]][2]
  g2 <- strsplit(colnames(pyracipe_params)[colID], "_")[[1]][3]
  
  edgeType <- topo[which(topo$Source == g1 & topo$Target == g2), "Type"]
  
  if(edgeType == 2) {
    pyracipe_params[,colID] <- 1/pyracipe_params[,colID]
  }
  
}


########## IDENTIFY CLAMP VALUES ############
## Target: dataframe with model id, gene, expression, cluster (very long, maybe a wider format)
combined_df <- cbind(pyracipe_metadata, pyracipe_raw_states)
colnames(combined_df)[1:4] <- c("Model", "NO_STATES", "STATE_NO", "StateIndex")

clamp_df_fname <- file.path(dataDir,"clamp_values_2025.Rds")
clamp_df_full_fname <- file.path(dataDir,"clamp_values_all_2025.Rds")
if(!file.exists(clamp_df_fname)) {
  
  clamp_df <- pivot_longer(combined_df[keep_states,], cols = all_of(genes), 
                           names_to = "Gene", values_to = "Expression")
  clamp_df$ModelIndex <- as.numeric(factor(clamp_df$Model))
  
  saveRDS(clamp_df, clamp_df_fname)
  
  
  clamp_df_full <- pivot_longer(combined_df[,], cols = all_of(genes), 
                           names_to = "Gene", values_to = "Expression")
  clamp_df_full$ModelIndex <- as.numeric(factor(clamp_df_full$Model))
  saveRDS(clamp_df_full, clamp_df_full_fname)
  
} else {
  clamp_df <- readRDS(clamp_df_fname)
  clamp_df_full <- readRDS(clamp_df_full_fname)
}





########## PREP SIGNALING SIMULATIONS ############

# Using only bistable models (and beginning in E state), we will arrange the signaling simulations now
# necessary inputs:
# racipe object (using E states of bistable models as ICs, and matched parameters)
# pca (only of E/M bistable states)
# cluster labels
# clamp df
# noise level, tcorr
# num genes
# sim time & relaxation time

signal_simTime <- 500
signal_relaxTime <- 50
signal_nGenes <- c(1,2)
signal_noise <- c(0, 0.02, 0.04) 
signal_tcorr <- 10
initClust <- 1
tgtClust <- 2


numModels_signaling <- length(keep_states_e)
ics_signaling <- t(pyracipe_raw_states[keep_states_e,genes])
params_signaling <- pyracipe_params[keep_models,]

# placeholder simulation to initialize RACIPE object
placeholder_fname <- file.path(dataDir,"racipe_signaling_placeholder.Rds") 
if (!file.exists(placeholder_fname)) {
  racipe_signaling <- sracipeSimulate(topo, numModels = numModels_signaling, nIC = 1, 
                                      genIC = T, genParams = T, integrate = F,
                                      ouNoise_t = signal_tcorr, initialNoise = 0, 
                                      simulationTime = 1, simDet = T)
  sracipeIC(racipe_signaling) <- ics_signaling
  
  
  #assay(racipe_signaling, withDimnames = F) <- ics_signaling
  sracipeParams(racipe_signaling) <- params_signaling
  racipe_signaling@metadata$config$options[["genIC"]] <- F
  racipe_signaling@metadata$config$options[["genParams"]] <- F
  
  racipe_signaling_relaxed <- sracipeSimulate(racipe_signaling, numModels = numModels_signaling, nIC = 1, 
                                      genIC = F, genParams = F, integrate = T,
                                      ouNoise_t = signal_tcorr, initialNoise = 0, 
                                      simulationTime = 50, simDet = T)
  
  saveRDS(racipe_signaling_relaxed, placeholder_fname)
  racipe_signaling <- racipe_signaling_relaxed
} else {
  racipe_signaling <- readRDS(placeholder_fname)
}


pca_signaling <- pca
pca_signaling$x <- pca_signaling$x[keep_states_e,]

clust_signaling <- clust_full[keep_states_e]




######## NOISE-ONLY CONTROL SIMULATIONS #####
## Simulate trials with various noise levels & no (effective) signal
## Track number of transitions over time

ctrl_simTime <- 300
ctrl_relaxTime <- 50
ctrl_tcorr <- 10
initClust <- 1
tgtClust <- 2



ctrl_data_dir <- file.path(dataDir, "noise_only_controls")
if(!dir.exists(ctrl_data_dir)) {
  dir.create(ctrl_data_dir)
}


ctrl_noise_levels <- c(0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.5, 1, 2)
ctrl_trial_expr <- as.data.frame(t(assay(racipe_signaling)))
ctrl_trials_per_noise <- 5

ctrl_indices <- c()
for(model in keep_models) {
  clust_pick <- sample(c(1,2),1)
  idx_add <- pyracipe_metadata[which(pyracipe_metadata$MODEL_NO == model), "Index"]
  ctrl_indices <- c(ctrl_indices, idx_add)
}

# create RACIPE object with equal split of E and M steady states to begin
racipe_ctrl_placeholder <- sracipeSimulate(topo, numModels = (2*length(keep_models)), nIC = 1, integrate = F,
                                           genIC = T, genParams = T, 
                                           ouNoise_t = signal_tcorr, initialNoise = 0, 
                                           simulationTime = 1, simDet = T)
sracipeParams(racipe_ctrl_placeholder) <- pyracipe_params[rep(as.numeric(unname(keep_models)), each=2),]
sracipeIC(racipe_ctrl_placeholder) <- t(pyracipe_raw_states[ctrl_indices,genes])
racipe_ctrl_placeholder@metadata$config$simParams[["numModels"]] <- length(ctrl_indices)
racipe_ctrl_placeholder@metadata$config$simParams["nIC"] <- 1
ctrl_init_data_norm <- pyracipe_raw_states[ctrl_indices,genes]
ctrl_init_data_norm[,genes] <- log2(1+ctrl_init_data_norm[,genes]) # Log transform
ctrl_init_data_norm[,genes] <- sweep(ctrl_init_data_norm[,genes], 2, tmpMeans, FUN = "-") # scale
ctrl_init_data_norm[,genes] <- sweep(ctrl_init_data_norm[,genes], 2, tmpSds, FUN = "/") # scale
ctrl_init_pca <- scale(ctrl_init_data_norm, pca$center, pca$scale) %*% pca$rotation
#ctrl_init_clusts_knn <- knn_classifier(ctrl_init_pca[,1:num_pcs_clustering], pca_df[,1:num_pcs_clustering], clust_full, k=25)
ctrl_init_clusts = predict(gmm, newdata = ctrl_init_pca[,1:num_pcs_clustering])

# compare ctrl cluster labels to original labels
ggplot(ctrl_init_pca, aes(x=PC1, y=PC2, color=ctrl_init_clusts)) + 
  geom_point() + xlim(-5, 14) + ylim(-6, 7) + labs(color="Cluster")
ggplot(ctrl_init_pca, 
       aes(x=PC1, y=PC2, 
           color=factor(pyracipe_metadata[which(pyracipe_metadata$MODEL_NO %in% keep_models),"Cluster"]))) + 
  geom_point() + xlim(-5, 14) + ylim(-6, 7) + labs(color="Cluster")
# look at an example model which has its cluster label change
orig_clusts <- pyracipe_metadata[which(pyracipe_metadata$MODEL_NO %in% keep_models),"Cluster"]
consistent_models <- which((ctrl_init_clusts == 1 & orig_clusts == 1) |(ctrl_init_clusts == 2 & orig_clusts == 2))
e_to_m_models <- which(ctrl_init_clusts == 2 & orig_clusts == 1)
m_to_e_models <- which(ctrl_init_clusts == 1 & orig_clusts == 2)
m_to_e_models_indices <- ctrl_indices[m_to_e_models]
m_to_e_models_modelNos <- pyracipe_metadata[which(pyracipe_metadata$Index %in% m_to_e_models_indices),"MODEL_NO"]
m_to_e_models_other_states <- pyracipe_metadata[which(pyracipe_metadata$MODEL_NO %in% m_to_e_models_modelNos),"Index"]
m_to_e_models_other_states <- m_to_e_models_other_states[which(!m_to_e_models_other_states %in% m_to_e_models_indices)]
ggplot(ctrl_init_pca[m_to_e_models,], aes(x=PC1, y=PC2, color=ctrl_init_clusts[m_to_e_models])) + 
  geom_point() + xlim(-5, 14) + ylim(-6, 7) + labs(color="Cluster")
ggplot(pca_df[m_to_e_models_other_states,], aes(x=PC1, y=PC2)) + 
  geom_point() + xlim(-5, 14) + ylim(-6, 7) + 
  geom_point(data=ctrl_init_pca[m_to_e_models,], color="red")


ctrl_results_df_fname <- file.path(ctrl_data_dir,
                                   paste0("ctrl_summary_trials=",ctrl_trials_per_noise,
                                          "_noise=",paste0(ctrl_noise_levels, collapse = ","),
                                          "_tcorr=",ctrl_tcorr))
if(!file.exists(ctrl_results_df_fname)) {
  
  ctrl_results_df <- data.frame(Noise=rep(ctrl_noise_levels, each=ctrl_trials_per_noise),
                                Trial=rep(1:ctrl_trials_per_noise, length(ctrl_noise_levels)),
                                Init_E=length(which(ctrl_init_clusts==1)),
                                Final_E=NA,
                                Init_M=length(which(ctrl_init_clusts==2)),
                                Final_M=NA,
                                NumConverting=NA,
                                NumRebellious=NA,
                                PctConverting=NA,
                                PctRebellious=NA)
  
  for(ctrl_current_noise in ctrl_noise_levels) {
    print(paste0("Beginning control trials on noise level ", ctrl_current_noise,"..."))
    for(ctrl_current_trial in 1:ctrl_trials_per_noise) {
      print(paste0("Beginning trial  ", ctrl_current_trial,"..."))
      
      # Set up placeholder racipe object
      racipe_ctrl <- racipe_ctrl_placeholder
      
      # Set ICs to WT steady states
      sracipeIC(racipe_ctrl) <- sracipeIC(racipe_ctrl_placeholder)[,]
      
      # Simulate with noise for a duration
      fname_racipe_ctrl <- file.path(ctrl_data_dir,
                                     paste0("racipe_ctrl_noise=",ctrl_current_noise,
                                            "_trial=",ctrl_current_trial,".Rds"))
      if(!file.exists(fname_racipe_ctrl)) {
        racipe_ctrl <- sracipeSimulate(racipe_ctrl, genIC = F, genParams = F, simulationTime = ctrl_simTime,
                                       integrate = T,
                                       initialNoise=ctrl_current_noise, nNoise=1, scaledNoise = T,
                                       integrateStepSize = 0.2, simDet=F, stepper = "EM_OU", ouNoise_t=ctrl_tcorr,
                                       nCores = 1)
        saveRDS(racipe_ctrl, fname_racipe_ctrl)
        
      } else {
        racipe_ctrl <- readRDS(fname_racipe_ctrl)
      }
      
      # Remove noise and simulate relaxation phase
      racipe_ctrl_relax_fname <- file.path(ctrl_data_dir,
                                           paste0("racipe_ctrl_noise=",ctrl_current_noise,
                                                  "_trial=",ctrl_current_trial,"_relaxed.Rds"))
      if(!file.exists(racipe_ctrl_relax_fname)) {
        # Remove signal
        racipe_ctrl_relax <- racipe_ctrl
        sracipeParams(racipe_ctrl_relax) <- sracipeParams(racipe_ctrl)
        # Set ICs to perturbed steady states
        sracipeIC(racipe_ctrl_relax) <- assay(racipe_ctrl)[,]  #, as.character(ctrl_current_noise)
        
        
        racipe_ctrl_relax <- sracipeSimulate(racipe_ctrl_relax, genIC = F, genParams = F, simulationTime = ctrl_relaxTime,
                                             integrate=T,
                                             integrateStepSize = 0.2, nNoise = 0, iNoise = 0, nCores = 1, simDet = T)
        saveRDS(racipe_ctrl_relax, racipe_ctrl_relax_fname)
        
        
      } else {
        racipe_ctrl_relax <- readRDS(racipe_ctrl_relax_fname)
      }
      
      # Rescale data & assign clusters to perturbed data
      racipe_ctrl_final <- as.data.frame(t(assay(racipe_ctrl_relax)))
      nCols <- ncol(racipe_ctrl_final)
      racipe_ctrl_final[,1:nCols] <- log2(1+racipe_ctrl_final[,1:nCols]) # Log transform
      racipe_ctrl_final[,1:nCols] <- sweep(racipe_ctrl_final[,1:nCols], 2, tmpMeans, FUN = "-") # scale
      racipe_ctrl_final[,1:nCols] <- sweep(racipe_ctrl_final[,1:nCols], 2, tmpSds, FUN = "/") # scale
      newpca <- scale(racipe_ctrl_final, pca$center, pca$scale) %*% pca$rotation
      
      # assign clusters to final states
      #newlabels <- knn_classifier(newpca[,1:num_pcs_clustering], pca_df[,1:num_pcs_clustering], clust_full, k=25)
      newlabels = predict(gmm, newdata = newpca[,1:num_pcs_clustering])
      
      
      # Compute & store # of transitions
      numTransitions <- length(which(newlabels == tgtClust & ctrl_init_clusts == initClust))
      numRebellious <- length(which(newlabels == initClust & ctrl_init_clusts == tgtClust))
      
      
      ctrl_results_df[which(ctrl_results_df$Noise == ctrl_current_noise &
                              ctrl_results_df$Trial == ctrl_current_trial), "Final_E"] <- length(which(newlabels == 1))
      ctrl_results_df[which(ctrl_results_df$Noise == ctrl_current_noise &
                              ctrl_results_df$Trial == ctrl_current_trial), "Final_M"] <- length(which(newlabels == 2))
      ctrl_results_df[which(ctrl_results_df$Noise == ctrl_current_noise &
                              ctrl_results_df$Trial == ctrl_current_trial), "NumConverting"] <- numTransitions
      ctrl_results_df[which(ctrl_results_df$Noise == ctrl_current_noise &
                              ctrl_results_df$Trial == ctrl_current_trial), "NumRebellious"] <- numRebellious
      ctrl_results_df[which(ctrl_results_df$Noise == ctrl_current_noise &
                              ctrl_results_df$Trial == ctrl_current_trial), "PctConverting"] <- numTransitions / length(which(ctrl_init_clusts == 1))
      ctrl_results_df[which(ctrl_results_df$Noise == ctrl_current_noise &
                              ctrl_results_df$Trial == ctrl_current_trial), "PctRebellious"] <- numRebellious / length(which(ctrl_init_clusts == 2))
      
    }
  }
  
  saveRDS(ctrl_results_df, ctrl_results_df_fname)
  
} else {
  ctrl_results_df <- readRDS(ctrl_results_df_fname)
}



########## SIMULATIONS WITH CLAMPING ############

expName <- paste0("ffctopo_t=",signal_simTime,"_relax_OUnoise=",paste0(signal_noise, collapse = "."),
                  "_tau=",signal_tcorr,"_genes=",paste0(signal_nGenes,collapse = "."),"_CLAMPS_2025")


undebug(optimizeST_Clamp)
paramSets <- optimizeST_Clamp(racipe_signaling, # NON NORMALIZED
                             pca_signaling, # same size as racipe
                             clust_signaling, # vector of length nrow(pca$x)
                             initialClust = 1, # should correspond to clust
                             targetClust = 2, # should correspond to clust
                             nSigGenes = signal_nGenes, # single value or vector of choices
                             clamp_df = clamp_df,
                             outDir = file.path(topoDir),
                             expName = expName,
                             plot=F,
                             noise = signal_noise, # single value or vector of choices
                             forceRerun = F, # whether to rerun signal simulations
                             forceRecompute = F, # whether to recompute final scores
                             checkpointSize=5, # how often to report progress to user
                             totalPerturbations = 500, # Only considered if randomParams is TRUE
                             nPerturbations=10,
                             relax=T, # whether to continue simulations for some time without a signal
                             simTime = signal_simTime,
                             simTimeRelax = signal_relaxTime,
                             onlyParams = F,
                             noise_tcorr = signal_tcorr)

debug(optimizeST_parallel)
paramSets <- optimizeST_parallel(racipe_signaling,
                                 pca_signaling,
                                 clust_signaling,
                                 totalPerturbations = 500,
                                 initialClust=1,
                                 targetClust=2,
                                 nSigGenes=signal_nGenes,
                                 outDir = file.path(topoDir),
                                 expName = expName,
                                 paramType="G",
                                 nPerturbations = 200,
                                 plot=F,
                                 fcVals = 1,
                                 noise = signal_noise,
                                 checkpointSize = 5,
                                 forceRerun = F,
                                 forceRecompute=F,
                                 anneal = F,
                                 randomParams = F,
                                 relax=T,
                                 simTime = signal_simTime,
                                 simTimeRelax = signal_relaxTime,
                                 nCores = 8,
                                 saveTrajectory = F,
                                 printInterval = 10,
                                 noise_tcorr = signal_tcorr,
                                 clamp=T,
                                 clamp_df=clamp_df,
                                 tmpMeans=tmpMeans,
                                 tmpSds=tmpSds)



########## SANITY CHECK ############

# "racipe_SNAI1_ERK_noise=0_relaxed.Rds" - top by outdegree
# "racipe_LIF_noise=0_relaxed.Rds" - input node, should have no effect
sig1_rs <- readRDS(file.path(topoDir, expName, "racipe_JAK_noise=0.04_relaxed.Rds"))
sig1_rs_clamped <- readRDS(file.path(topoDir, expName, "racipe_JAK_noise=0.04.Rds"))
#sig2_rs <- readRDS(file.path(topoDir, expName, "racipe_SNAI1_ERK_noise=0.04_relaxed.Rds"))


sig_use <- sig1_rs
sig_init <- t(sracipeIC(sig1_rs_clamped))
sig_init_norm <- sig_init
sig_init_norm[,1:nCols] <- log2(1+sig_init_norm[,1:nCols]) # Log transform
sig_init_norm[,1:nCols] <- sweep(sig_init_norm[,1:nCols], 2, tmpMeans, FUN = "-") # scale
sig_init_norm[,1:nCols] <- sweep(sig_init_norm[,1:nCols], 2, tmpSds, FUN = "/") # scale
sig_init_pca <- scale(sig_init_norm, pca$center, pca$scale) %*% pca$rotation 
sig_init_labels <- knn_classifier(sig_init_pca[,1:num_pcs_clustering], pca$x[,1:num_pcs_clustering], clust_full, k=25)


sig_clamped <- t(assay(sig1_rs_clamped))
sig_clamped_norm <- sig_clamped
sig_clamped_norm[,1:nCols] <- log2(1+sig_clamped_norm[,1:nCols]) # Log transform
sig_clamped_norm[,1:nCols] <- sweep(sig_clamped_norm[,1:nCols], 2, tmpMeans, FUN = "-") # scale
sig_clamped_norm[,1:nCols] <- sweep(sig_clamped_norm[,1:nCols], 2, tmpSds, FUN = "/") # scale
sig_clamped_pca <- scale(sig_clamped_norm, pca$center, pca$scale) %*% pca$rotation 
sig_clamped_labels <- knn_classifier(sig_clamped_pca[,1:num_pcs_clustering], pca$x[,1:num_pcs_clustering], clust_full, k=25)



sig_final <- t(assay(sig_use))
sig_final_norm <- sig_final
sig_final_norm[,1:nCols] <- log2(1+sig_final_norm[,1:nCols]) # Log transform
sig_final_norm[,1:nCols] <- sweep(sig_final_norm[,1:nCols], 2, tmpMeans, FUN = "-") # scale
sig_final_norm[,1:nCols] <- sweep(sig_final_norm[,1:nCols], 2, tmpSds, FUN = "/") # scale
sig_final_pca <- scale(sig_final_norm, pca$center, pca$scale) %*% pca$rotation 
sig_final_labels <- knn_classifier(sig_final_pca[,1:num_pcs_clustering], pca$x[,1:num_pcs_clustering], clust_full, k=25)

keep_states_m <- keep_states[which(!keep_states %in% keep_states_e)]
m_states <- pyracipe_raw_states[keep_states_m,genes]
m_states_norm <- m_states
m_states_norm[,1:nCols] <- log2(1+m_states_norm[,1:nCols]) # Log transform
m_states_norm[,1:nCols] <- sweep(m_states_norm[,1:nCols], 2, tmpMeans, FUN = "-") # scale
m_states_norm[,1:nCols] <- sweep(m_states_norm[,1:nCols], 2, tmpSds, FUN = "/") # scale
m_states_pca <- scale(m_states_norm, pca$center, pca$scale) %*% pca$rotation 



ggplot() +
  geom_point(data=sig_init_pca, aes(x=PC1, y=PC2), color="black", size=3) +
  geom_point(data=m_states_pca, aes(x=PC1, y=PC2), color="grey", size=3) +
  geom_point(data=sig_clamped_pca, aes(x=PC1, y=PC2), color="pink") +
  geom_point(data=sig_final_pca, aes(x=PC1, y=PC2), color="red", alpha=0.5) 
  
  

all.equal(unname(as.matrix(ics_signaling))[1:72,1:100], as.matrix(sracipeIC(sig1_rs_clamped))[1:72,1:100])

sig_init[1:5,20:25]
t(assay(sig1_rs_clamped))[1:5,20:25]
sig_final[1:5,20:25]

########## ASYMMETRY ANALYSIS - TOPOLOGY ############
# We observe that EMT seems to be easier to trigger than MET in the 26-node network
# Here, we check network topology metrics of E vs M genes to try to understand why



## Assign E and M genes
assign_cluster_genes <- function(expr, cluster_labels, p_thresh = 0.05) {
  group1 <- which(cluster_labels == 1)
  group2 <- which(cluster_labels == 2)
  
  res <- apply(expr, 1, function(x) {
    t_out <- t.test(x[group1], x[group2])
    c(p = t_out$p.value, mean_E = mean(x[group1]), mean_M = mean(x[group2]))
  })
  
  res <- as.data.frame(t(res))
  res$gene <- rownames(res)
  res <- res %>%
    mutate(cluster = case_when(
      p < p_thresh & mean_E > mean_M ~ "E",
      p < p_thresh & mean_E < mean_M ~ "M",
      TRUE ~ "Unassigned"
    ))
  
  return(res)
}

genes_assigned <- assign_cluster_genes(t(pyracipeNorm), clust_full)
table(genes_assigned$cluster)

## Prevalence of self-activation
count_self_regulation <- function(topo, cluster_genes, clust, reg_type) {
  sub_genes <- cluster_genes$gene[cluster_genes$cluster == clust]
  self_edges <- topo %>%
    filter(Source == Target, Type == reg_type, Source %in% sub_genes)
  return(nrow(self_edges))
}

# E and M self-activation
count_self_regulation(topo, genes_assigned, "E", 1)
count_self_regulation(topo, genes_assigned, "M", 1)

# E and M self-inhibition
count_self_regulation(topo, genes_assigned, "E", 2)
count_self_regulation(topo, genes_assigned, "M", 2)


## Centrality and out-degree
summarize_topology <- function(topo, cluster_genes, clust) {
  g <- graph_from_data_frame(topo, directed = TRUE)
  sub_genes <- cluster_genes$gene[cluster_genes$cluster == clust]
  
  bc <- betweenness(g, directed = TRUE, normalized = TRUE)
  od <- degree(g, mode = "out")
  
  data.frame(
    avg_out_degree = mean(od[sub_genes], na.rm = TRUE),
    avg_betweenness = mean(bc[sub_genes], na.rm = TRUE)
  )
}

# E and M betweenness centrality and out-degree
summarize_topology(topo, genes_assigned, "E")
summarize_topology(topo, genes_assigned, "M")

# Overall network betweenness & out-degree
g <- graph_from_data_frame(topo, directed = TRUE)
bc <- betweenness(g, directed = TRUE, normalized = TRUE)
od <- degree(g, mode = "out")

mean(bc)
mean(od)


## Interaction bias - how often do E/M genes activate/inhibit E/M genes?
count_directional_interactions <- function(topo, cluster_genes) {
  pheno_list <- split(cluster_genes$gene, cluster_genes$cluster)
  combos <- expand.grid(from = names(pheno_list), to = names(pheno_list), stringsAsFactors = FALSE)
  
  results <- lapply(1:nrow(combos), function(i) {
    from_genes <- pheno_list[[combos$from[i]]]
    to_genes <- pheno_list[[combos$to[i]]]
    
    subset <- topo %>%
      filter(Source %in% from_genes, Target %in% to_genes)
    
    tibble::tibble(
      from = combos$from[i],
      to = combos$to[i],
      activation = sum(subset$Type == 1),
      inhibition = sum(subset$Type == 2)
    )
  })
  
  do.call(rbind, results)
}


count_directional_interactions(topo, genes_assigned[which(genes_assigned$cluster %in% c("E","M")),])


## Network clustering coefficient - are E or M genes better connected?
calculate_clustering <- function(topo, phenotype_genes, clust) {
  g <- graph_from_data_frame(topo, directed = TRUE)
  sub_genes <- phenotype_genes$gene[phenotype_genes$cluster == clust]
  subg <- induced_subgraph(g, vids = sub_genes[sub_genes %in% V(g)$name])
  
  transitivity(subg, type = "global")
}

calculate_clustering(topo, genes_assigned, "E")
calculate_clustering(topo, genes_assigned, "M")


## Shortest paths - mean shortest path length between E and M genes
calculate_path_asymmetry <- function(topo, cluster_genes) {
  g <- graph_from_data_frame(topo, directed = TRUE)
  E_genes <- cluster_genes$gene[cluster_genes$cluster == "E"]
  M_genes <- cluster_genes$gene[cluster_genes$cluster == "M"]
  
  # Only consider genes present in the graph
  E_genes <- E_genes[E_genes %in% V(g)$name]
  M_genes <- M_genes[M_genes %in% V(g)$name]
  
  # Get shortest paths
  E_to_M <- unlist(lapply(E_genes, function(e) mean(suppressWarnings(distances(g, v = e, to = M_genes, mode = "out")), na.rm = TRUE)))
  M_to_E <- unlist(lapply(M_genes, function(m) mean(suppressWarnings(distances(g, v = m, to = E_genes, mode = "out")), na.rm = TRUE)))
  
  tibble::tibble(
    direction = c("E_to_M", "M_to_E"),
    mean_shortest_path = c(mean(E_to_M, na.rm = TRUE), mean(M_to_E, na.rm = TRUE))
  )
}


calculate_path_asymmetry(topo, genes_assigned)


########## ASYMMETRY ANALYSIS - DYNAMICS ############
## Some more simulations to understand EMT/MET symmetry
## All will be done on a given selected model
## First, simulate the model from many random ICs & estimate attractor boundaries
## Next, simulate the model with noise for a long duration, then plot trajectory as a heatmap
## Simulate model after many small perturbations from each state - estimate local stability 

# Using model 3, the first E/M bistable model in the list
# Also try model q, which is easily inducible to MET
# Also try model qq, which is easily inducible to EMT
asym_model_use <- 3 # Using 
asym_model_params <- pyracipe_params[asym_model_use,]
asym_model_states <- pyracipe_raw_states[which(pyracipe_metadata$MODEL_NO == asym_model_use),genes]


## Multiple ICs (estimate attractor boundary)
asym_multIC_fname <- file.path(dataDir, paste0("multIC_asym_model_",asym_model_use,".Rds"))
if(!file.exists(asym_multIC_fname)) {
  asym_multIC_placeholder <- sracipeSimulate(topo, numModels = 1, nIC = 500, integrate = F,
                                             genIC = T, genParams = T, simDet = T,
                                             initialNoise = 0, stepper = "EM_OU", ouNoise_t = 10,
                                             integrateStepSize = 0.2)
  sracipeParams(asym_multIC_placeholder) <- asym_model_params
  # sample ICs from random states?
  #sracipeIC(asym_multIC_placeholder) <- t(pyracipe_raw_states[sample(1:nrow(pyracipe_raw_states), 500),genes])
  asym_multIC <- sracipeSimulate(asym_multIC_placeholder, simulationTime = 100, integrate = T, 
                                 genParams = F, genIC = F, numModels = 500, nIC = 1)
  saveRDS(asym_multIC, asym_multIC_fname)
} else {
  asym_multIC <- readRDS(asym_multIC_fname)
}

## Long duration with noise
asym_noise_levels <- c(0.02, 0.04, 0.08, 0.1, 0.2)
for(asym_noise_level in asym_noise_levels) {
  asym_longStoch_fname <- file.path(dataDir, 
                                    paste0("longStoch_asym_model_",asym_model_use,
                                           "_noise=",asym_noise_level,".Rds"))
  if(!file.exists(asym_longStoch_fname)) {
    asym_longStoch_placeholder <- sracipeSimulate(topo, numModels = 1, nIC = 1, integrate = F,
                                                  genIC = T, genParams = T, simDet = T,
                                                  initialNoise = 0.1, stepper = "EM_OU", ouNoise_t = 10,
                                                  timeSeries = T, integrateStepSize = 0.2, printInterval = 0.2)
    sracipeParams(asym_longStoch_placeholder) <- asym_model_params
    
    system.time({
      asym_longStoch <- sracipeSimulate(asym_longStoch_placeholder, simulationTime = 10000, integrate = T, 
                                        genParams = F, genIC = F, numModels = 1, nIC = 1,
                                        simDet = F,
                                        timeSeries = T, integrateStepSize = 0.2, printInterval = 2, printStart = 0.2,
                                        initialNoise = asym_noise_level, stepper = "EM_OU", ouNoise_t = 10)
    })
    
    saveRDS(asym_longStoch, asym_longStoch_fname)
  } else {
    asym_longStoch <- readRDS(asym_longStoch_fname)
  }
}



## Small perturbation responses
asym_num_perturbations <- 500
asym_perturb_state_idx <- c()
asym_perturb_state <- c() # select state to begin from
asym_perturb_placeholder <- sracipeSimulate(topo, numModels = 1, nIC = 1, integrate = F,
                                            genIC = T, genParams = T, simDet = T,
                                            initialNoise = 0.1, stepper = "EM_OU", ouNoise_t = 10)
sracipeParams(asym_multIC_placeholder) <- asym_model_params

asym_perturb_summary_fname <- file.path(dataDir, 
                                  paste0("perturb_asym_model_",asym_model_use,
                                         "_noise=",asym_noise_level,
                                         "_fromStateIndex=",asym_perturb_state_idx,".Rds"))
if(!file.exists(asym_perturb_summary_fname)) {
  for(asym_perturb_num in 1:asym_num_perturbations) {
    # Generate small perturbation to initial state
    
    # Set IC
    
    # Simulate
    
    # Analyze results & add to summary
  }
  
  # Save summary table
  saveRDS(asym_perturb_summary, asym_perturb_summary_fname)
  
}





######## AGGREGATE ANALYSIS #####
resultSet_fname <- file.path(topoDir,expName,"result_summary.Rds")
genes_x_transitions_matrix_fname <- file.path(topoDir, expName, "genes_x_transitions_matrix.Rds")
resultSet_full <- readRDS(resultSet_fname)
if("Species 2" %in% colnames(resultSet_full)) {
  resultSet_full$SigName_Short <- paste0(resultSet_full[,c("Species 1")], "_",
                                         resultSet_full[,c("Species 2")])  
}

resultSet <- resultSet_full
rs_det <- resultSet_full[which(resultSet_full$Noise == 0),]

if(!"ConversionPct" %in% colnames(resultSet_full)) {
  rs_full_list <- genes_x_transitions(resultSet_full, # dataframe w/ columns: ModelNo, SetName
                                      topoName = topoName, # string
                                      collectionName = expName, # string
                                      initialClust = 1, # int
                                      targetClust = 2, # int
                                      wt_data = pyracipeNorm[,genes], # matrix/dataframe of original data
                                      clust = clust_signaling, # vector of length numSamples containing integers
                                      clust_all = clust_full, # full cluster labels matching wt_data
                                      tmpMeans = tmpMeans,
                                      tmpSds = tmpSds
  )
  
  rsMatrix_full <- rs_full_list[[1]]
  rownames(rsMatrix_full) <- resultSet$SetName[which(!is.na(resultSet$ConversionPct))]
  
  resultSet_full <- rs_full_list[[2]]
  resultSet_full$ConversionPct <- resultSet_full$ConvertingModels / resultSet_full$startingInitPopulation
  
  saveRDS(resultSet_full, resultSet_fname)
  saveRDS(rsMatrix_full, genes_x_transitions_matrix_fname)
} 


## Basic stats
# overall conversion pcts (deterministic)
summary(resultSet_full[which(resultSet_full$Noise == 0),"ConvertingModels"])
summary(resultSet_full[which(resultSet_full$Noise == 0),"ConversionPct"])
summary(resultSet_full[which(resultSet_full$Noise == 0 & resultSet_full$NumGenes == 1),"ConversionPct"])

# signals above 50% eff
length(which(resultSet_full[which(resultSet_full$Noise == 0),"ConversionPct"] >= 0.5))
summary(resultSet_full[which(resultSet_full$Noise == 0 & 
                               resultSet_full$ConversionPct > 0.5),"ConversionPct"])
resultSet_full[which(resultSet_full$Noise == 0 & resultSet_full$ConversionPct >= 0.5),"SetName"]


# signals above 5% eff
length(which(resultSet_full[which(resultSet_full$Noise == 0),"ConversionPct"] >= 0.05))

# 1-gene best signals
View(resultSet_full[which(resultSet_full$Noise == 0 & resultSet_full$NumGenes == 1),
                    c("SetName", "ConversionPct")])

# 2-gene best signals
View(resultSet_full[which(resultSet_full$Noise == 0 & resultSet_full$NumGenes == 2),
                    c("SetName", "ConversionPct")])

## 0.04 noise level:
# overall conversion pcts
summary(resultSet_full[which(resultSet_full$Noise == 0.04),"ConvertingModels"])
summary(resultSet_full[which(resultSet_full$Noise == 0.04),"ConversionPct"])
summary(resultSet_full[which(resultSet_full$Noise == 0.04 & resultSet_full$NumGenes == 1),"ConversionPct"])

# signals above 50% eff
length(which(resultSet_full[which(resultSet_full$Noise == 0.04),"ConversionPct"] >= 0.5))
resultSet_full[which(resultSet_full$Noise == 0.04 & resultSet_full$ConversionPct > 0.5),"SetName"]

# in earlier workflows, I subset the data here - now it just takes all noise levels by default
selectedNoise <- signal_noise
resultSet <- resultSet_full[which(resultSet_full$Noise %in% selectedNoise),]



######### Signal Characteristics #########
# Import network into igraph
g <- igraph::graph_from_edgelist(as.matrix(topo[,c("Source","Target")]))

w <- as.matrix(as_adjacency_matrix(g))
gene_id_list <- rownames(w)
colnames(w) <- 1:ncol(w)
rownames(w) <- 1:nrow(w)

if(!"GroupBetweenCentrality" %in% colnames(resultSet)) {
  for(i in rownames(resultSet)) {
    
    # Identify nodes involved in signal
    genesInSignal <-  c(resultSet[i,"Species 1"], resultSet[i,"Species 2"])
    genesInSignal <- genesInSignal[which(!is.na(genesInSignal))]
    
    genesInSignalIDs <- which(gene_id_list %in% genesInSignal)
    
    # Compute group betweenness centrality
    bet_cent <- kpcent(w, genesInSignalIDs, type="betweenness")
    
    
    # Compute group closeness centrality
    close_cent <- kpcent(w, genesInSignalIDs, type="closeness")
    
    
    # Add to resultSet
    resultSet[i,"GroupBetweenCentrality"] <- bet_cent
    resultSet[i,"GroupClosenessCentrality"] <- close_cent
    
  }
  
  saveRDS(resultSet, resultSet_fname)
} 






######## EFFICACY ANALYSIS ############
#debug(cells_x_signals)
cell_signal_df_fname <- file.path(topoDir,expName,
                                  paste0("cell_signal_df_noise=", 
                                         paste0(selectedNoise, collapse = ","),".Rds")) 
if(!file.exists(cell_signal_df_fname)) {
  undebug(cells_x_signals)
  cell_signal_df <- cells_x_signals(paramSets = resultSet, # dataframe w/ columns: ModelNo, SetName
                                    topoName = topoName, # string
                                    collectionName = expName, # string
                                    initialClust = 1, # int
                                    targetClust = 2, # int
                                    wt_data = pyracipeNorm[,1:nGenes], # matrix/dataframe of size numSamples x numFeatures
                                    clust = clust_signaling, # vector of length numSamples containing integers
                                    clust_all = clust_full,
                                    tmpMeans = tmpMeans,
                                    tmpSds = tmpSds,
                                    InitRawStates = pyracipe_raw_states[keep_states_e,genes],
                                    useDiffs = F)
  saveRDS(cell_signal_df, cell_signal_df_fname)
  
  
} else {
  cell_signal_df <- readRDS(cell_signal_df_fname)
  #cell_signal_df <- cell_signal_df[,which(resultSet_full$Noise %in% selectedNoise)]
}








