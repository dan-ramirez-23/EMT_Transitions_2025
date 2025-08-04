########## SETUP ############

rm(list=ls())
library(sRACIPE)#, lib.loc = "/Users/danramirez/localR/4.2.2-arm")
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
library(ComplexHeatmap)
library(igraph)
library(keyplayer)
library(mclust)
library(cluster)
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

########## UNTREATED CONDITION ############

numModels <- 1000
numICs <- 100
simTime <- 200


# simulate WT with RACIPE
racipe_fname <- file.path(topoDir, "racipe_100IC.Rds")
if(!file.exists(racipe_fname)) {
  racipe <- sRACIPE::sracipeSimulate(circuit = topo, numModels = numModels,
                                     nIC = numICs, simulationTime = simTime, initialNoise = 0, nNoise = 0)

  saveRDS(racipe, racipe_fname)
} else {
  racipe <- readRDS(racipe_fname)
  racipe_old <- racipe
}

# continuation to reach fuller convergence
racipe_cont_fname <- file.path(topoDir, "racipe_100IC_continued.Rds")
if(!file.exists(racipe_cont_fname)) {
  racipe <- readRDS(racipe_fname)
  racipe_cont <- racipe
  sracipeIC(racipe_cont) <- assay(racipe)
  racipe_cont@metadata$config$simParams[["nIC"]] <- 1
  racipe_cont <- sracipeSimulate(racipe_cont, nIC = 1, numModels = 100000,
                                 genIC = F, genParams = F,
                                 simDet = T, initialNoise = 0, simulationTime = 200,
                                 nCores = 1, integrateStepSize = 0.2, stepper = "EM_OU")
  
  saveRDS(racipe_cont, racipe_cont_fname)
  racipe <- racipe_cont
} else {
  racipe <- readRDS(racipe_cont_fname)
}


# get moments for normalizing other cases later
genes <- rownames(racipe)
unnormData <- t(assay(racipe))
simExp <- assay(racipe, 1)
simExp <- log2(1+simExp)
tmpMeans <- rowMeans(simExp)
tmpSds <- apply(simExp,1,sd)

racipeNorm <- sracipeNormalize(racipe)
racipeData <- as.data.frame(t(assay(racipeNorm)))
exprMat_norm <- racipeData


########## IDENTIFY UNIQUE STATES ############
## For each model, id unique states (doesn't need much precision)
ss_unique_fname <- file.path(dataDir,"ss_unique_df.Rds")

if(!file.exists(ss_unique_fname)) {
  # Find unique steady states per model
  ss_rounded <- round(as.data.frame(unnormData), 1)
  ss_rounded$Model <- rep(1:numModels, numICs)
  unique_state_idx <- which(!duplicated(ss_rounded[,]))
  unique_state_models <- ss_rounded$Model[unique_state_idx]
  ss_unique <- ss_rounded[,] %>%
    distinct()
  all.equal(ss_unique$Model, unique_state_models)
  ss_unique$StateIndex <- unique_state_idx


  saveRDS(ss_unique, ss_unique_fname)
} else {
  ss_unique <- readRDS(ss_unique_fname)
}




########## PCA ON UNIQUE STATES ############
pca_fname <- file.path(dataDir,"pca_2025.Rds")
if(!file.exists(pca_fname)) {
  # PCA on full data
  pca <- prcomp(exprMat_norm[ss_unique$StateIndex,1:nGenes])
  pca_df_full <- as.data.frame(pca$x)

  # We will flip PC2 for consistency with other visualizations
  #pca_df_full$PC2 <- -1 * pca_df_full$PC2
  #pca$rotation[,2] <- -1 * pca$rotation[,2]
  #pca$x[,2] <- -1 * pca$x[,2]

  # sanity check: plot colored by Cdh1
  ggplot(pca_df_full, aes(x=PC1, y=PC2, color=ss_unique$Cdh1)) + geom_point()

  saveRDS(pca, pca_fname)
} else {
  pca <- readRDS(pca_fname)
  pca_df_full <- as.data.frame(pca$x)
}


########## CLUSTERING ############
sil_df_fname <- file.path(dataDir, "silhouette_df.Rds")

# Helper function to select k for clustering unique states
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

# Identify optimal k
if(!file.exists(sil_df_fname)) {
  k_range <- 2:8

  sil_df <- calc_silhouette(exprMat_norm[ss_unique$StateIndex,1:nGenes], k_range)
  saveRDS(sil_df, sil_df_fname)
}



clust_all_fname <- file.path(dataDir,"clust_all_2025.Rds")
num_pcs_clustering <- 15 # First 15 PCs cover ~93% of variance
if(!file.exists(clust_all_fname)) {
  # fit GMM model
  gmm = GMM(pca_df_full[,1:num_pcs_clustering], 2, dist_mode = "eucl_dist", seed_mode = "random_subset", km_iter = 10,
            em_iter = 10, verbose = F)

  # predict centroids, covariance matrix and weights
  clust_full = predict(gmm, newdata = pca_df_full[,1:num_pcs_clustering])

  # Im reversing it so the initial (left-side) state is cluster 1
  #revClust <- clust_full
  #revClust <- ifelse(clust_full == 1, 2, ifelse(clust_full == 2, 1, clust_full))
  #clust_full <- revClust

  ggplot(pca_df_full, aes(x=PC1, y=PC2, color=as.factor(clust_full))) + geom_point()
  ggplot(pca_df_full, aes(x=PC1, y=PC2, color=exprMat_norm[ss_unique$StateIndex,"Cdh1"])) + geom_point()

  saveRDS(clust_full, file = clust_all_fname)
  ss_unique$Cluster <- clust_full
  saveRDS(ss_unique, ss_unique_fname)
} else {
  clust_full <- readRDS(clust_all_fname)
}


########## SELECT E/M BISTABLE MODELS ############
racipe_bistable_fname <- file.path(dataDir,"racipe_bistable_2025.Rds")
racipe_bistable_indices_fname <- file.path(dataDir,"racipe_bistable_indices_2025.Rds")
unique_state_idx_list_fname <- file.path(dataDir,"unique_state_bistable_e_indices_2025")
unique_state_idx_list_all_fname <- file.path(dataDir,"unique_state_bistable_indices_all_2025")
summary_df_fname <- file.path(dataDir,"state_summary_df.Rds")
models_selected_fname <- file.path(dataDir,"racipe_bistable_indexMap_2025.Rds")
if(!file.exists(racipe_bistable_fname)) {
  # update state summary df
  summary_df <- ss_unique %>%
    select(all_of(c("Model", "Cluster"))) %>%
    group_by(Model) %>%
    summarise(
      NumStates = n(),
      Stability = case_when(
        all(Cluster == 1) ~ '1',
        all(Cluster == 2) ~ '2',
        any(Cluster == 1) & any(Cluster == 2) ~ 'bistable'
      )
    )

  # filter for models with <10 states, only bistable
  # save new racipe object w/ only bistable
  models_selected <- unlist(summary_df[which(summary_df$NumStates == 2 &
                                               summary_df$Stability == "bistable"),"Model"])[1:500]
  keepIdx <- c()
  unique_state_idx_list <- c()
  unique_state_idx_list_all <- c()
  for(model in models_selected) {
    # Select only epithelial state, and first epithelial state if multiple
    modelStates <- ss_unique[which(ss_unique$Model == model & ss_unique$Cluster == 1),]
    model_idx_list <- seq(model, numModels*numICs, numModels)
    addIdx <- sample(which(row.match(
      as.data.frame(t(round(assay(racipe), 1)))[model_idx_list,],
      modelStates[,genes],
      nomatch = NA) == 1), 1)
    addIdx <- model_idx_list[addIdx] # bring back index for original racipe object
    keepIdx <- c(keepIdx, addIdx)

    unique_state_idx <- which(ss_unique$Model == model & ss_unique$Cluster == 1)
    unique_state_idx_list <- c(unique_state_idx_list, unique_state_idx)
    unique_state_idx_all <- which(ss_unique$Model == model)
    unique_state_idx_list_all <- c(unique_state_idx_list_all, unique_state_idx_all)
  }

  # subset racipe object for bistable models, update parameters
  racipe_bistable <- racipe[,keepIdx]
  sracipeParams(racipe_bistable) <- sracipeParams(racipe)[as.numeric(unname(models_selected)),] # parameter numbering is messssssed up in sRACIPE, not my fault
  racipe_bistable@metadata$config$simParams[["numModels"]] <- length(keepIdx)

  saveRDS(racipe_bistable, racipe_bistable_fname)
  saveRDS(models_selected, models_selected_fname)
  saveRDS(keepIdx, racipe_bistable_indices_fname)
  saveRDS(unique_state_idx_list, unique_state_idx_list_fname)
  saveRDS(unique_state_idx_list_all, unique_state_idx_list_all_fname)
  saveRDS(summary_df, summary_df_fname)
  racipe_bistable_indices <- keepIdx
  pca_df <- as.data.frame(pca$x)[as.character(unique_state_idx_list),]
  clust <- clust_full[unique_state_idx_list]
} else {
  racipe_bistable <- readRDS(racipe_bistable_fname)
  models_selected <- readRDS(models_selected_fname)
  racipe_bistable_indices <- readRDS(racipe_bistable_indices_fname)
  unique_state_idx_list <- readRDS(unique_state_idx_list_fname)
  unique_state_idx_list_all <- readRDS(unique_state_idx_list_all_fname)
  summary_df <- readRDS(summary_df_fname)
  pca_df <- as.data.frame(pca$x)[as.character(unique_state_idx_list),]
  clust <- clust_full[unique_state_idx_list]
}

racipe_bistable_raw <- racipe_bistable
racipe_bistable_raw@metadata$config$simParams["nIC"] <- 1
unnormData <- t(assay(racipe_bistable_raw))
racipe_bistable <- sracipeNormalize(racipe_bistable)
exprMat <- as.data.frame(t(assay(racipe_bistable)))
exprMat$Cluster <- clust
exprMat_norm <- as.data.frame(t(assay(racipeNorm)))

########## IDENTIFY CLAMP VALUES ############
## Target: dataframe with model id, gene, expression, cluster (very long, maybe a wider format)
clamp_df_fname <- file.path(dataDir,"clamp_values_2025.Rds")
clamp_df_full_fname <- file.path(dataDir,"clamp_values_all_2025.Rds")
if(!file.exists(clamp_df_fname)) {
  keepIdx <- c()
  for(model in models_selected) {
    # add steady states for cluster 1 and 2
    addIdx <- which(ss_unique$Model == model)
    keepIdx <- c(keepIdx, addIdx)

  }
  clamp_df <- pivot_longer(ss_unique[keepIdx,], cols = all_of(genes),
                           names_to = "Gene", values_to = "Expression")
  clamp_df$ModelIndex <- as.numeric(factor(clamp_df$Model))
  saveRDS(clamp_df, clamp_df_fname)

  clamp_df_full <- pivot_longer(ss_unique[,], cols = all_of(genes),
                                names_to = "Gene", values_to = "Expression")
  clamp_df_full$ModelIndex <- as.numeric(factor(clamp_df_full$Model))
  saveRDS(clamp_df_full, clamp_df_full_fname)

} else {
  clamp_df <- readRDS(clamp_df_fname)
  clamp_df_full <- readRDS(clamp_df_full_fname)
}

######## NOISE-ONLY CONTROL SIMULATIONS #####
## Simulate trials with various noise levels & no (effective) signal
## Track number of transitions over time

ctrl_simTime <- 200
ctrl_relaxTime <- 50
ctrl_tcorr <- 10
initClust <- 1
tgtClust <- 2



ctrl_data_dir <- file.path(dataDir, "noise_only_controls")
if(!dir.exists(ctrl_data_dir)) {
  dir.create(ctrl_data_dir)
}


ctrl_noise_levels <- c(0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.5, 1, 2)
ctrl_trial_expr <- as.data.frame(t(assay(racipe_bistable)))
ctrl_trials_per_noise <- 10

ctrl_indices <- c()
for(model in models_selected) {
  clust_pick <- sample(c(1,2),1)
  idx_add <- ss_unique[which(ss_unique$Model == model), "StateIndex"]
  ctrl_indices <- c(ctrl_indices, idx_add)
}

# create RACIPE object with equal split of E and M steady states to begin
racipe_ctrl_placeholder <- racipe[, ctrl_indices]
sracipeParams(racipe_ctrl_placeholder) <- sracipeParams(racipe)[rep(as.numeric(unname(models_selected)), each=2),]
sracipeIC(racipe_ctrl_placeholder) <- assay(racipe_ctrl_placeholder)[,]
racipe_ctrl_placeholder@metadata$config$simParams[["numModels"]] <- length(ctrl_indices)
racipe_ctrl_placeholder@metadata$config$simParams["nIC"] <- 1
ctrl_init_data_norm <- as.data.frame(t(assay(racipe_ctrl_placeholder)))
ctrl_init_data_norm[,genes] <- log2(1+ctrl_init_data_norm[,genes]) # Log transform
ctrl_init_data_norm[,genes] <- sweep(ctrl_init_data_norm[,genes], 2, tmpMeans, FUN = "-") # scale
ctrl_init_data_norm[,genes] <- sweep(ctrl_init_data_norm[,genes], 2, tmpSds, FUN = "/") # scale
ctrl_init_pca <- scale(ctrl_init_data_norm, pca$center, pca$scale) %*% pca$rotation
ctrl_init_clusts <- knn_classifier(ctrl_init_pca[,1:num_pcs_clustering], pca_df_full[,1:num_pcs_clustering], clust_full, k=25)




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
    for(ctrl_current_trial in 1:ctrl_trials_per_noise) {

      # Set up placeholder racipe object
      racipe_ctrl <- racipe_ctrl_placeholder

      # Set ICs to WT steady states
      sracipeIC(racipe_ctrl) <- assay(racipe_ctrl_placeholder)[,]

      # Simulate with noise for a duration
      fname_racipe_ctrl <- file.path(ctrl_data_dir,
                                     paste0("racipe_ctrl_noise=",ctrl_current_noise,
                                            "_trial=",ctrl_current_trial,".Rds"))
      if(!file.exists(fname_racipe_ctrl)) {
        racipe_ctrl <- sracipeSimulate(racipe_ctrl, genIC = F, genParams = F, simulationTime = ctrl_simTime,
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
        sracipeIC(racipe_ctrl_relax) <- assay(racipe_ctrl, as.character(ctrl_current_noise))[,]


        racipe_ctrl_relax <- sracipeSimulate(racipe_ctrl_relax, genIC = F, genParams = F, simulationTime = ctrl_relaxTime,
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
      newlabels <- knn_classifier(newpca[,1:num_pcs_clustering], pca_df_full[,1:num_pcs_clustering], clust_full, k=25)


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


########## ASYMMETRY ANALYSIS - DYNAMICS ############
## Some more simulations to understand EMT/MET symmetry
## All will be done on a given selected model
## First, simulate the model from many random ICs & estimate attractor boundaries
## Next, simulate the model with noise for a long duration, then plot trajectory as a heatmap
## Simulate model after many small perturbations from each state - estimate local stability 

# Identify an EMT-inducible and MET-inducible model
ctrl_summary_df <- data.frame(Model=rep(models_selected, length(ctrl_noise_levels)), 
                              Noise=rep(ctrl_noise_levels, each=length(models_selected)), 
                              EMT_Count=NA, MET_Count=NA)
ctrl_summary_df_fname <- file.path(ctrl_data_dir,"ctrl_summary_df_modelwise.Rds")
if(!file.exists(ctrl_summary_df_fname)) {
  for(ctrl_current_noise in ctrl_noise_levels) {
    modelwise_met_counts <- rep(0, length(models_selected))
    modelwise_emt_counts <- rep(0, length(models_selected))
    
    for(ctrl_current_trial in 1:ctrl_trials_per_noise) { 
      racipe_ctrl_relax_fname <- file.path(ctrl_data_dir,
                                           paste0("racipe_ctrl_noise=",ctrl_current_noise,
                                                  "_trial=",ctrl_current_trial,"_relaxed.Rds"))
      racipe_ctrl_relax <- readRDS(racipe_ctrl_relax_fname)
      
      # Rescale data & assign clusters to perturbed data
      racipe_ctrl_final <- as.data.frame(t(assay(racipe_ctrl_relax)))
      nCols <- ncol(racipe_ctrl_final)
      racipe_ctrl_final[,1:nCols] <- log2(1+racipe_ctrl_final[,1:nCols]) # Log transform
      racipe_ctrl_final[,1:nCols] <- sweep(racipe_ctrl_final[,1:nCols], 2, tmpMeans, FUN = "-") # scale
      racipe_ctrl_final[,1:nCols] <- sweep(racipe_ctrl_final[,1:nCols], 2, tmpSds, FUN = "/") # scale
      newpca <- scale(racipe_ctrl_final, pca$center, pca$scale) %*% pca$rotation
      
      # assign clusters to final states
      newlabels <- knn_classifier(newpca[,1:num_pcs_clustering], pca_df_full[,1:num_pcs_clustering], clust_full, k=25)
      
      # Compute & store modelwise transitions
      old_mat <- matrix(ctrl_init_clusts, nrow = 2)
      new_mat <- matrix(newlabels, nrow = 2)
      
      # EMT: started in 1 (E), ended in 2 (M)
      EMT <- colSums(old_mat == 1 & new_mat == 2) > 0
      
      # MET: started in 2 (M), ended in 1 (E)
      MET <- colSums(old_mat == 2 & new_mat == 1) > 0
      
      # Convert to integer vectors (0 or 1)
      EMT <- as.integer(EMT)
      MET <- as.integer(MET)
      
      
      modelwise_emt_counts <- modelwise_emt_counts + EMT
      modelwise_met_counts <- modelwise_met_counts + MET
      
    }
    
    # Add results to summary df
    for(modelNo in seq_along(models_selected)) {
      ctrl_summary_model <- models_selected[modelNo]
      ctrl_summary_df[which(ctrl_summary_df$Model == ctrl_summary_model &
                              ctrl_summary_df$Noise == ctrl_current_noise),"EMT_Count"] <- modelwise_emt_counts[modelNo]
      ctrl_summary_df[which(ctrl_summary_df$Model == ctrl_summary_model &
                              ctrl_summary_df$Noise == ctrl_current_noise),"MET_Count"] <- modelwise_met_counts[modelNo]
    }
  }
  
  saveRDS(ctrl_summary_df, file = ctrl_summary_df_fname)
} else {
  ctrl_summary_df <- readRDS(ctrl_summary_df_fname)
}

# now we will summarize the bias per model
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

# try model 144, which is easily inducible to MET
# Also try model 589, which is easily inducible to EMT
# Maybe model 601, which is nearly perfectly neutral
asym_models_use <- c(144, 589, 601)

for(asym_model_use in asym_models_use) {
  asym_model_params <- sracipeParams(racipe)[asym_model_use,]
  asym_model_states <- ss_unique[which(ss_unique$Model == asym_model_use),genes]
  
  
  ## Multiple ICs (estimate attractor boundary)
  asym_multIC_fname <- file.path(dataDir, paste0("multIC_asym_model_",asym_model_use,".Rds"))
  if(!file.exists(asym_multIC_fname)) {
    print(paste0("Beginning multi-IC simulations for model ",asym_model_use))
    asym_multIC_placeholder <- sracipeSimulate(topo, numModels = 1, nIC = 500, integrate = F,
                                               genIC = T, genParams = T, simDet = T,
                                               initialNoise = 0, stepper = "EM_OU", ouNoise_t = 10,
                                               integrateStepSize = 0.2)
    sracipeParams(asym_multIC_placeholder) <- asym_model_params
    
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
      print(paste0("Beginning long noisy simulations for model ",asym_model_use," with noise level ",asym_noise_level))
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
  
  
}





########## SIMULATIONS WITH CLAMPING ############
signal_simTime <- 500
signal_relaxTime <- 50
signal_nGenes <- c(1,2)
signal_noise <- c(0, 0.04, 0.2)
signal_tcorr <- 10
initClust <- 1
tgtClust <- 2


expName <- paste0("bhtopo_t=",signal_simTime,"_relax_OUnoise=",paste0(signal_noise, collapse = "."),
                  "_tau=",signal_tcorr,"_genes=",paste0(signal_nGenes,collapse = "."),"_CLAMPS_2025")
resultSet_fname <- file.path(topoDir,expName,"result_summary.Rds")

# For signaling, use prcomp object but replace embeddings with subset of models
pca_st <- pca
pca_st$x <- pca_df

# can use single-threaded mode as well alternatively
# undebug(optimizeST_Clamp)
# paramSets <- optimizeST_Clamp(racipe_bistable_raw, # NON NORMALIZED
#            pca_st, # same size as racipe
#            clust, # vector of length nrow(pca$x)
#            initialClust = 1, # should correspond to clust
#            targetClust = 2, # should correspond to clust
#            nSigGenes = signal_nGenes, # single value or vector of choices
#            clamp_df = clamp_df,
#            outDir = file.path(topoDir),
#            expName = expName,
#            plot=F,
#            noise = signal_noise, # single value or vector of choices
#            forceRerun = F, # whether to rerun signal simulations
#            forceRecompute = F, # whether to recompute final scores
#            checkpointSize=25, # how often to report progress to user
#            totalPerturbations = 500, # Only considered if randomParams is TRUE
#            nPerturbations=30,
#            relax=T, # whether to continue simulations for some time without a signal
#            simTime = signal_simTime,
#            simTimeRelax = signal_relaxTime,
#            onlyParams = F,
#            noise_tcorr = signal_tcorr)

if(!file.exists(resultSet_fname)) {
  undebug(optimizeST_parallel)
  paramSets <- optimizeST_parallel(racipe_bistable_raw,
                                   pca_st,
                                   clust,
                                   totalPerturbations = 500,
                                   initialClust=1,
                                   targetClust=2,
                                   nSigGenes=signal_nGenes,
                                   outDir = file.path(topoDir),
                                   expName = expName,
                                   paramType="G",
                                   nPerturbations = 351,
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
                                   nCores = 6,
                                   saveTrajectory = F,
                                   printInterval = 10,
                                   noise_tcorr = signal_tcorr,
                                   clamp=T,
                                   clamp_df=clamp_df,
                                   tmpMeans=tmpMeans,
                                   tmpSds=tmpSds)
} else {
  resultSet_full <- readRDS(resultSet_fname)
}



######## AGGREGATE ANALYSIS #####
genes_x_transitions_matrix_fname <- file.path(topoDir, expName, "genes_x_transitions_matrix.Rds")
resultSet_full <- readRDS(resultSet_fname)
if("Species 2" %in% colnames(resultSet_full)) {
  resultSet_full$SigName_Short <- paste0(resultSet_full[,c("Species 1")], "_",
                                         resultSet_full[,c("Species 2")])
}

rs_det <- resultSet_full[which(resultSet_full$Noise == 0),]

if(!"ConversionPct" %in% colnames(resultSet_full)) {
  rs_full_list <- genes_x_transitions(resultSet_full, # dataframe w/ columns: ModelNo, SetName
                                      topoName = topoName, # string
                                      collectionName = expName, # string
                                      initialClust = 1, # int
                                      targetClust = 2, # int
                                      wt_data = exprMat_norm[ss_unique$StateIndex,genes], # matrix/dataframe of original data
                                      clust = clust, # vector of length numSamples containing integers
                                      clust_all = clust_full, # full cluster labels matching wt_data
                                      tmpMeans = tmpMeans,
                                      tmpSds = tmpSds
  )

  rsMatrix_full <- rs_full_list[[1]]
  rownames(rsMatrix_full) <- resultSet_full$SetName[which(!is.na(resultSet_full$ConversionPct))]

  resultSet_full <- rs_full_list[[2]]
  resultSet_full$ConversionPct <- resultSet_full$ConvertingModels / resultSet_full$startingInitPopulation

  saveRDS(resultSet_full, resultSet_fname)
  saveRDS(rsMatrix_full, genes_x_transitions_matrix_fname)
}


# 1-gene efficacy, deterministic and 0.04 noise
summary(resultSet_full[which(resultSet_full$NumGenes == 1 & resultSet_full$Noise == 0),"ConversionPct"])
summary(resultSet_full[which(resultSet_full$NumGenes == 1 & resultSet_full$Noise == 0.04),"ConversionPct"])
summary(resultSet_full[which(resultSet_full$NumGenes == 1 & resultSet_full$Noise == 0.2),"ConversionPct"])

# 1- and 2-gene efficacy
summary(resultSet_full[which(resultSet_full$Noise == 0),"ConversionPct"])
summary(resultSet_full[which(resultSet_full$Noise == 0.04),"ConversionPct"])
summary(resultSet_full[which(resultSet_full$Noise == 0.2),"ConversionPct"])

# Zeb1 efficacy
summary(resultSet_full[which(resultSet_full$Noise == 0 & 
                               (resultSet_full$`Species 1` == "Zeb1" | 
                                  resultSet_full$`Species 2` == "Zeb1")),"ConversionPct"])
summary(resultSet_full[which(resultSet_full$Noise == 0.04 & 
                               (resultSet_full$`Species 1` == "Zeb1" | 
                                  resultSet_full$`Species 2` == "Zeb1")),"ConversionPct"])
summary(resultSet_full[which(resultSet_full$Noise == 0.2 & 
                               (resultSet_full$`Species 1` == "Zeb1" | 
                                  resultSet_full$`Species 2` == "Zeb1")),"ConversionPct"])
# Snai1 efficacy
summary(resultSet_full[which(resultSet_full$Noise == 0.04 & 
                               (resultSet_full$`Species 1` == "Snai1" | 
                                  resultSet_full$`Species 2` == "Snai1")),"ConversionPct"])


######### Signal Characteristics #########
# Import network into igraph
g <- igraph::graph_from_edgelist(as.matrix(topo[,c("Source","Target")]))

w <- as.matrix(as_adjacency_matrix(g))
gene_id_list <- rownames(w)
colnames(w) <- 1:ncol(w)
rownames(w) <- 1:nrow(w)

if(!"GroupBetweenCentrality" %in% colnames(resultSet_full)) {
  for(i in rownames(resultSet_full)) {

    # Identify nodes involved in signal
    genesInSignal <-  c(resultSet_full[i,"Species 1"], resultSet_full[i,"Species 2"])
    genesInSignal <- genesInSignal[which(!is.na(genesInSignal))]

    genesInSignalIDs <- which(gene_id_list %in% genesInSignal)

    # Compute group betweenness centrality
    bet_cent <- kpcent(w, genesInSignalIDs, type="betweenness")


    # Compute group closeness centrality
    close_cent <- kpcent(w, genesInSignalIDs, type="closeness")


    # Add to resultSet
    resultSet_full[i,"GroupBetweenCentrality"] <- bet_cent
    resultSet_full[i,"GroupClosenessCentrality"] <- close_cent

  }

  saveRDS(resultSet_full, resultSet_fname)
}


######### Effect of Noise by Signal #########
signal_summary_df_fname <- file.path(dataDir,"signal_summary_df.Rds")

if(!file.exists(signal_summary_df_fname)) {
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
    #signal_summary_df[which(signal_summary_df$Signal==sig),"OutDegree"] <- 
    #  unique(resultSet_full[which(resultSet_full$SigName_Short == sig),"TotalOutDegree"])
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
  
  saveRDS(signal_summary_df, signal_summary_df_fname)

  
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
  
} else{
  signal_summary_df <- readRDS(signal_summary_df_fname)
}


######## EFFICACY ANALYSIS ############
#debug(cells_x_signals)
# in earlier workflows, I subset the data here - now it just takes all noise levels by default
selectedNoise <- signal_noise
cell_signal_df_fname <- file.path(topoDir,expName,
                                  paste0("cell_signal_df_noise=",
                                         paste0(selectedNoise, collapse = ","),".Rds"))
if(!file.exists(cell_signal_df_fname)) {
  #undebug(cells_x_signals)
  cell_signal_df <- cells_x_signals(paramSets = resultSet_full, # dataframe w/ columns: ModelNo, SetName
                                    topoName = topoName, # string
                                    collectionName = expName, # string
                                    initialClust = 1, # int
                                    targetClust = 2, # int
                                    wt_data = exprMat_norm[ss_unique$StateIndex,1:nGenes], # matrix/dataframe of size numSamples x numFeatures
                                    clust = clust, # vector of length numSamples containing integers
                                    clust_all = clust_full,
                                    tmpMeans = tmpMeans,
                                    tmpSds = tmpSds,
                                    InitRawStates = unnormData,
                                    useDiffs = F)
  saveRDS(cell_signal_df, cell_signal_df_fname)


} else {
  cell_signal_df <- readRDS(cell_signal_df_fname)
}




########## RANKING ANALYSIS ############

boolean_compare_noise <- 0.04
sigEffs_1gene_racipe_fname <- file.path(dataDir, "sigEffs_1gene_comparison.Rds")
sigEffs_2gene_racipe_fname <- file.path(dataDir, "sigEffs_2gene_comparison.Rds")

if(!file.exists(sigEffs_1gene_racipe_fname)) {
  ## Correlation of 1-gene signals
  sigEffs_1gene_racipe <- resultSet_full[which(resultSet_full$NumGenes == 1 & resultSet_full$Noise == boolean_compare_noise), c("Species 1", "ConversionPct")]
  # manually putting in signal data from spin model
  sigEffs_1gene_spin <- data.frame(Species=sigEffs_1gene_racipe$`Species 1`, ConversionPct=NA)
  sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "Cdh1"),"ConversionPct"] <- 0.65
  sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "miR200b"),"ConversionPct"] <- 0.45
  sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "miR200c"),"ConversionPct"] <- 0.25
  sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "miR34a"),"ConversionPct"] <- 0.30
  sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "Ovol2"),"ConversionPct"] <- 0.45
  sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "Grhl2"),"ConversionPct"] <- 0.05
  sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "Foxc2"),"ConversionPct"] <- 0.35
  sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "Zeb1"),"ConversionPct"] <- 1
  sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "Zeb2"),"ConversionPct"] <- 1
  sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "Snai1"),"ConversionPct"] <- 1
  sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "Vim"),"ConversionPct"] <- 0.30
  sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "Twist1"),"ConversionPct"] <- 0.85

  sigEffs_1gene_racipe$ConversionPct_Spin <- sigEffs_1gene_spin$ConversionPct




  ## Correlation of 2-gene signals
  sigEffs_2gene_racipe <- resultSet_full[which(resultSet_full$NumGenes == 2 & resultSet_full$Noise == boolean_compare_noise), c("Species 1", "Species 2", "ConversionPct")]
  sigEffs_2gene_spin <- read.table(file.path(dataDir, "spin_2node_effs.dat"), sep = " ")
  rownames(sigEffs_2gene_spin) <- genes_reordered
  colnames(sigEffs_2gene_spin) <- genes_reordered

  sigEffs_2gene_racipe$ConversionPct_Spin <- NA

  find_index <- function(df, set) {
    for(row in rownames(df)) {
      rowGeneSet <- c(df[row,"Species 1"], df[row, "Species 2"])
      if(length(setdiff(rowGeneSet, set)) == 0) {
        return(row)
      }
    }
  }

  # populate matrix
  for(gene1 in genes_reordered) {
    for(gene2 in genes_reordered) {
      spinEff <- sigEffs_2gene_spin[gene1, gene2]
      idx <- find_index(sigEffs_2gene_racipe, c(gene1, gene2))

      sigEffs_2gene_racipe[idx, "ConversionPct_Spin"] <- spinEff
    }
  }

  # Save results
  saveRDS(sigEffs_2gene_racipe, file=sigEffs_2gene_racipe_fname)
  saveRDS(sigEffs_1gene_racipe, file=sigEffs_1gene_racipe_fname)
} else {
  sigEffs_1gene_racipe <- readRDS(sigEffs_1gene_racipe_fname)
  sigEffs_2gene_racipe <- readRDS(sigEffs_2gene_racipe_fname)
}



######## NONLINEAR EFFECTS #####







######## TRANSITION RATE VS NOISE SIMULATIONS #####

# Set up signals to simulate
time_trial_sig_gene <- c("Zeb1")
time_trial_noise_levels <- c(0, 0.02, 0.04, 0.08, 0.1, 0.2, 0.5)
time_trial_expr <- as.data.frame(t(assay(racipe_bistable)))
time_trial_sig_names <- paste0(time_trial_sig_gene,"_noise=",time_trial_noise_levels)
num_trials <- 10
#time_trial_resultSet[time_trial_setIDList,"SetName"]

# Prepare parameter sets
time_trial_resultSet_fname <- file.path(dataDir,"Zeb1_timeTrial_paramSets.Rds")
if(!file.exists(time_trial_resultSet_fname)) {
  time_trial_geneinfo <- idGenes(clust=clust, initialClust=1, targetClust=2, expr=time_trial_expr)
  time_trial_geneinfo <- topoOutDegrees(geneinfo = time_trial_geneinfo,
                                  topo = sracipeCircuit(racipe))
  time_trial_resultSet <- optimizeSTParams_Clamp(targetClust = 2,
                                           nSigGenes = 1,
                                           geneinfo = time_trial_geneinfo,
                                           noise = time_trial_noise_levels
  )
  time_trial_setIDList <- c(which(time_trial_resultSet$SetName %in% time_trial_sig_names))
  time_trial_setIDList <- rep(time_trial_setIDList, num_trials)

  # Additional parameters to vary in time_trial simulations
  time_trial_resultSet$Tau <- rep(signal_tcorr, nrow(time_trial_resultSet))
  time_trial_resultSet$ScaledNoise <- rep(T, nrow(time_trial_resultSet))
  time_trial_resultSet <- time_trial_resultSet[time_trial_setIDList,]
  rownames(time_trial_resultSet) <- 1:nrow(time_trial_resultSet)
  time_trial_setIDList <- rownames(time_trial_resultSet)

  saveRDS(time_trial_resultSet, file = time_trial_resultSet_fname)

} else {
  time_trial_resultSet <- readRDS(time_trial_resultSet_fname)
  time_trial_setIDList <- rownames(time_trial_resultSet)
}


# Set up new experiment to hold simulation data
time_trial_simTime <- 300
time_trial_expName_list <- paste0("bhtopo_timeTrial_t=",time_trial_simTime,
                           "_relax_OUnoise=",time_trial_resultSet[time_trial_setIDList, "Noise"],
                           "_tau=",signal_tcorr,
                           "_SIG=",time_trial_resultSet[time_trial_setIDList, "SetName"],
                           "_runNo=",rep(seq(num_trials),each=length(time_trial_noise_levels)))

#times <- seq(2, 10, 2)
times <- c(seq(2, 30, 2), seq(35, 100, 5), seq(120, time_trial_simTime, 20))
for(setIDNo in 1:length(time_trial_setIDList)) {
  print(paste0("Beginning set ",setIDNo," of ",length(time_trial_setIDList)))
  setID = time_trial_setIDList[setIDNo]
  setName = time_trial_resultSet[setID, "SetName"]
  expName_new = time_trial_expName_list[setIDNo]
  current_noise <- time_trial_resultSet[setID,"Noise"]
  #debug(calcTransitionRate)
  sampleSet_times <- calcTransitionRate(paramSets = time_trial_resultSet,
                                        setID = setID,
                                        racipe = racipe_bistable_raw,
                                        pca = pca_st,
                                        wt_data = exprMat_norm[ss_unique$StateIndex,genes], # matrix/dataframe of original data
                                        clust = clust, # vector of length numSamples containing integers
                                        clust_all = clust_full, # full cluster labels matching wt_data
                                        tmpMeans = tmpMeans,
                                        tmpSds = tmpSds,
                                        initialClust = 1,
                                        targetClust = 2,
                                        sigName = setName,
                                        outDir = file.path(topoDir),
                                        expName = expName_new,
                                        plot=F,
                                        noise = current_noise,
                                        forceRerun = F,
                                        forceRecompute = F,
                                        anneal=F,
                                        relax=T,
                                        simTimes = times,
                                        simTimeRelax = signal_relaxTime,
                                        save=T,
                                        clamp=T)
}

## Aggregate results
multiSet_fname <- file.path(dataDir,
                            paste0("Zeb1_transition_time_trials=",num_trials,
                                   "_noise=",paste0(time_trial_noise_levels, collapse = ","),"_",".Rds"))
if(!file.exists(multiSet_fname)) {
  multiSet_times <- list()
  for(setIDNo in 1:length(time_trial_setIDList)) {
    setID <- time_trial_setIDList[setIDNo]
    expDir <- file.path(topoDir, time_trial_expName_list[setIDNo])
    sampleSet_times <- readRDS(file.path(expDir,
                                         paste0("transitionTimes_",
                                                time_trial_resultSet[setID,"SetName"],".Rds")))
    multiSet_times[[setIDNo]] <- sampleSet_times
  }
  saveRDS(multiSet_times, multiSet_fname)
} else {
  multiSet_times <- readRDS(multiSet_fname)
}


######## TIME SERIES, ONE MODEL ############
## Simulate a selected model under signaling and save the full trajectory
## We will use Zeb1 clamping as the signal 
## Models: 9, 205, (legacy) 144MET, 589EMT, 601N (new asym models)
hd_traj_sim_time <- 10
hd_traj_step <- 2
hd_traj_noise <- 0.04
hd_traj_tcorr <- 10
hd_traj_sig <- "Zeb1_noise=0.04"
hd_traj_sig <- gsub("0.04", hd_traj_noise, hd_traj_sig)
hd_traj_sigID <- which(resultSet_full$SetName == hd_traj_sig)

# # ?
# setName <- hd_traj_sig
# resultSet$Tau <- signal_tcorr
# resultSet$ScaledNoise <- T
# setIDList <- which(resultSet$`Species 1` == "Zeb1" & resultSet$NumGenes == 1)
# setID <- setIDList[1]
#trajPlots <- T


hd_traj_model_list <- c(10, 297, 144, 589, 601)
hd_traj_no_attempts <- 10


for(hd_traj_model in hd_traj_model_list) {
  for(attemptNo in seq(1, hd_traj_no_attempts)) {
    
    # set up data storage
    modelDataDir <- file.path(dataDir,paste0("trajectories_model",hd_traj_model))
    if(!dir.exists(modelDataDir)) {
      dir.create(modelDataDir)
    }
    
    # Setup
    hd_traj_model_idx <- as.numeric(gsub("Model","",names(models_selected)[which(models_selected == hd_traj_model)])) # for retrieving relevant state
    model_states <- ss_unique[which(ss_unique$Model == hd_traj_model), ]
    
    # Sim parameters
    hd_traj_startTime <- 2
    hd_traj_simTime <- 100
    hd_traj_resolution <- 0.2
    
    
    # may need to insert this after hd_traj_sig: "_tau=",hd_traj_tcorr,
    traj_fname <- file.path(modelDataDir,paste0("trajectory_model",hd_traj_model,
                                                "_simTime=",hd_traj_simTime,
                                                "_SIG=",hd_traj_sig,"_v",attemptNo,".Rds"))
    
    traj_pca_fname <- file.path(modelDataDir,paste0("trajectoryPCA_model",hd_traj_model,
                                                    "_simTime=",hd_traj_simTime,
                                                    "_SIG=",hd_traj_sig,"_v",attemptNo,".Rds"))
    
    
    
    # run simulations if not already done
    if(!file.exists(traj_pca_fname) | !file.exists(traj_fname)) {
      # Set up simulation from original steady states
      hd_traj_checkpoint <- racipe_bistable_raw 
      #dim(hd_traj_checkpoint)
      
      # Filter selected model & set ICs (to WT steady state)
      hd_traj_p2 <- hd_traj_checkpoint[,hd_traj_model_idx]
      hd_traj_p2@metadata$config$simParams["numModels"] <- 1
      sracipeIC(hd_traj_p2) <- assay(racipe_bistable_raw)[,as.numeric(hd_traj_model_idx)] 
      
      # confirm parameters match
      orig_params <- as.numeric(sracipeParams(racipe_bistable_raw)[as.numeric(hd_traj_model_idx),])
      new_params <- as.numeric(sracipeParams(hd_traj_p2))
      # all.equal(new_params, orig_params)
      
      # confirm ICs match previous SS
      orig_SS <- as.numeric(unname(assay(racipe_bistable_raw)[,as.numeric(hd_traj_model_idx)]))
      new_IC <- as.numeric(sracipeIC(hd_traj_p2))
      # all.equal(orig_SS, new_IC)
      
      # QC simulation (ensure steady state is the same, simulate for 10 @ 0.1)
      hd_traj_p2@metadata$config$simParams["nIC"] <- 1
      hd_traj_p2 <- sracipeSimulate(hd_traj_p2, timeSeries = T, printStart = 0, 
                                    printInterval = 0.1, simulationTime = 10,
                                    genParams = F, genIC = F, integrateStepSize = 0.01, initialNoise=0.0, scaledNoise = T, 
                                    stepper = "EM", ouNoise_t = signal_tcorr)
      
      # # confirm parameters still match
      # all.equal(as.numeric(sracipeParams(hd_traj_p2)), 
      #           orig_params)
      # 
      # # confirm ICs still match previous SS
      # all.equal(as.numeric(sracipeIC(hd_traj_p2)), 
      #           orig_SS )
      # 
      # # check if new SS matches ICs and previous SS
      # all.equal(rownames(sracipeIC(hd_traj_p2)), rownames(hd_traj_p2@metadata$timeSeriesDet))
      # 
      # # early
      # qc_ts <- as.numeric(unlist(as.data.frame(hd_traj_p2@metadata$timeSeriesDet)[,1]))
      # all.equal(orig_SS, as.numeric(qc_ts))
      # # mid
      # qc_ts <- as.numeric(unlist(as.data.frame(hd_traj_p2@metadata$timeSeriesDet)[,10]))
      # all.equal(orig_SS, as.numeric(qc_ts))
      # # late
      # qc_ts <- as.numeric(unlist(as.data.frame(hd_traj_p2@metadata$timeSeriesDet)[,101]))
      # all.equal(orig_SS, as.numeric(qc_ts))
      # 
      # # does it approach its other state?
      # model_states <- ss_unique[which(ss_unique$Model == hd_traj_model_idx), ]
      # 
      # model_states
      # round(orig_SS, 1)
      # round(qc_ts, 1)
      
      
      # Implement signal
      sig_clamp_genes <- getClampGenes(resultSet_full, hd_traj_sigID)
      sig_clamp_gene_ids <- unlist(lapply(sig_clamp_genes, function(x) which(rownames(racipe) == x)))
      sig_clamp_df <- getClampDF(clamp_df, sig_clamp_genes, 2)
      colnames(sig_clamp_df) <- as.numeric(sig_clamp_gene_ids)
      sig_clamp_df <- as.matrix(sig_clamp_df)
      
      sig_1step_tgts <- topo[which(topo$Source %in% sig_clamp_genes), "Target"]
      sig_2step_tgts <- unique(topo[which(topo$Source %in% sig_1step_tgts), "Target"])
      sig_3step_tgts <- unique(topo[which(topo$Source %in% sig_2step_tgts), "Target"])
      
      ## Simulate with signal
      hd_traj_p2@metadata$config$simParams["nIC"] <- 1
      hd_traj_p2 <- sracipeSimulate(hd_traj_p2, timeSeries = T, printStart = 0, 
                                    printInterval = hd_traj_resolution, simulationTime = hd_traj_simTime,
                                    genParams = F, genIC = F, integrateStepSize = 0.1, initialNoise=hd_traj_noise, scaledNoise = T, 
                                    stepper = "EM_Clamp", ouNoise_t = signal_tcorr, clampGenes=sig_clamp_gene_ids,
                                    clampValues=sig_clamp_df[hd_traj_model_idx,])
      
      
      hd_traj_df <- as.data.frame(t(hd_traj_p2@metadata$timeSeries)) # for stoch
      #hd_traj_df <- as.data.frame(t(hd_traj_p2@metadata$timeSeriesDet)) # for det
      

      # relaxation simulation 
      hd_traj_relaxation <- hd_traj_p2
      sracipeIC(hd_traj_relaxation) <- as.numeric(t(as.data.frame(t(hd_traj_p2@metadata$timeSeries))[501,]))
      hd_traj_relaxation <- sracipeSimulate(hd_traj_relaxation, timeSeries = T, printStart = 0, 
                                            printInterval = hd_traj_resolution, simulationTime = hd_traj_simTime,
                                            genParams = F, genIC = F, integrateStepSize = 0.1, initialNoise=0.0, scaledNoise = T, 
                                            stepper = "EM_Clamp", ouNoise_t = hd_traj_tcorr)
      
      # Append signaling & relaxation trajectories
      hd_traj_df_relaxation <- as.data.frame(t(hd_traj_relaxation@metadata$timeSeries))
      hd_traj_df_relaxation$Phase <- "Relax"
      hd_traj_df_relaxation$Time <- as.numeric(gsub("X", "", rownames(hd_traj_df_relaxation)))
      
      hd_traj_df$Phase <- "Signal"
      hd_traj_df$Time <- as.numeric(gsub("X", "", rownames(hd_traj_df)))
      
      hd_traj_df_full <- rbind(hd_traj_df, hd_traj_df_relaxation)
      hd_traj_df_full <- hd_traj_df_full[,c(names(tmpMeans), "Phase","Time")]

      # normalize trajectory
      for(gene in genes) {
        hd_traj_df_full[,gene] <- as.numeric(hd_traj_df_full[,gene])
      }
      hd_traj_df_full[,1:length(genes)] <- log2(1+hd_traj_df_full[,1:length(genes)]) # Log transform
      hd_traj_df_full[,1:length(genes)] <- sweep(hd_traj_df_full[,1:length(genes)], 2, tmpMeans, FUN = "-") # scale
      hd_traj_df_full[,1:length(genes)] <- sweep(hd_traj_df_full[,1:length(genes)], 2, tmpSds, FUN = "/") # scale
      
      # Clean up trajectory format
      hd_traj_df_full[which(hd_traj_df_full$Phase == "Relax"),"Time"] <- hd_traj_df_full[which(hd_traj_df_full$Phase == "Relax"),"Time"] + max(hd_traj_df_full[which(hd_traj_df_full$Phase == "Signal"),"Time"])
      rownames(hd_traj_df_full) <- hd_traj_df_full$Time
      hd_traj_df_full <- hd_traj_df_full[,c("Time","Phase",genes_reordered)]
      
      
      # Add ICs to trajectory
      modelICs <- as.data.frame(t(sracipeIC(hd_traj_p2)))[,genes_reordered] # use checkpoint SS
      modelICs[,names(tmpMeans)] <- log2(1+modelICs[,names(tmpMeans)]) # Log transform
      modelICs[,names(tmpMeans)] <- sweep(modelICs[,names(tmpMeans)], 2, tmpMeans, FUN = "-") # scale
      modelICs[,names(tmpMeans)] <- sweep(modelICs[,names(tmpMeans)], 2, tmpSds, FUN = "/") # scale
      
      modelICs$Time <- 0
      modelICs$Phase <- "IC"
      modelICs <- modelICs[,c("Time","Phase",genes_reordered)]
      
      hd_traj_df_full <- as.data.frame(rbind(modelICs, hd_traj_df_full))
      rownames(hd_traj_df_full)[1] <- "0"
      rownames(hd_traj_df_full) <- hd_traj_df_full$Time
      
      
      # Save trajectory & pca projection for later
      hd_traj_df_long <- pivot_longer(hd_traj_df_full, cols=all_of(genes), 
                                      names_to = "Gene", values_to = "Expression")
      
      hd_traj_pca <- as.data.frame(predict(pca, hd_traj_df_full[,names(tmpMeans)]))
      hd_traj_pca$Time <- hd_traj_df_full$Time
      hd_traj_pca$Model <- hd_traj_model
      
      modelDataDir <- file.path(dataDir,paste0("trajectories_model",hd_traj_model))
      if(!dir.exists(modelDataDir)) {
        dir.create(modelDataDir)
      }
      
      traj_fname <- file.path(modelDataDir,paste0("trajectory_model",hd_traj_model,
                                                  "_simTime=",hd_traj_simTime,
                                                  "_SIG=",hd_traj_sig,"_tau=",hd_traj_tcorr,"_v",attemptNo,".Rds"))
      saveRDS(hd_traj_df_long, traj_fname)
      
      traj_pca_fname <- file.path(modelDataDir,paste0("trajectoryPCA_model",hd_traj_model,
                                                      "_simTime=",hd_traj_simTime,
                                                      "_SIG=",hd_traj_sig,"_tau=",hd_traj_tcorr,"_v",attemptNo,".Rds"))
      saveRDS(hd_traj_pca, traj_pca_fname)
      
      
    } else {
      hd_traj_pca <- readRDS(traj_pca_fname)
      hd_traj_df_long <- readRDS(traj_fname)
      hd_traj_df_full <- pivot_wider(
        hd_traj_df_long,
        names_from = "Gene",      # Column with names to become new columns
        values_from = "Expression" # Column with values to fill the new columns
      )
    }
  }

}



