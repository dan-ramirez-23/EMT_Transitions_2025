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
  ss_rounded$Model <- rep(1:numModels, each=numICs)
  unique_state_idx <- which(!duplicated(ss_rounded))
  ss_unique <- ss_rounded %>%
    distinct()
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
  pca_df_full$PC2 <- -1 * pca_df_full$PC2
  pca$rotation[,2] <- -1 * pca$rotation[,2]
  pca$x[,2] <- -1 * pca$x[,2]

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
  revClust <- clust_full
  revClust <- ifelse(clust_full == 1, 2, ifelse(clust_full == 2, 1, clust_full))
  clust_full <- revClust

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
  for(model in models_selected) {
    # Select only epithelial state, and first epithelial state if multiple
    modelStates <- ss_unique[which(ss_unique$Model == model & ss_unique$Cluster == 1),]
    addIdx <- sample(which(row.match(
      as.data.frame(t(round(assay(racipe), 1)))[(numICs*(model-1)+1):(numICs*model),],
      modelStates[,genes],
      nomatch = NA) == 1), 1)
    addIdx <- numICs*(model-1)+addIdx # bring back index for original racipe object
    keepIdx <- c(keepIdx, addIdx)

    unique_state_idx <- which(ss_unique$Model == model & ss_unique$Cluster == 1)
    unique_state_idx_list <- c(unique_state_idx_list, unique_state_idx)
  }

  # subset racipe object for bistable models, update parameters
  racipe_bistable <- racipe[,keepIdx]
  sracipeParams(racipe_bistable) <- sracipeParams(racipe)[as.numeric(unname(models_selected)),] # parameter numbering is messssssed up in sRACIPE, not my fault
  racipe_bistable@metadata$config$simParams[["numModels"]] <- length(keepIdx)

  saveRDS(racipe_bistable, racipe_bistable_fname)
  saveRDS(models_selected, models_selected_fname)
  saveRDS(keepIdx, racipe_bistable_indices_fname)
  saveRDS(unique_state_idx_list, unique_state_idx_list_fname)
  saveRDS(summary_df, summary_df_fname)
  racipe_bistable_indices <- keepIdx
  pca_df <- as.data.frame(pca$x)[as.character(racipe_bistable_indices),]
  clust <- clust_full[unique_state_idx_list]
} else {
  racipe_bistable <- readRDS(racipe_bistable_fname)
  models_selected <- readRDS(models_selected_fname)
  racipe_bistable_indices <- readRDS(racipe_bistable_indices_fname)
  unique_state_idx_list <- readRDS(unique_state_idx_list_fname)
  summary_df <- readRDS(summary_df_fname)
  pca_df <- as.data.frame(pca$x)[as.character(racipe_bistable_indices),]
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

  for(ctrl_current_noise in ctrl_noise_levels[1:4]) {
    for(ctrl_current_trial in 1) {

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

# For signaling, use prcomp object but replace embeddings with subset of models
pca_st <- pca
pca_st$x <- pca_df

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


######## AGGREGATE ANALYSIS #####
resultSet_fname <- file.path(topoDir,expName,"result_summary.Rds")
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
time_trial_noise_levels <- c(0, 0.02, 0.04, 0.08, 0.1, 0.15, 0.2)
time_trial_expr <- as.data.frame(t(assay(racipe_bistable)))
time_trial_sig_names <- paste0(time_trial_sig_gene,"_noise=",time_trial_noise_levels)
#time_trial_resultSet[time_trial_setIDList,"SetName"]

# Prepare parameter sets
time_trial_resultSet_fname <- file.path(dataDir,"Zeb1_timeTrial_paramSets.Rds")
if(!file.exists(time_trial_resultSet_fname)) {
  num_trials <- 10
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
  setID = time_trial_setIDList[setIDNo]
  setName = time_trial_resultSet[setID, "SetName"]
  expName_new = time_trial_expName_list[setIDNo]
  current_noise <- time_trial_resultSet[setID,"Noise"]
  debug(calcTransitionRate)
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










######## NOISE THRESHOLD SIMULATIONS #####

