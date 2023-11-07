################################################################################
#                    functions needed to impute transcriptome                  #
################################################################################
# rm(list = ls())
gc()
# .libPaths("/Users/msmuhammad/workbench/miniconda3/envs/tximpute2/lib/R/library")
# source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(tidyverse, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/tximpute/lib/R/library")
library(doMC, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/tximpute/lib/R/library")
library(readr, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/tximpute/lib/R/library")
library(data.table, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/tximpute/lib/R/library")
registerDoMC(cores = 6)
################################################################################
################################################################################
# functions you need
# ifelse(verbose, print(paste0()), NULL)
################################################################################
#####
# saving files in rds.pxz format and reading them
pdssave <- function(...,file){  
  #' @description this function saves files in rds format and compresses it using pixz
  #' @param file character. takes the file name of where to save the file. Do NOT include ".pxz" in the file name here
  con = pipe(paste("/Dedicated/jmichaelson-wdata/msmuhammad/workbench/pixz -2 -q 80 -f 3 > ",file,".pxz",sep=""),"wb") 
  saveRDS(...,file=con)
}
pdsload <- function(fname,envir=.GlobalEnv){
  con <- pipe(paste("/Dedicated/jmichaelson-wdata/msmuhammad/workbench/pixz -d <",fname),"rb")
  return(readRDS(con))
}
#####
################################################################################
#####
# subset genotypes to the ones with eQTL weights
subset_genotypes <- function(genotypes_path_base_name, 
                             project_dir, 
                             tissue = T, 
                             celltype = F, 
                             tissue_type = "brain",
                             celltype_type = "all", 
                             build = "hg19",
                             with_chr = T,
                             verbose = T) {
  #' @description This function subsets the huge genotypes bed file
  #' to only keeps genotypes of interest for the selected tissues or celltypes
  #'
  #' @param genotypes_path_base_name character. give this the path to your genotypes bed file and only include the base name of the bed
  #' @param project_dir character. give this the main project directory that you want to save files at. 
  #' the function will build data/derivatives/genotypes-subset directory inside this directory by default
  #' @param tissue_type logical. This to identify if you want to subset genotypes for a tissue or not. default = T.
  #' @param celltype_type character. this is to identify if you just want certain celltype or all tissues from GTEx or a specific tissue. either "brain" or "all". default = "all"
  #' @param celltype logical. This to identify if you want to subset genotypes for a celltype or not. default = F.
  #' @param tissue character. this is to identify if you just want brain tissues or all tissues from GTEx or a specific tissue. either "brain" or "all". default = "brain"
  #' @param build character. this identifies if your genotypes build is in "hg19" or not. 
  #' @param with_chr logical. this indicates the naming format for the variants in the genotypes BED file. TRUE means the genotype is names in this format: chr4:101668718:A:G, and FALSE means it's named in this format 4:101668718:A:G
  #' @param verbose logical. print messages or not. 
  #' 
  #' @return the function does not return anything. it saves the subsetted genotypes by default
  
  if (tissue == T) {
    tissues <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Liver", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")
    if (length(tissue_type) == 1) {
      if (tissue_type == "brain") {
        tissues <- tissues[7:16]
      } else if (tissue_type == "all") {
        tissues <- tissues
      }
    } else {
      tissues <- tissue_type
    }
    registerDoMC(cores = 6)
    foreach::foreach(i = 1:length(tissues)) %dopar% {
      tissue <- tissues[i]
      weights.path <- "/Dedicated/jmichaelson-wdata/msmuhammad/data/UTMOST-GTEx-model-weights/tmp/"
      if (build == "hg19") {
        if (with_chr == T) {
          ids <- paste0(weights.path, "rsid-ID02-UTMOST-for-", tissue)
        } else {
          ids <- paste0(weights.path, "rsid-ID02-UTMOST-no-chr-for-", tissue)
        }
        
        gcta_command <- paste(
          "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/gcta/gcta-1.94.1", # GCTA executable
          "--bfile", genotypes_path_base_name, #PLINK files
          "--extract", ids, #list of SNPs
          "--thread-num", 2, # threads
          "--recode",
          "--out", paste0(project_dir, "/data/derivatives/genotypes-subset/", tissue),
          sep = " ")
        system(paste0("mkdir -p ", project_dir, "/data/derivatives/genotypes-subset/"))
        ifelse(verbose, print(paste0("this is the GCTA command for tissue: ", 
                                     tissue, "\n",
                                     gcta_command)),
               NULL)
        system(gcta_command)
        ifelse(verbose, print(paste0("done subsetting for tissue: ", tissue)), NULL)
      }
    }
  }
  if (celltype == T) {
    # print("pipeline not ready")
    celltypes <- c("Excitatory", "Astrocytes", "Endothelial", "Inhibitory", "Microglia", "Oligodendrocytes", "OPCs", "Pericytes", "pb")
    if (length(celltype_type) ==1) {
      if (celltype_type == "all") {
        celltypes <- celltypes
      } 
    } else {
      celltypes <- celltype_type
    }
    registerDoMC(cores = 6)
    foreach::foreach(i = 1:length(celltypes)) %dopar% {
      celltype <- celltypes[i]
      weights.path <- "/Dedicated/jmichaelson-wdata/msmuhammad/data/celltypes-cis-eQTLs/data/derivatives/"
      if (build == "hg19") {
        if (with_chr == T) {
          ids <- paste0(weights.path, celltype, "-ID_37-FDR")
        } else {
          ids <- paste0(weights.path, celltype, "-ID_37-FDR-no-chr")
        }
        gcta_command <- paste(
          "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/gcta/gcta-1.94.1", # GCTA executable
          "--bfile", genotypes_path_base_name, #PLINK files
          "--extract", ids, #list of SNPs
          "--thread-num", 2, # threads
          "--recode",
          "--out", paste0(project_dir, "/data/derivatives/genotypes-subset/", celltype),
          sep = " ")
        system(paste0("mkdir -p ", project_dir, "/data/derivatives/genotypes-subset/"))
        ifelse(verbose, print(paste0("this is the GCTA command for celltype: ", 
                                     celltype, "\n",
                                     gcta_command)),
               NULL)
        system(gcta_command)
        ifelse(verbose, print(paste0("done subsetting for celltype: ", celltype)), NULL)
      }
      
  }
  }
}
#####
################################################################################
################################################################################
#####
# impute tx for selected tissues 
impute.tx <- function(project_dir, 
                      tissue = T, 
                      celltype = F, 
                      tissue_type = "brain",
                      celltype_type = "all", 
                      build = "hg19",
                      verbose = T, 
                      with_chr = T,
                      threads = 6) {
  #' @description This function imputes transcriptome for selected tissue using genotypes subset
  #'
  #' @param project_dir character. give this the main project directory that you want to save files at. 
  #' the function will build data/derivatives/imputed-tx directory inside this directory by default. 
  #' It assumes having the subsetted directory created by "subset_genotypes" function.
  #' @param tissue logical. This to identify if you want to subset genotypes for a tissue or not. default = T.
  #' @param celltype logical. This to identify if you want to subset genotypes for a celltype or not. default = F.
  #' @param tissue_type character. this is to identify if you just want brain tissues or all tissues from GTEx or a specific tissue. either "brain" or "all". default = "brain"
  #' @param celltype_type character. this is to identify if you just want certain celltype or all tissues from GTEx or a specific tissue. either "brain" or "all". default = "all"
  #' @param build character. this identifies if your genotypes build is in "hg19" or not. 
  #' @param verbose logical. print messages or not. 
  #' @param with_chr logical. this indicates the naming format for the variants in the genotypes BED file. TRUE means the genotype is names in this format: chr4:101668718:A:G, and FALSE means it's named in this format 4:101668718:A:G
  #' @param threads integer. number of threads to use
  #' 
  #' @return the function does not return anything. it saves the imputed-tx by default
  
  if (tissue == T) {
    tissues <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Liver", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")
    if (length(tissue_type) == 1) {
      if (tissue_type == "brain") {
        tissues <- tissues[7:16]
      } else if (tissue_type == "all") {
        tissues <- tissues
      }
    } else {
      tissues <- tissue_type
    }
    registerDoMC(cores = 6)
    foreach::foreach(i = 1:length(tissues)) %dopar% {
      tissue <- tissues[i]
      weights.path <- "/Dedicated/jmichaelson-wdata/msmuhammad/data/UTMOST-GTEx-model-weights/tmp/"
      if (build == "hg19") {
        genotypes <- fread(file = paste0(project_dir, 
                                         "/data/derivatives/genotypes-subset/", 
                                         tissue, ".xmat.gz"), 
                           header = T, 
                           nThread = 2)
        genotypes <- genotypes[-1,-1]
        gc()
        ifelse(verbose, print(paste0("Done reading genotypes file for tissue: ", tissue)), NULL)
        
        # get tissue weights
        tissue.weights <- read.table(paste0(weights.path, "rsid-for-", tissue), row.names = 1)
        if (with_chr == F) {
          ready.weights <- tissue.weights %>% 
            mutate(ID_02_UTMOST = sub("chr", "", ID_02_UTMOST)) %>%
            dplyr::select(variant = ID_02_UTMOST, gene, weight) %>%
            distinct(variant, gene, .keep_all = T) %>%
            filter(variant %in% colnames(genotypes))
        } else {
          ready.weights <- tissue.weights %>% 
            dplyr::select(variant = ID_02_UTMOST, gene, weight) %>%
            distinct(variant, gene, .keep_all = T) %>%
            filter(variant %in% colnames(genotypes))
        }
        gc()
        ifelse(verbose, print(paste0("Done reading weights file for tissue: ", tissue)), NULL)
        system(paste0("mkdir -p ", project_dir, "/data/derivatives/imputed-tx/"))
        
        ge <- intersect(colnames(genotypes), ready.weights$variant)
        if (length(ge)==0) {
          print("no genotypes found to have weight")
          return(NULL)
        }
        filt.genotypes <- data.frame(lapply(genotypes, function(x) as.numeric(x)))[,-1]
        rownames(filt.genotypes) <- genotypes$IID
        colnames(filt.genotypes) <- colnames(genotypes)[-1]
        filt.genotypes <- filt.genotypes[,ge]
        gc()
        
        filt.weights <- ready.weights %>%
          filter(variant %in% ge)
        gc()
        imputed <- foreach(j=1:length(unique(filt.weights$gene)), .combine = cbind) %dopar% {
          # j=1
          gen <- unique(filt.weights$gene)[j]
          gene.weights <- filt.weights %>% filter(gene%in%gen) %>%
            pivot_wider(names_from = "variant", values_from = "weight") %>%
            column_to_rownames("gene")
          imp <- as.matrix(filt.genotypes[,colnames(gene.weights)]) %*% t(gene.weights)
          return(imp)
        }
        pdssave(imputed, file = paste0(project_dir, "/data/derivatives/imputed-tx/", 
                                       tissue, ".rds"))
        ifelse(verbose, print(paste0("done imputing tx for tissue: ", tissue)), NULL)
      }
    }
  }
  if (celltype == T) {
    # print("pipeline not ready")
    celltypes <- c("Excitatory", "Astrocytes", "Endothelial", "Inhibitory", "Microglia", "Oligodendrocytes", "OPCs", "Pericytes", "pb")
    if (length(celltype_type) ==1) {
      if (celltype_type == "all") {
        celltypes <- celltypes
      } 
    } else {
      celltypes <- celltype_type
    }
    
    registerDoMC(cores = 6)
    foreach::foreach(i = 1:length(celltypes)) %dopar% {
      celltype <- celltypes[i]
      weights.path <- "/Dedicated/jmichaelson-wdata/msmuhammad/data/celltypes-cis-eQTLs/data/derivatives/"
      if (build == "hg19") {
        genotypes <- fread(file = paste0(project_dir, 
                                         "/data/derivatives/genotypes-subset/", 
                                         celltype, ".xmat.gz"), 
                           header = T, 
                           nThread = 2)
        genotypes <- genotypes[-1,-1]
        gc()
        ifelse(verbose, print(paste0("Done reading genotypes file for celltype: ", celltype)), NULL)
        
        # get celltype weights
        celltype.weights <- read_rds(paste0(weights.path, celltype, "-weights-fdr-sig.rds"))
        if (with_chr == F) {
          ready.weights <- celltype.weights %>% 
            filter(FDR<0.05) %>%
            mutate(ID_37 = sub("chr", "", ID_37)) %>%
            dplyr::select(variant=ID_37, gene, weight=beta) %>%
            distinct(variant, gene, .keep_all = T) %>%
            filter(variant %in% colnames(genotypes))
        } else {
          ready.weights <- celltype.weights %>% 
            filter(FDR<0.05) %>%
            dplyr::select(variant=ID_37, gene, weight=beta) %>%
            distinct(variant, gene, .keep_all = T) %>%
            filter(variant %in% colnames(genotypes))
        }
        gc()
        ifelse(verbose, print(paste0("Done reading weights file for celltype: ", celltype)), NULL)
        system(paste0("mkdir -p ", project_dir, "/data/derivatives/imputed-tx/"))
        
        ge <- intersect(colnames(genotypes), ready.weights$variant)
        if (length(ge)==0) {
          print("no genotypes found to have weight")
          return(NULL)
        }
        filt.genotypes <- data.frame(lapply(genotypes, function(x) as.numeric(x)))[,-1]
        rownames(filt.genotypes) <- genotypes$IID
        colnames(filt.genotypes) <- colnames(genotypes)[-1]
        filt.genotypes <- filt.genotypes[,ge]
        gc()
        
        filt.weights <- ready.weights %>%
          filter(variant %in% ge)
        gc()
        imputed <- foreach(j=1:length(unique(filt.weights$gene)), .combine = cbind) %dopar% {
          # j=1
          gen <- unique(filt.weights$gene)[j]
          gene.weights <- filt.weights %>% filter(gene%in%gen) %>%
            pivot_wider(names_from = "variant", values_from = "weight") %>%
            column_to_rownames("gene")
          imp <- as.matrix(filt.genotypes[,colnames(filt.genotypes) %in% colnames(gene.weights)]) %*% t(gene.weights)
          return(imp)
        }
        pdssave(imputed, file = paste0(project_dir, "/data/derivatives/imputed-tx/", 
                                       celltype, ".rds"))
        ifelse(verbose, print(paste0("done imputing tx for celltype: ", celltype)), NULL)
      }
    }
  }
}
#####
################################################################################
################################################################################
