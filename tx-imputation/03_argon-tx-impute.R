################################################################################
################################################################################
# tempdir <- function() "/Dedicated/jmichaelson-wdata/msmuhammad/scratch/tmp"
# unlockBinding("tempdir", baseenv())
# assignInNamespace("tempdir", tempdir, ns="base", envir=baseenv())
# assign("tempdir", tempdir, baseenv())
# lockBinding("tempdir", baseenv())
# base::tmp
# 
# qlogin -q JM-GPU -pe 80cpn 160 -l h_rt=60:00:00 -N tx-imputation
################################################################################
################################################################################
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
.libPaths("/old_Users/msmuhammad/workbench/miniconda3/envs/tximpute2/lib/R/library")
library(tidyverse);library(data.table);library(doMC);library(readr)
psave <- function(...,file){  
  con = pipe(paste("/Dedicated/jmichaelson-wdata/msmuhammad/workbench/pixz -2 -q 80 -f 3 > ",file,".pxz",sep=""),"wb") 
  save(...,file=con,envir=.GlobalEnv); close(con) 
} 
################################################################################
################################################################################
################################################################################
################################################################################
args <- commandArgs(trailingOnly = T)
# args <- c("tissue",
#           "Brain_Anterior_cingulate_cortex_BA24",
#           "hg19",
#           "/Dedicated/jmichaelson-wdata/msmuhammad/projects/tx-imputation/2024",
#           30)


type <- as.character(args[1]) # either tissue or cell
if (type == "tissue") {
  tissue =T
}
if (type == "celltype") {
  celltype = T
}

type_2 = as.character(args[2]) # type of tissue or celltype
if (!type_2 %in% c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", 
                   "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", 
                   "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Liver", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood",
                   "Excitatory", "Astrocytes", "Endothelial", "Inhibitory", "Microglia", "Oligodendrocytes", "OPCs", "Pericytes", "pb")) {
  message("entered celltype or tissue is not available")
  q(save = "no", status = 1)
}


hg <- as.character(args[3])
if (hg != "hg19") {
  message("lift over your genotypes to hg19")
  q(save = "no", status = 1)
}


project.dir <- as.character(args[4])
setwd(project.dir)


threads <- as.numeric(args[5])
registerDoMC(cores = 10)

with_chr <- T

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

if (type == "tissue") {
  ## define tissue weights
  tissue <- type_2
  weights.path <- "/Dedicated/jmichaelson-wdata/msmuhammad/data/UTMOST-GTEx-model-weights/tmp/"
  
  
  
  ## read genotypes matrix
  genotypes <- fread(file = paste0(project.dir, 
                                   "/data/derivatives/genotypes-subset/", 
                                   tissue, ".xmat.gz"), 
                     header = T, 
                     nThread = 80, verbose = T, 
                     tmpdir = "/Dedicated/jmichaelson-wdata/msmuhammad/scratch/tmp")
  gc()
  # drop the first column and row (family ID and description row for variants)
  # genotypes <- genotypes[-1,-1]
  # gc()
  
  message(paste0("Done reading genotypes file for tissue: ", tissue))
  
  
  
  # get tissue weights
  tissue.weights <- read.table(paste0(weights.path, "rsid-for-", tissue), row.names = 1)
  
  
  
  ## clean variant names if needed
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
  rm(tissue.weights)
  gc()
  
  
  message(paste0("Done reading weights file for tissue: ", tissue))
  
  
  ## make dir for outputs
  system(paste0("mkdir -p ", project.dir, "/data/derivatives/imputed-tx/"))
  
  
  # ## identify intersection between weights and genotypes matrix
  # ge <- intersect(colnames(genotypes), ready.weights$variant)
  # ge <- colnames(genotypes)[(colnames(genotypes) %in% ready.weights$variant)]
  # if (length(ge)==0) {
  #   message("no genotypes found to have weight")
  #   q(save = "no", status = 1)
  # }
  # gc()
  
  ## filter genotypes to the ones in tissue weights
  ## reorder genotypes matrix to match weights variants order
  # genotypes <- data.frame(lapply(genotypes, function(x) as.numeric(x)))[,-1]
  # rownames(genotypes) <- genotypes$IID
  # colnames(genotypes) <- colnames(genotypes)[-1]
  # genotypes <- genotypes[,ge]
  # gc()
  # 
  # genotypes.filt <- genotypes[-1,-1] %>%
  #   column_to_rownames("IID") %>%
  #   select()
  # 
  
  
  ## filter weights to keep intersecting variants
  # filt.weights <- ready.weights %>%
  #   filter(variant %in% ge)
  # # rm(genotypes)
  # gc()
  
  
  ## run imputation
  ## basically multiply variants' weights by the genotypes matrix
  ## do that iteratively by gene
  # filt.weights <- ready.weights
  # chunks <- c(2,seq.int(from = 1000, to = 160000, by = 10000), dim(genotypes)[1])
  
  # imputed.2 <- foreach::foreach(k = 1:length(chunks)-1, .combine = rbind) %dopar% {
  # # imputed.2 <- foreach::foreach(k = 1:5, .combine = rbind) %dopar% {
  #   tmp <- genotypes[chunks[k]:(chunks[k+1]-1), ]
    
    imputed <- foreach(j=1:length(unique(ready.weights$gene)), .combine = cbind) %dopar% {
    # imputed <- foreach(j=1:10, .combine = cbind) %dopar% {
      # j=1
      gen <- unique(ready.weights$gene)[j]
      gene.weights <- ready.weights %>% filter(gene%in%gen) %>%
        pivot_wider(names_from = "variant", values_from = "weight") %>%
        column_to_rownames("gene")
      # imp <- as.matrix(genotypes[,colnames(gene.weights)]) %*% t(gene.weights)
      imp <- as.matrix(genotypes %>% 
                         select(colnames(gene.weights)) %>%
                         mutate_all(.funs = function(x) as.numeric(x))) %*% t(gene.weights)
      return(imp)
    }
  #   return(imputed)
  # }
  gc()
  rownames(imputed) <- genotypes$IID
  psave(imputed, file = paste0(project.dir, "/data/derivatives/imputed-tx/", 
                               tissue, ".rda"))
  message(paste0("done imputing tx for tissue: ", tissue))
}

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

if (type == "celltype") {
  ## define tissue weights
  celltype <- type_2
  weights.path <- "/Dedicated/jmichaelson-wdata/msmuhammad/data/celltypes-cis-eQTLs/data/derivatives/"
  
  
  
  ## read genotypes matrix
  genotypes <- fread(file = paste0(project.dir, 
                                   "/data/derivatives/genotypes-subset/", 
                                   celltype, ".xmat.gz"), 
                     header = T, 
                     nThread = 120)
  # drop the first column and row (family ID and description row for variants)
  genotypes <- genotypes[-1,-1]
  gc()
  
  message(paste0("Done reading genotypes file for celltype: ", celltype))
  
  
  # get celltype weights
  celltype.weights <- read_rds(paste0(weights.path, celltype, "-weights-fdr-sig.rds"))
  
  
  ## clean variant names if needed
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
  message(paste0("Done reading weights file for celltype: ", celltype))
  
  
  ## mkdir for outputs
  system(paste0("mkdir -p ", project.dir, "/data/derivatives/imputed-tx/"))
  
  
  
  # ## identify intersection between weights and genotypes matrix
  # ge <- intersect(colnames(genotypes), ready.weights$variant)
  # if (length(ge)==0) {
  #   message("no genotypes found to have weight")
  #   q(save = "no", status = 1)
  # }
  # 
  # 
  # ## filter genotypes to the ones in celltype weights
  # ## reorder genotypes matrix to match weights variants order
  # filt.genotypes <- data.frame(lapply(genotypes, function(x) as.numeric(x)))[,-1]
  # rownames(filt.genotypes) <- genotypes$IID
  # colnames(filt.genotypes) <- colnames(genotypes)[-1]
  # filt.genotypes <- filt.genotypes[,ge]
  # gc()
  # 
  # 
  # ## filter weights to keep intersecting variants
  # filt.weights <- ready.weights %>%
  #   filter(variant %in% ge)
  # rm(genotypes)
  # gc()
  
  
  
  ## run imputation
  ## basically multiply variants' weights by the genotypes matrix
  ## do that iteratively by gene
  imputed <- foreach(j=1:length(unique(ready.weights$gene)), .combine = cbind) %dopar% {
    # j=1
    gen <- unique(ready.weights$gene)[j]
    gene.weights <- ready.weights %>% filter(gene%in%gen) %>%
      pivot_wider(names_from = "variant", values_from = "weight") %>%
      column_to_rownames("gene")
    # imp <- as.matrix(filt.genotypes[,colnames(filt.genotypes) %in% colnames(gene.weights)]) %*% t(gene.weights)
    imp <- as.matrix(genotypes %>% 
                       select(colnames(gene.weights)) %>%
                       mutate_all(.funs = function(x) as.numeric(x))) %*% t(gene.weights)
    return(imp)
  }
  rownames(imputed) <- genotypes$IID
  
  psave(imputed, file = paste0(project.dir, "/data/derivatives/imputed-tx/", 
                               celltype, ".rda"))
  
  message(paste0("done imputing tx for celltype: ", celltype))
}

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
