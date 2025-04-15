################################################################################
#                       pipeline for drug response prediction                  #
################################################################################
# rm(list = ls())
gc()
.libPaths("/old_Users/msmuhammad/workbench/miniconda3/envs/tximpute2/lib/R/library")
# source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
# library(tidyverse, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/theone/lib/R/library")
library(tidyverse)
# library(doMC, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/theone/lib/R/library")
library(doMC)
# library(readr, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/theone/lib/R/library")
library(readr)
library(data.table, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/theone/lib/R/library")
# library(bigsnpr, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/theone/lib/R/library")
library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE);options(default.nproc.blas = NULL)
################################################################################
################################################################################
# functions you need
# you need input of: 
#     run_type: either local or as a job script (LOCAL RUN IS NOT READY)
#     source name (either tissue name or celltype)
#     type: either tissue or celltype
#     genotypes file path
#     project.dir to save files
#     # threads to use
#     where to save: folder of where to save imputed tx
################################################################################
################################################################################
################################################################################
################################################################################
args <- commandArgs(trailingOnly = T)
# args <- c("job",
#           "tissue",
#           "Brain_Frontal_Cortex_BA9",
#           "/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/genetics/data/derivatives/FINAL_NDVR-merged-w-1KG-hmp3-bigsnpr.rds",
#           "/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/genetics",
#           "15",
#           "/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/genetics/data/derivatives/tx-imputation")

if (args[1] == "job") {
  source <- args[2]
  type <- args[3]
  genotypes.f.path <- args[4]
  project.dir <- args[5]
  threads <- as.numeric(args[6])
  dir.save <- args[7]
  setwd(project.dir)
  system(paste0("mkdir -p ", dir.save))
  
  if (!file.exists(genotypes.f.path)) {
    # Reading the bedfile and storing the data in temporary directory
    rds <- snp_readBed2(bedfile = sub("-bigsnpr\\.rds", "\\.bed",
                                      genotypes.f.path), 
                        backingfile = sub("\\.rds", "",genotypes.f.path),
                        ncores = 5)
  }
  # Loading the data from backing files
  bed.f <- snp_attach(genotypes.f.path)
} 
################################################################################
############################ tx imputation function ############################
################################################################################
impute.tx <- function(genotypes, weights, nThread) {
  # the genotypes matrix is expected to have participants as rownames and genotypes as colnames
  # the weights matrix is expected to have 3 columns: variant, gene, weight
  ge <- intersect(colnames(genotypes), weights$variant)
  if (length(ge)==0) {
    print("no genotypes found to have weight")
    return(NULL)
  }
  filt.genotypes <- genotypes[,ge]
  gc()
  
  filt.weights <- weights %>%
    filter(variant %in% ge)
  gc()
  registerDoMC(cores = nThread)
  
  imputed <- foreach(j=1:length(unique(filt.weights$gene)), .combine = cbind) %dopar% {
    # j=1
    gen <- unique(filt.weights$gene)[j]
    gene.weights <- filt.weights %>% filter(gene%in%gen) %>%
      pivot_wider(names_from = "variant", values_from = "weight") %>%
      column_to_rownames("gene")
    imp <- as.matrix(filt.genotypes[,colnames(gene.weights)]) %*% t(gene.weights)
    return(imp)
  }
  return(imputed)
}

################################################################################
################################################################################
################################################################################
################################################################################
if (args[1] == "job") {
  genotypes <- bed.f$genotypes[]
  rownames(genotypes) <- bed.f$fam$sample.ID
  colnames(genotypes) <- bed.f$map$marker.ID
  gc()
  print(paste0("Done with: ", "reading genotypes file for tissue"))
  if (source == "tissue") {
    # get tissue weights
    tissue.weights <- read.table(paste0("/Dedicated/jmichaelson-wdata/msmuhammad/projects/tx-imputation/UTMOST-GTEx-model-weights/",
                                        "tmp/rsid-for-", type), 
                                 row.names = 1)
    ready.weights <- tissue.weights %>% 
      dplyr::select(variant = KG_ID, gene, weight) %>%
      distinct(variant, gene, .keep_all = T) %>%
      filter(variant %in% colnames(genotypes))
    gc()
    print(paste0("Done with: ", "reading weights file for tissue"))
  } else if(source == "celltype") {
    # get celltype weights
    celltype.weights <- read_rds(paste0("/Dedicated/jmichaelson-wdata/msmuhammad/data/celltypes-cis-eQTLs/data/derivatives/", 
                                        type,"-weights-fdr-sig.rds"))
    ready.weights <- celltype.weights %>%
      filter(FDR<0.05) %>%
      mutate(KG_ID = sub("chr","", ID_37)) %>%
      dplyr::select(variant=KG_ID, gene, weight=beta) %>%
      distinct(variant, gene, .keep_all = T) %>%
      filter(variant %in% colnames(genotypes))
    rm(celltype.weights)
    gc()
    print(paste0("Done with: ", "reading weights file for celltype"))
  }
  
  # impute tx for all genes
  imputed.tx <- impute.tx(genotypes = genotypes, 
                          weights = ready.weights, 
                          nThread = threads) %>% 
    as.data.frame()
  gc()
  print(paste0("Done with: ", "imputing tx"))
  
  # save imputed tx
  write_rds(imputed.tx, 
            file = paste0(dir.save, "/",
                          source, "_", type, ".rds"),
            compress = "gz")
  gc()
  print(paste0("Done with: ", "saving imputed tx file"))
}

################################################################################
################################################################################
################################################################################
