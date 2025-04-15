################################################################################
#                            PGS calculation pipeline                          #
################################################################################
rm(list = ls()); gc()
# .libPaths(c("/old_Users/msmuhammad/workbench/miniconda3/envs/tximpute2/lib/R/library"))
# library(data.table, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/theone/lib/R/library")
# library(tidyverse, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/theone/lib/R/library")
library(tidyverse);library(doMC);library(foreach);set.seed(123);library(data.table);library(runonce);library(remotes);library(bigsnpr);library(magrittr);library(bigreadr)
options(bigstatsr.check.parallel.blas = FALSE);options(default.nproc.blas = NULL)
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
# argon <- T
if (device == "argon") {
  args <- commandArgs(trailingOnly = T)
  id <- as.numeric(args[1])
  cores <- as.numeric(args[2])
  bed.path <- args[3]
  output.dir <- args[4]
  project.dir <- args[5]
  message(paste0("here's the id: ", id, " and here's the cores' number: ", cores, " and here's the bed file path :",
                 bed.path, " and here's the output dir path: ", output.dir))
} else {
  device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
  project.dir <- paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
                        "/jmichaelson-wdata/msmuhammad/projects/RPOE/genetics")
  bed.path<-paste0(project.dir, "/data/derivatives/FINAL_NDVR-merged-w-1KG-hmp3.bed")
  output.dir <- paste0(project.dir, "/data/derivatives/pgs")
}
################################################################################
################################################################################
setwd(project.dir)
################################################################################
################################################################################
# this script is to have a standard pipeline for calculation PGS
# the script will use the clean/reformatted GWAS sum stats from:
# /wdata/msmuhammad/data/gwas-sumstats/ALL/LC
# the script will use ldpred to compute these PGS
# all you need as an input is the bed file path for your samples
# you MUST have tidyverse, bigsnpr, and other packages loaded above
# if you want to calculate PGS for certain categories only, make sure to filter the ss.meta file at first
################################################################################
################################################################################
# load the GWAS sum stats metafile
# if you want fewer PGS, filter this file
ss.meta.r <- read_csv(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
                           "/jmichaelson-wdata/msmuhammad/data/gwas-sumstats/summary-stats-metadata.csv"))
# ss.meta.r <- ss.meta.r[119:nrow(ss.meta.r),]
if (device=="argon") {
  ss.meta <- ss.meta.r[id,]
  # load my PGS calc function
  source("/Dedicated/jmichaelson-wdata/msmuhammad/pipelines/PGS/calc_pgs_ME_V.R")
  
  calc_pgs_ME_V(bed_filepath = bed.path,
                device = "argon",
                n_cores = 3, 
                build = "hg19", 
                output_directory = output.dir, # will save a tsv file for each PGS
                ss.meta = ss.meta, # give it the list of PGS scores you want
                sd = ifelse(is.na(ss.meta$sd[1]), 0, ss.meta$sd[1]),
                combine = F) # this is to combine the PGS at the end or not. it will make a wide dataframe for combining all PGS files in the output directory
  # clean and remove tmp files
  system(paste0("rm -rf ", project.dir, "/tmp-", ss.meta$file[1]))
  
} else {
  source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
                "/jmichaelson-wdata/msmuhammad/pipelines/PGS/calc_pgs_ME_V.R"))
  registerDoMC(cores = 2)
  foreach(i = 1:nrow(ss.meta.r)) %dopar% {
    calc_pgs_ME_V(bed_filepath = bed.path,
                  device = device,
                  n_cores = 3,
                  build = "hg19",
                  output_directory = output.dir, # will save a tsv file for each PGS
                  ss.meta = ss.meta.r[i,] %>% 
                    mutate(full_path = sub("/Dedicated", 
                                           ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
                                           full_path)), # give it the list of PGS scores you want
                  sd = ifelse(is.na(ss.meta.r$sd[i]), 0, ss.meta.r$sd[i]),
                  combine = F) # this is to combine the PGS at the end or not. it will make a wide dataframe for combining all PGS files in the output directory
    # clean and remove tmp files
    system(paste0("rm -rf ", project.dir, "/tmp-", ss.meta.r$file[i]))
    gc()
  }
  
}
####
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
