################################################################################
#                            PGS calculation pipeline                          #
################################################################################
rm(list = ls()); gc()
library(tidyverse);library(doMC);library(foreach);set.seed(123);library(data.table);library(runonce);library(remotes);library(bigsnpr);library(magrittr);library(bigreadr)
options(bigstatsr.check.parallel.blas = FALSE);options(default.nproc.blas = NULL)
argon <- T
if (argon == T) {
  args <- commandArgs(trailingOnly = T)
  id <- as.numeric(args[1])
  cores <- as.numeric(args[2])
  bed.path <- args[3]
  output.dir <- args[4]
  message(paste0("here's the id: ", id, " and here's the cores' number: ", cores, " and here's the bed file path :",
                 bed.path, " and here's the output dir path: ", output.dir))
}
################################################################################
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/genetics"
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
ss.meta <- read_csv("/Dedicated/jmichaelson-wdata/msmuhammad/data/gwas-sumstats/summary-stats-metadata.csv")
if (argon==T) {
  ss.meta <- ss.meta[id,]
}
####
# load my PGS calc function
source("/Dedicated/jmichaelson-wdata/msmuhammad/workbench/customized-functions/calc_pgs_ME_V.R")

calc_pgs_ME_V(bed_filepath = bed.path,
              n_cores = cores, 
              build = "hg19", 
              output_directory = output.dir, # will save a tsv file for each PGS
              ss.meta = ss.meta, # give it the list of PGS scores you want
              combine = F) # this is to combine the PGS at the end or not. it will make a wide dataframe for combining all PGS files in the output directory
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
