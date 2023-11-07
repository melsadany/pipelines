########################################################################
#                         make GRMs trial 1                            #
########################################################################

muhammad.wd <- "/wdata/msmuhammad/cratch/GWAS/"
library(tidyverse)
my_wd <- "/wdata/trthomas/spark_sleep/2022-06/"
setwd(my_wd)

pheno <- read_tsv("data/sleep_pheno_ABCD_SPARK_combined.tsv")
master <- read_tsv("/wdata/trthomas/array/merged_2022_ABCD_iWES1_WGS_2-4/master.tsv", guess_max = 90000)
pcs <- read_tsv("/Dedicated/jmichaelson-wdata/trthomas/array/merged_2022_ABCD_iWES1_WGS_2-4/PCA/all/PCs.tsv")
pcs <- pcs[,1:21]
cca_covars <- read_tsv("/wdata/jmichaelson/SPARK/ResearchMatch/eat_sleep/sleep/2022/CC_covariates.txt", col_names = c("IID", "CCA1", "CCA2", "CCA3", "CCA4", "CCA5"))

sum_covars <- read_tsv("EFA/sum_scores.tsv")
