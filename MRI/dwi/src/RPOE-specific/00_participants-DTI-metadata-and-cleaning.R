################################################################################
#                    building metadata and cleaning DTI data                   #
################################################################################
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(oro.nifti, lib.loc = sub("tximpute", "ENA", lib.location))
################################################################################
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/dwi"
setwd(project.dir)
################################################################################
################################################################################
# get a list of participants with data collected
ids <- readxl::read_xlsx("../shared_data/data/RPOE_meta.xlsx",
                         sheet = "MRI-processing") %>%
  select(1,2)
# filter to new participants
ids <- ids[c(77:80),]

# list data files and make a sym link in my directory
registerDoMC(4)
meta <- foreach(i = 1:nrow(ids), .combine = rbind) %dopar% {
  id <- ids$te_id[i]
  df <- data.frame(file = list.files(paste0("/Dedicated/jmichaelson-sdata-new/private-data/RPOE_MR/rawdata/sub-", id),
                                     recursive = T, pattern = "DTI")) %>%
    mutate(session = as.numeric(str_sub(sub("_.+", "", sub("/dwi.*", "", sub("ses-", "", file))), start = 1, end = 8))) %>%
    filter(session == max(session)) %>%
    mutate(te_id = id, devGenes_id = ids$devGenes_id[i],
           type = sub("_", "", sub("[0-9]+", "", sub(".+_ses-", "", file))),
           ext = sub(".+\\.", "", type),
           ext = ifelse(ext == "gz", "nii.gz", ext),
           type = sub("\\..*", "", type),
           type = case_when(type == "DTI_32_DIR" ~ "DTI_32_DIR",
                            type == "DTI_32_DIR_5" ~ "DTI_32_DIR",
                            type == "DTI_32_DIR_6" ~ "DTI_32_DIR",
                            type == "DTI_32_DIR_7" ~ "DTI_32_DIR",
                            type == "DTI_32_DIR_8" ~ "DTI_32_DIR",
                            type == "1_DTI_32_DIR_7" ~ "DTI_32_DIR",
                            type == "2_DTI_32_DIR_7" ~ "DTI_32_DIR",
                            type == "DTI_b__32_Dir_6" ~ "DTI_32_DIR",
                            type == "DTI__DIR_5" ~ "DTI_32_DIR",
                            type == "DTI__DIR_6" ~ "DTI_32_DIR",
                            type == "DTI__DIR_7" ~ "DTI_32_DIR",
                            type == "DTI__DIR_8" ~ "DTI_32_DIR",
                            type == "DTI__DIR_3" ~ "DTI_32_DIR",
                            type == "DTI__DIR_4" ~ "DTI_32_DIR",
                            type == "DTI_Rev_PE" ~ "DTI_Rev_PE",
                            type == "DTI_Rev_PE_8" ~ "DTI_Rev_PE",
                            type == "1_DTI_Rev_PE_8" ~ "DTI_Rev_PE",
                            type == "2_DTI_Rev_PE_8" ~ "DTI_Rev_PE",
                            type == "DTI_Rev_PE_9" ~ "DTI_Rev_PE",
                            type == "DTI_Rev_PE_5" ~ "DTI_Rev_PE",
                            type == "DTI_Rev_PE_7" ~ "DTI_Rev_PE",
                            type == "DTI_Rev_PE_" ~ "DTI_Rev_PE",
                            type == "DTI_Flip_Phase_" ~ "DTI_Rev_PE",
                            type == "ORIG_DTI_32_DIR" ~ "ORIG_DTI_32_DIR",
                            type == "ORIG_DTI_32_DIR_710" ~ "ORIG_DTI_32_DIR",
                            type == "ORIG_DTI_32_DIR_810" ~ "ORIG_DTI_32_DIR",
                            type == "1_ORIG_DTI_32_DIR_710" ~ "ORIG_DTI_32_DIR",
                            type == "2_ORIG_DTI_32_DIR_710" ~ "ORIG_DTI_32_DIR",
                            type == "ORIG_DTI__DIR_710" ~ "ORIG_DTI_32_DIR",
                            type == "ORIG_DTI__DIR_810" ~ "ORIG_DTI_32_DIR",
                            type == "ORIG_DTI__DIR_610" ~ "ORIG_DTI_32_DIR",
                            type == "ORIG_DTI__DIR_310" ~ "ORIG_DTI_32_DIR",
                            type == "ORIG_DTI__DIR_510" ~ "ORIG_DTI_32_DIR",
                            type == "ORIG_DTI__DIR_410" ~ "ORIG_DTI_32_DIR"),
           new_name = paste0("sub-", te_id, "_", type, ".", ext))
  if(nrow(df)<1){
    return()
  }
  system(paste0("rm -rf ", "/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/dwi/data/raw/clean/", id))
  system(paste0("mkdir -p ", "/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/dwi/data/raw/clean/", id))
  for (i2 in 1:nrow(df)) {
    system(paste0("ln -s ", "/Dedicated/jmichaelson-sdata-new/private-data/RPOE_MR/rawdata/sub-", id, "/",
                  df$file[i2], " ", 
                  "/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/dwi/data/raw/clean/", id, "/",
                  df$new_name[i2]))
    
  }
  return(df)
}
# read old
meta.o <- read_csv("data/collected-dti-list.csv")
meta.n <- full_join(meta.o, meta)

write_csv(meta.n, "data/collected-dti-list.csv")

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
