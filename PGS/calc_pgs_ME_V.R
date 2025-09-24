################################################################################
#                            PGS calculation function                          #
################################################################################
################################################################################
# this script is to have a standard pipeline for calculation PGS
# the script will use the clean/reformatted GWAS sum stats from:
# /wdata/msmuhammad/data/gwas-sumstats/ALL/LC
# the script will use ldpred to compute these PGS
# all you need as an input is the bed file path for your samples
# you MUST have tidyverse, bigsnpr
# if you want to calculate PGS for certain categories only, make sure to filter the ss.meta file at first
################################################################################
################################################################################
# load the GWAS sum stats metafile
# if you want fewer PGS, filter this file
# ss.meta <- read_csv("/Dedicated/jmichaelson-wdata/msmuhammad/data/gwas-sumstats/sum-stats-metadata.csv")
####
# my function to run ldpred and out pgs by trait
calc_pgs_ME_V <- function(bed_filepath,
                          device,
                          build = "hg19",
                          output_directory = "PGS-output",
                          n_cores = 1,
                          ss.meta = ss.meta, 
                          combine = T,
                          sd) {
  # build the output directory, if it doesn't exist
  system(paste0("mkdir -p ", output_directory))
  
  #### get LDPred2 function ####
  source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
                "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/LDPred2_argon_hm3plus_ME.R"))
  map_ld = read_rds(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
                           '/jmichaelson-wdata/lcasten/tools/LDPred2/map_hm3_plus.rds'))
  map_ldref = read_rds(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
                              '/jmichaelson-wdata/lcasten/tools/LDPred2/map_hm3_plus.rds'))
  ####
  
  #### next chunk is from Lucas Casten ####
  ## subset to overlapping SNPs
  bim_snp <- read_tsv(file = str_replace_all(bed_filepath, pattern = '[.]bed', replacement = '.bim'), col_names = F) %>%
    rename(chr = X1, pos = X4) %>% select(chr, pos, major = X6, minor = X5) %>%
    filter(nchar(major) == 1 & nchar(minor) == 1)
  map_ld <- map_ld %>%
    inner_join(bim_snp) %>%
    filter(major == a0 | major == a1) %>% filter(minor == a0 | minor == a1) %>%
    select(-c(major, minor))
  map_ldref <- map_ldref %>%
    inner_join(bim_snp) %>%
    filter(major == a0 | major == a1) %>% filter(minor == a0 | minor == a1) %>%
    select(-c(major, minor))
  ##### STEP 0: get geno data ready ####
  rds.path <- str_replace_all(bed_filepath, pattern = '[.]bed', replacement = '.rds')
  if (!file.exists(rds.path)) {
    # run this command, just once before anything:
    # this is a function from bigsnpr, that will make an rds file for the bed file of your samples.
    # it's important to have this file before running any ldpred or bigsnpr-related scripts
    snp_readBed(bedfile = bed_filepath)
    print(paste0("Done making the rds file from bigsnpr for: ", bed_filepath))
  } else {
    print(paste0("The rds file from bigsnpr for: ", bed_filepath, " was found to be done already"))
  }
  
  ##### loop over these summary stats and calculate the PGS #####
  # registerDoMC(cores=n_cores)
  # foreach(i = 1:nrow(ss.meta)) %dopar% {
  for(i in 1:nrow(ss.meta)) {
    phenotype = ss.meta$phenotype[i]
    print(paste0("Started working on phenotype: ", phenotype, 
                 " from: ", ss.meta$source[i], " done at: ", ss.meta$year[i]))
    # load the summary stats file
    ss <- bigreadr::fread2(ss.meta$full_path[i]) 
    # Filter out hapmap SNPs
    ss.f <- ss %>%
      inner_join(select(bim_snp, chr, pos)) %>% inner_join(select(map_ldref, chr, pos))
    rm(ss);gc()
    ###
    print(paste0("PGS calc mark for: ", phenotype))
    calculate_ldpred_pgs(sumstats = ss.f,
                         device = device,
                         plink_rds = rds.path,
                         sd_y = sd,
                         pheno_name = paste0(ss.meta$source[i], "_", phenotype, "_", ss.meta$year[i]),
                         n_core = n_cores,
                         build = build,
                         output_path = output_directory)
    print(paste0("Done calculating PGS for: ", phenotype))
  }
  if (combine == T) {
    df <- data.frame(file = list.files(output_directory, full.names = T, pattern = "_PGS"))
    all <- foreach(j = 1:nrow(df), .combine = inner_join) %dopar% {
      r <- read_tsv(df$file[j])
      return(r)
    }
    write_rds(all, file = paste0(output_directory, "/ALL-PGS.rds"))
  }
}


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
