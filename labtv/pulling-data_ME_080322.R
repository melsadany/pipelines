####################################################################
#      pulling recent publications from pubmed using keywords      #
####################################################################
#                       Muhammad Elsadany                          #
####################################################################

# I used "test2" conda env

library(easyPubMed)
library(RISmed)

setwd("/Dedicated/jmichaelson-wdata/msmuhammad/extra/labtv/")

# might need it. keep it here in case. That's mine
api.key <- "fdf77f5b0f5aa282c19167f0d91bf649de08"

# keywords <- c("psychiatry", "genetics", "neuroscience", "autism", "bipolar", "computational",
              # "deep learning", "bioinformatics", "neuroimaging", "GWAS", "transcriptome",
              # "genome", "biobank", "ABCD", "human brain", "genome-wide", "transgender")
keywords <- c("computational AND psychiatry", "genetics AND neuroimaging", "psychiatry AND statistics")

my.query <- keywords[1]
for (i in 2:length(keywords)) {
  my.query <- paste0(my.query, " OR ",keywords[i])
}

my.query

entrez.id <- get_pubmed_ids(my.query, api_key = api.key)
# abstracts.txt <- fetch_pubmed_data(entrez.id, format = "abstract")
# head(abstracts.txt)
abstracts.xml <- fetch_pubmed_data(pubmed_id_list = entrez.id)
# class(abstracts.xml)

titles <- custom_grep(abstracts.xml, "ArticleTitle", "char")
# head(titles)
TTM <- nchar(titles) > 75
titles[TTM] <- paste(substr(titles[TTM], 1, 70), "...", sep = "")
# head(titles)



# extract affiliation data 
PM.list <- articles_to_list(pubmed_data = abstracts.xml)
# class(PM.list[1])

# substr(PM.list[4], 1, 510)


xx <- lapply(PM.list, article_to_df, autofill = TRUE, max_chars = -1)
full.df <- do.call(rbind, xx)

# next step is to filter the df of publications list you have, because any 
# pub title will be repeated n times (n is number of authors on it)

refined.df <- unique(full.df[which(nchar(full.df$title) > 30 & nchar(full.df$abstract) > 500), c("title", "jabbrv", "month", "year", "abstract", "pmid")]) 
refined.df.2 <- na.omit(refined.df)
over.abstracts <- nchar(refined.df.2$abstract) > 500
refined.df.2$abstract[over.abstracts] <- paste(substr(refined.df.2$abstract, 1, 500), "...", sep = "")
refined.df.2$pmid <- as.numeric(refined.df.2$pmid)



# connecting to googlesheets and uploading the datframe
library(googlesheets4)
library(googlesheets)

# gs4_auth()

sheet.URL <- "https://docs.google.com/spreadsheets/d/1wlCxgHltOvwaBRYB_uwzOpJs_mIp-EOwBmdmAawb79E/edit#gid=0"
sheet <- read_sheet("https://docs.google.com/spreadsheets/d/1wlCxgHltOvwaBRYB_uwzOpJs_mIp-EOwBmdmAawb79E/edit#gid=0")
sheet_write(refined.df.2, ss = sheet.URL, sheet = "Data")


