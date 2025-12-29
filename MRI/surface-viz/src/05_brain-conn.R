source("/wdata/msmuhammad/msmuhammad-source.R")
library(brainconn)

custom_atlas <- read.csv("/wdata/msmuhammad/refs/Schaefer2018/Parcellations/MNI/Centroid_coordinates/Schaefer2018_100Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv")%>%
  select(ROI.Name=2, x.mni=3,y.mni=4,z.mni=5) %>% mutate(ROI.Name=sub("17Networks","",ROI.Name))
check_atlas(custom_atlas)

brainconn(atlas = custom_atlas, conmat = example_unweighted_undirected,node.color = "black")