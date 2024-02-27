################################################################################
#                                     DL lecture 1                             #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(keras)
library(reticulate)
library(tensorflow)
use_condaenv("EDL")
#################################################################################

project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/extra/DL-lab"
setwd(project.dir)


# setwd("/Dedicated/jmichaelson-wdata/jmichaelson/NN/course/DL_in_biomedicine/session_1")
source("/Dedicated/jmichaelson-wdata/jmichaelson/NN/AA/AA_NN.R"); k_clear_session(); gc()
source("/Dedicated/jmichaelson-wdata/jmichaelson/NN/course/DL_in_biomedicine/utils/functions.R")

pload("perturbation_matrix.Rdata.pxz") ## add details on provenance

### this is what the data look like
img(pert[1:100,1:100])

### read in the DD evidence list (from https://pubmed.ncbi.nlm.nih.gov/33932580/)
dd = read.table("DD_genes.txt",sep="\t",header=T,stringsAsFactors=F)
dd0 = dd[,5:ncol(dd)]=="True"
rownames(dd0) = dd[,2]
mode(dd0) = "integer"

### what genes are common between the perturbation data and the DD data?
nn = rownames(dd0)[rownames(dd0)%in%rownames(pert)]

### arrange X and Y in the same row order
X = pert[nn,]
range(X)
X[is.infinite(X)] = 125 ## check your data to make sure there aren't NAs, etc.
X[X>5] = 5
X = scale(X)

Y = dd0[nn,]
range(Y)

###########################################################################
### subset data into train/test sets
###########################################################################
set.seed(8762)
tst = sample(nrow(X),1000)

Xtr = X[-tst,]
Xts = X[tst,]

Ytr = Y[-tst,]
Yts = Y[tst,]


###########################################################################
### input layer
###########################################################################
rm(ipt,hidden,out,m); gc(); k_clear_session()
ipt = layer_input(ncol(X),name="input")


###########################################################################
### intermediate layers
###########################################################################
hidden = ipt %>%
  layer_dense(64,activation="relu") %>%
  layer_dense(32,activation="relu") %>%
  #layer_batch_normalization() %>%
  layer_dropout(0.5) %>%
  layer_dense(2,activation="linear") #%>%
#layer_gaussian_dropout(0.1)

###########################################################################
### output layer
###########################################################################
out = hidden %>%
  layer_dense(4,activation="relu") %>%
  layer_dense(8,activation="relu") %>%
  layer_dense(ncol(Y),activation="sigmoid",name="out")


###########################################################################
### compile the model
###########################################################################
m = keras_model(ipt,out)

m %>% compile(
  optimizer = optimizer_adam(1e-4),
  loss = "binary_crossentropy",
  metrics = "mae",
)

###########################################################################
### look at the model's structure
###########################################################################
summary(m)

###########################################################################
### train the model
###########################################################################
history = m %>% keras::fit(
  x = Xtr,
  y = Ytr,
  epochs = 50, 
  batch_size = 64,
  validation_data=list(Xts,Yts),
  verbose=T,
)


###########################################################################
###########################################################################
### view performance over the training process
###########################################################################
plot(history)


###########################################################################
### use the model to make predictions
###########################################################################

### training set predictions
prd = predict(m,Xtr,batch_size=128)
img(cor(prd,Ytr))
plot(jitter(Ytr[,3]),prd[,3],xlab="truth",ylab="prediction")

### test set predictions
prd = predict(m,Xts,batch_size=128)
img(cor(prd,Yts))
plot(jitter(Yts[,3]),prd[,3],xlab="truth",ylab="prediction")

### all predictions
prd = predict(m,X,batch_size=128)
rownames(prd) = nn
g = rowMeans(prd)
head(sort(g,dec=T),20)
plot(g,rowSums(Y))


### access the latent space (bottleneck) and take a look
m2 = keras_model(ipt,hidden)
prd = predict(m2,X,batch_size=128)
rownames(prd) = nn

### are there outliers in the latent space? would removing them
### from training improve things?
hplot = function(x){
  plot(x,col=ifelse(Y[,3]==1,'red','grey'),
       xlab="latent dim 1",ylab="latent dim 2")
  OL = rowSums(apply(x,2,ol,6))>0
  text(x[OL,1],x[OL,2],rownames(x)[OL],cex=0.6,col='blue',pos=4)
}

hplot(prd)

### mark outliers
OL = rownames(prd)[rowSums(apply(prd,2,ol,6))>0]
img(X[OL,])


###########################################################################
### tweaks and discussion questions:
### * add units and/or hidden layers (i.e., "capacity") to the model until
### the training data is fit perfectly in 25 epochs
### * how can we "regularize" the model using only number of layers
###	and number of units?
### * how does changing the batch size and learning rate affect
###	training and performance?
### * inspect the output at the penultimate layer - how do those 
###	intermediate outputs correlate with the final predictions?
### * access the kernel weights of the model; which weights might
###	be most informative? what would you want to pay attention to?
### * kernel constraints as a way to enforce certain "desirable" 
###     properties
###########################################################################

















