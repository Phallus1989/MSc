#DIABLO: Data intergration analyiss for biomarker discovery using latent components
setwd("C:/")

library(knitr)
library(mixOmics)
library(ggplot2)
library(ggsci)
library(lme4)
library(glmnetUtils)
library(limma)
library(edgeR)
library(survival)
library(Glimma)
library(tidyverse)
library(Biobase)

#data normalization transform
relt <- read.csv("n_counts.csv", row.names = 1)
reln <- read.csv("t_counts.csv", row.names = 1)
imt <- read.csv("Tumour_immune.csv", row.names = 1)
imn <- read.csv("Normal_immune.csv", row.names = 1)
ngc <- read.csv("bacn.csv", row.names = 1)
tgc <- read.csv("bact.csv", row.names = 1)

#add +1 offset for log transform
tgcdat <- tgc + 1
ngcdat <- ngc + 1
rnadat <- relt + 1
rnadan <- reln + 1
imdat <- imt + 1
imdan <- imn + 1

#log transfrom CLR (centered log ratio)
gct <- logratio.transfo(as.matrix(tgcdat), logratio = 'CLR', offset = 0)
gcn <- logratio.transfo(as.matrix(ngcdat), logratio = 'CLR', offset = 0)
rmt <- logratio.transfo(as.matrix(rnadat), logratio = 'CLR', offset = 0)
rmn <- logratio.transfo(as.matrix(rnadan), logratio = 'CLR', offset = 0)
imt <- logratio.transfo(as.matrix(imdat), logratio = 'CLR', offset = 0)
imn <- logratio.transfo(as.matrix(imdan), logratio = 'CLR', offset = 0)

#check
Meta <- read.csv("Rec_meta.csv", row.names = 1)
nimm <- imn
timm <- imt
nrna <- t(rmn)
trna <- t(rmt)
ntax <- gcn
ttax <- gct


#All
x <- list(NTaxa = ntax, TTaxa = ttax, NImmune = nimm, TImmune = timm, NGenes =  nrna, TGenes = trna)
dim(imn)
dim(imt)
dim(nrna)
dim(trna)
dim(ntax)
dim(ttax)
summary(x)

#running diablo
y <- Meta$Cohort
Diablo_rectum <- block.splsda(x, y)
indvplot <- plotIndiv(Diablo_rectum, 
                      ind.names = TRUE, 
                      legend=TRUE
                      )

#Removing batch effect
#design
c.et <- Meta$Cohort
t.et <- Meta$Dworak

#treatment effects for preservation
rmt.mod <- model.matrix(~0 + t.et)

#removing batch effect
rmt.fix <- t(removeBatchEffect(t(ntax), batch = c.et, 
                               design = rmt.mod))
rmn.fix <- t(removeBatchEffect(t(ttax), batch = c.et, 
                               design = rmt.mod))
rit.fix <- t(removeBatchEffect(t(timm), batch = c.et, 
                               design = rmt.mod))
rin.fix <- t(removeBatchEffect(t(nimm), batch = c.et, 
                               design = rmt.mod))
rgt.fix <- t(removeBatchEffect(t(trna), batch = c.et, 
                               design = rmt.mod))
rgn.fix <- t(removeBatchEffect(t(nrna), batch = c.et, 
                               design = rmt.mod))


#check
Meta <- read.csv("Rec_meta.csv", row.names = 1)
nim <- rin.fix
tim <- rit.fix
ntax <- rmn.fix
ttax <- rmt.fix
tgc <- rgt.fix
ngc <- rgn.fix

#all
x <- list(NTaxa = ntax, TTaxa = ttax, NImmune = nim, TImmune = tim, NGenes =  ngc, TGenes = tgc)

#Checking batch removal
y <- Meta$Cohort
Diablo_rectum <- block.splsda(x, y, ncomp = 6)
indvplot <- plotIndiv(Diablo_rectum, 
                      ind.names = TRUE, 
                      legend=TRUE)
#raw result
y <- Meta$CIP
Diablo_rectum <- block.splsda(x, y, ncomp = 6)
cimDiablo(Diablo_rectum,
          margins = c(15, 15),
          row.names = TRUE,
          col.names = TRUE,
          size.legend = 0.8,
          legend.position = "right")


set.seed(123)
MyPerf.diablo <- perf(Diablo_rectum, validation = 'loo', folds = 5, progressBar = TRUE,
                      dist = 'centroids.dist')
MyPerf.diablo$WeightedVote.error.rate$centroids.dist
#plot(MyPerf.diablo)

#$X.test <- list(NormGene =  ngc, TumGene = tgc)
#X.test <- list(NormImm = nImm, TumImm = tImm)
#X.test <- list(NormTax = ntax, TumourTax = ttax)
X.test <- list(NTaxa = ntax, TTaxa = ttax, NImmune = nim, TImmune = tim, NGenes =  ngc, TGenes = tgc)
Mypredict.diablo <- predict(Diablo_rectum, newdata = X.test)
confusion.mat <- get.confusion_matrix(
  truth = Meta$CIP, 
  predicted = Mypredict.diablo$WeightedVote$centroids.dist[,2])
confusion.mat

#
#variable selection
y <- Meta$CIP
#Normal tax
x <- ntax
res <- spls(x,y, mode = "regression", ncomp = 4)
var <-selectVar(res, comp =4)
write.csv(var$X, "ntax_filtd.csv")
perf.pls <- perf(res, validation = "Mfold", folds = 5,
                 progressBar = TRUE, nrepeat = 50)
plot(perf.pls$Q2.total)
abline(h = 0.0975)
#Tumour tax
x <- ttax
res <- spls(x,y, mode = "regression", ncomp = 4)
var <-selectVar(res, comp =4)
write.csv(var$X, "ttax_filtd.csv")
perf.pls <- perf(res, validation = "Mfold", folds = 5,
                 progressBar = TRUE, nrepeat = 50)
plot(perf.pls$Q2.total)
abline(h = 0.0975)
#Normal immune
x <- nim
res <- spls(x,y, mode = "regression", ncomp = 4)
var <-selectVar(res, comp =4)
write.csv(var$X, "nimm_filtd.csv")
perf.pls <- perf(res, validation = "Mfold", folds = 5,
                 progressBar = TRUE, nrepeat = 50)
plot(perf.pls$Q2.total)
abline(h = 0.0975)

#Tumour immune
x <- tim
res <- spls(x,y, mode = "regression", ncomp = 4)
var <-selectVar(res, comp =4)
write.csv(var$X, "timm_filtd.csv")
perf.pls <- perf(res, validation = "Mfold", folds = 5,
                 progressBar = TRUE, nrepeat = 50)
plot(perf.pls$Q2.total)
abline(h = 0.0975)

#Normal genes
x <- ngc
res <- spls(x,y, mode = "regression", ncomp = 4)
var <-selectVar(res, comp =4)
write.csv(var$X, "ngene_filtd.csv")
perf.pls <- perf(res, validation = "Mfold", folds = 5,
                 progressBar = TRUE, nrepeat = 50)
plot(perf.pls$Q2.total)
abline(h = 0.0975)

#Tumour genes
x <- tgc
res <- spls(x,y, mode = "regression", ncomp = 4)
var <-selectVar(res, comp =4)
write.csv(var$X, "tgene_filtd.csv")
perf.pls <- perf(res, validation = "Mfold", folds = 5,
                 progressBar = TRUE, nrepeat = 50)
plot(perf.pls$Q2.total)
abline(h = 0.0975)

#filtering datasets to most powerful predictors

#bacterial transcription
t_tax <- read.csv(as.matrix("t_tax.csv"))
n_tax <- read.csv(as.matrix("n_tax.csv"))
totax <- (t(ttax))
notax <- (t(ntax))
write.csv(totax, "totax.csv")
write.csv(notax, "notax.csv")
totax <- read.csv("totax.csv")
notax <- read.csv("notax.csv")
newttax <- totax[totax$taxa %in% t_tax$name, ]
dim(newttax)
newntax <- notax[notax$tax %in% n_tax$name, ]
dim(newntax)
write.csv(t(newntax), "mnf_b.csv")
write.csv(t(newttax), "mtf_b.csv")

#immune abundance
t_imm <- read.csv(as.matrix("t_imm.csv"))
n_imm <- read.csv(as.matrix("n_imm.csv"))
toimm <- (t(tim))
noimm <- (t(nim))
write.csv(toimm, "toimm.csv")
write.csv(noimm, "noimm.csv")
toimm <- read.csv("toimm.csv")
noimm <- read.csv("noimm.csv")
newtim <- toimm[toimm$cell %in% t_imm$name, ]
dim(newtim)
newnim <- noimm[noimm$cell %in% n_imm$name, ]
dim(newnim)
write.csv(t(newnim), "inf_b.csv")
write.csv(t(newtim), "itf_b.csv")

#gene expression
t_gc <- read.csv(as.matrix("t_genes.csv"))
n_gc <- read.csv(as.matrix("n_genes.csv"))
togc <- (t(tgc))
nogc <- (t(ngc))
write.csv(togc, "togc.csv")
write.csv(nogc, "nogc.csv")
togc <- read.csv("togc.csv")
nogc <- read.csv("nogc.csv")
newtgc <- togc[togc$gene %in% t_gc$name, ]
dim(newtgc)
newngc <- nogc[nogc$gene %in% n_gc$name, ]
dim(newngc)
write.csv(t(newngc), "gnf_b.csv")
write.csv(t(newtgc), "gtf_b.csv")

#applying blocksplda with refined variable set
ttax2 <- read.csv("mtf_b.csv", row.names = 1)
ntax2 <- read.csv("mnf_b.csv", row.names = 1)
tim2 <- read.csv("itf_b.csv", row.names = 1)
nim2 <- read.csv("inf_b.csv", row.names = 1)
tgc2 <- read.csv("gtf_b.csv", row.names = 1)
ngc2 <- read.csv("gnf_b.csv", row.names = 1)

#all
x <- list(NTaxa = ntax2, TTaxa = ttax2, NImmune = nim2, TImmune = tim2, NGenes =  ngc2, TGenes = tgc2)
#Running diablo
y <- Meta$CIP
Diablo_rectum <- block.splsda(x, y, ncomp = 13)
indvplot <- plotIndiv(Diablo_rectum, 
                      ind.names = TRUE, 
                      legend=TRUE,
                      ellipse = TRUE)
cimDiablo(Diablo_rectum,
          margins = c(15, 15),
          row.names = TRUE,
          col.names = FALSE,
          size.legend = 0.8,
          legend.position = "right")
set.seed(123)
MyPerf.diablo <- perf(Diablo_rectum, validation = 'loo', folds = 5, progressBar = TRUE,
                      dist = 'centroids.dist')
MyPerf.diablo$WeightedVote.error.rate$centroids.dist
plot(MyPerf.diablo)

#breaking into test and train
test <- read.csv(as.matrix("test_meta.csv"))
train <- read.csv(as.matrix("train_meta.csv"))
ttax2 <- read.csv("mtf_cbt.csv")
ntax2 <- read.csv("mnf_cbt.csv")
tim2 <- read.csv("itf_cbt.csv")
nim2 <- read.csv("inf_cbt.csv")
tgc2 <- read.csv("tumour_test.csv")
ngc2 <- read.csv("normal_test.csv")

#gene expression
test_tgc <- tgc2[tgc2$gene %in% test$Tumour_number, ]
dim(test_tgc)
train_tgc <- tgc2[tgc2$gene %in% train$ï..Tumour_number, ]
dim(train_tgc)

test_ngc <- ngc2[tgc2$gene %in% test$Tumour_number, ]
dim(test_ngc)
train_ngc <- ngc2[ngc2$gene %in% train$ï..Tumour_number, ]
dim(train_ngc)

#immune
test_tic <- tim2[tim2$cell %in% test$Tumour_number, ]
dim(test_tic)
train_tic <- tim2[tim2$cell %in% train$ï..Tumour_number, ]
dim(train_tic)
test_nic <- nim2[nim2$cell %in% test$Tumour_number, ]
dim(test_nic)
train_nic <- nim2[nim2$cell %in% train$ï..Tumour_number, ]
dim(train_nic)

#taxa
test_ttc <- ttax2[ttax2$taxa %in% test$Tumour_number, ]
dim(test_ttc)
train_ttc <- ttax2[ttax2$taxa %in% train$ï..Tumour_number, ]
dim(train_ttc)
test_ntc <- ntax2[ntax2$taxa %in% test$Tumour_number, ]
dim(test_ntc)
train_ntc <- ntax2[ntax2$taxa %in% train$ï..Tumour_number, ]
dim(train_ntc)
#fixing data
library(textshape)
library(tidyverse)

test <- test[order(test$Tumour_number),]
train <- train[order(train$ï..Tumour_number),]
test_ngc <- test_ngc[order(test_ngc$gene),]
test_tgc <- test_tgc[order(test_tgc$gene),]
test_ntc <- test_ntc[order(test_ntc$taxa),]
test_ttc <- test_ttc[order(test_ttc$taxa),]
test_tic <- test_tic[order(test_tic$cell),]
test_nic <- test_nic[order(test_nic$cell),]
train_tgc <- train_tgc[order(train_tgc$gene),]
train_ngc <- train_ngc[order(train_ngc$gene),]
train_ntc <- train_ntc[order(train_ntc$taxa),]
train_ttc <- train_ttc[order(train_ttc$taxa),]
train_tic <- train_tic[order(train_tic$cell),]
train_nic <- train_nic[order(train_nic$cell),]
library(textshape)
write.csv(test, "text.csv")
write.csv(train, "train.csv")
write.csv(test_ngc , "test_ngc.csv")
write.csv(test_tgc , "test_tgc.csv")
write.csv(test_ntc , "test_ntc.csv")
write.csv(test_ttc , "test_ttc.csv")
write.csv(test_tic , "test_tic.csv")
write.csv(test_nic , "test_nic.csv")
write.csv(train_ngc, "train_ngc.csv")
write.csv(train_tgc, "train_tgc.csv")
write.csv(train_ntc, "train_ntc.csv")
write.csv(train_ttc, "train_ttc.csv")
write.csv(train_tic, "train_tic.csv")
write.csv(train_nic, "train_nic.csv")

test <- read.csv("text.csv", row.names = 1)
train <- read.csv("train.csv", row.names = 1)
test_ngc  <- read.csv("test_ngc.csv", row.names = 1)
test_tgc  <- read.csv("test_tgc.csv", row.names = 1)
test_ntc  <- read.csv("test_ntc.csv", row.names = 1)
test_ttc  <- read.csv("test_ttc.csv", row.names = 1)
test_tic  <- read.csv("test_tic.csv", row.names = 1)
test_nic  <- read.csv("test_nic.csv", row.names = 1)
train_ngc <- read.csv("train_ngc.csv", row.names = 1)
train_tgc <- read.csv("train_tgc.csv", row.names = 1)
train_ntc <- read.csv("train_ntc.csv", row.names = 1)
train_ttc <- read.csv("train_ttc.csv", row.names = 1)
train_tic <- read.csv("train_tic.csv", row.names = 1)
train_nic <- read.csv("train_nic.csv", row.names = 1)


#Running DIABLO
x <- list(NTaxa = train_ntc, TTaxa = train_ttc, NImmune = train_nic, TImmune = train_tic, NGenes =  train_ngc, TGenes = train_tgc)
y <- train$CIP
summary(x)

Diablo_rectum <- block.splsda(x, y, ncomp = 5)
indvplot <- plotIndiv(Diablo_rectum, 
                      ind.names = TRUE, 
                      legend=TRUE,
                      ellipse = TRUE)
cimDiablo(Diablo_rectum,
          margins = c(15, 15),
          row.names = TRUE,
          col.names = TRUE,
          size.legend = 0.8,
          legend.position = "right")
set.seed(123)

MyPerf.diablo <- perf(Diablo_rectum, validation = 'loo', folds = 5,  progressBar = TRUE,
                      dist = 'centroids.dist')
MyPerf.diablo$WeightedVote.error.rate$centroids.dist
plot(MyPerf.diablo)

x.test <- list(NTaxa = test_ntc, TTaxa = test_ttc, NImmune = test_nic, TImmune = test_tic, NGenes =  test_ngc, TGenes = test_tgc)
Mypredict.diablo <- predict(Diablo_rectum, newdata = x.test)
confusion.mat <- get.confusion_matrix(
  truth = test$CIP, 
  predicted = Mypredict.diablo$WeightedVote$centroids.dist[,5])
confusion.mat


auroc(Diablo_rectum, roc.block = "TTaxa", roc.comp = 5)
auroc(Diablo_rectum, roc.block = "NTaxa", roc.comp = 5)
auroc(Diablo_rectum, roc.block = "TImmune", roc.comp = 5)
auroc(Diablo_rectum, roc.block = "NImmune", roc.comp = 5)
auroc(Diablo_rectum, roc.block = "TGenes", roc.comp = 5)
auroc(Diablo_rectum, roc.block = "NGenes", roc.comp = 5)


plotLoadings(Diablo_rectum, comp =5, contrib = "max", block = "NGenes", method = "mean",
             size.name = 0.8, xlim = c(-1, 1), size.title = 1, size.legend = 1)

plotLoadings(Diablo_rectum, comp =5, contrib = "max", block = "TGenes", method = "mean",
             size.name = 0.8, xlim = c(-1, 1), size.title = 1, size.legend = 1)
plotLoadings(Diablo_rectum, comp =5, contrib = "max", block = "TImmune", method = "mean",
             size.name = 0.8, xlim = c(-1, 1), size.title = 1, size.legend = 1)


Diablo_rectum$loadings$TImmune


