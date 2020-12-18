#Immune cell analysis and correlations with response with response and bacterial transcription

#William Taylor 18/12/2020
#wil_tr@live.co.uk

setwd("C:/")

library(ltm)
library(corrplot)
library(textshape)
library(tidyverse)
library(ggsci)
library(rstatix)
library(psych)
library(ggpubr)
library(ggplot2)
library(lemon)


#signifance of differences
norm <- read.csv("immune_tumo.csv", row.names = 1)
tumo <- read.csv("immune_norm.csv", row.names = 1)
imm <- read.csv("immune_box.csv", row.names = 1)


d4t <- read.csv("d4t.csv", row.names = 1)
d4n <- read.csv("d4n.csv", row.names = 1)
d3t <- read.csv("d3t.csv", row.names = 1)
d3n <- read.csv("d3n.csv", row.names = 1)
d2t <- read.csv("d2t.csv", row.names = 1)
d2n <- read.csv("d2n.csv", row.names = 1)
d1t <- read.csv("d1t.csv", row.names = 1)
d1n <- read.csv("d1n.csv", row.names = 1)

ht <- read.csv("hight.csv", row.names = 1)
hn <- read.csv("highn.csv", row.names = 1)
lt <- read.csv("lowt.csv", row.names = 1)
ln <- read.csv("lown.csv", row.names = 1)

ncrn <- read.csv("ncrn.csv", row.names = 1)
ncrt <- read.csv("ncrt.csv", row.names = 1)

#normality test
lapply(d4t, shapiro.test)
lapply(d3t, shapiro.test)
lapply(d2t, shapiro.test)
lapply(d1t, shapiro.test)
lapply(d4n, shapiro.test)
lapply(d3n, shapiro.test)
lapply(d2n, shapiro.test)
lapply(d1n, shapiro.test)

#non-parametric recommended

library(EMA)
library(data.table)


#tumour vs normal
out <- lapply(1:22, function(x) pairwise.wilcox.test(imm[[x]], imm$Tissue, p.adjust = 'fdr', paired = T))
names(out) <- names(imm)[1:22]
out
outp <- sapply(out, function(x) {
  p <- x$p.value
  n <- outer(rownames(p), colnames(p), paste, sep='v')
  p <- as.vector(p)
  names(p) <- n
  p
})

write.csv(outp, "outp.csv")

t.test_results <- mapply(wilcox.test, x =norm, y =tumo, paired = T, SIMPLIFY = F, conf.int = T, conf.level = 0.95)
write.matrix(t.test_results, "tvn_ttest.csv")

outp <- sapply(t.test_results, function(x) {
  p <- x$p.value
  n <- outer(rownames(p), colnames(p), paste, sep='v')
  p <- as.vector(p)
  names(p) <- n
  p
})
outp2 <- p.adjust(outp,method='fdr')
write.csv(outp2, "outp.csv")
t.test_results

#plot
sigtvn <- read.csv("sig_tvn.csv")
sigtvn$ï..Cell.type <- factor(sigtvn$ï..Cell.type, levels=unique(as.character(sigtvn$ï..Cell.type)))

fd_plot <- ggplot(sigtvn,aes(x=ï..Cell.type,y=Fold.Difference,fill=Fold.Difference>0))+
  geom_col() + coord_flip() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(),
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "right",) +
  xlab(NULL) + ylab("Greater in incomplete") +
  scale_y_continuous(breaks = seq(-5, 9, 1))

  
fd <- fd_plot + scale_fill_npg()
fd

#method
#use full comparison list to generate significant associations, then include only those with a p > 0.1 in the published figure

imb <- read.csv("immune_profile_box.csv")

imb$Tissue <- factor(imb$Tissue, levels = c("Tumour", "Normal"))
comp <- list(c("Tumour", "Normal"))
TN <- imb$Tumour_number


bcm <- ggboxplot(imb, x = "Tissue", y = "B.cells.memory", color = "Tissue", 
                 add = "jitter", legend = "right") +stat_compare_means(comparisons = comp, label = "p.signif", paired = TRUE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())
bcm <- bcm + ggtitle("a) Memory B cells")

cd8 <- ggboxplot(imb, x = "Tissue", y = "T.cells.CD8", color = "Tissue", 
                add = "jitter", legend = "none") +
  stat_compare_means(comparisons = comp, label = "p.signif", paired = TRUE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())
cd8 <- cd8 + ggtitle("c) CD8 T Cells")

pc <- ggboxplot(imb, x = "Tissue", y = "Plasma.cells", color = "Tissue", 
                 add = "jitter", legend = "none") +
  stat_compare_means(comparisons = comp, label = "p.signif", paired = TRUE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())
pc <- pc + ggtitle("b) Plasma cells")

mca <- ggboxplot(imb, x = "Tissue", y = "Mast.cells.activated", color = "Tissue", 
                add = "jitter", legend = "none") +
  stat_compare_means(comparisons = comp, label = "p.signif", paired = TRUE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())
mca <- mca + ggtitle("d) Activated mast cells")

m0 <- ggboxplot(imb, x = "Tissue", y = "Macrophages.M0", color = "Tissue", 
                add = "jitter", legend = "none") +
  stat_compare_means(comparisons = comp, label = "p.signif", paired = TRUE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())
m0 <- m0 + ggtitle("e) M0 macrophages")

nc <- ggboxplot(imb, x = "Tissue", y = "Neutrophils", color = "Tissue", 
                add = "jitter", legend = "none") +
  stat_compare_means(comparisons = comp, label = "p.signif", paired = TRUE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())
nc <- nc + ggtitle("f) Neutrophils")

m1 <- ggboxplot(imb, x = "Tissue", y = "Macrophages.M1", color = "Tissue", 
                   add = "jitter", legend = "none") +
  stat_compare_means(comparisons = comp, label = "p.signif", paired = TRUE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())
m1 <- m1 + ggtitle("g) M1 macrophages")

ggarrange(bcm, pc, cd8, mca, m0, nc, m1, nrow = 2, ncol = 4, legend = "right")

legend_b <- get_legend(bcm + theme(legend.position="right"))
plot2 <- ggarrange(
  bcm + theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank(),),
  pc + theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank(),),
  cd8 + theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank(),),
  mca + theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank(),),
  m0 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  nc + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  m1 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  legend_b,
  align = "hv", ncol = 4, nrow = 2)
plot2
dca <- annotate_figure(plot2, left = "Estimated abundance")
dca




#D4
t.test_results <- mapply(wilcox.test, x= d4t, y = d4n, paired = T, SIMPLIFY = F, conf.int = T, conf.level = 0.95)
t.test_results
write.matrix(t.test_results, "d4_ttest.csv")
outp <- sapply(t.test_results, function(x) {
  p <- x$p.value
  n <- outer(rownames(p), colnames(p), paste, sep='v')
  p <- as.vector(p)
  names(p) <- n
  p
})
outp2 <- p.adjust(outp,method='fdr')
write.csv(outp2, "d4_outp.csv")

#D3
t.test_results <- mapply(wilcox.test, x= d3t, y = d3n, paired = T, SIMPLIFY = F, conf.int = T, conf.level = 0.95)
write.matrix(t.test_results, "d3_ttest.csv")
outp <- sapply(t.test_results, function(x) {
  p <- x$p.value
  n <- outer(rownames(p), colnames(p), paste, sep='v')
  p <- as.vector(p)
  names(p) <- n
  p
})
outp2 <- p.adjust(outp,method='fdr')
write.csv(outp2, "d3_outp.csv")

#D2
t.test_results <- mapply(wilcox.test, x= d2t, y = d2n, paired = T, SIMPLIFY = F, conf.int = T, conf.level = 0.95)
write.matrix(t.test_results, "d2_ttest.csv")
outp <- sapply(t.test_results, function(x) {
  p <- x$p.value
  n <- outer(rownames(p), colnames(p), paste, sep='v')
  p <- as.vector(p)
  names(p) <- n
  p
})
outp2 <- p.adjust(outp,method='fdr')
write.csv(outp2, "d2_outp.csv")

#D1
t.test_results <- mapply(wilcox.test, x= d1t, y = d1n, paired = T, SIMPLIFY = F, conf.int = T, conf.level = 0.95)
write.matrix(t.test_results, "d1_ttest.csv")
outp <- sapply(t.test_results, function(x) {
  p <- x$p.value
  n <- outer(rownames(p), colnames(p), paste, sep='v')
  p <- as.vector(p)
  names(p) <- n
  p
})
outp2 <- p.adjust(outp,method='fdr')
write.csv(outp2, "d1_outp.csv")

#high
t.test_results <- mapply(wilcox.test, x = ht, y = hn, paired = T, SIMPLIFY = F, conf.int = T, conf.level = 0.95)
write.matrix(t.test_results, "h_ttest.csv")
outp <- sapply(t.test_results, function(x) {
  p <- x$p.value
  n <- outer(rownames(p), colnames(p), paste, sep='v')
  p <- as.vector(p)
  names(p) <- n
  p
})
outp2 <- p.adjust(outp,method='fdr')
write.csv(outp2, "h_outp.csv")

#low
t.test_results <- mapply(wilcox.test, x = lt, y = ln, paired = T, SIMPLIFY = F, conf.int = T, conf.level = 0.95)
write.matrix(t.test_results, "l_ttest.csv")
outp <- sapply(t.test_results, function(x) {
  p <- x$p.value
  n <- outer(rownames(p), colnames(p), paste, sep='v')
  p <- as.vector(p)
  names(p) <- n
  p
})
outp2 <- p.adjust(outp,method='fdr')
write.csv(outp2, "l_outp.csv")

#high v low t
t.test_results <- mapply(wilcox.test, x = ht, y = lt, SIMPLIFY = F, conf.int = T, conf.level = 0.95)
write.matrix(t.test_results, "hlt_ttest.csv")
outp <- sapply(t.test_results, function(x) {
  p <- x$p.value
  n <- outer(rownames(p), colnames(p), paste, sep='v')
  p <- as.vector(p)
  names(p) <- n
  p
})
outp2 <- p.adjust(outp,method='fdr')
write.csv(outp2, "hlt_outp.csv")

#high v low n
t.test_results <- mapply(wilcox.test, x = hn, y = ln, SIMPLIFY = F, conf.int = T, conf.level = 0.95)
write.matrix(t.test_results, "hln_ttest.csv")
outp <- sapply(t.test_results, function(x) {
  p <- x$p.value
  n <- outer(rownames(p), colnames(p), paste, sep='v')
  p <- as.vector(p)
  names(p) <- n
  p
})
outp2 <- p.adjust(outp,method='fdr')
write.csv(outp2, "hln_outp.csv")

#D4 v NCR t
t.test_results <- mapply(wilcox.test, x = d4t, y = ncrt, SIMPLIFY = F, conf.int = T, conf.level = 0.95)
write.matrix(t.test_results, "cr_ncrt_ttest.csv")
outp <- sapply(t.test_results, function(x) {
  p <- x$p.value
  n <- outer(rownames(p), colnames(p), paste, sep='v')
  p <- as.vector(p)
  names(p) <- n
  p
})
outp2 <- p.adjust(outp,method='fdr')
write.csv(outp2, "cr_ncrt_outp.csv")

#D4 v NCR n
t.test_results <- mapply(wilcox.test, x = d4n, y = ncrn, SIMPLIFY = F, conf.int = T, conf.level = 0.95)
write.matrix(t.test_results, "cr_ncrn_ttest.csv")
outp <- sapply(t.test_results, function(x) {
  p <- x$p.value
  n <- outer(rownames(p), colnames(p), paste, sep='v')
  p <- as.vector(p)
  names(p) <- n
  p
})
outp2 <- p.adjust(outp,method='fdr')
write.csv(outp2, "cr_ncrn_outp.csv")

#plot
signcr <- read.csv("cr_ncr_sig.csv")

signcr$ï..Cell.Type <- factor(signcr$ï..Cell.Type, levels=unique(as.character(signcr$ï..Cell.Type)))

fd_plot <- ggplot(signcr,aes(x=ï..Cell.Type,y=FC.Complete,fill=FC.Complete>0))+
  geom_col() + coord_flip() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  xlab(NULL) + ylab("Fold Difference") +
  scale_y_continuous(breaks = seq(-6, 35, 5))
fd_plot

fd <- fd_plot + scale_fill_npg(name = NULL, labels=NULL)
fd

#method
#use full comparison list to generate significant associations, then include only those with a p > 0.1 in the published figure

imb <- read.csv("ncr_profile_box.csv")

imb$Response <- factor(imb$Response, levels = c("Complete", "Incomplete"))
comp <- list(c("Complete", "Incomplete"))
TN <- imb$Tumour_number

mca <- ggboxplot(imb, x = "Response", y = "Mast.cells.activated", color = "Response", 
                 add = "jitter", legend = "bottom") +
  stat_compare_means(comparisons = comp, label = "p.signif", method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())
mca <- mca + ggtitle("a) Activated mast cells")

mcr <- ggboxplot(imb, x = "Response", y = "Mast.cells.resting", color = "Response", 
                 add = "jitter", legend = "none") +
  stat_compare_means(comparisons = comp, label = "p.signif", method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())
mcr <- mcr + ggtitle("b) Resting mast cells")

treg <- ggboxplot(imb, x = "Response", y = "T.cells.regulatory..Tregs.", color = "Response", 
                add = "jitter", legend = "none") +
  stat_compare_means(comparisons = comp, label = "p.signif",  method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())
treg <- treg + ggtitle("c) Regulatory T cells")

m1 <- ggboxplot(imb, x = "Response", y = "Macrophages.M1", color = "Response", 
                add = "jitter", legend = "none") +
  stat_compare_means(comparisons = comp, label = "p.signif", method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())
m1 <- m1 + ggtitle("d) M1 macrophages")


ggarrange(mca,mcr,treg, m1, legend = TRUE)

legend_b <- get_legend(mca + theme(legend.position="bottom"))
plot2 <- ggarrange(
  mca + theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank(),),
  mcr + theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank(),),
  treg + theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank(),),
  m1 + theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank(),),
  align = "v", ncol = 2, nrow = 2, legend.grob = legend_b , legend = "bottom")
plot2
dca <- annotate_figure(plot2, left = "Estimated abundance")
dca




#correlations
tumo <- read.csv("Tumour_t.csv", row.names = 1)
norm <- read.csv("normal_t.csv", row.names = 1)

#use Spearman correlation to establish relationship between estimated immune cell presence and response
tcor <- cor(tumo, method = "spearman")
ncor <- cor(norm, method = "spearman")

RT <- tumo$Dworak

#function for pvalue estimation
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
tumotest <- cor.mtest(tumo, method = "spearman")
normtest <- cor.mtest(norm, method = "spearman")

pAdj <- p.adjust(tumotest, method = "fdr")
resAdj <- matrix(pAdj, ncol = dim(tumotest)[1])
colnames(resAdj) <- colnames(tumotest)
rownames(resAdj) <- rownames(tumotest)
corrplot(tcor, type="upper", method = "number",
         p.mat = resAdj, sig.level = 0.05, insig = "blank", tl.col = "black", diag = FALSE, addCoef.col = "black", main = "Tumour Tissue")
write.csv(tcor, "tcord.csv")
write.csv(tumotest, "tcordp.csv")
write.csv(resAdj, "tcorp.csv")


pAdj <- p.adjust(normtest, method = "fdr")
resAdj <- matrix(pAdj, ncol = dim(ncor)[1])
colnames(resAdj) <- colnames(ncor)
rownames(resAdj) <- rownames(ncor)

corrplot(ncor, type="upper",  method = "number",
         p.mat = resAdj, sig.level = 0.05, insig = "blank", tl.col = "black", diag = FALSE, addCoef.col = "black", main = "Normal Tissue")

write.csv(ncor, "ncord.csv")
write.csv(normtest, "ncordp.csv")
write.csv(resAdj, "ncorp.csv")
#corrleation scatter
tumo
p1 <- ggplot(tumo, aes(x=Dworak, y=Mast.cells.resting)) + 
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  geom_smooth(method=lm, level = 0.95) +
  labs(title = "Resting mast cells",
       subtitle = expression("Dworak: r = 0.511, " ~italic("p")~" = 0.0007, FDR = 0.009")) +
  theme(plot.subtitle=element_text(size=8))
p1
p2 <- ggplot(tumo, aes(x=Dworak, y=Macrophages.M1)) + 
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  geom_smooth(method=lm, level = 0.95) +
  labs(title = "M1 Macrophages",
       subtitle = expression("Dworak: r = 0.462, " ~italic("p")~" = 0.0027, FDR = 0.024")) +
  theme(plot.subtitle=element_text(size=8))
p2

scatfig <- ggarrange(p2, p1, labels = c("a)", "b)"))
scatfig <- annotate_figure(scatfig, left = "Immune cell abundance", bottom = "Dworak score")
scatfig

#estimate
library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)

library(estimate)

#test
setwd("C:/Users/William/Google Drive/Purcell studies/PhD/rectal_rna-seq/Txi/Unstranded/quants")
RectalGeneExp <- system.file("extdata", "counting.txt", package ="estimate")
#Intersecting inputdata with 10,412 common genes
filterCommonGenes(input.f = RectalGeneExp, output.f = "Rec_genes.gct", id = "GeneSymbol")
estimateScore("Rec_genes.gct", "rectal_estimate_score.gct", "platform" = c("illumina"))
RectalGeneExp <- system.file("extdata", "counts.txt", package="estimate")
out.file <- tempfile(pattern="estimate", fileext=".gct")
estimateScore(RectalGeneExp, out.file)

plotPurity(RectalGeneExp, output.dir = "C:/", platform = "illumina")
est <- read.csv("est_score.csv", row.names = 1)
c <- cor(est, method = "spearman", )

# matrix of the p-value of the correlation
hltest <- corr.test(est, method = 'spearman', adjust = "fdr", ci = TRUE, alpha = 0.05)
print(hltest, short = FALSE)
#signifance of differences
norm <- read.csv("normal_cells_comp.csv", row.names = 1)
tumo <- read.csv("tumour_cells_comp.csv", row.names = 1)
normC <- read.csv("complete_normal.csv", row.names = 1)
normIN <- read.csv("incomplete_normal.csv", row.names = 1)
tumC <- read.csv("complete_tumour.csv", row.names = 1)
tumIN <- read.csv("incomplete_tumour.csv", row.names = 1)
normC <- read.csv("complete_normal.csv", row.names = 1)
normIN <- read.csv("incomplete_normal.csv", row.names = 1)

#tumour vs normal
t.test_results <- mapply(t.test, x= norm, y = tumo, paired = T, SIMPLIFY = F)
write.matrix(t.test_results, "tvn_ttest.csv")

#incnormal vs cnormal
t.test_results <- mapply(t.test, x= normC, y = normIN, SIMPLIFY = F)
write.matrix(t.test_results, "n_cvin_ttest.csv")

#inctumour vs ctumour
t.test_results <- mapply(t.test, x= tumC, y = tumIN, SIMPLIFY = F)
write.matrix(t.test_results, "t_cvin_ttest.csv")

#ctum vs cnorm
t.test_results <- mapply(t.test, x= tumC, y = normC, paired = T, SIMPLIFY = F)
write.matrix(t.test_results, "tvn_c_ttest.csv")
#inctum vs incnorm
t.test_results <- mapply(t.test, x= tumIN, y = normIN, paired = T, SIMPLIFY = F)
write.matrix(t.test_results, "tvn_inc_ttest.csv")

#signficance of estimate scores
est <- read.csv("est_score.csv", row.names = 1)

shapiro.test(est$StromalScore)
shapiro.test(est$ImmuneScore)
shapiro.test(est$ESTIMATEScore)
shapiro.test(est$tumourPurity)

et_cor <- cor(est, method = "spearman")

et_cor_test <- cor.mtest(est, method = "spearman")

write.csv(et_cor, "et_cor.csv")
write.csv(et_cor_test, "et_cor_test.csv")

corrplot(et_cor, type = "full", method = "number",
         p.mat = et_cor_test, sig.level = 0.05, insig = "blank", tl.col = "black", addCoef.col = "black", diag = FALSE, main = "Estimate")



#box plots for "differntially expressed immune cells"

imb <- read.csv("immune_profile_box.csv")

imb$Macrophages.M1

#method
#use full comparison list to generate significant associations, then include only those with a p > 0.1 in the published figure

imb$RT <- factor(imb$ResT, levels = c("Normal One", "Tumour One", "Normal Two","Tumour Two","Normal Three", "Tumour Three", "Normal Four", "Tumour Four"))

my_comparisons <- list(c("Normal Four", "Tumour Four"), c("Normal Three", "Tumour Three"), c("Normal Two", "Tumour Two"), c("Normal One", "Tumour One"), 
                       c("Normal Four", "Normal Three"), c("Normal Four", "Normal Two"), c("Normal Four", "Normal One"), 
                       c("Normal Three", "Normal Two"), c("Normal Three", "Normal One"), 
                       c("Normal Two", "Normal One"),
                       c("Tumour Four", "Tumour Three"), c("Tumour Four", "Tumour Two"), c("Tumour Four", "Tumour One"), 
                       c("Tumour Three", "Tumour Two"), c("Tumour Three", "Tumour One"), 
                       c("Tumour Two", "Tumour One"))

comp <- list(c("Normal Two", "Tumour Two"))
nu <- ggboxplot(imb, x = "RT", y = "Neutrophils", color = "RT", 
                add = "jitter", legend = "bottom") +
  rotate_x_text(angle = 0) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank())+
  stat_compare_means(method = "kruskal.test", label.y = 0.22, label.x = 2) +        
  stat_compare_means(comparisons = comp, label = "p.signif", method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE, p.adjust.methods = "fdr")
nu <- nu + scale_color_npg(name = "Tissue + Dworak") + ggtitle("j) Neutrophils")
nu


comp <- list(c("Normal Four", "Tumour Four"), c("Normal Two", "Tumour Two"),
             c("Tumour Four", "Tumour Three"), c("Tumour Four", "Tumour Two"), c("Tumour Four", "Tumour One"))
             

m1 <- ggboxplot(imb, x = "RT", y = "Macrophages.M1", color = "RT", 
                   add = "jitter", legend = "bottom") +
  rotate_x_text(angle = 0) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank())+
  stat_compare_means(method = "kruskal.test", label.y = 0.41, label.x = 2) +        
  stat_compare_means(comparisons = comp, label = "p.signif", method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE, p.adjust.methods = "fdr")
m1 <- m1 + scale_color_npg(name = "Tissue + Dworak") + ggtitle("a) M1 macrophages")
m1

comp <- list(c("Tumour One", "Tumour Four"), c("Tumour Four", "Tumour Two"),
             c("Tumour One", "Normal One"), c("Tumour Two", "Normal Two"))

mr <- ggboxplot(imb, x = "RT", y = "Mast.cells.resting", color = "RT", 
                   add = "jitter", legend = "none") +
  rotate_x_text(angle = 0) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank())+
  stat_compare_means(method = "kruskal.test", label.y = 0.4, label.x = 2) +
  stat_compare_means(comparisons = comp, label = "p.signif", method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE)

mr <- mr + scale_color_npg() + ggtitle("b) Resting mast cells")
mr

comp <- list(c("Tumour Four", "Tumour Three"), c("Tumour Four", "Tumour Two"),c("Tumour Four", "Tumour One"),
               c("Tumour One", "Normal One"), c("Tumour Two", "Normal Two"))

ma <- ggboxplot(imb, x = "RT", y = "Mast.cells.activated", color = "RT", 
                      add = "jitter", legend = "bottom") +
  rotate_x_text(angle = 0) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank())+
  stat_compare_means(method = "kruskal.test", label.y = 0.3, label.x = 2) +
  stat_compare_means(comparisons = comp, label = "p.signif", method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE)
ma <- ma + scale_color_npg() + ggtitle("c) Activated mast cells")

adc <- ggboxplot(imb, x = "RT", y = "Dendritic.cells.activated", color = "RT", 
                  add = "jitter", legend = "none") +
  rotate_x_text(angle = 0) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank())+
  stat_compare_means(method = "kruskal.test", label.y = 0.22, label.x = 2) 
adc <- adc + scale_color_npg() + ggtitle("e) Activated dendritic cells")
adc

rdc <- ggboxplot(imb, x = "RT", y = "Dendritic.cells.resting", color = "RT", 
                  add = "jitter", legend = "none") +
  rotate_x_text(angle = 0) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank())+
  stat_compare_means(method = "kruskal.test", label.y = 0.07, label.x = 2) 
rdc <- rdc + scale_color_npg() + ggtitle("d) Resting dendritic cells")
rdc

nbc <- ggboxplot(imb, x = "RT", y = "B.cells.naive", color = "RT", 
                 add = "jitter", legend = "none") +
  rotate_x_text(angle = 0) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank())+
  stat_compare_means(method = "kruskal.test", label.y = 0.33, label.x = 2)
nbc <- nbc + scale_color_npg() + ggtitle("f) Naive B cells")
nbc

comp <- list(c("Tumour One", "Normal One"))

cd8 <- ggboxplot(imb, x = "RT", y = "T.cells.CD8", color = "RT", 
                 add = "jitter", legend = "none") +
  rotate_x_text(angle = 0) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank())+
  stat_compare_means(method = "kruskal.test", label.y = 0.24, label.x = 2) +
  stat_compare_means(comparisons = comp, label = "p.signif", method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE)
cd8 <- cd8 + scale_color_npg() + ggtitle("g) CD8 T cells")
cd8

comp <- list(c("Normal One", "Normal Three"), c("Tumour Three", "Normal Three"), c("Normal Four", "Normal Three"))
eos <- ggboxplot(imb, x = "RT", y = "Eosinophils", color = "RT", 
                 add = "jitter", legend = "none") +
  rotate_x_text(angle = 0) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank())+
  stat_compare_means(method = "kruskal.test", label.y = 0.11, label.x = 2) +
  stat_compare_means(comparisons = comp, label = "p.signif", method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE)
eos <- eos + scale_color_npg() + ggtitle("h) Eosinophils")
eos

nkr <- ggboxplot(imb, x = "RT", y = "NK.cells.resting", color = "RT", 
                 add = "jitter", legend = "none") +
  rotate_x_text(angle = 0) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank())+
  stat_compare_means(method = "kruskal.test", label.y = 0.2, label.x = 2)
nkr <- nkr + scale_color_npg() + ggtitle("i) Resting Natural killer cells")
nkr

fht <- ggboxplot(imb, x = "RT", y = "T.cells.follicular.helper", color = "RT", 
                 add = "jitter", legend = "none") +
  rotate_x_text(angle = 0) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank())+
  stat_compare_means(method = "kruskal.test", label.y = 0.22, label.x = 2) 
fht <- fht + scale_color_npg() + ggtitle("i) Follicular T helper cells")
fht


legend_b <- get_legend(m1 + theme(legend.position="bottom"))

corfig <- ggarrange(m1+ theme(legend.position="none", axis.text.y = element_text(size = 10)),
                     mr+ theme(legend.position="none", axis.text.y = element_text(size = 10)), 
                     ma+ theme(legend.position="none", axis.text.y = element_text(size = 10)),
                     rdc+ theme(legend.position="none", axis.text.y = element_text(size = 10)),
                     adc+ theme(legend.position="none", axis.text.y = element_text(size = 10)),
                     nbc+ theme(legend.position="none", axis.text.y = element_text(size = 10)),
                     cd8+ theme(legend.position="none", axis.text.y = element_text(size = 10)),
                     eos+ theme(legend.position="none", axis.text.y = element_text(size = 10)),
                     fht+ theme(legend.position="none", axis.text.y = element_text(size = 10)),
                     nu+ theme(legend.position="none", axis.text.y = element_text(size = 10)),
                    align = "hv",
                     legend_b,
                     ncol = 4, nrow = 3,
                    hjust = 0.55)

tcorfig <- annotate_figure(corfig, left = "Estimated abundance", bottom = "Tissue + Dworak Score")
tcorfig

endfig<- ggarrange(tcorfig, legend_b, ncol =1, heights = c(1, 0.1))
endfig


corfig <- ggarrange(m1+ theme(legend.position="none", axis.text.y = element_text(size = 10)),
                    mr+ theme(legend.position="none", axis.text.y = element_text(size = 10)), 
                    align = "hv",
                    ncol = 2, nrow = 1)
corfig
tcorfig <- annotate_figure(corfig, left = "Estimated abundance")
tcorfig
endfig<- ggarrange(tcorfig, legend_b, ncol =1, heights = c(1, 0.1))
endfig



#microbiome correlations
library(corrplot)
library(Hmisc)
micor <-read.csv("micro_cor.csv", row.names = 1)
micort <-read.csv("cor_t.csv", row.names = 1)
micorn <-read.csv("cor_n.csv", row.names = 1)

#function for calculating P values across all variables
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

#function for dealing with unweildy correlation tables
ComCor <- function(cormat, pmat, padj) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut],
    fdr = padj[ut]
  )
}

cored <- cor(micor, method = "spearman")
coredp <- cor.mtest(micor, method = "spearman")
corpadj <- as.matrix(p.adjust(coredp, method = "fdr"))
pAdj <- p.adjust(coredp, method = "fdr")
resAdj <- matrix(pAdj, ncol = dim(coredp)[1])
colnames(resAdj) <- colnames(coredp)
rownames(resAdj) <- rownames(coredp)
core <- ComCor(cored, coredp, resAdj)
write.csv(core, "core.csv")

tcored <- cor(micort, method = "spearman")
tcoredp <- cor.mtest(micort, method = "spearman")
tcorpadj <- as.matrix(p.adjust(tcoredp, method = "fdr"))
pAdj <- p.adjust(tcoredp, method = "fdr")
resAdj <- matrix(pAdj, ncol = dim(tcoredp)[1])
colnames(resAdj) <- colnames(tcoredp)
rownames(resAdj) <- rownames(tcoredp)
Tcore <- ComCor(tcored, tcoredp, resAdj)
write.csv(Tcore, "Tcore.csv")

ncored <- cor(micorn, method = "spearman")
ncoredp <- cor.mtest(micorn, method = "spearman")
ncorpadj <- as.matrix(p.adjust(ncoredp, method = "fdr"))
pAdj <- p.adjust(ncoredp, method = "fdr")
resAdj <- matrix(pAdj, ncol = dim(ncoredp)[1])
colnames(resAdj) <- colnames(ncoredp)
rownames(resAdj) <- rownames(ncoredp)
Ncore <- ComCor(ncored, ncoredp, resAdj)
write.csv(Ncore, "Ncore.csv")

#scatter plots
library(ggpmisc)
scatdat <- read.csv("scatdat.csv")

#tumours
#not significant M1 macrophage bacterial correlations
#pos
p1 <- ggplot(data = scatdat, aes(x = Macrophages.M1, y = s_Campylobacter.ureolyticus, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Campylobacter ureolyticus")~ "vs M1 Macrophages"),
       subtitle = expression("Normal: r = 0.063, " ~italic("p")~" = 0.696, FDR = 0.136; Tumour: r = 0.341, " ~italic("p")~" = 0.031, FDR = 0.114"))+
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p1
p2 <- ggplot(data = scatdat, aes(x = Macrophages.M1, y = f_Campylobacteraceae, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Campylobacteraceae")~ "vs M1 Macrophages"),
       subtitle = expression("Normal: r = -0.150, " ~italic("p")~" = 0.354, FDR = 0.516; Tumour: r = 0.326, " ~italic("p")~" = 0.039, FDR = 0.137"))+
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p2


#neg
p3 <- ggplot(data = scatdat, aes(x = Macrophages.M1, y = p_Cyanobacteria, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Cyanobacteria")~ "vs M1 Macrophages"),
       subtitle = expression("Normal: r = -0.089, " ~italic("p")~" = 0.583, FDR = 0.718; Tumour: r = -0.365, " ~italic("p")~" = 0.020, FDR = 0.081"))+
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p3

p4 <- ggplot(data = scatdat, aes(x = Macrophages.M1, y = s_Bacteroides.vulgatus, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Bacteroides vulgatus")~ "vs M1 Macrophages"),
       subtitle = expression("Normal: r = -0.059, " ~italic("p")~" = 0.713, FDR = 0.818; Tumour: r = -0.360, " ~italic("p")~" = 0.022, FDR = 0.88"))+
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p4

p5 <- ggplot(data = scatdat, aes(x = Macrophages.M1, y = g_Desulfovibrio, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Desulfovibrio")~ "vs M1 Macrophages"),
       subtitle = expression("Normal: r = -0.200, " ~italic("p")~" = 0.213, FDR = 0.361; Tumour: r = -0.346, " ~italic("p")~" = 0.028, FDR = 0.106"))+
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p5

p6 <- ggplot(data = scatdat, aes(x = Macrophages.M1, y = s_Alistipes.finegoldii, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Alistipes finegoldii")~ "vs M1 Macrophages"),
       subtitle = expression("Normal: r = -0.065, " ~italic("p")~" = 0.690, FDR = 0.802; Tumour: r = -0.319, " ~italic("p")~" = 0.044, FDR = 0.149"))+
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p6
p7 <- ggplot(data = scatdat, aes(x = Macrophages.M1, y = g_Alistipes, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Alistipes")~ "vs M1 Macrophages"),
       subtitle = expression("Normal: r = -0.103, " ~italic("p")~" = 0.524, FDR = 0.670; Tumour: r = -0.317, " ~italic("p")~" = 0.046, FDR = 0.153"))+
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p7

ggarrange(p1,p2,p3,p4,p5,p6,p7, ncol = 2, nrow = 4)


plot2 <- ggarrange(
  p1 + theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p2 + theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p3 + theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p4 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p5 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p6 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p7 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  legend_b,
  align = "hv", ncol = 2, nrow = 4, labels = c("a)","b)","c)", "d)", "e)", "f)", "g)"))
plot2
dca <- annotate_figure(plot2, left = "Relative transcription", bottom = "Immune cell abundance")
dca




#mast cells
p2 <- ggplot(data = scatdat, aes(x = Mast.cells.resting, y = s_Bacteroides.dorei, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Bacteroides dorei")~ "vs Resting mast cells"),
  subtitle = expression("Normal: r = 0.301, " ~italic("p")~" = 0.058, FDR = 0.136; Tumour: r = -0.445, " ~italic("p")~" = 0.004, FDR = 0.021"))+
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p2

p1 <- ggplot(data = scatdat, aes(x = Mast.cells.activated, y = s_Bacteroides.dorei, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Bacteroides dorei")~ "vs Activated mast cells"),
  subtitle = expression("Normal: r = -0.157, " ~italic("p")~" = 0.332, FDR = 0.494; Tumour: r = 0.452, " ~italic("p")~" = 0.003, FDR = 0.018"))+
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p1

p3 <- ggplot(data = scatdat, aes(x = Mast.cells.resting, y = f_Bifidobacteriaceae, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = "Bifidobacteriaceae vs Resting mast cells",
       subtitle = expression("Normal: r = -0.39, " ~italic("p")~" = 0.010, FDR = 0.035; Tumour: r = 0.037, " ~italic("p")~" = 0.820, FDR = 0.919"))+
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p3


legend_b <- get_legend(p1 + theme(legend.position="right"))
plot2 <- ggarrange(
  p1 + theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p2 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p3 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  align = "v", ncol = 1, labels = c("a)","b)","c)"), legend.grob = legend_b , legend = "right")
mastcells <- annotate_figure(plot2, left = "Relative transcription", bottom = "Immune cell abundance")
mastcells

#m1 macrophages

m1 <- ggplot(data = scatdat, aes(x = Macrophages.M1, y = f_Pasteurellaceae, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  scale_x_continuous(breaks = seq(0, 0.24, 0.05), limits = c(0,0.24))+
  labs(title = "Pasteurellaceae vs M1 Macrophages",
       subtitle = expression("Normal: r = -0.439, " ~italic("p")~" = 0.004, FDR = 0.018; Tumour: r = 0.069, " ~italic("p")~" = 0.671, FDR = 0.835"))+
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
m1

#dendritic resting
dr <- ggplot(data = scatdat, aes(x = Dendritic.cells.activated, y = s_Bacteroides.vulgatus, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Bacteroides vulgatus")~ "vs Resting dendritic cells" ),
       subtitle = expression("Normal: r =  -0.229" ~italic("p")~" = 0.155, FDR = 0.287; Tumour: r = -0.453, " ~italic("p")~" = 0.003, FDR = 0.017"))+  
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
dr


#dendritic cells activated
p1 <- ggplot(data = scatdat, aes(x = Dendritic.cells.activated, y = s_Clostridium.saccharobutylicum, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Clostridium saccharobutylicum")~"vs Activated dendritic cells"),
       subtitle = expression("Normal: r =  -0.157" ~italic("p")~" = 0.331, FDR = 0.493; Tumour: r = 0.403, " ~italic("p")~" = 0.010, FDR = 0.047"))+  
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p1

p2 <- ggplot(data = scatdat, aes(x = Dendritic.cells.activated, y = s_Salmonella.enterica, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Salmonella enterica")~ "vs Activated dendritic cells"),
       subtitle = expression("Normal: r =  -0.164" ~italic("p")~" = 0.309, FDR = 0.470; Tumour: r = 0.431, " ~italic("p")~" = 0.005, FDR = 0.029"))+  
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p2

p3 <- ggplot(data = scatdat, aes(x = Dendritic.cells.activated, y = g_Escherichia, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Escherichia")~ "vs Activated dendritic cells"),
       subtitle = expression("Normal: r = -0.266, " ~italic("p")~" = 0.09, FDR = 0.200; Tumour: r = 0.450, " ~italic("p")~" = 0.003, FDR = 0.020"))+
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p3
p4 <- ggplot(data = scatdat, aes(x = Dendritic.cells.activated, y = s_Escherichia.coli, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Escherichia coli")~ "vs Activated dendritic cells"),
       subtitle = expression("Normal: r =  -0.176" ~italic("p")~" = 0.257, FDR = 0.433; Tumour: r = 0.476, " ~italic("p")~" = 0.002, FDR = 0.012"))+  
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p4

plot2 <- ggarrange(
  m1 + theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank(),),
  dr + theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p1 + theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p2 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p3 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p4 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  align = "v", ncol = 2, nrow = 3, labels = c("a)","b)","c)", "d)", "e)", "f)"), legend.grob = legend_b , legend = "right")
plot2
dca <- annotate_figure(plot2, left = "Relative transcription", bottom = "Immune cell abundance")
dca

#B frag
p1 <- ggplot(data = scatdat, aes(x = B.cells.naive, y = s_Bacteroides.fragilis, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Bacteroides fragilis")~ "vs Naive B cells"),
       subtitle = expression("Normal: r = -0.077, " ~italic("p")~" = 0.633, FDR = 0.757; Tumour: r = 0.431, " ~italic("p")~" = 0.005, FDR = 0.028"))+  
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p1

p2 <- ggplot(data = scatdat, aes(x = T.cells.CD8, y = s_Bacteroides.fragilis, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Bacteroides fragilis")~ "vs CD8 T cells"),
       subtitle = expression("Normal: r = 0.126, " ~italic("p")~" = 0.437, FDR = 0.593; Tumour: r = 0.490, " ~italic("p")~" = 0.001, FDR = 0.008"))+  
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p2

#eubacterium rectale

p3 <- ggplot(data = scatdat, aes(x = Eosinophils, y = s_rectale, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Eubacterium rectale")~ "vs Eosinophils"),
       subtitle = expression("Normal: r = 0.184, " ~italic("p")~" = 0.255, FDR = 0.410; Tumour: r = 0.427, " ~italic("p")~" = 0.006, FDR = 0.029"))+  
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()

p4 <- ggplot(data = scatdat, aes(x = NK.cells.resting, y = s_rectale, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Eubacterium rectale")~ "vs Resting natural killer (NK) cells"),
       subtitle = expression("Normal: r = 0.145, " ~italic("p")~" = 0.369, FDR = 0.531; Tumour: r = 0.428, " ~italic("p")~" = 0.005, FDR = 0.029"))+  
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()


#follicular t helper

p5 <- ggplot(data = scatdat, aes(x = T.cells.follicular.helper, y = s_Streptococcus.pyogenes, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Streptococcus pyogenes")~ "vs Follicular T helper cells"),
       subtitle = expression("Normal: r = 0.422, " ~italic("p")~" = 0.0066, FDR = 0.024; Tumour: r = -0.235, " ~italic("p")~" = 0.143, FDR = 0.334"))+  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p6 <- ggplot(data = scatdat, aes(x = T.cells.follicular.helper, y = f_Flavobacteriaceae, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = "Flavobacteriaceae vs Follicular T helper cells",
       subtitle = expression("Normal: r = -0.205, " ~italic("p")~" = 0.204, FDR = 0.350; Tumour: r = 0.413, " ~italic("p")~" = 0.007, FDR = 0.038"))+
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()


plot2 <- ggarrange(
  p1 + theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p2 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p3 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p4 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p5 + theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p6 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  align = "h", ncol = 2, nrow = 3, labels = c("a)","b)","c)", "d)", "e)", "f)"), legend.grob = legend_b , legend = "right")
plot2
combo1 <- annotate_figure(plot2, left = "Relative transcription", bottom = "Immune cell abundance")
combo1

#neutrophils
p1 <- ggplot(data = scatdat, aes(x = Neutrophils, y = p_Proteobacteria, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = "Proteobacteria vs Neutrophils",
       subtitle = expression("Normal: r = -0.468, " ~italic("p")~" = 0.002, FDR = 0.010; Tumour: r = -0.002, " ~italic("p")~" = 0.986, FDR = 0.994"))+
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p1

p2 <- ggplot(data = scatdat, aes(x = Neutrophils, y = p_Bacteroidetes, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = "Bacteroidetes vs Neutrophils",
       subtitle = expression("Normal: r = 0.396, " ~italic("p")~" = 0.011, FDR = 0.037; Tumour: r = 0.084, " ~italic("p")~" = 0.605, FDR = 0.791"))+
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p2

p3 <- ggplot(data = scatdat, aes(x = Neutrophils, y = g_Escherichia, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Escherichia")~ "vs Neutrophils"),
       subtitle = expression("Normal: r = -0.461, " ~italic("p")~" = 0.002, FDR = 0.012; Tumour: r = -0.022, " ~italic("p")~" = 0.892, FDR = 0.955"))+  
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p3

p4 <- ggplot(data = scatdat, aes(x = Neutrophils, y = s_Escherichia.coli, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Escherichia coli")~ "vs Neutrophils"),
       subtitle = expression("Normal: r = -0.437, " ~italic("p")~" = 0.004, FDR = 0.018; Tumour: r = -0.029, " ~italic("p")~" = 0.854, FDR = 0.937"))+
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p4

p5 <- ggplot(data = scatdat, aes(x = Neutrophils, y = g_Pseudomonas, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Pseudomonas")~ "vs Neutrophils"),
       subtitle = expression("Normal: r = -0.461, " ~italic("p")~" = 0.002, FDR = 0.011; Tumour: r = -0.057, " ~italic("p")~" = 0.723, FDR = 0.866"))+  
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p5

p6 <- ggplot(data = scatdat, aes(x = Neutrophils, y = s_Pseudomonas.sp..NC02, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Pseudomonas")~ "sp. NC02 vs Neutrophils"),
       subtitle = expression("Normal: r = -0.431, " ~italic("p")~" = 0.005, FDR = 0.020; Tumour: r = -0.134, " ~italic("p")~" = 0.408, FDR = 0.635"))+
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p6
p7 <- ggplot(data = scatdat, aes(x = Neutrophils, y = g_Hungatella, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Hungatella")~ "vs Neutrophils"),
       subtitle = expression("Normal: r = 0.444, " ~italic("p")~" = 0.004, FDR = 0.016; Tumour: r = 0.148, " ~italic("p")~" = 0.361, FDR = 0.594"))+  
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p7

p8 <- ggplot(data = scatdat, aes(x = Neutrophils, y = s_Hungatella.hathewayi, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Hungatella hathewayi")~ "vs Neutrophils"),
       subtitle = expression("Normal: r = 0.439, " ~italic("p")~" = 0.004, FDR = 0.017; Tumour: r = 0.140, " ~italic("p")~" = 0.388, FDR = 0.619"))+
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p8

p9 <- ggplot(data = scatdat, aes(x = Neutrophils, y = g_Butyricimonas, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Butyricimonas")~ "vs Neutrophils" ),
       subtitle = expression("Normal: r = 0.503, " ~italic("p")~" < 0.001, FDR = 0.005; Tumour: r = 0.189, " ~italic("p")~" = 0.241, FDR = 0.465"))+  
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p9

p10 <- ggplot(data = scatdat, aes(x = Neutrophils, y = s_Butyricimonas.faecalis, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Butyricimonas faecalis")~ "vs Neutrophils"),
       subtitle = expression("Normal: r = 0.504, " ~italic("p")~" < 0.001, FDR = 0.004; Tumour: r = 0.174, " ~italic("p")~" = 0.280, FDR = 0.509"))+
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p10

plot2 <- ggarrange(
  p1 + theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p2 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p3 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p4 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p5 + theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p6 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  align = "h", ncol = 2, nrow = 3, labels = c("a)","b)","c)", "d)","e)", "f)"), legend.grob = legend_b , legend = "right")
plot2
neut1 <- annotate_figure(plot2, left = "Relative transcription", bottom = "Immune cell abundance")
neut1

plot2 <- ggarrange(
  p7 + theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p8 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p9 + theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p10 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  align = "h", ncol = 2, nrow = 2, labels = c("a)","b)","c)", "d)","e)", "f)"), legend.grob = legend_b , legend = "right")
plot2
neut2 <- annotate_figure(plot2, left = "Relative transcription", bottom = "Immune cell abundance")
neut2

#rouge neutrophils
p4 <- ggplot(data = scatdat, aes(x = Neutrophils, y = s_Salmonella.enterica, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Salmonella enterica")~ "vs Neutrophils" ),
       subtitle = expression("Normal: r = -0.431, " ~italic("p")~" = 0.005, FDR = 0.020; Tumour: r = -0.055, " ~italic("p")~" = 0.732, FDR = 0.871"))+  
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p4

p2 <- ggplot(data = scatdat, aes(x = Neutrophils, y = s_Clostridium.saccharobutylicum, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Clostridium saccharobutylicum")~ "vs Neutrophils"),
       subtitle = expression("Normal: r = -0.401, " ~italic("p")~" = 0.010, FDR = 0.034; Tumour: r = -0.056, " ~italic("p")~" = 0.729, FDR = 0.869"))+
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p2

p3 <- ggplot(data = scatdat, aes(x = Neutrophils, y = s_Flavonifractor.plautii, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Flavonifractor plautii")~ "vs Neutrophils"),
       subtitle = expression("Normal: r = 0.389, " ~italic("p")~" = 0.013, FDR = 0.041; Tumour: r = 0.284, " ~italic("p")~" = 0.074, FDR = 0.215"))+  
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p3

p1 <- ggplot(data = scatdat, aes(x = Neutrophils, y = g_Alistipes, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Alistipes")~ "vs Neutrophils"),
       subtitle = expression("Normal: r = 0.392, " ~italic("p")~" = 0.012, FDR = 0.040; Tumour: r = 0.086, " ~italic("p")~" = 0.595, FDR = 0.785"))+
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p1
p5 <- ggplot(data = scatdat, aes(x = Neutrophils, y = s_Bacteroides.fragilis, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Bacteroides fragilis")~ "vs Neutrophils"),
       subtitle = expression("Normal: r = 0.393, " ~italic("p")~" = 0.012, FDR = 0.039; Tumour: r = -0.343, " ~italic("p")~" = 0.030, FDR = 0.110"))+  
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p5

p6 <- ggplot(data = scatdat, aes(x = Neutrophils, y = s_Odoribacter.splanchnicus, color = Tissue)) +
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method = lm, level = 0.95) +
  labs(title = expression(~italic("Odoribacter splanchnicus")~ "vs Neutrophils"),
       subtitle = expression("Normal: r = 0.436, " ~italic("p")~" = 0.005, FDR = 0.018; Tumour: r = 0.023, " ~italic("p")~" = 0.884, FDR = 0.951"))+
  theme(plot.subtitle=element_text(size=8))  +
  scale_color_npg()
p6

plot2 <- ggarrange(
  p1 + theme(legend.position="none",axis.title.x = ecallement_blank(), axis.title.y = element_blank(),),
  p2 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p3 + theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p4 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p5 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  p6 + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(),),
  align = "h", ncol = 2, nrow = 3, labels = c("a)","b)","c)", "d)","e)", "f)"), legend.grob = legend_b , legend = "right")
plot2
neut3 <- annotate_figure(plot2, left = "Relative transcription", bottom = "Immune cell abundance")
neut3

