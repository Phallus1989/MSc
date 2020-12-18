#Rectal cancer microbiome analysis

#William Taylor 18/12/2020
#wil_tr@live.co.uk


library("phyloseq")
library("vegan")
library("ggplot2")
library("plyr")
library("tidyverse")
library("ape")
library("ggsci")
library("ggpubr")
library("reshape2")
library("RColorBrewer")
library("PCAtools")
library("mixOmics")
library("scales")
library("dplyr")
library("magrittr")
library("reshape2")
library("knitr")
library("grid")
library("rfUtilities")
library("nlme")
library("compositions")
library("limma")
library("edgeR")

scripts <- c("graphical_methods.R",
             "tree_methods.R",
             "plot_merged_trees.R",
             "specificity_methods.R",
             "ternary_plot.R",
             "richness.R",
             "edgePCA.R",
             "copy_number_correction.R",
             "import_frogs.R",
             "prevalence.R",
             "compute_niche.R")
urls <- paste0("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/R/", scripts)

for (url in urls) {
  source(url)
}


setwd("C:/")

#Importing Biom files
rna <- import_biom("RNA.unmapped.biom")
ont <- import_biom("ONT.biom")
amp <- import_biom("amp.biom")


#getting sample names for metadata construction
write.table(sample_names(rna), "rna-seq_biom_names.txt", row.names = F)
write.table(sample_names(ont), "ont-seq_biom_names.txt", row.names = F)
write.table(sample_names(amp), "amp-seq_biom_names.txt", row.names = F)

#attaching metadata
sample_data(rna) <- read.delim(("meta_data.txt"), row.names = 1)
sample_data(ont) <- read.delim(("ont_meta.txt"), row.names = 1)
sample_data(amp) <- read.delim(("amp_meta.txt"), row.names = 1)

sample_names(amp) <- sample_data(amp)$Pairwise
sample_names(ont) <- sample_data(ont)$Pairwise
sample_names(rna) <- sample_data(rna)$Identifier

#fixing ranking issue
new_rank <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(tax_table(rna)) <- new_rank
colnames(tax_table(ont)) <- new_rank
colnames(tax_table(amp)) <- new_rank
rna
ont
amp

#filtering Non-bacterial/non-informative taxa
rna <- subset_taxa(rna, Kingdom == "k__Bacteria")
ont <- subset_taxa(ont, Kingdom == "k__Bacteria")
amp <- subset_taxa(amp, Kingdom == "k__Bacteria")
rna
ont
amp

#rarefaction curves
rna.r = rarefy_even_depth(rna, rngseed=1, sample.size=0.9*min(sample_sums(rna)), replace=F)
amp.r = rarefy_even_depth(amp, rngseed=1, sample.size=0.9*min(sample_sums(amp)), replace=F)
ont.r = rarefy_even_depth(ont, rngseed=1, sample.size=0.9*min(sample_sums(ont)), replace=F)

write.csv(sample_sums(rna.r), 'sums_rnar.csv')
write.csv(sample_sums(amp.r), 'sums_ampr.csv')
write.csv(sample_sums(ont.r), 'sums_ontr.csv')
write.csv(sample_sums(rna), 'sums_rna.csv')
write.csv(sample_sums(amp), 'sums_amp.csv')
write.csv(sample_sums(ont), 'sums_ont.csv')

col <- c("black", "darkred", "forestgreen", "orange", "blue", "brown", "hotpink")
lty <- c("solid", "dashed", "longdash", "dotdash")
pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)

ontcurve <- with(pars[1:26, ],
            rarecurve(otu_table(ont), step = 5, col = col,
                      lty = lty, label = FALSE))

rnacurve <- with(pars[1:26, ],
                 rarecurve(otu_table(rna), step = 5, col = col,
                           lty = lty, label = FALSE))

ampcurve <- with(pars[1:26, ],
                 rarecurve(otu_table(amp), step = 5, col = col,
                           lty = lty, label = FALSE))


#looking at the impact of rarefaction
par(mfcol = c(1,2))
plot(sort(taxa_sums(rna), TRUE), type="h", ylim=c(0, 3000), main = "a)", xlab = "Number of taxa", ylab = "Taxa sums")
plot(sort(taxa_sums(rna.r), TRUE), type="h", ylim=c(0, 3000), xlim = c(0, 2000), main = "b)", xlab = "Number of taxa", ylab = "Taxa sums")

dev.off()

#alpha diversity
#Tissue plot
sample_data(rna.r)$Tissue <- factor(sample_data(rna.r)$Tissue, levels = c("Tumour", "Normal"))

comp <- list(c("Tumour", "Normal"))

rna_diversity_plot_av <- plot_richness(rna.r, x="Tissue", color = "Tissue", measures=c("Observed", "Simpson", "Shannon")) + geom_boxplot(aes(fill = Tissue), alpha = 0.2)  + scale_fill_npg() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.text.x= element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  stat_compare_means(comparisons = comp, label = "p.signif", method = "wilcox.test", hide.ns = FALSE, p.adjust.methods = "fdr")
rna_diversity_plot_av

rtfig <- annotate_figure(rna_diversity_plot_av, left = text_grob("Alpha diversity score", rot = 90))
rtfig

#Dworak Response

sample_data(rna.r)$RT <- factor(sample_data(rna.r)$RT, levels =c("Tumour One", "Normal One", "Tumour Two", "Normal Two", "Tumour Three", "Normal Three", "Tumour Four", "Normal Four"))

comp <- list(c("Tumour Four", "Normal Four"), c("Tumour Three", "Normal Three"), c("Tumour Two", "Normal Two"),c("Tumour One", "Normal One"),
             c("Tumour Four", "Tumour Three"), c("Tumour Four", "Tumour Two"), c("Tumour Four", "Tumour One"),
             c("Normal Four", "Normal Three"), c("Normal Four", "Normal Two"), c("Normal Four", "Normal One"),
             c("Tumour Three", "Tumour Two"), c("Tumour Three", "Tumour One"),
             c("Normal Three", "Normal Two"), c("Normal Three", "Normal One"),
             c("Tumour Two", "Tumour One"),
             c("Normal Two", "Normal One"))

rna_comp <- list( c("Normal Two", "Normal One"), c("Tumour Two", "Normal Two"))

rna_diversity_response <- plot_richness(rna.r, x="RT", color = "Tissue", measures=c("Observed", "Simpson", "Shannon")) + geom_boxplot(aes(fill = Tissue), alpha = 0.2)  + scale_fill_npg() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.2), axis.title.x = element_blank(), axis.title.y = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  stat_compare_means(comparisons = rna_comp, label = "p.format", method = "wilcox.test", hide.ns = FALSE, p.adjust.methods = "fdr")
rna_diversity_response

rdfig <- annotate_figure(rna_diversity_response, bottom = text_grob("Tissue and Dworak score"))
rdfig

#High and low response plot
sample_data(rna.r)$HT <- factor(sample_data(rna.r)$HT, levels = c("Tumour Low", "Normal Low", "Tumour High", "Normal High"))

comp <- list(c("Tumour High", "Normal High"), c("Tumour Low", "Normal Low"), 
             c("Tumour High", "Tumour Low"), c("Normal High", "Normal Low"))

tcomp <- list(c("Tumour Low", "Normal Low"),c("Tumour High", "Normal High"))

rna_diversity_response <- plot_richness(rna.r, x="HT", color = "Tissue", measures=c("Observed", "Simpson", "Shannon")) + geom_boxplot(aes(fill = Tissue), alpha = 0.2)   + scale_fill_npg() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.2), axis.title.y = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Response") +
  stat_compare_means(comparisons = tcomp, label = "p.signif", method = "wilcox.test", hide.ns = FALSE, p.adjust.methods = "fdr")
rna_diversity_response

hlfig<- annotate_figure(rna_diversity_response, bottom = text_grob("Tissue and response grade"))
hlfig

ggarrange(rtfig, rdfig, hlfig, ncol = 3, common.legend = TRUE, labels = c("a)","b)", "c)"))



#Complete response plot
sample_data(rna.r)$ResT <- factor(sample_data(rna.r)$ResT, levels = c("Tumour Incomplete", "Normal Incomplete", "Tumour Complete", "Normal Complete"))

comp <- list(c("Tumour Incomplete", "Normal Incomplete"), c("Tumour Complete", "Normal Complete"),
             c("Tumour Incomplete","Tumour Complete"),c("Normal Incomplete","Normal Complete"))
       

rna_diversity_response <- plot_richness(rna.r, x="ResT", color = "Tissue", measures=c("Observed", "Simpson", "Shannon")) + geom_boxplot(aes(fill = Tissue), alpha = 0.2)   + scale_fill_npg() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.2), axis.title.x = element_blank(), axis.title.y = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  stat_compare_means(comparisons = comp, label = "p.signif", method = "wilcox.test", hide.ns = FALSE, p.adjust.methods = "fdr")
rna_diversity_response
annotate_figure(rna_diversity_response, bottom = text_grob("Tissue and response"), left = text_grob("Alpha diversity score", rot = 90))

#Filter unclassified and low abundance species
rna1 <- filter_taxa(rna.r, function(x) sum(x > 30) > (0.2*length(x)), TRUE)

#building plots of top 9 + other bacterial genera
#RNA

#Phyla
rna_p <- tax_glom(rna.r, "Phylum")
rna_p <- subset_taxa(rna_p, !Phylum=="p__")
rna_p <- transform_sample_counts(rna_p, function(x) x / sum(x))
#removing kraken2 naming scheme suffix (p__)
write.csv(otu_table(rna_p),"rna_p_otu1.csv")
write.csv(tax_table(rna_p),"rna_p_tax1.csv")
rop <- as.matrix(read.csv("rna_p_otu2.csv", row.names = 1))
rtp <- as.matrix(read.csv("rna_p_tax1.csv", row.names = 1))
TAX <- tax_table(rtp)
OTU <- otu_table(rop, taxa_are_rows = TRUE)
RNAP <- phyloseq(TAX, OTU)
sample_data(RNAP) <- sample_data(rna.r)
sample_data(RNAP)$Dworak <- factor(sample_data(RNAP)$Dworak, levels =c("One", "Two", "Three", "Four"))
rna_HL <- plot_bar(RNAP, fill = "Phylum") +
  geom_bar(stat = "identity", position = "stack") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(vjust = 0.2))
rna_HL <- rna_HL + facet_wrap(~Tissue + Dworak, scales = "free_x", nrow = 2)
rnap_plot <- rna_HL + scale_fill_npg()
rnap_plot

#Genera
rna_g <- tax_glom(rna.r, "Genus")
rna_g <- subset_taxa(rna_g, !Genus == "g__")
rna_g <- transform_sample_counts(rna_g, function(x) x / sum(x))
#removing kraken2 naming scheme suffix (g__)
write.csv(otu_table(rna_g),"rna_genus_otu.csv")
write.csv(tax_table(rna_g),"rna_genus_tax.csv")
rog <- as.matrix(read.csv("rna_g_otu.csv", row.names = 1))
rtg <- as.matrix(read.csv("rna_g_tax.csv", row.names = 1))
TAX <- tax_table(rtg)
OTU <- otu_table(rog, taxa_are_rows = TRUE)
RNAG <- phyloseq(TAX, OTU)
sample_data(RNAG) <- sample_data(rna1)
sample_names(RNAG) <- sample_data(rna1)$Identifier
sample_data(RNAG)$Dworak <- factor(sample_data(RNAG)$Dworak, levels =c("One", "Two", "Three", "Four"))
rna_HL <- plot_bar(RNAG, fill = "Genus") +
  geom_bar(stat = "identity", position = "stack") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(vjust = 0.2))
rna_HL <- rna_HL + facet_wrap(~Tissue + Dworak, scales = "free_x", nrow = 2)
rnag_plot <- rna_HL  + scale_fill_npg()
rnag_plot

#Differential expression
#RNA-seq
#Phyla
#adjust naming scheme suffix for nicer plots
rna_p <- tax_glom(rna, "Phylum")
rna_p <- subset_taxa(rna_p, !Phylum == "p__")
filter <- phyloseq::genefilter_sample(rna_p, filterfun_sample(function(x) x >= 30), 
                                      A = 0.2*nsamples(rna_p))
rna_pf <- prune_taxa(filter, rna_p)
rna_pf
write.csv(otu_table(rna_pf),"rna_diff_otup.csv")
write.csv(tax_table(rna_pf),"rna_diff_taxp.csv")

counts <- read.csv('rna_diff_otup.csv', row.names = 1)

counts <- counts + 1 
des_T <- read.delim("meta_data.txt", row.names = 1)
#assigning parameters
Subject <- des_T$Patient
Tissue <- factor(des_T$Tissue, levels=c("Tumour", "Normal"))
Dworak <- factor(des_T$Dworak, levels= c("Four", "Three", "Two", "One"))
HL <- factor(des_T$HL, levels= c("High", "Low"))
#bulding DGE object
xb <- DGEList(counts = counts, lib.size = colSums(counts), norm.factors =calcNormFactors(counts), samples = des_T, group = des_T$Patient, remove.zeros = TRUE)
xb$samples$lib.size <- colSums(xb$counts)
xb <- calcNormFactors(xb, method="RLE")
xb <- estimateCommonDisp(xb, verbose=F)
xb <- estimateTagwiseDisp(xb)
dim(xb)
#first design for tumour vs normal, with patient block
design <- model.matrix(~0+Subject + Tissue)
design
rownames(design) <- samplenames
colnames(design) <- make.names(colnames(design))
xb <- estimateDisp(xb, design, robust = TRUE)
sqrt(xb$common.dispersion)
dim(xb)
contr.matrix <- makeContrasts(
  TVN = TissueTumour-TissueNormal,
  levels = colnames(design))
contr.matrix
fit <- glmFit(xb, design)
TVN <- glmLRT(fit, contrast=contr.matrix[,"TVN"])

summary(decideTests(TVN, adjust.method = "fdr"))
write.csv(topTags(TVN, adjust.method = "fdr"), "tvn_p.csv")
topTags(TVN, adjust.method = "fdr")

#plotting DET
taxa_tvn <- read.csv("tvn_p.csv")
Ptvnplot <- ggplot(taxa_tvn,aes(x=reorder(taxa, logFC), y=logFC,fill=logFC>0))+
  geom_col() + coord_flip()+ xlab(NULL) + ylab(NULL) + theme(legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_line(colour = "black"), panel.background = element_blank())+ 
  scale_y_continuous(breaks = seq(0, 2.5, 0.5), limits = c(0,2.5)) + scale_fill_manual(values = pal_npg("nrc")(2)[2])
Ptvnplot

#Genera
#adjust naming scheme suffix for nicer plots
rna_g <- tax_glom(rna, "Genus")
rna_g <- subset_taxa(rna_g, !Genus == "g__")
filter <- phyloseq::genefilter_sample(rna_g, filterfun_sample(function(x) x >= 30), 
                                      A = 0.2*nsamples(rna_g))
rna_gf <- prune_taxa(filter, rna_g)
rna_gf
write.csv(otu_table(rna_gf),"rna_diff_otug.csv")
write.csv(tax_table(rna_gf),"rna_diff_taxg.csv")

counts <- read.csv('rna_diff_otug.csv', row.names = 1)

counts <- counts + 1 
des_T <- read.delim("meta_data.txt", row.names = 1)
#assigning parameters
Subject <- des_T$Patient
Tissue <- factor(des_T$Tissue, levels=c("Tumour", "Normal"))
HL <- factor(des_T$HL, levels= c("High", "Low"))
#bulding DGE object
xb <- DGEList(counts = counts, lib.size = colSums(counts), norm.factors =calcNormFactors(counts), samples = des_T, group = des_T$Patient, remove.zeros = TRUE)
xb$samples$lib.size <- colSums(xb$counts)
xb <- calcNormFactors(xb, method="RLE")
xb <- estimateCommonDisp(xb, verbose=F)
xb <- estimateTagwiseDisp(xb)
dim(xb)
#first design for tumour vs normal, with patient block
design <- model.matrix(~0+Subject + Tissue)
design
xb <- estimateDisp(xb, design, robust = TRUE)
sqrt(xb$common.dispersion)
dim(xb)
contr.matrix <- makeContrasts(
  TVN = TissueTumour-TissueNormal,
  levels = colnames(design))
contr.matrix
fit <- glmFit(xb, design)
TVN <- glmLRT(fit, contrast=contr.matrix[,"TVN"])

summary(decideTests(TVN, adjust.method = "fdr"))
write.csv(topTags(TVN, adjust.method = "fdr"), "tvn_g.csv")
topTags(TVN, adjust.method = "fdr")
#plotting DET
taxa_tvn <- read.csv("tvn_g.csv")
Gplot <- ggplot(taxa_tvn,aes(x=reorder(taxa, logFC), y=logFC,fill=logFC>0))+
  geom_col() + coord_flip()+ xlab("Genus") + ylab(NULL) + theme(axis.text.y=element_text(face = "italic"), legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_line(colour = "black"), panel.background = element_blank())+
  scale_y_continuous(breaks = seq(-0.5, 7.5, 0.5), limits = c(-0.5,7.5)) + scale_fill_npg()
Gplot 

#high low comparisons
#phylum
counts <- read.csv('rna_diff_otup.csv', row.names = 1)
counts <- counts + 1 
des_T <- read.delim("meta_data.txt", row.names = 1)
#assigning parameters
Subject <- des_T$Patient
Tissue <- factor(des_T$Tissue, levels=c("Tumour", "Normal"))
Dworak <- factor(des_T$Dworak, levels= c("Four", "Three", "Two", "One"))
HL <- factor(des_T$HL, levels= c("High", "Low"))
#bulding DGE object
xb <- DGEList(counts = counts, lib.size = colSums(counts), norm.factors =calcNormFactors(counts), samples = des_T, group = des_T$Patient, remove.zeros = TRUE)
xb$samples$lib.size <- colSums(xb$counts)
xb <- calcNormFactors(xb, method="RLE")
xb <- estimateCommonDisp(xb, verbose=F)
xb <- estimateTagwiseDisp(xb)
dim(xb)
#assigning parameters
Subject <- des_T$Patient
Tissue <- factor(des_T$Tissue, levels=c("Tumour", "Normal"))
HL <- factor(des_T$HL, levels= c("High", "Low"))
#bulding DGE object
xb <- DGEList(counts = counts, lib.size = colSums(counts), norm.factors =calcNormFactors(counts), samples = des_T, group = des_T$Patient, remove.zeros = TRUE)
xb$samples$lib.size <- colSums(xb$counts)
xb <- calcNormFactors(xb, method="RLE")
xb <- estimateCommonDisp(xb, verbose=F)
xb <- estimateTagwiseDisp(xb)
dim(xb)
#first design for tumour vs normal, with patient block
design <- model.matrix(~0 + Subject + Tissue:HL)
colnames(design) <- c("Subject", "TumourHigh", "NormalHigh", "TumourLow", "NormalLow")
design
xb <- estimateDisp(xb, design, robust = TRUE)
sqrt(xb$common.dispersion)
dim(xb)
contr.matrix <- makeContrasts(
  TumourHigh - NormalHigh,
  TumourLow - NormalLow,
  TumourHigh - TumourLow,
  NormalHigh - NormalLow,
  (TumourHigh - NormalHigh) - (TumourLow - NormalLow),
  levels = colnames(design))
contr.matrix
fit <- glmFit(xb, design)
H_TVN <- glmLRT(fit, contrast=contr.matrix[,"TumourHigh - NormalHigh"])
L_TVN <- glmLRT(fit, contrast=contr.matrix[,"TumourLow - NormalLow"])
HvL_TVN <- glmLRT(fit, contrast=contr.matrix[,"(TumourHigh - NormalHigh) - (TumourLow - NormalLow)"])
tvt_hl <- glmLRT(fit, contrast=contr.matrix[,"TumourHigh - TumourLow"])
nvn_hl <- glmLRT(fit, contrast=contr.matrix[,"NormalHigh - NormalLow"])

summary(decideTests(H_TVN, adjust.method = "fdr"))
summary(decideTests(L_TVN, adjust.method = "fdr"))
summary(decideTests(HvL_TVN, adjust.method = "fdr"))
summary(decideTests(tvt_hl, adjust.method = "fdr"))
summary(decideTests(nvn_hl, adjust.method = "fdr"))

topTags(H_TVN, adjust.method = "fdr")
topTags(L_TVN, adjust.method = "fdr")
topTags(HvL_TVN, adjust.method = "fdr")
topTags(tvt_hl, adjust.method = "fdr")
topTags(nvn_hl, adjust.method = "fdr")

write.csv(topTags(H_TVN, adjust.method = "fdr"), "H_TVNp.csv")
write.csv(topTags(L_TVN, adjust.method = "fdr"), "L_TVNp.csv")
write.csv(topTags(HvL_TVN, adjust.method = "fdr"), "HvL_TVNp.csv")
write.csv(topTags(tvt_hl, adjust.method = "fdr"), "tvt_hlp.csv")
write.csv(topTags(nvn_hl, adjust.method = "fdr"), "nvn_hlp.csv")

#plotting DET
htaxa_tvn <- read.csv("H_TVNp.csv")
ltaxa_tvn <- read.csv("L_TVNp.csv")

Phplot <- ggplot(htaxa_tvn,aes(x=reorder(taxa, logFC), y=logFC,fill=logFC>0))+
  geom_col() + coord_flip()+ xlab(NULL) + ylab(NULL) + theme(legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_line(colour = "black"), panel.background = element_blank())+
  scale_y_continuous(breaks = seq(0, 2.5, 0.5), limits = c(0,2.5)) + scale_fill_manual(values = pal_npg("nrc")(2)[2])
Phplot
Plplot <- ggplot(ltaxa_tvn,aes(x=reorder(taxa, logFC), y=logFC,fill=logFC>0))+
  geom_col() + coord_flip()+ xlab(NULL) + ylab(NULL) + theme(legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_line(colour = "black"), panel.background = element_blank())+
  scale_y_continuous(breaks = seq(0, 2.5, 0.5), limits = c(0,2.5)) + scale_fill_manual(values = pal_npg("nrc")(2)[2])
Plplot
hl_tvn_plot <- ggarrange(Ptvnplot, Phplot, Plplot, ncol = 1, labels = c("a)", "b)","c)"))
hl_tvn_plot
hl_p_fig <- annotate_figure(hl_tvn_plot, bottom = text_grob("log2 fold change"))
hl_p_fig


#genus
counts <- read.csv('rna_diff_otug.csv', row.names = 1)
counts <- counts + 1 
des_T <- read.delim("meta_data.txt", row.names = 1)
#assigning parameters
Subject <- des_T$Patient
Tissue <- factor(des_T$Tissue, levels=c("Tumour", "Normal"))
Dworak <- factor(des_T$Dworak, levels= c("Four", "Three", "Two", "One"))
HL <- factor(des_T$HL, levels= c("High", "Low"))
#bulding DGE object
xb <- DGEList(counts = counts, lib.size = colSums(counts), norm.factors =calcNormFactors(counts), samples = des_T, group = des_T$Patient, remove.zeros = TRUE)
xb$samples$lib.size <- colSums(xb$counts)
xb <- calcNormFactors(xb, method="RLE")
xb <- estimateCommonDisp(xb, verbose=F)
xb <- estimateTagwiseDisp(xb)
dim(xb)
#assigning parameters
Subject <- des_T$Patient
Tissue <- factor(des_T$Tissue, levels=c("Tumour", "Normal"))
HL <- factor(des_T$HL, levels= c("High", "Low"))
#bulding DGE object
xb <- DGEList(counts = counts, lib.size = colSums(counts), norm.factors =calcNormFactors(counts), samples = des_T, group = des_T$Patient, remove.zeros = TRUE)
xb$samples$lib.size <- colSums(xb$counts)
xb <- calcNormFactors(xb, method="RLE")
xb <- estimateCommonDisp(xb, verbose=F)
xb <- estimateTagwiseDisp(xb)
dim(xb)
#first design for tumour vs normal, with patient block
design <- model.matrix(~0 + Subject + Tissue:HL)
colnames(design) <- c("Subject", "TumourHigh", "NormalHigh", "TumourLow", "NormalLow")
design
xb <- estimateDisp(xb, design, robust = TRUE)
sqrt(xb$common.dispersion)
dim(xb)
contr.matrix <- makeContrasts(
  TumourHigh - NormalHigh,
  TumourLow - NormalLow,
  TumourHigh - TumourLow,
  NormalHigh - NormalLow,
  (TumourHigh - NormalHigh) - (TumourLow - NormalLow),
  levels = colnames(design))
contr.matrix
fit <- glmFit(xb, design)
H_TVN <- glmLRT(fit, contrast=contr.matrix[,"TumourHigh - NormalHigh"])
L_TVN <- glmLRT(fit, contrast=contr.matrix[,"TumourLow - NormalLow"])
HvL_TVN <- glmLRT(fit, contrast=contr.matrix[,"(TumourHigh - NormalHigh) - (TumourLow - NormalLow)"])
tvt_hl <- glmLRT(fit, contrast=contr.matrix[,"TumourHigh - TumourLow"])
nvn_hl <- glmLRT(fit, contrast=contr.matrix[,"NormalHigh - NormalLow"])

summary(decideTests(H_TVN, adjust.method = "fdr"))
summary(decideTests(L_TVN, adjust.method = "fdr"))
summary(decideTests(HvL_TVN, adjust.method = "fdr"))
summary(decideTests(tvt_hl, adjust.method = "fdr"))
summary(decideTests(nvn_hl, adjust.method = "fdr"))

topTags(H_TVN, adjust.method = "fdr")
topTags(L_TVN, adjust.method = "fdr")
topTags(HvL_TVN, adjust.method = "fdr")
topTags(tvt_hl, adjust.method = "fdr")
topTags(nvn_hl, adjust.method = "fdr")

write.csv(topTags(H_TVN, adjust.method = "fdr", n = Inf), "H_TVNg.csv")
write.csv(topTags(L_TVN, adjust.method = "fdr"), "L_TVNg.csv")
write.csv(topTags(HvL_TVN, adjust.method = "fdr"), "HvL_TVNg.csv")
write.csv(topTags(tvt_hl, adjust.method = "fdr", n = Inf), "tvt_hlg.csv")
write.csv(topTags(nvn_hl, adjust.method = "fdr"), "nvn_hlg.csv")

#plotting DET
htaxa_tvn <- read.csv("H_TVNg.csv")
ltaxa_tvn <- read.csv("L_TVNg.csv")
hvltaxa_tvn <- read.csv("HvL_TVNg.csv")
tvttaxa_tvn <- read.csv("tvt_hlg.csv")
nvntaxa_tvn <- read.csv("nvn_hlg.csv")

hgplot <- ggplot(htaxa_tvn,aes(x=reorder(taxa, logFC), y=logFC,fill=logFC>0))+
  geom_col() + coord_flip()+ xlab("Genus") + ylab(NULL) + theme(axis.text.y=element_text(face = "italic"), legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_line(colour = "black"), panel.background = element_blank())+
  scale_y_continuous(breaks = seq(-0.5, 7.5, 0.5), limits = c(-0.5,7.2)) + scale_fill_manual(values = pal_npg("nrc")(2)[2]) 
hgplot

lgplot <- ggplot(ltaxa_tvn,aes(x=reorder(taxa, logFC), y=logFC,fill=logFC>0))+
  geom_col() + coord_flip()+ xlab("Genus") + ylab(NULL) + theme(axis.text.y=element_text(face = "italic"), legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_line(colour = "black"), panel.background = element_blank())+
  scale_y_continuous(breaks = seq(-0.5, 7.5, 0.5), limits = c(-0.5,7.2))+ scale_fill_manual(values = pal_npg("nrc")(2)[2]) 
lgplot 

hvlplot <- ggplot(hvltaxa_tvn,aes(x=reorder(taxa, logFC), y=logFC,fill=logFC>0))+
  geom_col() + coord_flip()+ xlab("Genus") + ylab(NULL) + theme(axis.text.y=element_text(face = "italic"), legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_line(colour = "black"), panel.background = element_blank())+
  scale_y_continuous(breaks = seq(-0.5, 7.5, 0.5), limits = c(-0.5,7.2)) + scale_fill_manual(values = pal_npg("nrc")(2)[2]) 
hvlplot 

hl_tvn_gplot <- ggarrange(Gplot, hgplot, lgplot, hvlplot, ncol = 1, labels = c("a)", "b)", "c)","d)"), heights = c(1,1,0.5,0.5))
hl_tvn_gplot
hl_g_fig <- annotate_figure(hl_tvn_gplot, bottom = text_grob("log2 fold change"))
hl_g_fig

tvtplot <- ggplot(tvttaxa_tvn,aes(x=reorder(taxa, logFC), y=logFC,fill=logFC>0))+
  geom_col() + coord_flip()+ xlab(NULL) + ylab(NULL) + theme(axis.text.y=element_text(face = "italic"), legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_line(colour = "black"), panel.background = element_blank())+
  xlab("Genus") +
  scale_y_continuous(breaks = seq(-5.5, 5, 0.5), limits =c(-5.5,5)) +  scale_fill_npg()
tvtplot

nvnplot <- ggplot(nvntaxa_tvn,aes(x=reorder(taxa, logFC), y=logFC,fill=logFC>0))+
  geom_col() + coord_flip()+ xlab(NULL) + ylab(NULL) + theme(axis.text.y=element_text(face = "italic"), legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_line(colour = "black"), panel.background = element_blank())+
  xlab("Genus") +
  scale_y_continuous(breaks = seq(-5.5, 5, 0.5), limits =c(-5.5,5))+ scale_fill_npg()
nvnplot

hl_tissuecomp_gplot <- ggarrange(tvtplot, nvnplot, ncol = 1, labels = c("a)", "b)"), heights = c(1, 0.5))
hl_tissuecomp_gplot
hl_compg_fig <- annotate_figure(hl_tissuecomp_gplot, bottom = text_grob("log2 fold change"))
hl_compg_fig

hl_tvn_plot <- ggarrange(Ptvnplot, Phplot, Plplot, ncol = 1, labels = c("a)", "b)", "c)"))
hl_tvn_plot
hl_p_fig <- annotate_figure(hl_tvn_plot, bottom = text_grob("log2 fold change"))
hl_p_fig

#species level
rna_s <- tax_glom(rna, "Species")
rna_s <- subset_taxa(rna_s, !Species == "s__")
filter <- phyloseq::genefilter_sample(rna_s, filterfun_sample(function(x) x >= 30), 
                                      A = 0.2*nsamples(rna_s))
rna_pf <- prune_taxa(filter, rna_s)
rna_pf
write.csv(otu_table(rna_pf),"rna_diff_otus.csv")
write.csv(tax_table(rna_pf),"rna_diff_taxs.csv")

counts <- read.csv('rna_diff_otus2.csv', row.names = 1)

counts <- counts + 1 
des_T <- read.delim("meta_data.txt", row.names = 1)
#assigning parameters
Subject <- des_T$Patient
Tissue <- factor(des_T$Tissue, levels=c("Tumour", "Normal"))
Dworak <- factor(des_T$Dworak, levels= c("Four", "Three", "Two", "One"))
HL <- factor(des_T$HL, levels= c("High", "Low"))
#bulding DGE object
xb <- DGEList(counts = counts, lib.size = colSums(counts), norm.factors =calcNormFactors(counts), samples = des_T, group = des_T$Patient, remove.zeros = TRUE)
xb$samples$lib.size <- colSums(xb$counts)
xb <- calcNormFactors(xb, method="RLE")
xb <- estimateCommonDisp(xb, verbose=F)
xb <- estimateTagwiseDisp(xb)
dim(xb)
#first design for tumour vs normal, with patient block
design <- model.matrix(~0+Subject + Tissue)
design
rownames(design) <- samplenames
colnames(design) <- make.names(colnames(design))
xb <- estimateDisp(xb, design, robust = TRUE)
sqrt(xb$common.dispersion)
dim(xb)
contr.matrix <- makeContrasts(
  TVN = TissueTumour-TissueNormal,
  levels = colnames(design))
contr.matrix
fit <- glmFit(xb, design)
TVN <- glmLRT(fit, contrast=contr.matrix[,"TVN"])

summary(decideTests(TVN, adjust.method = "fdr"))
write.csv(topTags(TVN, adjust.method = "fdr"), "tvn_s.csv")
topTags(TVN, adjust.method = "fdr")

#plotting DET
taxa_tvn <- read.csv("tvn_s.csv")
stvnplot <- ggplot(taxa_tvn,aes(x=reorder(taxa, logFC), y=logFC,fill=logFC>0))+
  geom_col() + coord_flip()+ xlab(NULL) + ylab(NULL) + theme(axis.text.y=element_text(face = "italic"), legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_line(colour = "black"), panel.background = element_blank())+
  xlab("Species") +
  scale_y_continuous(breaks = seq(-3.5, 6.5, 0.5), limits = c(-3.5,6.5)) + scale_fill_npg()
stvnplot

#genus
counts <- read.csv('rna_diff_otus2.csv', row.names = 1)
counts <- counts + 1 
des_T <- read.delim("meta_data.txt", row.names = 1)
#assigning parameters
Subject <- des_T$Patient
Tissue <- factor(des_T$Tissue, levels=c("Tumour", "Normal"))
Dworak <- factor(des_T$Dworak, levels= c("Four", "Three", "Two", "One"))
HL <- factor(des_T$HL, levels= c("High", "Low"))
#bulding DGE object
xb <- DGEList(counts = counts, lib.size = colSums(counts), norm.factors =calcNormFactors(counts), samples = des_T, group = des_T$Patient, remove.zeros = TRUE)
xb$samples$lib.size <- colSums(xb$counts)
xb <- calcNormFactors(xb, method="RLE")
xb <- estimateCommonDisp(xb, verbose=F)
xb <- estimateTagwiseDisp(xb)
dim(xb)
#assigning parameters
Subject <- des_T$Patient
Tissue <- factor(des_T$Tissue, levels=c("Tumour", "Normal"))
HL <- factor(des_T$HL, levels= c("High", "Low"))
#bulding DGE object
xb <- DGEList(counts = counts, lib.size = colSums(counts), norm.factors =calcNormFactors(counts), samples = des_T, group = des_T$Patient, remove.zeros = TRUE)
xb$samples$lib.size <- colSums(xb$counts)
xb <- calcNormFactors(xb, method="RLE")
xb <- estimateCommonDisp(xb, verbose=F)
xb <- estimateTagwiseDisp(xb)
dim(xb)
#first design for tumour vs normal, with patient block
design <- model.matrix(~0 + Subject + Tissue:HL)
colnames(design) <- c("Subject", "TumourHigh", "NormalHigh", "TumourLow", "NormalLow")
design
xb <- estimateDisp(xb, design, robust = TRUE)
sqrt(xb$common.dispersion)
dim(xb)
contr.matrix <- makeContrasts(
  TumourHigh - NormalHigh,
  TumourLow - NormalLow,
  TumourHigh - TumourLow,
  NormalHigh - NormalLow,
  (TumourHigh - NormalHigh) - (TumourLow - NormalLow),
  levels = colnames(design))
contr.matrix
fit <- glmFit(xb, design)
H_TVN <- glmLRT(fit, contrast=contr.matrix[,"TumourHigh - NormalHigh"])
L_TVN <- glmLRT(fit, contrast=contr.matrix[,"TumourLow - NormalLow"])
HvL_TVN <- glmLRT(fit, contrast=contr.matrix[,"(TumourHigh - NormalHigh) - (TumourLow - NormalLow)"])
tvt_hl <- glmLRT(fit, contrast=contr.matrix[,"TumourHigh - TumourLow"])
nvn_hl <- glmLRT(fit, contrast=contr.matrix[,"NormalHigh - NormalLow"])

summary(decideTests(H_TVN, adjust.method = "fdr"))
summary(decideTests(L_TVN, adjust.method = "fdr"))
summary(decideTests(HvL_TVN, adjust.method = "fdr"))
summary(decideTests(tvt_hl, adjust.method = "fdr"))
summary(decideTests(nvn_hl, adjust.method = "fdr"))

topTags(H_TVN, adjust.method = "fdr")
topTags(L_TVN, adjust.method = "fdr")
topTags(HvL_TVN, adjust.method = "fdr")
topTags(tvt_hl, adjust.method = "fdr")
topTags(nvn_hl, adjust.method = "fdr")

write.csv(topTags(H_TVN, adjust.method = "fdr", n = Inf), "H_TVNs.csv")
write.csv(topTags(L_TVN, adjust.method = "fdr"), "L_TVNs.csv")
write.csv(topTags(HvL_TVN, adjust.method = "fdr"), "HvL_TVNs.csv")
write.csv(topTags(tvt_hl, adjust.method = "fdr", n = Inf), "tvt_hls.csv")
write.csv(topTags(nvn_hl, adjust.method = "fdr"), "nvn_hls.csv")

#plotting DET
htaxa_tvn <- read.csv("H_TVNs.csv")
ltaxa_tvn <- read.csv("L_TVNs.csv")
hvltaxa_tvn <- read.csv("HvL_TVNs.csv")
tvttaxa_tvn <- read.csv("tvt_hls.csv")
nvntaxa_tvn <- read.csv("nvn_hls.csv")

hsplot <- ggplot(htaxa_tvn,aes(x=reorder(taxa, logFC), y=logFC,fill=logFC>0))+
  geom_col() + coord_flip()+ xlab("Species") + ylab(NULL) + theme(axis.text.y=element_text(face = "italic"), legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_line(colour = "black"), panel.background = element_blank())+
  scale_y_continuous(breaks = seq(-3.5, 6.5, 0.5), limits = c(-3.5,6.5)) + scale_fill_manual(values = pal_npg("nrc")(2)[2]) 
hsplot

lsplot <- ggplot(ltaxa_tvn,aes(x=reorder(taxa, logFC), y=logFC,fill=logFC>0))+
  geom_col() + coord_flip()+ xlab("Species") + ylab(NULL) + theme(axis.text.y=element_text(face = "italic"), legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_line(colour = "black"), panel.background = element_blank())+
  scale_y_continuous(breaks = seq(-3.5, 6.5, 0.5), limits = c(-3.5,6.5)) + scale_fill_npg() 
lsplot 

hvlplot <- ggplot(hvltaxa_tvn,aes(x=reorder(taxa, logFC), y=logFC,fill=logFC>0))+
  geom_col() + coord_flip()+ xlab("Species") + ylab(NULL) + theme(axis.text.y=element_text(face = "italic"), legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_line(colour = "black"), panel.background = element_blank())+
  scale_y_continuous(breaks = seq(-3.5, 6.5, 0.5), limits = c(-3.5,6.5)) + scale_fill_manual(values = pal_npg("nrc")(2)[2]) 
hvlplot 

hl_tvn_splot <- ggarrange(stvnplot, hsplot, lsplot, hvlplot, ncol = 1, labels = c("a)", "b)", "c)","d)"))
hl_tvn_splot
hl_s_fig <- annotate_figure(hl_tvn_splot, bottom = text_grob("log2 fold change"))
hl_s_fig


tvtplot <- ggplot(tvttaxa_tvn,aes(x=reorder(taxa, logFC), y=logFC,fill=logFC>0))+
  geom_col() + coord_flip()+ xlab(NULL) + ylab(NULL) + theme(axis.text.y=element_text(face = "italic"), legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_line(colour = "black"), panel.background = element_blank())+
  xlab("Species") +
  scale_y_continuous(breaks = seq(-6, 5, 0.5), limits =c(-6,5)) +  scale_fill_npg()
tvtplot

nvnplot <- ggplot(nvntaxa_tvn,aes(x=reorder(taxa, logFC), y=logFC,fill=logFC>0))+
  geom_col() + coord_flip()+ xlab(NULL) + ylab(NULL) + theme(axis.text.y=element_text(face = "italic"), legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_line(colour = "black"), panel.background = element_blank())+
  xlab("Species") +
  scale_y_continuous(breaks = seq(-6, 5, 0.5), limits =c(-6,5))+ scale_fill_npg()
nvnplot

hl_tissuecomp_gplot <- ggarrange(tvtplot, nvnplot, ncol = 1, labels = c("a)", "b)"), heights = c(1.7, 0.3))
hl_tissuecomp_gplot
hl_compg_fig <- annotate_figure(hl_tissuecomp_gplot, bottom = text_grob("log2 fold change"))
hl_compg_fig

#sanity boxplots

rna_s <- tax_glom(rna.r, "Species")
rna_s <- subset_taxa(rna_s, !Species == "s__")
filter <- phyloseq::genefilter_sample(rna_s, filterfun_sample(function(x) x >= 30), A = 0.2*nsamples(rna_s))
rna_sd <- prune_taxa(filter, rna_s)
rna_sd <- transform_sample_counts(rna_sd, function(x) x / sum(x))
rna_sd
write.csv(otu_table(rna_sd), "sd_otu.csv")
write.csv(tax_table(rna_sd), "sd_tax.csv")
write.csv(sample_data(rna_sd), "sd_sam.csv")
imb <- read.csv("sd_dat.csv", row.names = 1)

imb$RT <- factor(imb$RT, levels = c("Tumour One", "Normal One", "Tumour Two", "Normal Two", "Tumour Three", "Normal Three", "Tumour Four", "Normal Four"))
imb$Dworak <- factor(imb$Dworak, levels = c("One", "Two", "Three","Four"))
imb$Tissue <- factor(imb$Tissue, levels = c("Tumour", "Normal"))
comp <- list(c("Tumour Four", "Normal Four"), c("Tumour Three", "Normal Three"), c("Tumour Two", "Normal Two"),c("Tumour One", "Normal One"),
             c("Tumour Four", "Tumour Three"), c("Tumour Four", "Tumour Two"), c("Tumour Four", "Tumour One"),
             c("Normal Four", "Normal Three"), c("Normal Four", "Normal Two"), c("Normal Four", "Normal One"),
             c("Tumour Three", "Tumour Two"), c("Tumour Three", "Tumour One"),
             c("Normal Three", "Normal Two"), c("Normal Three", "Normal One"),
             c("Tumour Two", "Tumour One"),
             c("Normal Two", "Normal One"))

ccomp <- list(c("Tumour Four", "Normal Four"))
cu <- ggboxplot(imb, x = "RT", y = "Campylobacter.ureolyticus", color = "RT", 
                 add = "jitter") +
  stat_compare_means(comparisons = ccomp, label = "p.format", paired = FALSE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE, p.adjust.methods = "fdr") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "right", axis.title.y = element_blank(), axis.text.x = element_blank()) + 
  scale_color_npg(name = "Tissue + Dworak", labels = c("Tumour One", "Normal One", "Tumour Two", "Normal Two", "Tumour Three", "Normal Three", "Tumour Four", "Normal Four"))
cu

#not significant
fn <- ggboxplot(imb, x = "RT", y = "Butyricimonas.faecalis", color = "RT", 
                add = "jitter") +
  stat_compare_means(comparisons = comp, label = "p.format", paired = FALSE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE, p.adjust.methods = "fdr") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "right", axis.title.y = element_blank(), axis.text.x = element_blank()) + 
  scale_color_npg(name = "Tissue + Dworak", labels = c("Tumour One", "Normal One", "Tumour Two", "Normal Two", "Tumour Three", "Normal Three", "Tumour Four", "Normal Four"))
fn

ccomp <- list(c("Tumour Two", "Normal Two"), c("Normal Three", "Normal Four"))
fn <- ggboxplot(imb, x = "RT", y = "Fusobacterium.nucleatum", color = "RT", 
                add = "jitter") +
  stat_compare_means(comparisons = ccomp, label = "p.format", paired = FALSE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE, p.adjust.methods = "fdr") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "right", axis.title.y = element_blank(), axis.text.x = element_blank()) + 
  scale_color_npg(name = "Tissue + Dworak", labels = c("Tumour One", "Normal One", "Tumour Two", "Normal Two", "Tumour Three", "Normal Three", "Tumour Four", "Normal Four"))
fn

ccomp <- list(c("Tumour Four", "Normal Four"), c("Tumour Three", "Tumour Four"), c("Normal Two", "Normal Four"), c("Normal One", "Normal Four"))
hh <- ggboxplot(imb, x = "RT", y = "Hungatella.hathewayi", color = "RT", 
                add = "jitter") +
  stat_compare_means(comparisons = ccomp, label = "p.format", paired = FALSE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE, p.adjust.methods = "fdr") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "right", axis.title.y = element_blank(), axis.text.x = element_blank()) + 
  scale_color_npg(name = "Tissue + Dworak", labels = c("Tumour One", "Normal One", "Tumour Two", "Normal Two", "Tumour Three", "Normal Three", "Tumour Four", "Normal Four"))
hh

#no significant differences
gam <- ggboxplot(imb, x = "RT", y = "Lachnospiraceae.bacterium.GAM79", color = "RT", 
                add = "jitter") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "none") + scale_color_npg()
gam

#no significant differences
ccomp <- list(c("Tumour Four", "Normal Four"), c("Tumour Three", "Tumour Four"), c("Normal Two", "Normal Four"), c("Normal One", "Normal Four"))
bt <- ggboxplot(imb, x = "RT", y = "Bacteroides.thetaiotaomicron", color = "RT", 
                add = "jitter") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "none") + scale_color_npg()
bt
#no significant differences
ccomp <- list(c("Tumour Four", "Normal Four"), c("Tumour Three", "Tumour Four"), c("Normal Two", "Normal Four"), c("Normal One", "Normal Four"))
bf <- ggboxplot(imb, x = "RT", y = "Butyricimonas.faecalis", color = "RT", 
                add = "jitter") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "none") + scale_color_npg()
bf

#no significant differences
ccomp <- list(c("Tumour Four", "Normal Four"), c("Tumour Three", "Tumour Four"), c("Normal Two", "Normal Four"), c("Normal One", "Normal Four"))
er <- ggboxplot(imb, x = "RT", y = "Eubacterium..rectale", color = "RT", 
                add = "jitter") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "none") + scale_color_npg()
er
#no significant differences
ccomp <- list(c("Tumour Four", "Normal Four"), c("Tumour Three", "Tumour Four"), c("Normal Two", "Normal Four"), c("Normal One", "Normal Four"))
af <- ggboxplot(imb, x = "RT", y = "Alistipes.finegoldii", color = "RT", 
                add = "jitter") +
  stat_compare_means(comparisons = comp, label = "p.signif", paired = FALSE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE, p.adjust.methods = "fdr") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "none") + scale_color_npg()
af

ccomp <- list(c("Tumour Four", "Normal Four"), c("Tumour Three", "Tumour Four"), c("Normal Two", "Normal Four"), c("Normal One", "Normal Four"))
sp <- ggboxplot(imb, x = "RT", y = "Streptococcus.pyogenes", color = "RT", 
                add = "jitter") +
  stat_compare_means(comparisons = comp, label = "p.signif", paired = FALSE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE, p.adjust.methods = "fdr") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "none") + scale_color_npg()
sp

ccomp <- list(c("Normal Two", "Normal One"))
cs <- ggboxplot(imb, x = "RT", y = "Clostridium.saccharobutylicum", color = "RT", 
                add = "jitter") +
  stat_compare_means(comparisons = ccomp, label = "p.signif", paired = FALSE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE, p.adjust.methods = "fdr") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "right", axis.title.y = element_blank(), axis.text.x = element_blank()) + 
  scale_color_npg(name = "Tissue + Dworak", labels = c("Tumour One", "Normal One", "Tumour Two", "Normal Two", "Tumour Three", "Normal Three", "Tumour Four", "Normal Four"))
cs

ccomp <- list(c("Normal Three", "Normal One"))
bl <- ggboxplot(imb, x = "RT", y = "Bifidobacterium.longum", color = "RT", 
                add = "jitter") +
  stat_compare_means(comparisons = ccomp, label = "p.signif", paired = FALSE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE, p.adjust.methods = "fdr") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "right", axis.title.y = element_blank(), axis.text.x = element_blank()) + 
  scale_color_npg(name = "Tissue + Dworak", labels = c("Tumour One", "Normal One", "Tumour Two", "Normal Two", "Tumour Three", "Normal Three", "Tumour Four", "Normal Four"))
bl

#not statistically significant
pa<- ggboxplot(imb, x = "RT", y = "Porphyromonas.asaccharolytica", color = "RT", 
               add = "jitter") +
  stat_compare_means(comparisons = comp, label = "p.signif", paired = FALSE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE, p.adjust.methods = "fdr") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "none") + scale_color_npg()
pa
#not statistically significant
ca<- ggboxplot(imb, x = "RT", y = "Cutibacterium.acnes", color = "RT", 
               add = "jitter") +
  stat_compare_means(comparisons = comp, label = "p.signif", paired = FALSE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE, p.adjust.methods = "fdr") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "none") + scale_color_npg()
ca

ccomp <- list(c("Tumour Two", "Normal Two"))
sm<- ggboxplot(imb, x = "RT", y = "Stenotrophomonas.maltophilia", color = "RT", 
               add = "jitter") +
  stat_compare_means(comparisons = ccomp, label = "p.signif", paired = FALSE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE, p.adjust.methods = "fdr") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "right", axis.title.y = element_blank(), axis.text.x = element_blank()) + 
  scale_color_npg(name = "Tissue + Dworak", labels = c("Tumour One", "Normal One", "Tumour Two", "Normal Two", "Tumour Three", "Normal Three", "Tumour Four", "Normal Four"))
sm

ccomp <- list(c("Tumour Two", "Normal Two"), c("Normal One", "Normal Two"))
pf<- ggboxplot(imb, x = "RT", y = "Paraburkholderia.fungorum", color = "RT", 
               add = "jitter") +
  stat_compare_means(comparisons = ccomp, label = "p.format", paired = FALSE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE, p.adjust.methods = "fdr") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "right", axis.title.y = element_blank(), axis.text.x = element_blank()) + 
  scale_color_npg(name = "Tissue + Dworak", labels = c("Tumour One", "Normal One", "Tumour Two", "Normal Two", "Tumour Three", "Normal Three", "Tumour Four", "Normal Four"))
pf
#not statistically significant
rp<- ggboxplot(imb, x = "RT", y = "Ralstonia.pickettii", color = "RT", 
               add = "jitter") +
  stat_compare_means(comparisons = comp, label = "p.signif", paired = FALSE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE, p.adjust.methods = "fdr") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "none") + scale_color_npg()
rp

ccomp <- list(c("Tumour Two", "Normal Two"), c("Normal One", "Normal Two"), c("Normal Three", "Normal Two"))
nc02<- ggboxplot(imb, x = "RT", y = "Pseudomonas.sp..NC02", color = "RT", 
                 add = "jitter") +
  stat_compare_means(comparisons = ccomp, label = "p.format", paired = FALSE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE, p.adjust.methods = "fdr") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "right", axis.title.y = element_blank(), axis.text.x = element_blank()) + 
  scale_color_npg(name = "Tissue + Dworak", labels = c("Tumour One", "Normal One", "Tumour Two", "Normal Two", "Tumour Three", "Normal Three", "Tumour Four", "Normal Four"))
nc02
#not statistically significant
ao<- ggboxplot(imb, x = "RT", y = "Actinomyces.oris", color = "RT", 
               add = "jitter") +
  stat_compare_means(comparisons = comp, label = "p.signif", paired = FALSE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE, p.adjust.methods = "fdr") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "none") + scale_color_npg()
ao

ccomp <- list(c("Tumour Two", "Normal Two"), c("Tumour One", "Normal One"), c("Normal One", "Normal Two"))
ri<- ggboxplot(imb, x = "RT", y = "Ralstonia.insidiosa", color = "RT", 
               add = "jitter") +
  stat_compare_means(comparisons = ccomp, label = "p.signif", paired = FALSE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE, p.adjust.methods = "fdr") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "right", axis.title.y = element_blank(), axis.text.x = element_blank()) + 
  scale_color_npg(name = "Tissue + Dworak", labels = c("Tumour One", "Normal One", "Tumour Two", "Normal Two", "Tumour Three", "Normal Three", "Tumour Four", "Normal Four"))
ri
#not statistically significant
cm<- ggboxplot(imb, x = "RT", y = "Cupriavidus.metallidurans", color = "RT", 
               add = "jitter") +
  stat_compare_means(comparisons = comp, label = "p.signif", paired = FALSE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE, p.adjust.methods = "fdr") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "none") + scale_color_npg()
cm
#not statistically significant
aj<- ggboxplot(imb, x = "RT", y = "Acinetobacter.johnsonii", color = "RT", 
               add = "jitter") +
  stat_compare_means(comparisons = comp, label = "p.signif", paired = FALSE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE, p.adjust.methods = "fdr") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "none") + scale_color_npg()
aj
#not statistically significant
r13<- ggboxplot(imb, x = "RT", y = "Paracoccus.sp..Arc7.R13", color = "RT", 
               add = "jitter") +
  stat_compare_means(comparisons = comp, label = "p.signif", paired = FALSE, method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE, p.adjust.methods = "fdr") +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank(), legend.position = "none") + scale_color_npg()
r13

#signficant plot
library(gridExtra)
library(patchwork)

legend_b <- get_legend(cu + theme(legend.position="bottom"))

plot2 <- grid.arrange(arrangeGrob)
plot2 <- ggarrange(
  cu + theme(legend.position="none", axis.text.y = element_text(size = 10)),
  fn + theme(legend.position="none", axis.text.y = element_text(size = 10)),
  hh + theme(legend.position="none", axis.text.y = element_text(size = 10)),
  cs + theme(legend.position="none", axis.text.y = element_text(size = 10)),
  bl + theme(legend.position="none", axis.text.y = element_text(size = 10)),
  sm + theme(legend.position="none", axis.text.y = element_text(size = 10)),
  pf + theme(legend.position="none", axis.text.y = element_text(size = 10)),
  nc02 + theme(legend.position="none", axis.text.y = element_text(size = 10)),
  ri + theme(legend.position="none", axis.text.y = element_text(size = 10)),
  align = "v")
plot2

plot(legend_b)


#Ordination / beta-diversity
#permanova testing for centroid differences/effect size:

#filtering taxa
wh0 <- genefilter_sample(rna.r, filterfun_sample(function(x) x > 5), A=0.2*nsamples(rna.r))
rna2 <- prune_taxa(wh0, rna.r)
#assign sample data to dataframe
rnadf <- data.frame(sample_data(rna.r))

#betadisp
#testing
bray <- phyloseq::distance(rna2, method = "bray")
adonis(bray ~ Ager , data = rnadf, strata = rnadf$Tissue)
adonis(bray ~ Gender , data = rnadf, strata = rnadf$Tissue)
adonis(bray ~ Cohort , data = rnadf, strata = rnadf$Tissue)
adonis(bray ~ Cohort , data = rnadf, strata = rnadf$Tissue)
adonis(bray ~ Dworak , data = rnadf, strata = rnadf$Tissue)
adonis(bray ~ Tissue , data = rnadf, strata = rnadf$Dworak)
adonis(bray ~ Cohort + Dworak + Tissue + Cohort/Dworak/Tissue, data = rnadf, strata = rnadf$Tissue)
adonis(bray ~ Cohort + Dworak + Tissue + Dworak/Tissue, data = rnadf, strata = rnadf$Tissue)

beta <- betadisper(bray, rnadf$Dworak, bias.adjust = TRUE, type = "centroid")
permutest(beta)
beta <- betadisper(bray, rnadf$Ager, bias.adjust = TRUE, type = "centroid")
permutest(beta)
beta <- betadisper(bray, rnadf$Cohort * rnadf$Dworak, bias.adjust = TRUE, type = "centroid")
permutest(beta)
beta <- betadisper(bray, rnadf$Gender, bias.adjust = TRUE, type = "centroid")
permutest(beta)
beta <- betadisper(bray, rnadf$Dworak, bias.adjust = TRUE, type = "centroid")
permutest(beta)
beta <- betadisper(bray, rnadf$Tissue, bias.adjust = TRUE, type = "centroid")
permutest(beta)

jsd <- phyloseq::distance(rna2, method = "jsd")
adonis(jsd ~  Ager , data = rnadf, strata = rnadf$Tissue)
adonis(jsd ~  Gender , data = rnadf, strata = rnadf$Tissue)
adonis(jsd ~  Cohort , data = rnadf, strata = rnadf$Tissue)
adonis(jsd ~  Dworak , data = rnadf, strata = rnadf$Tissue)
adonis(jsd ~  Tissue , data = rnadf, strata = rnadf$Dworak)
adonis(jsd ~ Cohort + Dworak + Tissue+  Cohort/Dworak/Tissue , data = rnadf, strata = rnadf$Tissue)
adonis(jsd ~ Cohort + Dworak + Tissue + Dworak/Tissue, data = rnadf, strata = rnadf$Tissue)
beta <- betadisper(jsd, rnadf$Ager, bias.adjust = TRUE, type = "centroid")
permutest(beta)
beta <- betadisper(jsd, rnadf$Cohort, bias.adjust = TRUE, type = "centroid")
permutest(beta)
beta <- betadisper(jsd, rnadf$Dworak, bias.adjust = TRUE, type = "centroid")
permutest(beta)
beta <- betadisper(jsd, rnadf$Tissue, bias.adjust = TRUE, type = "centroid")
permutest(beta)

jac <- phyloseq::distance(rna2, method = "jaccard")
adonis(jac ~  Ager , data = rnadf, strata = rnadf$Tissue)
adonis(jac ~  Gender , data = rnadf, strata = rnadf$Tissue)
adonis(jac ~  Cohort , data = rnadf, strata = rnadf$Tissue)
adonis(jac ~  Dworak , data = rnadf, strata = rnadf$Tissue)
adonis(jac ~  Tissue , data = rnadf, strata = rnadf$Dworak)
adonis(jac ~ Cohort + Dworak + Tissue+  Cohort/Dworak/Tissue , data = rnadf, strata = rnadf$Tissue)
adonis(jac ~ Cohort + Dworak + Tissue+  Dworak/Tissue , data = rnadf, strata = rnadf$Tissue)
beta <- betadisper(jac, rnadf$Ager, bias.adjust = TRUE, type = "centroid")
permutest(beta)
beta <- betadisper(jac, rnadf$Cohort, bias.adjust = TRUE, type = "centroid")
permutest(beta)
beta <- betadisper(jac, rnadf$Dworak, bias.adjust = TRUE, type = "centroid")
permutest(beta)
beta <- betadisper(jac, rnadf$Tissue, bias.adjust = TRUE, type = "centroid")
permutest(beta)



#Beta diversity plots
sample_data(rna2)$Dworak <- factor(sample_data(rna2)$Dworak ,levels=c("One","Two","Three", "Four"))
rnat <- subset_samples(rna2, Tissue == "Tumour")
rnan <- subset_samples(rna2, Tissue == "Normal")
#ordination and plots
rnao <- ordinate(rnat, "NMDS", "bray")
p1 = plot_ordination(rnat, rnao, type= "samples", color= "Dworak") + scale_fill_npg() + stat_conf_ellipse(size= 1.5) +
  geom_point(size = 2) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     plot.title = element_text(size= 20)) + scale_color_npg()
rnao <- ordinate(rnan, "NMDS", "bray")
p2 = plot_ordination(rnan, rnao, type= "samples", color= "Dworak") + scale_fill_npg() + stat_conf_ellipse(size= 1.5) +
  geom_point(size = 2) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     plot.title = element_text(size= 20)) + scale_color_npg()

rnao <- ordinate(rnat, "NMDS", "jsd")
p3 = plot_ordination(rnat, rnao, type= "samples", color= "Dworak") + scale_fill_npg() + stat_conf_ellipse(size= 1.5) +
  geom_point(size = 2) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     plot.title = element_text(size= 20)) + scale_color_npg()
rnao <- ordinate(rnan, "NMDS", "jsd")
p4 = plot_ordination(rnan, rnao, type= "samples", color= "Dworak") + scale_fill_npg() + stat_conf_ellipse(size= 1.5) +
  geom_point(size = 2) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     plot.title = element_text(size= 20)) + scale_color_npg()

rnao <- ordinate(rnat, "NMDS", "jaccard")
p5 = plot_ordination(rnat, rnao, type= "samples", color= "Dworak") + scale_fill_npg() + stat_conf_ellipse(size= 1.5) +
  geom_point(size = 2) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     plot.title = element_text(size= 20)) + scale_color_npg()
rnao <- ordinate(rnan, "NMDS", "jaccard")
p6 = plot_ordination(rnan, rnao, type= "samples", color= "Dworak") + scale_fill_npg() + stat_conf_ellipse(size= 1.5) +
  geom_point(size = 2) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     plot.title = element_text(size= 20)) + scale_color_npg()

rnaD <- ggarrange(p1,p2,p3,p4,p5,p6, labels = c("a)","b)","c)","d)","e)","f)"), common.legend = TRUE, legend = "right", ncol = 3, nrow = 2)
rnaD


#tissue comparisons

sample_data(rna2)$Patient <- factor(sample_data(rna2)$Patient)


rnao <- ordinate(rna2, "NMDS", "bray")
p7 = plot_ordination(rna2, rnao, type= "samples", color= "Tissue",shape = "Tissue")+  geom_line(aes(group = Patient, color = NULL)) + stat_conf_ellipse(size= 1.5)  + scale_color_npg() +
  geom_point(size = 2) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
rnao <- ordinate(rna2, "NMDS", "jsd")
p8 = plot_ordination(rna2, rnao, type= "samples", color= "Tissue",shape = "Tissue") +geom_line(aes(group = Patient, color = NULL))+  stat_conf_ellipse(size= 1.5) +scale_color_npg() +
  geom_point(size = 2) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
                     
rnao <- ordinate(rna2, "NMDS", "jaccard")
p9 = plot_ordination(rna2, rnao, type= "samples", color= "Tissue",shape = "Tissue") +geom_line(aes(group = Patient, color = NULL)) +  stat_conf_ellipse(size= 1.5) + scale_color_npg() +
  geom_point(size = 2) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

rnaD <- ggarrange(p1,p3,p5,p2,p4,p6, labels = c("a)","b)","c)","d)","e)","f)"), common.legend = TRUE, legend = "right", ncol = 3, nrow = 2)
rnaD
rnax <- ggarrange(p7,p8,p9, common.legend = TRUE, labels = c("g)", "h)", "i)"), legend = "right", ncol = 3, nrow = 1) 
rnax

rnaf <- ggarrange(rnaD, rnax, ncol = 1, heights = c(0.65, 0.35))
rnaf

#Correlation with response

#filter and agglomerate taxa for analysis in different tissues

filter <- phyloseq::genefilter_sample(rna.r, filterfun_sample(function(x) x >= 30), 
                                      A = 0.2*nsamples(rna.r))
rna.c <- prune_taxa(filter, rna.r)
rna.cg <- tax_glom(rna.c, taxrank = "Genus")
rna.cp <- tax_glom(rna.c, taxrank = "Phylum")
rna.cf <- tax_glom(rna.c, taxrank = "Family")
rna.cs <- tax_glom(rna.r, taxrank = "Species")

rna.cg <- subset_taxa(rna.cg, !Genus == "g__")
rna.cp <- subset_taxa(rna.cp, !Phylum == "p__")
rna.cf <- subset_taxa(rna.cf, !Family == "f__")
rna.cs <- subset_taxa(rna.cs, !Species == "s__")
rna.cs <- transform_sample_counts(rna.cs, function(x) x / sum(x))
filter <- phyloseq::genefilter_sample(rna.cs, filterfun_sample(function(x) x >= 0.02), 
                                      A = 0.1*nsamples(rna.cs))
rna.cs <- prune_taxa(filter, rna.cs)
rna.cs
rna.cg <- transform_sample_counts(rna.cg, function(x) x / sum(x))
rna.cp <- transform_sample_counts(rna.cp, function(x) x / sum(x))
rna.cf <- transform_sample_counts(rna.cf, function(x) x / sum(x))

rna.cg_t <- subset_samples(rna.cg, !Tissue == "Normal")
rna.cs_t <- subset_samples(rna.cs, !Tissue == "Normal")
rna.cp_t <- subset_samples(rna.cp, !Tissue == "Normal")
rna.cf_t <- subset_samples(rna.cf, !Tissue == "Normal")
rna.cg_n <- subset_samples(rna.cg, !Tissue == "Tumour")
rna.cs_n <- subset_samples(rna.cs, !Tissue == "Tumour")
rna.cp_n <- subset_samples(rna.cp, !Tissue == "Tumour")
rna.cf_n <- subset_samples(rna.cf, !Tissue == "Tumour")

write.csv(otu_table(rna.cg_t), "tcg_otu.csv")
write.csv(otu_table(rna.cs_t), "tcs_otu.csv")
write.csv(otu_table(rna.cp_t), "tcp_otu.csv")
write.csv(otu_table(rna.cf_t), "tcf_otu.csv")
write.csv(otu_table(rna.cg_n), "ncg_otu.csv")
write.csv(otu_table(rna.cs_n), "ncs_otu.csv")
write.csv(otu_table(rna.cp_n), "ncp_otu.csv")
write.csv(otu_table(rna.cf_n), "ncf_otu.csv")
write.csv(tax_table(rna.cg_t), "tcg_tax.csv")
write.csv(tax_table(rna.cs_t), "tcs_tax.csv")
write.csv(tax_table(rna.cp_t), "tcp_tax.csv")
write.csv(tax_table(rna.cf_t), "tcf_tax.csv")
write.csv(tax_table(rna.cg_n), "ncg_tax.csv")
write.csv(tax_table(rna.cs_n), "ncs_tax.csv")
write.csv(tax_table(rna.cp_n), "ncp_tax.csv")
write.csv(tax_table(rna.cf_n), "ncf_tax.csv")

write.csv(sample_data(rna.cg_t), "metat.csv")
write.csv(sample_data(rna.cg_n), "metan.csv")

#Bacterial correlations

#format naming schemes and adjust species abundance to reletive abudance
library(corrplot)

tcore <- read.csv("tcor.csv", row.names = 1)
ncore <- read.csv("ncor.csv", row.names = 1)

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

tcored <- cor(tcore, method = "spearman")
tcoredp <- cor.mtest(tcore, method = "spearman")
tcorpadj <- as.matrix(p.adjust(tcoredp, method = "fdr"))
pAdj <- p.adjust(tcoredp, method = "fdr")
resAdj <- matrix(pAdj, ncol = dim(tcoredp)[1])
colnames(resAdj) <- colnames(tcoredp)
rownames(resAdj) <- rownames(tcoredp)
write.csv(tcored, "tcord.csv")
write.csv(resAdj, "tcorp.csv")
write.csv(tcoredp, "tcorpunadj.csv")

corrplot(tcored, p.mat = resAdj, sig.level=0.05, insig= "pch", outline = TRUE, tl.cex = 0.7, method = "color",  tl.col = "black", cl.align="r", order="hclust", hclust.method = "centroid", diag = FALSE, type = "upper", pch.col = "black", pch.cex = 1)

ncored <- cor(ncore, method = "spearman")
ncoredp <- cor.mtest(ncore, method = "spearman")
ncorpadj <- as.matrix(p.adjust(ncoredp, method = "fdr"))
pAdj <- p.adjust(ncoredp, method = "fdr")
resAdj <- matrix(pAdj, ncol = dim(ncoredp)[1])
colnames(resAdj) <- colnames(ncoredp)
rownames(resAdj) <- rownames(ncoredp)
write.csv(ncored, "ncord.csv")
write.csv(resAdj, "ncorp.csv")
write.csv(ncoredp, "ncorpunadj.csv")

corrplot(ncored, p.mat = resAdj, sig.level=0.05, insig= "pch", outline = TRUE, tl.cex = 0.7, method = "color",  tl.col = "black", cl.align="r", order="hclust", hclust.method = "centroid", diag = FALSE, type = "upper", pch.col = "black", pch.cex = 1)

corrplot(ncore, method = "color",tl.col = "black",  cl.align="r", order="hclust", )

#scatter plots
library(ggpmisc)
#normal tissue 
#phyla
#positive
p1 <- ggplot(ncore, aes(x=Dworak, y=p_Proteobacteria)) + 
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  geom_smooth(method=lm, level = 0.95)+
labs(title = "Proteobacteria",
     subtitle = expression("Dworak: r = 0.256, " ~italic("p")~" = 0.115, FDR = 0.188; Dworak One: r = -0.397, " ~italic("p")~" = 0.0121, FDR = 0.030"))+
  theme(plot.subtitle=element_text(size=8))
p1

#negetive
p2 <- ggplot(ncore, aes(x=Dworak, y=p_Bacteroidetes)) + 
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  geom_smooth(method=lm, level = 0.95)+
  labs(title = "Bacteroidetes",
       subtitle = expression("Dworak: r = -0.295, " ~italic("p")~" = 0.068, FDR = 0.123; Dworak One: r = 0.391, " ~italic("p")~" = 0.013, FDR: 0.033"))+
  theme(plot.subtitle=element_text(size=8))
p2
ggarrange(p1, p2, ncol = 1)

#family
#positive
f1 <- ggplot(ncore, aes(x=Dworak, y=f_Pseudomonadaceae)) + 
geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  geom_smooth(method=lm, level = 0.95)+
  labs(title = "Pseudomonadaceae",
       subtitle = expression("Dworak: r = 0.075, " ~italic("p")~" = 0.64, FDR = 0.73; Dworak One: r = -0.366, " ~italic("p")~" = 0.021, FDR = 0.021"))+
  theme(plot.subtitle=element_text(size=8))
f1
#negetive
f2 <- ggplot(ncore, aes(x=Dworak, y=f_Bacteroidaceae)) + 
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  geom_smooth(method=lm, level = 0.95)+
  labs(title = "Bacteroidaceae",
       subtitle = expression("Dworak: r = -0.301 , " ~italic("p")~" = 0.062, FDR = 0.114; Dworak One: r = 0.391, " ~italic("p")~" = 0.013, FDR = 0.033"))+
  theme(plot.subtitle=element_text(size=8))
f2
#genera
#postive
g1 <- ggplot(ncore, aes(x=Dworak, y=g_Pseudomonas)) + 
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  geom_smooth(method=lm, level = 0.95)+
  labs(title = expression(~italic("Pseudomonas")),
       subtitle = expression("Dworak: r = 0.134, " ~italic("p")~" = 0.412, FDR = 0.511; Dworak One: r = -0.423, " ~italic("p")~" = 0.007, FDR = 0.020"))+
  theme(plot.subtitle=element_text(size=8))
g1

#negetive
g2 <- ggplot(ncore, aes(x=Dworak, y=g_Bacteroides)) + 
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  geom_smooth(method=lm, level = 0.95)+
  labs(title = expression(~italic("Bacteroides")),
       subtitle = expression("Dworak: r = -0.300, " ~italic("p")~" = 0.063, FDR = 0.115; Dworak One: r = 0.41, " ~italic("p")~" = 0.009, FDR = 0.025"))+
  theme(plot.subtitle=element_text(size=8))
g2

#Species
#postive
s1 <- ggplot(ncore, aes(x=Dworak, y=s_Klebsiella.pneumoniae)) + 
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  geom_smooth(method=lm, level = 0.95)+
  labs(title = expression(~italic("Klebsiella pneumoniae")),
       subtitle = expression("Dworak: r = 0.21, " ~italic("p")~" = 0.193, FDR = 0.286; Dworak One: r = -0.423, " ~italic("p")~" = 0.007, FDR = 0.020"))+
  theme(plot.subtitle=element_text(size=8))
s1
s2 <- ggplot(ncore, aes(x=Dworak, y=s_Escherichia.coli)) + 
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  geom_smooth(method=lm, level = 0.95)+
  labs(title = expression(~italic("Escherichia coli")),
       subtitle = expression("Dworak: r = 0.202, " ~italic("p")~" = 0.216, FDR = 0.313; Dworak One: r = -0.366, " ~italic("p")~" = 0.021, FDR = 0.049"))+
  theme(plot.subtitle=element_text(size=8))
s2
s3 <- ggplot(ncore, aes(x=Dworak, y=s_Salmonella.enterica)) + 
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  geom_smooth(method=lm, level = 0.95)+
  labs(title = expression(~italic("Salmonella enterica")),
       subtitle = expression("Dworak: r = 0.082, " ~italic("p")~" = 0.617, FDR = 0.708; Dworak One: r = -0.410, " ~italic("p")~" = 0.009, FDR = 0.025"))+
  theme(plot.subtitle=element_text(size=8))
s3

#negetive
s4 <- ggplot(ncore, aes(x=Dworak, y=s_Bacteroides.fragilis)) + 
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  geom_smooth(method=lm, level = 0.95)+
  labs(title = expression(~italic("Bacteroides fragilis")),
       subtitle = expression("Dworak: r = -0.305, " ~italic("p")~" = 0.058, FDR = 0.108; Dworak One: r = 0.366, " ~italic("p")~" = 0.021, FDR = 0.049"))+
  theme(plot.subtitle=element_text(size=8))
s4
s5 <- ggplot(ncore, aes(x=Dworak, y=s_Campylobacter.ureolyticus)) + 
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  geom_smooth(method=lm, level = 0.95)+
  labs(title = expression(~italic("Campylobacter ureolyticus")),
       subtitle = expression("Dworak: r = -0.318," ~italic("p")~" = 0.047, FDR = 0.091; Dworak One: r = 0.230, " ~italic("p")~" = 0.158, FDR = 0.243"))+
  theme(plot.subtitle=element_text(size=8))
s5
s6 <- ggplot(ncore, aes(x=Dworak, y=s_Odoribacter.splanchnicus)) + 
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  geom_smooth(method=lm, level = 0.95)+
  labs(title = expression(~italic("Odoribacter splanchnicus")),
       subtitle = expression("Dworak: r = -0.087," ~italic("p")~" = 0.595, FDR = 0.689; Dworak One: r = 0.353, " ~italic("p")~" = 0.027, FDR = 0.058"))+
  theme(plot.subtitle=element_text(size=8))
s6


normal_plot1 <- ggarrange(p1, p2,f1,f2,g1,g2,labels =c("a)", "b)","c)","d)","e)","f)"), ncol = 2, nrow = 3)
normal_plot1
normal_plot1 <- annotate_figure(normal_plot1, left = "Reletive transcription", bottom = "Dworak score")
normal_plot1

normal_plot2 <- ggarrange(s2,s1, s3,s4,s5,s6,labels =c("a)", "b)","c)","d)","e)","f)"), ncol = 2, nrow = 3)
normal_plot2
normal_plot2 <- annotate_figure(normal_plot2, left = "Reletive transcription", bottom = "Dworak score")
normal_plot2


#tumour tissue 
#phyla
#positive
p1 <- ggplot(tcore, aes(x=Dworak, y=p_Fusobacteria)) + 
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  geom_smooth(method=lm, level = 0.95)+
  labs(title = "Fusobacteria",
       subtitle = expression("Dworak: r = 0.116," ~italic("p")~" = 0.481, FDR = 0.723; Dworak One: r = -0.33, " ~italic("p")~" = 0.0391, FDR = 0.130"))+
  theme(plot.subtitle=element_text(size=8))
p1
#family
#positive
f1 <- ggplot(tcore, aes(x=Dworak, y=f_Fusobacteriaceae)) + 
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  geom_smooth(method=lm, level = 0.95)+
  labs(title = "Fusobacteriaceae",
       subtitle = expression("Dworak: r = 0.090," ~italic("p")~" = 0.584, FDR = 0.794; Dworak One: r = -0.331, " ~italic("p")~" = 0.0391, FDR = 0.130"))+
  theme(plot.subtitle=element_text(size=8))
f1
#genera
#negetive
g1 <- ggplot(tcore, aes(x=Dworak, y=g_Escherichia)) + 
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  geom_smooth(method=lm, level = 0.95)+
  labs(title = expression(~italic("Escherichia")),
       subtitle = expression("Dworak: r = -0.150," ~italic("p")~" = 0.361, FDR = 0.629; Dworak One: r = 0.328, " ~italic("p")~" = 0.041 , FDR = 0.135"))+
  theme(plot.subtitle=element_text(size=8))
g1
#Species
#postive

s1<- ggplot(tcore, aes(x=Dworak, y=s_Fusobacterium.nucleatum)) + 
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  geom_smooth(method=lm, level = 0.95)+
  labs(title = expression(~italic("Fusobacterium nucleatum")),
       subtitle = expression("Dworak: r = 0.082," ~italic("p")~" = 0.616, FDR = 0.809; Dworak One: r = -0.300, " ~italic("p")~" = 0.063, FDR = 0.188"))+
  theme(plot.subtitle=element_text(size=8))
s1
#negetive
s2<- ggplot(tcore, aes(x=Dworak, y=s_Bacteroides.caccae)) + 
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  geom_smooth(method=lm, level = 0.95)+
  labs(title = expression(~italic("Bacteroides caccae")),
       subtitle = expression("Dworak: r = -0.130," ~italic("p")~" = 0.430, FDR = 0.691; Dworak One: r = 0.461, " ~italic("p")~" = 0.003, FDR = 0.016"))+
  theme(plot.subtitle=element_text(size=8))

s2
s3<- ggplot(tcore, aes(x=Dworak, y=s_Bacteroides.vulgatus)) + 
  geom_jitter()+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.title = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  geom_smooth(method=lm, level = 0.95)+
  labs(title = expression(~italic("Bacteroides vulgatus")),
       subtitle = expression("Dworak: r = -0.197," ~italic("p")~" = 0.228, FDR = 0.475; Dworak One: r = 0.347, " ~italic("p")~" = 0.029, FDR = 0.107")) +
  theme(plot.subtitle=element_text(size=8))
s3

tumour_plot <- ggarrange(p1, f1, g1, s2, s3, s1, ncol = 2, nrow = 3,labels =c("a)", "b)","c)","d)","e)","f)","g)","h)","i)"))
tumour_plot
tumour_plot<- annotate_figure(tumour_plot, left = "Reletive transcription", bottom = "Dworak score")
tumour_plot
