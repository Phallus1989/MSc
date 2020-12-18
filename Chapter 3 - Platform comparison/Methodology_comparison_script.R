#Processing Microbiomes for Host Mapping Alone is Insufficient.

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
library("lubridate")
library("reshape2")
library("pavian")

setwd("C:/")

#CRC-RNA
CR_h_map <- import_biom("CRC-RNA_H_MAP.biom")
CR_h_nomap <- import_biom("CRC-RNA_H_NOMAP.biom")
CR_bac_map <- import_biom("RNA_Bac_Map.biom")

#Rec-ONT
O_h_map <- import_biom("ONT_H_MAP.biom")
O_h_nomap <- import_biom("ONT_H_NOMAP.biom")
O_bac_map <- import_biom("ONT_Bac_Map.biom")

#getting sample names for metadata construction
write.table(sample_names(CR_bac_map), "CRC_bac_map.txt", row.names = F)
write.table(sample_names(O_bac_map), "ONT_bac_map.txt", row.names = F)
write.table(sample_names(CR_h_map), "CRC_h_map.txt", row.names = F)
write.table(sample_names(CR_h_nomap), "CRC_h_nomap.txt", row.names = F)
write.table(sample_names(O_h_map), "O_h_map.txt", row.names = F)
write.table(sample_names(O_h_nomap), "O_h_nomap.txt", row.names = F)

#fixing ranking issue
new_rank <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(tax_table(CR_h_map)) <- new_rank
colnames(tax_table(CR_h_nomap)) <- new_rank
colnames(tax_table(CR_bac_map)) <- new_rank
colnames(tax_table(O_h_map)) <- new_rank
colnames(tax_table(O_h_nomap)) <- new_rank
colnames(tax_table(O_bac_map)) <- new_rank

#attaching metadata
sample_data(CR_h_map) <- read.table(("meta_crc_h_map.txt"), row.names = 1)
sample_data(CR_h_nomap) <- read.table(("meta_crc_h_nomap.txt"), row.names = 1)
sample_data(CR_bac_map) <- read.table(("meta_crc_bac_map.txt"), row.names = 1)
sample_data(O_h_map) <- read.table(("meta_ont_h_map.txt"), row.names = 1)
sample_data(O_h_nomap) <- read.table(("meta_ont_h_nomap.txt"), row.names = 1)
sample_data(O_bac_map) <- read.table(("meta_ont_bac_map.txt"), row.names = 1)

#Alpha Diversity
sample_data(O_h_map)$Tissue <- factor(sample_data(O_h_map)$Tissue, levels = c("tumour", "normal"))
sample_data(O_h_nomap)$Tissue <- factor(sample_data(O_h_nomap)$Tissue, levels = c("tumour", "normal"))
sample_data(O_bac_map)$Tissue <- factor(sample_data(O_bac_map)$Tissue, levels = c("tumour", "normal"))

comp <- list(c("tumour", "normal"))

ONT_plot1 <- plot_richness(O_h_map, x="Tissue", measures = c(measures=c("Observed", "Simpson", "Shannon"))) + 
  geom_boxplot(aes(fill = Tissue), alpha = 0.2, outlier.colour = NULL) + scale_fill_npg() + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x= element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())+
  stat_compare_means(comparisons = comp, label = "p.format", method = "wilcox.test",  hide.ns = FALSE)
ONT_plot1$layers <- ONT_plot1$layers[-1]
ONT_plot1

ONT_plot2 <- plot_richness(O_h_nomap, x="Tissue", measures = c(measures=c("Observed", "Simpson", "Shannon"))) + 
  geom_boxplot(aes(fill = Tissue), alpha = 0.2, outlier.colour = NULL) + scale_fill_npg() + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x= element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())+
  stat_compare_means(comparisons = comp, label = "p.format", method = "wilcox.test",  hide.ns = FALSE)
ONT_plot2$layers <- ONT_plot2$layers[-1]
ONT_plot2

ONT_plot3 <- plot_richness(O_bac_map, x="Tissue", measures = c(measures=c("Observed", "Simpson", "Shannon"))) + 
  geom_boxplot(aes(fill = Tissue), alpha = 0.2, outlier.colour = NULL) + scale_fill_npg() + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x= element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())+
  stat_compare_means(comparisons = comp, label = "p.format", method = "wilcox.test",  hide.ns = FALSE)
ONT_plot3$layers <- ONT_plot3$layers[-1]
ONT_plot3

ont_omega <- ggarrange(ONT_plot1, ONT_plot2, ONT_plot3, common.legend = TRUE, legend = "right", ncol = 3, labels = c("a)", "b)", "c)"))
annotate_figure(ont_omega, left = text_grob("Alpha Diversity Measure", rot = 90))


#CRC
sample_data(CR_h_map) <- read.table(("meta_crc_h_map.txt"), row.names = 1)
sample_data(CR_h_nomap) <- read.table(("meta_crc_h_nomap.txt"), row.names = 1)
sample_data(CR_bac_map) <- read.table(("meta_crc_bac_map.txt"), row.names = 1)

sample_data(CR_h_map)$Side <- factor(sample_data(CR_h_map)$Side, levels = c("Left", "Right"))
sample_data(CR_h_nomap)$Side <- factor(sample_data(CR_h_nomap)$Side, levels = c("Left", "Right"))
sample_data(CR_bac_map)$Side <- factor(sample_data(CR_bac_map)$Side, levels = c("Left", "Right"))

comp <- list(c("Left", "Right"))

CRC_plot1 <- plot_richness(CR_h_map, x="Side", measures = c(measures=c("Observed", "Simpson", "Shannon"))) + 
  geom_boxplot(aes(fill = Side), alpha = 0.2, outlier.colour = NULL) + scale_fill_npg() + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x= element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())+
  stat_compare_means(comparisons = comp, label = "p.signif", method = "wilcox.test",  hide.ns = FALSE)

CRC_plot1$layers <- CRC_plot1$layers[-1]
CRC_plot1

CRC_plot2 <- plot_richness(CR_h_nomap, x="Side", measures = c(measures=c("Observed", "Simpson", "Shannon"))) + 
  geom_boxplot(aes(fill = Side), alpha = 0.2, outlier.colour = NULL) + scale_fill_npg() + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x= element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())+
  stat_compare_means(comparisons = comp, label = "p.signif", method = "wilcox.test",  hide.ns = FALSE)

CRC_plot2$layers <- CRC_plot2$layers[-1]
CRC_plot2

CRC_plot3 <- plot_richness(CR_bac_map, x="Side", measures = c(measures=c("Observed", "Simpson", "Shannon"))) + 
  geom_boxplot(aes(fill = Side), alpha = 0.2, outlier.colour = NULL) + scale_fill_npg() + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x= element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())+
  stat_compare_means(comparisons = comp, label = "p.signif", method = "wilcox.test",  hide.ns = FALSE)

CRC_plot3$layers <- CRC_plot3$layers[-1]
CRC_plot3

crc_omega <- ggarrange(CRC_plot1, CRC_plot2, CRC_plot3, common.legend = TRUE, legend = "right", ncol = 3, labels = c("a)", "b)", "c)"))
annotate_figure(crc_omega, left = text_grob("Alpha Diversity Measure", rot = 90))



#Alpha diversity plots
RNA <- import_biom("CRC_total.biom")
ONT <- import_biom("ONT_total.biom")
write.table(sample_names(RNA), "RNA.txt", row.names = F)
write.table(sample_names(ONT), "ONT.txt", row.names = F)
#modify to fit metadata
sample_data(RNA) <- read.delim("C:/Users/William/Google Drive/Purcell studies/PhD/Mapping_paper/Second_version/meta_RNA.txt", row.names=1)
sample_data(ONT) <- read.delim("C:/Users/William/Google Drive/Purcell studies/PhD/Mapping_paper/Second_version/meta_ONT.txt", row.names=1)

new_rank <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(tax_table(RNA)) <- new_rank
colnames(tax_table(ONT)) <- new_rank
RNA <- subset_samples(RNA, !Side == "Rectum")
RNA <- subset_taxa(RNA, Kingdom == "k__Bacteria")
RNA <- prune_taxa(taxa_sums(RNA) > 2, RNA)
ONT <- subset_taxa(ONT, Kingdom == "k__Bacteria")
ONT <- prune_taxa(taxa_sums(ONT) > 2, ONT)


sample_data(RNA)$Method <- factor(sample_data(RNA)$Method, levels = c("All_DB_no_mapping","All_DB_prior mapping", "Bac_DB_prior_mapping"))

comp <- list(c("All_DB_no_mapping","All_DB_prior mapping"),
             c("Bac_DB_prior_mapping", "All_DB_prior mapping"),
             c("Bac_DB_prior_mapping", "All_DB_no_mapping"))

RNA_plot <- plot_richness(RNA, x="Method", measures = c(measures=c("Observed", "Simpson", "Shannon"))) + geom_boxplot(aes(fill = Side), alpha = 0.2, outlier.colour = NULL) + scale_fill_npg() + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x= element_text(size = 15), axis.title.x = element_blank())+
  stat_compare_means(comparisons = comp, label = "p.signif", method = "wilcox.test", ref.group = ".all.", hide.ns = FALSE) +
  scale_x_discrete(labels = c("a)", "b)", "c)")) 

RNA_plot$layers <- RNA_plot$layers[-1]
RNA_plot

sample_data(ONT)$Method <- factor(sample_data(ONT)$Method, levels = c("All_DB_no_mapping","All_DB_w_mapping", "Bac_DB_w_mapping"))

comp <- list(c("All_DB_no_mapping","All_DB_w_mapping"),
             c("Bac_DB_w_mapping", "All_DB_w_mapping"),
             c("Bac_DB_w_mapping", "All_DB_no_mapping"))

ONT_plot <- plot_richness(ONT, x="Method", measures = c(measures=c("Observed", "Simpson", "Shannon"))) + geom_boxplot(aes(fill = Tissue), alpha = 0.2, outlier.colour = NULL) + scale_fill_npg() + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x= element_text(size = 15), axis.title.x = element_blank())+
  stat_compare_means(comparisons = comp, label = "p.signif", method = "wilcox.test",  hide.ns = FALSE) +
  scale_x_discrete(labels = c("a)", "b)", "c)"))

ONT_plot$layers <- ONT_plot$layers[-1]
ONT_plot

#Statistical significance of diveristy differences
write.csv(O_h_map_Diversity, "O_h_map_Diversity.csv")
write.csv(O_h_nomap_Diversity, "O_h_nomap_Diversity.csv")
write.csv(O_bac_map_Diversity, "O_bac_map_Diversity.csv")
write.csv(CR_h_map_Diversity, "CR_h_map_Diversity.csv")
write.csv(CR_h_nomap_Diversity, "CR_h_nomap_Diversity.csv")
write.csv(CR_bac_map_Diversity, "CR_bac_map_Diversity.csv")

#compiled off screen...
CRC_diversity <- read.csv("CRC_diversity_2.csv")
ONT_diversity <- read.delim("ONT_diversity_tab.txt")

#CRC-RNA bac vs all_DB with prior mapping
wilcox.test(CRC_diversity$Observed_bac_m, CRC_diversity$Observed_All_m, conf.int = TRUE, conf.level = 0.95)
wilcox.test(CRC_diversity$Shannon_bac_m, CRC_diversity$Shannon_All_m, conf.int = TRUE, conf.level = 0.95)
wilcox.test(CRC_diversity$Simpson_bac_m, CRC_diversity$Simpson_All_m, conf.int = TRUE, conf.level = 0.95)

#CRC-RNA bac vs all_DB with no prior host mapping
wilcox.test(CRC_diversity$Observed_bac_m, CRC_diversity$Observed_All_nm, conf.int = TRUE, conf.level = 0.95)
wilcox.test(CRC_diversity$Shannon_bac_m, CRC_diversity$Shannon_All_nm, conf.int = TRUE, conf.level = 0.95)
wilcox.test(CRC_diversity$Simpson_bac_m, CRC_diversity$Simpson_All_nm, conf.int = TRUE, conf.level = 0.95)

#CRC-RNA H_DB with prior host mapping vs H_DB without prior host mapping
wilcox.test(CRC_diversity$Observed_All_m, CRC_diversity$Observed_All_nm, conf.int = TRUE, conf.level = 0.95)
wilcox.test(CRC_diversity$Shannon_All_m, CRC_diversity$Shannon_All_nm, conf.int = TRUE, conf.level = 0.95)
wilcox.test(CRC_diversity$Simpson_All_m, CRC_diversity$Simpson_All_nm, conf.int = TRUE, conf.level = 0.95)

#Rec-ONT NH_DB vs H_DB no mapping
wilcox.test(ONT_diversity$Observed_bac_m, ONT_diversity$Observed_All_m)
wilcox.test(ONT_diversity$Shannon_bac_m, ONT_diversity$Shannon_All_m, conf.int = TRUE, conf.level = 0.95)
wilcox.test(ONT_diversity$Simpson_bac_m, ONT_diversity$Simpson_All_m, conf.int = TRUE, conf.level = 0.95)

#ONT-RNA bac vs all_DB with no prior host mapping
wilcox.test(ONT_diversity$Observed_bac_m, ONT_diversity$Observed_All_nm, conf.int = TRUE, conf.level = 0.95)
wilcox.test(ONT_diversity$Shannon_bac_m, ONT_diversity$Shannon_All_nm, conf.int = TRUE, conf.level = 0.95)
wilcox.test(ONT_diversity$Simpson_bac_m, ONT_diversity$Simpson_All_nm, conf.int = TRUE, conf.level = 0.95)

#ONT-RNA H_DB with prior host mapping vs H_DB without prior host mapping
wilcox.test(ONT_diversity$Observed_All_m, ONT_diversity$Observed_All_nm, conf.int = TRUE, conf.level = 0.95)
wilcox.test(ONT_diversity$Shannon_All_m, ONT_diversity$Shannon_All_nm, conf.int = TRUE, conf.level = 0.95)
wilcox.test(ONT_diversity$Simpson_All_m, ONT_diversity$Simpson_All_nm, conf.int = TRUE, conf.level = 0.95)

#NMDS plots of microbiomes

#filtering so singletons done interfer
CR_h_map_p <- subset_samples(CR_h_map, !Side == "Rectum")
CR_h_map_p <- subset_taxa(CR_h_map_p, Kingdom == "k__Bacteria")
CR_h_map_p <- prune_taxa(taxa_sums(CR_h_map_p) > 2, CR_h_map_p)
CR_h_nomap_p <- subset_samples(CR_h_nomap, !Side == "Rectum")
CR_h_nomap_p <- subset_taxa(CR_h_nomap_p, Kingdom == "k__Bacteria")
CR_h_nomap_p <- prune_taxa(taxa_sums(CR_h_nomap_p) > 2, CR_h_nomap_p)
CR_bac_map_p <- subset_samples(CR_bac_map, !Side == "Rectum")
CR_bac_map_p <- subset_taxa(CR_bac_map_p, Kingdom == "k__Bacteria")
CR_bac_map_p <- prune_taxa(taxa_sums(CR_bac_map_p) > 2, CR_bac_map_p)

O_h_map_p <- subset_taxa(O_h_map, Kingdom == "k__Bacteria")
O_h_map_p <- prune_taxa(taxa_sums(O_h_map_p) > 2, O_h_map_p)
O_h_nomap_p <- subset_taxa(O_h_nomap, Kingdom == "k__Bacteria")
O_h_nomap_p <- prune_taxa(taxa_sums(O_h_nomap_p) > 2, O_h_nomap_p)
O_bac_map_p <- subset_taxa(O_bac_map, Kingdom == "k__Bacteria")
O_bac_map_p <- prune_taxa(taxa_sums(O_bac_map_p) > 2, O_bac_map_p)

#data transform
CR_h_map_p <- transform_sample_counts(CR_h_map_p, function(x) x / sum(x))
CR_h_nomap_p <- transform_sample_counts(CR_h_nomap_p, function(x) x / sum(x))
CR_bac_map_p <- transform_sample_counts(CR_bac_map_p, function(x) x / sum(x))

O_h_map_p <- transform_sample_counts(O_h_map_p, function(x) x / sum(x))
O_h_nomap_p <- transform_sample_counts(O_h_nomap_p, function(x) x / sum(x))
O_bac_map_p <- transform_sample_counts(O_bac_map_p, function(x) x / sum(x))

#CRC-RNA

h_map.ord <- ordinate(CR_h_map_p, "NMDS", "bray")
p2 = plot_ordination(CR_h_map_p, h_map.ord, type="samples", color="Side", title="a)") + scale_fill_npg() + stat_ellipse() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

h_nomap.ord <- ordinate(CR_h_nomap_p, "NMDS", "bray")
p4 = plot_ordination(CR_h_nomap_p, h_nomap.ord, type="samples", color="Side", title="b)") + scale_fill_npg() + stat_ellipse() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

bac_map.ord <- ordinate(CR_bac_map_p, "NMDS", "bray")
p3 = plot_ordination(CR_bac_map_p, bac_map.ord, type="samples", color="Side", title="c)") + scale_fill_npg() + stat_ellipse() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggarrange(p2, p4, p3, common.legend = TRUE, legend = "right", ncol = 3)

h_map_dist = phyloseq::distance(CR_h_map_p, method="bray")
h_nomap_dist = phyloseq::distance(CR_h_nomap_p, method="bray")
h_bac_dist = phyloseq::distance(CR_bac_map_p, method="bray")

adonis2(h_map_dist ~ sample_data(CR_h_map_p)$Side)
adonis2(h_nomap_dist ~ sample_data(CR_h_nomap_p)$Side)
adonis2(h_bac_dist ~ sample_data(CR_bac_map_p)$Side)

#Rec-ONT
O_h_map.ord <- ordinate(O_h_map_p, "NMDS", "bray")
op2 = plot_ordination(O_h_map_p, O_h_map.ord, type="samples", color="Tissue", title="a)") + scale_fill_npg() + stat_ellipse() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


O_h_nomap.ord <- ordinate(O_h_nomap_p, "NMDS", "bray")
op4 = plot_ordination(O_h_nomap_p, O_h_nomap.ord, type="samples", color="Tissue", title="b)") + scale_fill_npg() + stat_ellipse() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

O_bac_map.ord <- ordinate(O_bac_map_p, "NMDS", "bray")
op3 = plot_ordination(O_bac_map_p, O_bac_map.ord, type="samples", color="Tissue", title="c)") + scale_fill_npg() + stat_ellipse() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


ggarrange(op2, op4, op3, common.legend = TRUE, legend = "right", ncol =3)

h_map_dist = phyloseq::distance(O_h_map_p, method="bray")
h_nomap_dist = phyloseq::distance(O_h_nomap_p, method="bray")
h_bac_dist = phyloseq::distance(O_bac_map_p, method="bray")

adonis2(h_map_dist ~ sample_data(O_h_map_p)$Tissue)
adonis2(h_nomap_dist ~ sample_data(O_h_nomap_p)$Tissue)
adonis2(h_bac_dist ~ sample_data(O_bac_map_p)$Tissue)


#testing significance of differences betweeen mapping and databases

crc <- read.csv("crc_test.csv", row.names = 1)
ont <- read.csv("ont_test.csv", row.names = 1)

wilcox.test(crc$Bac, crc$All, conf.int = TRUE, conf.level = 0.95)
wilcox.test(ont$bac, ont$all, conf.int = TRUE, conf.level = 0.95)

test <- read.csv("CRC_mapstat.csv", row.names = 1)
head(test)
wilcox.test(test$nomaph, test$maph, conf.int = TRUE, conf.level = 0.95)
wilcox.test(test$Bac, test$nomaph, conf.int = TRUE, conf.level = 0.95)
wilcox.test(test$Bac, test$maph, conf.int = TRUE, conf.level = 0.95)


test <- read.csv("ont_mapstat.csv", row.names = 1)
wilcox.test(test$bac, test$map, conf.int = TRUE, conf.level = 0.95)
wilcox.test(test$bac, test$nomap, conf.int = TRUE, conf.level = 0.95)
wilcox.test(test$map, test$nomap, conf.int = TRUE, conf.level = 0.95)



#Correlation amp vs RNA

library(textshape)
library(dplyr)
library(compare)

#reimport cleaned and checked files
rnaP <- read.csv("rna_P.csv")
ampP <- read.csv("amp_P.csv")
ontP <- read.csv("ont_P.csv")
rnaG <- read.csv("rna_G.csv")
ampG <- read.csv("amp_G.csv")
ontG <- read.csv("ont_G.csv")
rnaS <- read.csv("rna_S.csv")
ampS <- read.csv("amp_S.csv")
ontS <- read.csv("ont_S.csv")

#Phylum
#rna vs 16s
rna_dat_p <- rnaP[rnaP$X %in% ampP$X, ]
amp_dat_p <- ampP[ampP$X %in% rna_dat_p$X, ]
amp_dat_p
rna_dat_p

#convert taxa names to row name
amp_dat_p <- amp_dat_p %>% arrange(desc(X))
rna_dat_p <- rna_dat_p %>% arrange(desc(X))
rna_dat_p <- column_to_rownames(rna_dat_p, "X")
amp_dat_p <- column_to_rownames(amp_dat_p, "X")
amp_dat_p <- amp_dat_p[,order(colnames(amp_dat_p))]
rna_dat_p <- rna_dat_p[,order(colnames(rna_dat_p))]

rownames(amp_dat_p)
rownames(rna_dat_p)
colnames(amp_dat_p)
colnames(rna_dat_p)

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

#perform correlation and export a barplot based on the diagonal mean
c <- cor(rna_dat_p, amp_dat_p, method='spearman')

cor.test(as.matrix(rna_dat_p), as.matrix(amp_dat_p), method='spearman')

ptest <- cor.mtest(c)
write.csv(c, "pc.csv")
write.csv(ptest, "pt.csv")

cd

cd <- diag(c)
mean(cd)
barplot(cd, 
        names.arg=colnames(rna_dat_p),
        main = "RNA x 16S rRNA Phylum Level",
        ylim = c(0,1),  
        col="#495b5b", 
        xlab="Samples", 
        ylab="Correlation",
        border = NA, 
        cex.axis=1.4, 
        cex.names=0.7, 
        cex.lab=1.4,
        las = 2,
        lwd = 3
)
abline(h=mean(cd),col='#000000',lty=2, lwd=3)

#ont vs rna
rna_dat_p <- rnaP[rnaP$X %in% ontP$X, ]
ont_dat_p <- ontP[ontP$X %in% rna_dat_p$X, ]

rna_dat_p
ont_dat_p
#convert taxa names to row name
ont_dat_p <- ont_dat_p %>% arrange(desc(X))
rna_dat_p <- rna_dat_p %>% arrange(desc(X))
rna_dat_p <- column_to_rownames(rna_dat_p, "X")
ont_dat_p <- column_to_rownames(ont_dat_p, "X")
ont_dat_p <- ont_dat_p[,order(colnames(ont_dat_p))]
rna_dat_p <- rna_dat_p[,order(colnames(rna_dat_p))]
#check if they match
rownames(ont_dat_p)
rownames(rna_dat_p)
colnames(ont_dat_p)
colnames(rna_dat_p)

#perform correlation and export a barplot based on the diagonal mean
c <- cor(rna_dat_p, ont_dat_p, method='spearman')
cor.test(as.matrix(rna_dat_p), as.matrix(ont_dat_p), method='spearman')

#testing significance of RT17N
write.csv(ont_dat_p, "op.csv")
write.csv(rna_dat_p, "rp.csv")

phy <- read.csv("phy.csv", row.names = 1)
cor.test(phy$ONT, phy$RNA, method = 'spearman')


cd <- diag(c)
mean(cd)
barplot(cd, 
        names.arg=colnames(rna_dat_p),
        main = "RNA x ONT Phylum Level",
        ylim = c(0,1),  
        col="#495b5b", 
        xlab="Samples", 
        ylab="Correlation",
        border = NA, 
        cex.axis=1.4, 
        cex.names=0.7, 
        cex.lab=1.4,
        las = 2,
        lwd = 3
)
abline(h=mean(cd),col='#000000',lty=2, lwd=3)


#ont vs 16s
amp_dat_p <- ampP[ampP$X %in% ontP$X, ]
ont_dat_p <- ontP[ontP$X %in% amp_dat_p$X, ]

amp_dat_p
ont_dat_p
#convert taxa names to row name
ont_dat_p <- ont_dat_p %>% arrange(desc(X))
amp_dat_p <- amp_dat_p %>% arrange(desc(X))
amp_dat_p <- column_to_rownames(amp_dat_p, "X")
ont_dat_p <- column_to_rownames(ont_dat_p, "X")
ont_dat_p <- ont_dat_p[,order(colnames(ont_dat_p))]
amp_dat_p <- amp_dat_p[,order(colnames(amp_dat_p))]
#check if they match
rownames(ont_dat_p)
rownames(amp_dat_p)
colnames(ont_dat_p)
colnames(amp_dat_p)

#perform correlation and export a barplot based on the diagonal mean
c <- cor(amp_dat_p, ont_dat_p, method='spearman')

cor.test(as.matrix(amp_dat_p), as.matrix(ont_dat_p), method='spearman')

cd <- diag(c)
mean(cd)
barplot(cd, 
        names.arg=colnames(amp_dat_p),
        main = "16S rRNA x ONT Phylum Level",
        ylim = c(0,1),  
        col="#495b5b", 
        xlab="Samples", 
        ylab="Correlation",
        border = NA, 
        cex.axis=1.4, 
        cex.names=0.7, 
        cex.lab=1.4,
        las = 2,
        lwd = 3
)
abline(h=mean(cd),col='#000000',lty=2, lwd=3)

#Genus correlation amp vs RNA

#Genus
#rna vs 16s
rna_dat_g <- rnaG[rnaG$X %in% ampG$Genus, ]
amp_dat_g <- ampG[ampG$Genus %in% rna_dat_g$X, ]
amp_dat_g
rna_dat_g

#convert taxa names to row name
amp_dat_g <- amp_dat_g %>% arrange(desc(Genus))
rna_dat_g <- rna_dat_g %>% arrange(desc(X))
rna_dat_g <- column_to_rownames(rna_dat_g, "X")
amp_dat_g <- column_to_rownames(amp_dat_g, "Genus")
amp_dat_g <- amp_dat_g[,order(colnames(amp_dat_g))]
rna_dat_g <- rna_dat_g[,order(colnames(rna_dat_g))]

rownames(amp_dat_g)
rownames(rna_dat_g)
colnames(amp_dat_g)
colnames(rna_dat_g)

#perform correlation and export a barplot based on the diagonal mean
c <- cor(rna_dat_g, amp_dat_g, method='spearman')

cor.test(as.matrix(rna_dat_g), as.matrix(amp_dat_g), method = "spearman")

cd <- diag(c)
mean(cd)
barplot(cd, 
        names.arg=colnames(rna_dat_p),
        main = "RNA x 16S rRNA Genus Level",
        ylim = c(0,1),  
        col="#495b5b", 
        xlab="Samples", 
        ylab="Correlation",
        border = NA, 
        cex.axis=1.4, 
        cex.names=0.7, 
        cex.lab=1.4,
        las = 2,
        lwd = 3
)
abline(h=mean(cd),col='#000000',lty=2, lwd=3)

#ont vs rna
rna_dat_g <- rnaG[rnaG$X %in% ontG$X, ]
ont_dat_g <- ontG[ontG$X %in% rna_dat_g$X, ]

rna_dat_g
ont_dat_g
#convert taxa names to row name
ont_dat_g <- ont_dat_g %>% arrange(desc(X))
rna_dat_g <- rna_dat_g %>% arrange(desc(X))
rna_dat_g <- column_to_rownames(rna_dat_g, "X")
ont_dat_g <- column_to_rownames(ont_dat_g, "X")
ont_dat_g <- ont_dat_g[,order(colnames(ont_dat_g))]
rna_dat_g <- rna_dat_g[,order(colnames(rna_dat_g))]
#check if they match
rownames(ont_dat_g)
rownames(rna_dat_g)
colnames(ont_dat_g)
colnames(rna_dat_g)

#perform correlation and export a barplot based on the diagonal mean
c <- cor(rna_dat_g, ont_dat_g, method='spearman')

cor.test(as.matrix(rna_dat_g), as.matrix(ont_dat_g), method = "spearman")

cd <- diag(c)
mean(cd)
barplot(cd, 
        names.arg=colnames(rna_dat_g),
        main = "RNA x ONT Genus Level",
        ylim = c(0,1),  
        col="#495b5b", 
        xlab="Samples", 
        ylab="Correlation",
        border = NA, 
        cex.axis=1.4, 
        cex.names=0.7, 
        cex.lab=1.4,
        las = 2,
        lwd = 3
)
abline(h=mean(cd),col='#000000',lty=2, lwd=3)


#ont vs 16s
amp_dat_g <- ampG[ampG$Genus %in% ontG$X, ]
ont_dat_g <- ontG[ontG$X %in% amp_dat_g$Genus, ]

amp_dat_g
ont_dat_g
#convert taxa names to row name
ont_dat_g <- ont_dat_g %>% arrange(desc(X))
amp_dat_g <- amp_dat_g %>% arrange(desc(Genus))
amp_dat_g <- column_to_rownames(amp_dat_g, "Genus")
ont_dat_g <- column_to_rownames(ont_dat_g, "X")
ont_dat_g <- ont_dat_g[,order(colnames(ont_dat_g))]
amp_dat_g <- amp_dat_g[,order(colnames(amp_dat_g))]
#check if they match
rownames(ont_dat_g)
rownames(amp_dat_g)
colnames(ont_dat_g)
colnames(amp_dat_g)

#perform correlation and export a barplot based on the diagonal mean
c <- cor(amp_dat_g, ont_dat_g, method='spearman')

cor.test(as.matrix(amp_dat_g), as.matrix(ont_dat_g), method = "spearman")

cd <- diag(c)
mean(cd)
barplot(cd, 
        names.arg=colnames(amp_dat_g),
        main = "16S rRNA x ONT Genus Level",
        ylim = c(0,1),  
        col="#495b5b", 
        xlab="Samples", 
        ylab="Correlation",
        border = NA, 
        cex.axis=1.4, 
        cex.names=0.7, 
        cex.lab=1.4,
        las = 2,
        lwd = 3
)
abline(h=mean(cd),col='#000000',lty=2, lwd=3)

#species

#ont vs 16s
amp_dat_s <- ampS[ampS$X %in% ontS$X, ]
ont_dat_s <- ontS[ontS$X %in% amp_dat_s$X, ]

amp_dat_s
ont_dat_s
#convert taxa names to row name
ont_dat_s <- ont_dat_s %>% arrange(desc(X))
amp_dat_s <- amp_dat_s %>% arrange(desc(X))
amp_dat_s <- column_to_rownames(amp_dat_s, "X")
ont_dat_s <- column_to_rownames(ont_dat_s, "X")
ont_dat_s <- ont_dat_s[,order(colnames(ont_dat_s))]
amp_dat_s <- amp_dat_s[,order(colnames(amp_dat_s))]
#check if they match
rownames(ont_dat_s)
rownames(amp_dat_s)
colnames(ont_dat_s)
colnames(amp_dat_s)

#perform correlation and export a barplot based on the diagonal mean
c <- cor(amp_dat_s, ont_dat_s, method='spearman')
#removing NA sample from correlation matrix
c <- subset(c, select = -c(RT7N))

cor.test(as.matrix(amp_dat_s), as.matrix(ont_dat_s), method = "spearman")

cd <- diag(c)
mean(cd)
barplot(cd, 
        names.arg=colnames(c),
        main = "16S rRNA x ONT Species Level",
        ylim = c(0,1),  
        col="#495b5b", 
        xlab="Samples", 
        ylab="Correlation",
        border = NA, 
        cex.axis=1.4, 
        cex.names=0.7, 
        cex.lab=1.4,
        las = 2,
        lwd = 3
)
abline(h=mean(cd),col='#000000',lty=2, lwd=3)

#ont vs rna
rna_dat_s <- rnaS[rnaS$X %in% ontS$X, ]
ont_dat_s <- ontS[ontS$X %in% rna_dat_s$X, ]

rna_dat_s
ont_dat_s
#convert taxa names to row name
ont_dat_s <- ont_dat_s %>% arrange(desc(X))
rna_dat_s <- rna_dat_s %>% arrange(desc(X))
rna_dat_s <- column_to_rownames(rna_dat_s, "X")
ont_dat_s <- column_to_rownames(ont_dat_s, "X")
ont_dat_s <- ont_dat_s[,order(colnames(ont_dat_s))]
rna_dat_s <- rna_dat_s[,order(colnames(rna_dat_s))]
#check if they match
rownames(ont_dat_s)
rownames(rna_dat_s)
colnames(ont_dat_s)
colnames(rna_dat_s)

#perform correlation and export a barplot based on the diagonal mean
c <- cor(rna_dat_s, ont_dat_s, method='spearman')

cor.mtest(c)

cor.test(as.matrix(rna_dat_s), as.matrix(ont_dat_s), method = "spearman")

cd <- diag(c)
mean(cd)
barplot(cd, 
        names.arg=colnames(rna_dat_s),
        main = "RNA x ONT Species Level",
        ylim = c(0,1),  
        col="#495b5b", 
        xlab="Samples", 
        ylab="Correlation",
        border = NA, 
        cex.axis=1.4, 
        cex.names=0.7, 
        cex.lab=1.4,
        las = 2,
        lwd = 3
)
abline(h=mean(cd),col='#000000',lty=2, lwd=3)

#16s vs rna
rna_dat_s <- rnaS[rnaS$X %in% ampS$X, ]
amp_dat_s <- ampS[ampS$X %in% rna_dat_s$X, ]

rna_dat_s
amp_dat_s
#convert taxa names to row name
amp_dat_s <- amp_dat_s %>% arrange(desc(X))
rna_dat_s <- rna_dat_s %>% arrange(desc(X))
rna_dat_s <- column_to_rownames(rna_dat_s, "X")
amp_dat_s <- column_to_rownames(amp_dat_s, "X")
amp_dat_s <- amp_dat_s[,order(colnames(amp_dat_s))]
rna_dat_s <- rna_dat_s[,order(colnames(rna_dat_s))]
#check if they match
rownames(amp_dat_s)
rownames(rna_dat_s)
colnames(amp_dat_s)
colnames(rna_dat_s)

#perform correlation and export a barplot based on the diagonal mean
c <- cor(rna_dat_s, amp_dat_s, method='spearman')
cor.test(as.matrix(rna_dat_s), as.matrix(amp_dat_s), method = "spearman")

#testing signfiicance of RT13N
write.csv(rnaS, "rspec.csv")
write.csv(ontS, "ospec.csv")
uro <- read.csv("uro.csv", row.names = 1)
cor.test(uro$RNA, uro$ONT, method = 'spearman')


cd <- diag(c)

mean(cd)
barplot(cd, 
        names.arg=colnames(rna_dat_s),
        main = "RNA x 16S rRNA Species Level",
        ylim = c(0,1),  
        col="#495b5b", 
        xlab="Samples", 
        ylab="Correlation",
        border = NA, 
        cex.axis=1.4, 
        cex.names=0.7, 
        cex.lab=1.4,
        las = 2,
        lwd = 3
)
abline(h=mean(cd),col='#000000',lty=2, lwd=3)

#plots for reletive abundance comparison

phyla <- read.csv("all_phyla.csv")
genus <- read.csv("all_genus.csv")

p2 <- melt(phyla)
p2
g2 <- melt(genus)
g2

order = c("RNA-seq", "ONT", "16S rRNA")

p_plot <- ggplot() + geom_bar(aes(y = value, x = variable, fill = ï..Taxa), data = p2,
                              stat="identity", position = "fill", colour = "black") +
  xlab("Platform") +
  ylab("Proportion") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(labels = c("RNA-seq", "ONT", "16S rRNA")) +
  ggtitle("Phyla Representation by Platform") +
  theme(legend.title = element_text(size=14)) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) 

p_plot
p_plot + scale_fill_npg() + labs(fill = "Phyla") + 
  theme(panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
  )

g_plot <- ggplot() + geom_bar(aes(y = value, x = variable, fill = Taxa), data = g2,
                              stat="identity", position = "fill", colour = "black") +
  xlab("Platform") +
  ylab("Proportion") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(labels = c("RNA-seq", "ONT", "16S rRNA")) +
  ggtitle("Genus Representation by Platform") +
  theme(legend.title = element_text(size=14)) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) 

g_plot
g_plot + scale_fill_npg() + labs(fill = "Genera") + 
  theme(panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
  )

#community standard
mbs <- read.csv("mbs.csv")
m2 <- melt(mbs)
m2

mbs_plot <- ggplot() + geom_bar(aes(y = value, x = variable, fill = Species), data = m2,
                                stat="identity", position = "fill", colour = "black") +
  xlab("Confidence Score") +
  ylab("Proportion") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(labels = c("Actual", "0.1", "0.2", "0.3", "None")) +
  ggtitle("Kraken2 Confidence Interval Effect on Mock Community Detection") +
  theme(legend.title = element_text(size=14)) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) 

mbs_plot
mbs_plot + scale_fill_npg() + labs(fill = "Species") + 
  theme(panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
  )


