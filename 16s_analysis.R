library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(DESeq2)
library(vegan)
library(phyloseq)
library(dplyr)
library(reshape2)
library(biomformat)
library(FSA)
library(lme4)

##### 1.validation experiment 16S data analysis ######

taxa.table<- read.csv2("taxon_table.csv",sep = ",",row.names=1)
taxa.table<-as.matrix(taxa.table)
md = read.csv("group_table_1567sample.csv",header = T,sep = ',')
asv.table<-read.csv2("otu_table_1567sample.csv",row.names=1,sep = ",")

colnames(asv.table) = paste(md$sample, md$genotype, md$rep, md$treatment, md$maizeline,md$Compartment, md$program, md$group1, md$group2, sep = '_')
md$sampleID = colnames(asv.table)

# create phyloseq object
OTU.tab = otu_table(asv.table, taxa_are_rows = TRUE)
TAX.tab = tax_table(taxa.table)
rownames(md) = md$sampleID
sample.tab = sample_data(md,errorIfNULL = TRUE)
ps = phyloseq(OTU.tab, TAX.tab, sample.tab)

# Associate the tree file with your phyloseq file
library(ape)
phy_tree <- read.tree("newtree.nwk")
ps<- merge_phyloseq(ps, phy_tree)

# keep only bacteria and remove Mitochondria and Chloroplast ASVs 
ps.prefiltered = subset_taxa(ps, Kingdom == 'Bacteria' & 
                               Family !='Mitochondria' & 
                               Order !='Chloroplast' )
# check if every sample has >1000 reads
range(colSums(otu_table(ps.prefiltered)) )
range(rowSums(otu_table(ps.prefiltered)) )

sum(taxa_sums(ps.prefiltered)) / sum(taxa_sums(ps))    # =0.96
total_reads <- sum(colSums(otu_table(ps.prefiltered)))
avg_reads_per_sample <- mean(colSums(otu_table(ps.prefiltered)))


# ra > 0.1% at least in 1% samples
ps.RA = transform_sample_counts(ps.prefiltered, function(x) x/sum(x))
ps.ra.abund = filter_taxa(ps.RA, function(x) sum(x > 0.001) >= 0.01*nsamples(ps.RA), TRUE)  #  745 taxa
ps.abund = subset_taxa(ps.prefiltered, taxa_names(ps.prefiltered)%in%taxa_names(ps.ra.abund))
ps.RA = transform_sample_counts(ps.abund, function(x) x/sum(x))
range(colSums(otu_table(ps.abund)) )  # all samples >1000 reads
sum(taxa_sums(ps.abund)) / sum(taxa_sums(ps.prefiltered))   # 0.95
saveRDS(ps.abund, 'ps.abund1567_0.001_0.01.RDS')
saveRDS(ps.RA, 'ps.RAabund1567_0.001_0.01.RDS')
saveRDS(ps.prefiltered, 'ps.prefiltered1567_0.001_0.01.RDS')


##### 2.beta diversity######
library("limma")
# DESeq2 Normalization 
ps.abund<-readRDS('ps.abund1567_0.001_0.01.RDS')
bac.deseq = phyloseq_to_deseq2(ps.abund, ~group1)
bac.deseq = estimateSizeFactors(bac.deseq, type="poscount")
bac.vst = varianceStabilizingTransformation(bac.deseq, blind=F)   
bac.vst.norm = ps.abund
otu_table(bac.vst.norm) = otu_table(assay(bac.vst), taxa_are_rows = T)
#Eliminate batch effect
batch <- colData(bac.deseq)$rep
vsd_batch_corrected <- removeBatchEffect(assay(bac.vst), batch=batch)
assay(bac.vst) <- vsd_batch_corrected

# Round negative values up to zeroes, to enable Bray-Curtis calculations
bac.vst.norm = transformSampleCounts(bac.vst.norm,function(x) ifelse(x<0,0,x)) 
# do PCoA analysis
PCoA.bray.vst = ordinate(bac.vst.norm, method = "PCoA", distance = "bray")
pcoa.df = cbind.data.frame(PCoA.bray.vst$vectors[, 1:5], sample_data(bac.vst.norm))
pct = PCoA.bray.vst$values$Relative_eig[1:2]
percentage = paste( c("PCoA1", "PCoA2"), " (", paste( as.character(round(pct*100,2)), "%", ")", sep=""), sep = "" )
pcoa.df$group1 = factor(pcoa.df$group1, levels = c('LI', 'LH','HI','HH', 'LBS','HBS'))
# PCoA plot
pal <- c("#ffd490","#f38d67", "#a4d0db","#014f91", "brown","black")
# define my plot theme
mytheme = theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"), plot.title = element_text(size = 12),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 12)) 
fig1b = ggplot(data= pcoa.df, aes(x = Axis.1, y = Axis.2, color = group1,fill = group1, shape = treatment)) +
  geom_point(alpha= 0.6,size=5) +
  scale_shape_manual(values = c(15:17)) +
  scale_color_manual(values = c(pal)) +
  scale_fill_manual(values = alpha(pal, alpha= 0.8)) + 
  stat_ellipse(data=pcoa.df,geom = "polygon",level=0.95,linetype = 2,aes(fill=group1),alpha=0.1,show.legend = T) +
  ggtitle("PCoA plot ") + 
  mytheme +
  xlab(percentage[1]) + ylab(percentage[2]) +
  coord_fixed()

fig1b
ggsave(fig1b, file = "1567bray_PCOA.pdf", width = 15, height = 10, units = "cm")

## PERMANOVA test ##
#test if treatment affect the beta diversity
bray.dist.part = phyloseq::distance(bac.vst.norm, method = "bray") 
permanova_treatment= adonis2(bray.dist.part~treatment,
                          data = as(sample_data(bac.vst.norm), "data.frame"), 
                          permutations = 999, parallel = 18)
permanova_treatment
