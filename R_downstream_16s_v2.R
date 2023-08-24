BiocManager::install("ashr")
library("phyloseq")
library("ggplot2")
library("ggbiplot")
library("ggalt")
library("ampvis2")
library("randomcoloR")
library("DESeq2")
library("EnhancedVolcano")
library("apeglm")
library("ashr")
library("Rmisc")
library("RColorBrewer")
install.packages("randomcoloR")
remotes::install_github("kasperskytte/ampvis2", Ncpus = 6)

setwd("C:/Users/dansa/Documents/Exjobb/summer_project")
# File paths to count & metadata:
path_to_metadata_file <- "metadata_chickens.csv"
path_metadata_lesions <- "metadata_chickens_lesion.csv"
# Read in files:
metadata <- read.delim(path_to_metadata_file, sep = ";", check.names=FALSE, stringsAsFactors=FALSE)
chicken_metadata <- metadata 
less_metadata <- cbind(metadata$`#SampleID`, metadata$TreatmentGroup, metadata$InfectionStatus, metadata$Timepoint, metadata$Description)
colnames(less_metadata) <- c("SampleID", "TreatmentGroup", "InfectionStatus", "Timepoint", "Description")
metadata_lesion <- read.delim(path_metadata_lesions, sep = ";", check.names=FALSE, stringsAsFactors=FALSE)
less_metadata_lesion <- cbind(metadata_lesion$`#SampleID`, metadata_lesion$TreatmentGroup, metadata_lesion$InfectionStatus, 
                              metadata_lesion$Timepoint, metadata_lesion$Description, metadata_lesion$LesionScore)
colnames(less_metadata_lesion) <- c("SampleID", "TreatmentGroup", "InfectionStatus", "Timepoint", "Description", "LesionScore")

otutable = "outputs_DADA/OTU_table.txt"
otus <- read.delim(otutable, sep = "\t", check.names=FALSE, stringsAsFactors=FALSE)
colnames(otus)[1] <- "OTU"
taxonomy = "outputs_DADA/taxonomy_table.txt"
taxtable <- read.delim(taxonomy, sep = "\t", check.names=FALSE, stringsAsFactors=FALSE, row.names = 1)
BIOMfilename <- "outputs_downstream_3/input_with_annotations.biom"
treefile <- "outputs_DADA/phylo/rooted_tree.nwk"
rare <- "outputs_DADA/graphs/rarecurve.txt"
rare <- read.delim(rare, sep = "\t", check.names=FALSE, stringsAsFactors=FALSE)

#Create phyloseq class:
phyloclass <- import_biom(BIOMfilename,
            treefilename=treefile, refseqfilename=NULL, refseqFunction=readDNAStringSet, refseqArgs=NULL,
            parseFunction=parse_taxonomy_default, parallel=FALSE, version=1.0)

#Set up the 8 categories:
metadata$Description <- paste(metadata$TreatmentGroup,metadata$InfectionStatus,metadata$Timepoint, sep="_")
phyloclass@sam_data$Description <- metadata$Description


#Normalize
microbiota <- transform_sample_counts(phyloclass, function(x) 1 * x / sum(x) )
#Transform into relative abundance
microbiota_r_a = transform_sample_counts(microbiota, function(x) x/sum(x)*100)

k <- length(levels(factor(microbiota@tax_table[,4])))
colors <- distinctColorPalette(k = 8)
colors <- c("#B4A7EB", "#9C746E", "#6EB1A7", "#9D4D81", "#9E95A1", "#62F29C", "#E7BF78",
                     "#8671BA", "#E7C7BC", "#A4E7CE", "#E867B1", "#DCEFE4", "#A34BDA",
                     "#4A6C9D", "#ED563F", "#92AA7D", "#D46EDE", "#E7B03B", "#AAEAAE",
                     "#E59AE8", "#DE85AA", "#7331ED", "#5E92E0", "#6DE867", "#E252DC",
                     "#DB33EB", "#53AF7D", "#A5D7DB", "#64E3E8", "#BFB4DC", "#B4E93E",
                     "#577AE5", "#E7E33E", "#39AFED", "#E9BEEA", "#5A9DB9", "#EAE4A6",
                     "#6CCEEC", "#87A9DC", "#51F144", "#9691EF", "#E6877D", "#E6AFC3",
                     "#60E6D0", "#E04C71", "#D4B289", "#6C3FA5", "#C086C4", "#D6884A",
                     "#A5CAEB", "#E935B3", "#63B658", "#A5E48E", "#97A24C", "#E3E688",
                     "#5CE9BA", "#E3DAEA", "#D4E6BD", "#C6E067", "#665BE5", "#4A7574",
                     "#A973DC")

colors_U_I <- c("#DD6A9A","#E04C71","#60E6D0","#6EB1A7","#E6877D","#ED563F", "#39AFED","#4A6C9D" ) 



#Plot each group
  #Plot each of the 8 categories
p <- plot_bar(microbiota, x="Description", y="Abundance", fill="Rank4")
p <- p + scale_fill_manual(name="Rank4", values=colors)+
  scale_color_manual(name="Rank4", values=colors)+
  geom_bar(aes(color=Rank4, fill=Rank4), stat="identity", position="stack")
  #Plot difference between algae/control groups
q <- plot_bar(microbiota, x="TreatmentGroup", y="Abundance", fill="Rank4")
q <- q +  scale_fill_manual(name="Rank4", values=colors)+
  scale_color_manual(name="Rank4", values=colors)+
  geom_bar(aes(color=Rank4, fill=Rank4), stat="identity", position="stack")
  #Plot difference infected/uninfected
r <- plot_bar(microbiota, x="InfectionStatus", y="Abundance", fill="Rank4")
r <- r +  scale_fill_manual(name="Rank4", values=colors)+
  scale_color_manual(name="Rank4", values=colors)+
  geom_bar(aes(color=Rank4, fill=Rank4), stat="identity", position="stack")
  #Plot difference D7/D10
s <- plot_bar(microbiota, x="Timepoint", y="Abundance", fill="Rank4")
s <- s +  scale_fill_manual(name="Rank4", values=colors)+
  scale_color_manual(name="Rank4", values=colors)+
  geom_bar(aes(color=Rank4, fill=Rank4), stat="identity", position="stack")
plot_list = list(p, q, r, s)
plot_path <- "plots/all_groups.png"
png(plot_path, height = 1200, width = 1800)
multiplot(plotlist = plot_list, layout = matrix(c(1,2,3,4), nrow = 2, byrow = TRUE))
dev.off()

t <- plot_bar(microbiota, x="SampleID", y="Abundance", fill="Rank4")+
  facet_grid(~TreatmentGroup, scales = "free", space = "free")
t <- t + scale_fill_manual(name="Rank4", values=colors)+
  scale_color_manual(name="Rank4", values=colors)+
  geom_bar(aes(color=Rank4, fill=Rank4), stat="identity", position="stack")
plot_path <- "plots/algae_control_bar_4.png"
png(plot_path, height = 1200, width = 1800)
t
dev.off()
u <- plot_bar(microbiota, x="SampleID", y="Abundance", fill="Rank4")+
  facet_grid(~InfectionStatus, scales = "free", space = "free")
u <- u + scale_fill_manual(name="Rank4", values=colors)+
  scale_color_manual(name="Rank4", values=colors)+
  geom_bar(aes(color=Rank4, fill=Rank4), stat="identity", position="stack")
plot_path <- "plots/uninfected_infected_bar_4.png"
png(plot_path, height = 1200, width = 1800)
u
dev.off()

v <- plot_bar(microbiota, x="SampleID", y="Abundance", fill="Rank4")+
  facet_grid(~factor(Description,levels=c("control_uninfected_D7", "control_uninfected_D10",
                                          "control_infected_D7", "control_infected_D10",
                                          "algae_uninfected_D7", "algae_uninfected_D10",
                                          "algae_infected_D7", "algae_infected_D10")), 
                     scales = "free", space = "free")
v <- v + scale_fill_manual(name="Rank4", values=colors)+
  scale_color_manual(name="Rank4", values=colors)+
  geom_bar(aes(color=Rank4, fill=Rank4), stat="identity", position="stack")+
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 18),legend.text = element_text(size = 15),
        plot.title = element_text(size = 20),strip.text.x = element_text(size = 12),
        legend.position = "right")+
        guides(fill=guide_legend(ncol=2))
  
plot_path <- "plots/all_samples_bar_4_text.png"
png(plot_path, height = 1250, width = 1800)
par(mar = c(5,5,5,5))
v
dev.off()

# Find top 20
top20otus = names(sort(taxa_sums(microbiota), TRUE)[1:2100])
taxtab20 = cbind(tax_table(microbiota), Rank4_2100 = NA)
taxtab20[top20otus, "Rank4_2100"] <- as(tax_table(microbiota)[top20otus, "Rank4"], 
                                      "character")
tax_table(microbiota) <- tax_table(taxtab20)
microbiota20 = prune_taxa(top20otus, microbiota)


#Plot only top 20

u <- plot_bar(microbiota20, x="SampleID", y="Abundance", fill="Rank4_2100")+
  facet_grid(~InfectionStatus, scales = "free", space = "free")
u <- u + scale_fill_manual(name="Rank4", values=colors)+
  scale_color_manual(name="Rank4", values=colors)+
  geom_bar(aes(color=Rank4_2100, fill=Rank4_2100), stat="identity", position="stack")
plot_path <- "plots/uninfected_infected_bar_4_top20.png"
png(plot_path, height = 1200, width = 1800)
u
dev.off()
t <- plot_bar(microbiota20, x="SampleID", y="Abundance", fill="Rank4_2100")+
  facet_grid(~TreatmentGroup, scales = "free", space = "free")
t <- t + scale_fill_manual(name="Rank4", values=colors)+
  scale_color_manual(name="Rank4", values=colors)+
  geom_bar(aes(color=Rank4_2100, fill=Rank4_2100), stat="identity", position="stack")
plot_path <- "plots/algae_control_bar_4_top20.png"
png(plot_path, height = 1200, width = 1800)
t
dev.off()
v <- plot_bar(microbiota20, x="SampleID", y="Abundance", fill="Rank4_2100")+
  facet_grid(~factor(Description,levels=c("control_uninfected_D7", "control_uninfected_D10",
                                          "control_infected_D7", "control_infected_D10",
                                          "algae_uninfected_D7", "algae_uninfected_D10",
                                          "algae_infected_D7", "algae_infected_D10")), 
             scales = "free", space = "free")
v <- v + scale_fill_manual(name="Rank4", values=colors)+
  scale_color_manual(name="Rank4", values=colors)+
  geom_bar(aes(color=Rank4_2100, fill=Rank4_2100), stat="identity", position="stack")+
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 18),legend.text = element_text(size = 15),
        plot.title = element_text(size = 20),strip.text.x = element_text(size = 14))
plot_path <- "plots/all_samples_bar_4_top20.png"
png(plot_path, height = 1200, width = 1800)
v
dev.off()



p <- plot_richness(phyloclass, x= "Description", measures=c("Chao1", "Shannon"), color = "Description")
p <- p + geom_boxplot(alpha=0.6)+ 
  aes(fill = Description, color = Description)+
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))+
  scale_fill_manual(name = "Description", values= c(colors_U_I))+
  scale_color_manual(name = "Description", values= c(colors_U_I))
plot_path <- "plots/diversity_category.png"
png(plot_path, height = 600, width = 1000)
p
dev.off()


#Ordination
pcoa_phy <- ordinate(phyloclass, method = "PCoA", "wunifrac")
p1 <- plot_ordination(phyloclass, pcoa_phy, type = "samples", color="Description")
p1 <- p1 + geom_encircle(aes(fill=Description, alpha = 0.1, s_shape = 1, expand = 0,)) + 
  geom_point(size=3)+
  scale_fill_manual(name="Description", values=c(colors_U_I))+
  scale_color_manual(name="Description", values=c(colors_U_I))+
  guides(alpha = FALSE)+
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 18),legend.text = element_text(size = 12),
        legend.title = element_text(size = 15))

p2 <-plot_ordination(phyloclass, pcoa_phy, type = "samples", color="TreatmentGroup")
p2 <- p2 + geom_encircle(aes(fill=TreatmentGroup, alpha = 0.1, s_shape = 1, expand = 0,)) + 
  geom_point(size=3)+
  guides(alpha = FALSE)+
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 18),legend.text = element_text(size = 12),
        legend.title = element_text(size = 15))

p3 <- plot_ordination(phyloclass, pcoa_phy, type = "samples", color="InfectionStatus")
p3 <- p3 + geom_encircle(aes(fill=InfectionStatus, alpha = 0.1, s_shape = 1, expand = 0,)) + 
  geom_point(size=3)+
  guides(alpha = FALSE)+
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 18),legend.text = element_text(size = 12),
        legend.title = element_text(size = 15))

p4 <- plot_ordination(phyloclass, pcoa_phy, type = "samples", color="Timepoint")
p4 <- p4 + geom_encircle(aes(fill=Timepoint, alpha = 0.1, s_shape = 1, expand = 0,)) + 
  geom_point(size=3)+
  guides(alpha = FALSE)+
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 18),legend.text = element_text(size = 12),
        legend.title = element_text(size = 15))

plot_list = list(p1, p2, p3, p4)
plot_path <- "plots/all_ord.png"
png(plot_path, height = 1200, width = 1800)
multiplot(plotlist = plot_list, layout = matrix(c(1,2,3,4), nrow = 2, byrow = TRUE))
dev.off()

plot_path <- "plots/pcoa_all_groups.png"
png(plot_path, height = 600, width = 800)
p1
dev.off()

plot_list = list(p2, p4)
plot_path <- "plots/pcoa_AC_D7D10.png"
png(plot_path, height = 800, width = 600)
multiplot(plotlist = plot_list)
dev.off()


# DEseq analysis

deseq_data<- phyloseq_to_deseq2(phyloclass, design= ~Description)
keep <- rowSums(counts(deseq_data)) >= 10
dds <- deseq_data[keep,]
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, type="poscounts")
#dds = estimateSizeFactors(dds, geoMeans = geoMeans)
dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized=TRUE)
res <- results(dds)

# Create contrasts
# Uninfected vs Infected
resA_UI_D10 <- results(dds, contrast=c("Description", "algae_infected_D10", "algae_uninfected_D10"))
resA_UI_D7 <- results(dds, contrast=c("Description", "algae_infected_D7", "algae_uninfected_D7"))
resC_UI_D10 <- results(dds, contrast=c("Description", "control_infected_D10", "control_uninfected_D10"))
resC_UI_D7 <- results(dds, contrast=c("Description", "control_infected_D7", "control_uninfected_D7"))

# Day 7 vs Day 10
resA_U_D710 <- results(dds, contrast=c("Description", "algae_uninfected_D10", "algae_uninfected_D7"))
resA_I_D710 <- results(dds, contrast=c("Description", "algae_infected_D10", "algae_infected_D7"))
resC_U_D710 <- results(dds, contrast=c("Description", "control_uninfected_D10", "control_uninfected_D7"))
resC_I_D710 <- results(dds, contrast=c("Description", "control_infected_D10", "control_infected_D7"))

# control vs algae
resAC_I_D7 <- results(dds, contrast=c("Description", "algae_infected_D7", "control_infected_D7"))
resAC_I_D10 <- results(dds, contrast=c("Description", "algae_infected_D10", "control_infected_D10"))
resAC_U_D7 <- results(dds, contrast=c("Description", "algae_uninfected_D7", "control_uninfected_D7"))
resAC_U_D10 <- results(dds, contrast=c("Description", "algae_uninfected_D10", "control_uninfected_D10"))

# Filter on adjusted p-value and add taxonomy

sigtabA_UI_D10 <- subset(resA_UI_D10, padj < 0.01)
sigtabA_UI_D10 = cbind(as(sigtabA_UI_D10, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabA_UI_D10), ], "matrix"))
sigtabA_UI_D7 <- subset(resA_UI_D7, padj < 0.01)
sigtabA_UI_D7 = cbind(as(sigtabA_UI_D7, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabA_UI_D7), ], "matrix"))
sigtabC_UI_D10 <- subset(resC_UI_D10, padj < 0.01)
sigtabC_UI_D10 = cbind(as(sigtabC_UI_D10, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabC_UI_D10), ], "matrix"))
sigtabC_UI_D7 <- subset(resC_UI_D7, padj < 0.01)
sigtabC_UI_D7 = cbind(as(sigtabC_UI_D7, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabC_UI_D7), ], "matrix"))

sigtabA_I_D710 <- subset(resA_I_D710, padj < 0.01)
sigtabA_I_D710 = cbind(as(sigtabA_I_D710, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabA_I_D710), ], "matrix"))
sigtabA_U_D710 <- subset(resA_U_D710, padj < 0.01)
sigtabA_U_D710 = cbind(as(sigtabA_U_D710, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabA_U_D710), ], "matrix"))
sigtabC_I_D710 <- subset(resC_I_D710, padj < 0.01)
sigtabC_I_D710 = cbind(as(sigtabC_I_D710, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabC_I_D710), ], "matrix"))
sigtabC_U_D710 <- subset(resC_U_D710, padj < 0.01)
sigtabC_U_D710 = cbind(as(sigtabC_U_D710, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabC_U_D710), ], "matrix"))

sigtabAC_I_D7 <- subset(resAC_I_D7, padj < 0.01)
sigtabAC_I_D7 = cbind(as(sigtabAC_I_D7, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabAC_I_D7), ], "matrix"))
sigtabAC_I_D10 <- subset(resAC_I_D10, padj < 0.01)
sigtabAC_I_D10 = cbind(as(sigtabAC_I_D10, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabAC_I_D10), ], "matrix"))
sigtabAC_U_D7 <- subset(resAC_U_D7, padj < 0.01)
sigtabAC_U_D7 = cbind(as(sigtabAC_U_D7, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabAC_U_D7), ], "matrix"))
sigtabAC_U_D10 <- subset(resAC_U_D10, padj < 0.01)
sigtabAC_U_D10 = cbind(as(sigtabAC_U_D10, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabAC_U_D10), ], "matrix"))



conditions <- c("A_UI_D10","A_UI_D7", "C_UI_D10", "C_UI_D7", 
                "A_U_D710", "A_I_D710", "C_U_D710", "C_I_D710",
                "AC_I_D7", "AC_I_D10", "AC_U_D7", "AC_U_D10")

reslist <- list(resA_UI_D10, resA_UI_D7, resC_UI_D10, resC_UI_D7, 
                resA_U_D710, resA_I_D710, resC_U_D710, resC_I_D710,
                resAC_I_D7, resAC_I_D10, resAC_U_D7, resAC_U_D10)

sigtabs <- list(sigtabA_UI_D10, sigtabA_UI_D7, sigtabC_UI_D10, sigtabC_UI_D7, 
                sigtabA_U_D710, sigtabA_I_D710, sigtabC_U_D710, sigtabC_I_D710,
                sigtabAC_I_D7, sigtabAC_I_D10, sigtabAC_U_D7, sigtabAC_U_D10)

#Plot volcanoes
logfc_threshold <- 1
fdr_threshold <- 0.05
plot_list <- list()
i <- 1
while (i <= length(sigtabs)){
  p <- EnhancedVolcano(sigtabs[[i]],
                     lab = rownames(sigtabs[[i]]),
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = conditions[i],
                     subtitle = "",
                     xlim = c(-25,25),
                     ylim = c(0,20),
                     ylab = bquote(~-Log[10]~italic(FDR)),
                     caption = "",
                     pointSize = 4,
                     labSize = 5,
                     pCutoff = fdr_threshold,
                     FCcutoff = logfc_threshold) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
        plot.subtitle = element_text(hjust = 0.5, size = 16),
        legend.text = element_text(size = 16))
  plot_list[[i]] <- p
  i = i + 1
} 
volcano_path <- "plots/all_chicken_plots_volcano.png"
png(volcano_path, height = 1200, width = 1800)
multiplot(plotlist = plot_list, layout = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), nrow = 3, byrow = TRUE))
dev.off()


# Plot lfc by taxonomy
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
i = 1
plot2_list <- list()
while (i <= length(sigtabs)){
  x = tapply(sigtabs[[i]]$log2FoldChange, sigtabs[[i]]$Rank3, function(x) max(x))
  x = sort(x, TRUE)
  sigtabs[[i]]$Rank3 = factor(as.character(sigtabs[[i]]$Rank3), levels=names(x))
  # Genus order
  x = tapply(sigtabs[[i]]$log2FoldChange, sigtabs[[i]]$Rank4, function(x) max(x))
  x = sort(x, TRUE)
  sigtabs[[i]]$Rank4 = factor(as.character(sigtabs[[i]]$Rank4), levels=names(x))

  p <- ggplot(sigtabs[[i]], aes(x=Rank4, y=log2FoldChange, color=Rank3)) + labs(title=conditions[i])+ geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  plot2_list[[i]] <- p
  i = i + 1
}
plot_path <- "plots/all_plots_lfc.png"
png(plot_path, height = 1200, width = 1800)
multiplot(plotlist = plot2_list, layout = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), nrow = 3, byrow = TRUE))
dev.off()


# Write out lists sorted on adjusted p value

j <- 1
while (j <= length(sigtabs)){
  sigtab <- sigtabs[[j]]
  ord_sigtab <- sigtab[order(sigtab$padj),]
  ord_sigtab["Taxonomy"] <- paste(ord_sigtab$Rank1, ord_sigtab$Rank2, ord_sigtab$Rank3,
                                  ord_sigtab$Rank4, ord_sigtab$Rank5, ord_sigtab$Rank6, sep = ", ")
  
  i <- 0
  while (i < nrow(ord_sigtab)) {
    i <- i + 1
   df_temp <- cbind(rownames(ord_sigtab[i,]), ord_sigtab[i,][,c("log2FoldChange", "pvalue", "padj", "Taxonomy")])
   if (i == 1) {
      all_de_genes <- df_temp
   } else {
     all_de_genes <- rbind(all_de_genes, df_temp)
   }
  }
  file_path <- paste("plots/all_de_genes_",c(conditions[j]), ".csv", sep="")
  write.table(all_de_genes, file_path, sep = ",", row.names = FALSE)
  j <- j + 1
}



# Correlation analysis 

trait_meta <- less_metadata_lesion[,c(1,2,3,4,6)]
rownames(trait_meta) <- trait_meta[,1]
trait_meta[trait_meta[,2]=="control",2] <- 0
trait_meta[trait_meta[,2]=="algae",2] <- 1
trait_meta[trait_meta[,3]=="uninfected",3] <- 0
trait_meta[trait_meta[,3]=="infected",3] <- 1
trait_meta[trait_meta[,4]=="D7",4] <- 0
trait_meta[trait_meta[,4]=="D10",4] <- 1
trait_meta <- trait_meta[,2:5]
names(trait_meta) <- c("TreatmentGroup", "InfectionStatus", "Timepoint", "LesionScore")
trait_df <- data.frame(trait_meta)
trait_df$TreatmentGroup <- as.numeric(trait_df$TreatmentGroup)
trait_df$InfectionStatus <- as.numeric(trait_df$InfectionStatus)
trait_df$Timepoint <- as.numeric(trait_df$Timepoint)
trait_df$LesionScore <- as.numeric(trait_df$LesionScore)


i = 1
j = 1
while (i <= nrow(normalized_counts)){
  countTraitCorP = cor.test(normalized_counts[i,], trait_df$TreatmentGroup, method = "pearson", use = "complete.obs")
  if ((is.na(countTraitCorP$p.value) == FALSE) & (countTraitCorP$p.value <= 0.05)){
    df_temp <- cbind(rownames(normalized_counts)[i], countTraitCorP$estimate, countTraitCorP$p.value, paste(c(taxtable[c(rownames(normalized_counts)[i]),])))
    if (j == 1) {
      all_sig_corr <- df_temp
    } else {
      all_sig_corr <- rbind(all_sig_corr, df_temp)
    }
    j = j + 1
  }
  i = i + 1
}
file_path <- paste("plots/tables/pearson_corr_algae.csv")
write.table(all_sig_corr, file_path, sep = ",", row.names = FALSE)

i = 1
j = 1
while (i <= nrow(normalized_counts)){
  countTraitCorP = cor.test(normalized_counts[i,], trait_df$TreatmentGroup, method = "spearman", use = "complete.obs")
  if ((is.na(countTraitCorP$p.value) == FALSE) & (countTraitCorP$p.value <= 0.05)){
    df_temp <- cbind(rownames(normalized_counts)[i], countTraitCorP$estimate, countTraitCorP$p.value, paste(c(taxtable[c(rownames(normalized_counts)[i]),])))
    if (j == 1) {
      all_sig_corr <- df_temp
    } else {
      all_sig_corr <- rbind(all_sig_corr, df_temp)
    }
    j = j + 1
  }
  i = i + 1
}
file_path <- paste("plots/tables/spearman_corr_algae.csv")
write.table(all_sig_corr, file_path, sep = ",", row.names = FALSE)




i = 1
j = 1
while (i <= nrow(normalized_counts)){
  countTraitCortest_lesion = cor.test(normalized_counts[i,], trait_df$LesionScore, method = "spearman", use = "complete.obs")
  if ((is.na(countTraitCortest_lesion$p.value) == FALSE) & (countTraitCortest_lesion$p.value <= 0.05)){
    df_temp <- cbind(rownames(normalized_counts)[i], countTraitCortest_lesion$estimate, countTraitCortest_lesion$p.value, paste(c(taxtable[c(rownames(normalized_counts)[i]),])))
    if (j == 1) {
      all_sig_corr <- df_temp
    } else {
      all_sig_corr <- rbind(all_sig_corr, df_temp)
    }
    j = j + 1
  }
  i = i + 1
}
file_path <- paste("plots/tables/spearman_corr_lesion.csv")
write.table(all_sig_corr, file_path, sep = ",", row.names = FALSE)


i = 1
j = 1

while (i <= nrow(normalized_counts)){
  countTraitCortest_lesion = cor.test(normalized_counts[i,], trait_df$LesionScore, method = "pearson", use = "complete.obs")
  if ((is.na(countTraitCortest_lesion$p.value) == FALSE) & (countTraitCortest_lesion$p.value <= 0.05)){
    df_temp <- cbind(rownames(normalized_counts)[i], countTraitCortest_lesion$estimate, countTraitCortest_lesion$p.value, paste(c(taxtable[c(rownames(normalized_counts)[i]),])))
    if (j == 1) {
      all_sig_corr <- df_temp
    } else {
      all_sig_corr <- rbind(all_sig_corr, df_temp)
    }
    j = j + 1
  }
  i = i + 1
}
file_path <- paste("plots/tables/pearson_corr_lesion.csv")
write.table(all_sig_corr, file_path, sep = ",", row.names = FALSE)
