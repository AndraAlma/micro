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
install.packages("randomcoloR")
remotes::install_github("kasperskytte/ampvis2", Ncpus = 6)


# File paths to count & metadata:
path_to_metadata_file <- "metadata_chickens.csv"
# Read in files:
metadata <- read.delim(path_to_metadata_file, sep = ";", check.names=FALSE, stringsAsFactors=FALSE)
chicken_metadata <- metadata 
less_metadata <- cbind(metadata$`#SampleID`, metadata$TreatmentGroup, metadata$InfectionStatus, metadata$Timepoint, metadata$Description)
colnames(less_metadata) <- c("SampleID", "TreatmentGroup", "InfectionStatus", "Timepoint", "Description")
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
colors <- distinctColorPalette(k = k)
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
  facet_grid(~Description, scales = "free", space = "free")
v <- v + scale_fill_manual(name="Rank4", values=colors)+
  scale_color_manual(name="Rank4", values=colors)+
  geom_bar(aes(color=Rank4, fill=Rank4), stat="identity", position="stack")
plot_path <- "plots/all_samples_bar_4.png"
png(plot_path, height = 1200, width = 1800)
v
dev.off()

# Find top 20
top20otus = names(sort(taxa_sums(microbiota), TRUE)[1:2100])
taxtab20 = cbind(tax_table(microbiota), Rank4_2100 = NA)
taxtab20[top20otus, "Rank4_2100"] <- as(tax_table(microbiota)[top20otus, "Rank4"], 
                                      "character")
tax_table(microbiota) <- tax_table(taxtab20)
microbiota20 = prune_taxa(top20otus, microbiota)

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
  facet_grid(~Description, scales = "free", space = "free")
v <- v + scale_fill_manual(name="Rank4", values=colors)+
  scale_color_manual(name="Rank4", values=colors)+
  geom_bar(aes(color=Rank4_2100, fill=Rank4_2100), stat="identity", position="stack")
plot_path <- "plots/all_samples_bar_4_top20.png"
png(plot_path, height = 1200, width = 1800)
v
dev.off()

myPalette=colors
names(myPalette) = paste("row", 1:length(colors))

# create data that references the palette
colorKey = data.frame(colorName=names(myPalette))
# plot with ggplot, referencing the palette
ggplot(data=colorKey, aes(x=1, y = 1:nrow(colorKey), fill=colorName, label=colorName)) +
  geom_tile() +
  scale_fill_manual(values = myPalette) +
  theme_void()+
  theme(legend.position="none") + 
  geom_text()


#Ordination
pcoa_phy <- ordinate(phyloclass, method = "PCoA", "wunifrac")
p1 <- plot_ordination(phyloclass, pcoa_phy, type = "samples", color="Description")
p1 <- p1 + geom_encircle(aes(fill=Description, alpha = 0.1, s_shape = 1, expand = 0,)) + geom_point(size=3)
p2 <-plot_ordination(phyloclass, pcoa_phy, type = "samples", color="TreatmentGroup")
p2 <- p2 + geom_encircle(aes(fill=TreatmentGroup, alpha = 0.1, s_shape = 1, expand = 0,)) + geom_point(size=3)
p3 <- plot_ordination(phyloclass, pcoa_phy, type = "samples", color="InfectionStatus")
p3 <- p3 + geom_encircle(aes(fill=InfectionStatus, alpha = 0.1, s_shape = 1, expand = 0,)) + geom_point(size=3)
p4 <- plot_ordination(phyloclass, pcoa_phy, type = "samples", color="Timepoint")
p4 <- p4 + geom_encircle(aes(fill=Timepoint, alpha = 0.1, s_shape = 1, expand = 0,)) + geom_point(size=3)

plot_list = list(p1, p2, p3, p4)
plot_path <- "plots/all_ord.png"
png(plot_path, height = 1200, width = 1800)
multiplot(plotlist = plot_list, layout = matrix(c(1,2,3,4), nrow = 2, byrow = TRUE))
dev.off()

#Create Amp object
data_amp <- amp_load(
  otutable = otus,
  metadata = less_metadata
)
amp_rarefaction_curve(
  data_amp,
  color_by = "TreatmentGroup")
amp_rarefaction_curve(
  data_amp,
  color_by = "InfectionStatus")


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
head(sigtabA_UI_D10)
sigtabA_UI_D7 <- subset(resA_UI_D7, padj < 0.01)
sigtabA_UI_D7 = cbind(as(sigtabA_UI_D7, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabA_UI_D7), ], "matrix"))
head(sigtabA_UI_D7)
sigtabC_UI_D10 <- subset(resC_UI_D10, padj < 0.01)
sigtabC_UI_D10 = cbind(as(sigtabC_UI_D10, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabC_UI_D10), ], "matrix"))
head(sigtabC_UI_D10)
sigtabC_UI_D7 <- subset(resC_UI_D7, padj < 0.01)
sigtabC_UI_D7 = cbind(as(sigtabC_UI_D7, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabC_UI_D7), ], "matrix"))
head(sigtabC_UI_D7)

sigtabA_I_D710 <- subset(resA_I_D710, padj < 0.01)
sigtabA_I_D710 = cbind(as(sigtabA_I_D710, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabA_I_D710), ], "matrix"))
head(sigtabA_I_D710)
sigtabA_U_D710 <- subset(resA_U_D710, padj < 0.01)
sigtabA_U_D710 = cbind(as(sigtabA_U_D710, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabA_U_D710), ], "matrix"))
head(sigtabA_U_D710)
sigtabC_I_D710 <- subset(resC_I_D710, padj < 0.01)
sigtabC_I_D710 = cbind(as(sigtabC_I_D710, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabC_I_D710), ], "matrix"))
head(sigtabC_I_D710)
sigtabC_U_D710 <- subset(resC_U_D710, padj < 0.01)
sigtabC_U_D710 = cbind(as(sigtabC_U_D710, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabC_U_D710), ], "matrix"))
head(sigtabC_U_D710)

sigtabAC_I_D7 <- subset(resAC_I_D7, padj < 0.01)
sigtabAC_I_D7 = cbind(as(sigtabAC_I_D7, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabAC_I_D7), ], "matrix"))
head(sigtabAC_I_D7)
sigtabAC_I_D10 <- subset(resAC_I_D10, padj < 0.01)
sigtabAC_I_D10 = cbind(as(sigtabAC_I_D10, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabAC_I_D10), ], "matrix"))
head(sigtabAC_I_D10)
sigtabAC_U_D7 <- subset(resAC_U_D7, padj < 0.01)
sigtabAC_U_D7 = cbind(as(sigtabAC_U_D7, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabAC_U_D7), ], "matrix"))
head(sigtabAC_U_D7)
sigtabAC_U_D10 <- subset(resAC_U_D10, padj < 0.01)
sigtabAC_U_D10 = cbind(as(sigtabAC_U_D10, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabAC_U_D10), ], "matrix"))
head(sigtabAC_U_D10)



#Shrink LFC
resshrink_A_UI_D10 <- lfcShrink(dds, contrast=c("Description", "algae_uninfected_D10", "algae_infected_D10"),type="ashr")
sigtabshrink_A_UI_D10 <- subset(resshrink_A_UI_D10, padj < 0.01)
sigtabshrink_A_UI_D10 = cbind(as(sigtabshrink_A_UI_D10, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabshrink_A_UI_D10), ], "matrix"))
resshrink_A_UI_D7 <- lfcShrink(dds, contrast=c("Description", "algae_uninfected_D7", "algae_infected_D7"),type="ashr")
sigtabshrink_A_UI_D7 <- subset(resshrink_A_UI_D7, padj < 0.01)
sigtabshrink_A_UI_D7 = cbind(as(sigtabshrink_A_UI_D7, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabshrink_A_UI_D7), ], "matrix"))
resshrink_C_UI_D10 <- lfcShrink(dds, contrast=c("Description", "control_uninfected_D10", "control_infected_D10"),type="ashr")
sigtabshrink_C_UI_D10 <- subset(resshrink_C_UI_D10, padj < 0.01)
sigtabshrink_C_UI_D10 = cbind(as(sigtabshrink_C_UI_D10, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabshrink_C_UI_D10), ], "matrix"))
resshrink_C_UI_D7 <- lfcShrink(dds, contrast=c("Description", "control_uninfected_D7", "control_infected_D7"),type="ashr")
sigtabshrink_C_UI_D7 <- subset(resshrink_C_UI_D7, padj < 0.01)
sigtabshrink_C_UI_D7 = cbind(as(sigtabshrink_C_UI_D7, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabshrink_C_UI_D7), ], "matrix"))

resshrink_A_U_D710 <- lfcShrink(dds, contrast=c("Description", "algae_uninfected_D7", "algae_uninfected_D10"),type="ashr")
sigtabshrink_A_U_D710 <- subset(resshrink_A_U_D710, padj < 0.01)
sigtabshrink_A_U_D710 = cbind(as(sigtabshrink_A_U_D710, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabshrink_A_U_D710), ], "matrix"))
resshrink_A_I_D710 <- lfcShrink(dds, contrast=c("Description", "algae_infected_D7", "algae_infected_D10"),type="ashr")
sigtabshrink_A_I_D710 <- subset(resshrink_A_I_D710, padj < 0.01)
sigtabshrink_A_I_D710 = cbind(as(sigtabshrink_A_I_D710, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabshrink_A_I_D710), ], "matrix"))
resshrink_C_U_D710 <- lfcShrink(dds, contrast=c("Description", "control_uninfected_D7", "control_uninfected_D10"),type="ashr")
sigtabshrink_C_U_D710 <- subset(resshrink_C_U_D710, padj < 0.01)
sigtabshrink_C_U_D710 = cbind(as(sigtabshrink_C_U_D710, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabshrink_C_U_D710), ], "matrix"))
resshrink_C_I_D710 <- lfcShrink(dds, contrast=c("Description", "control_infected_D7", "control_infected_D10"),type="ashr")
sigtabshrink_C_I_D710 <- subset(resshrink_C_I_D710, padj < 0.01)
sigtabshrink_C_I_D710 = cbind(as(sigtabshrink_C_I_D710, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabshrink_C_I_D710), ], "matrix"))

resshrink_AC_I_D7 <- lfcShrink(dds, contrast=c("Description", "algae_infected_D7", "control_infected_D7"),type="ashr")
sigtabshrink_AC_I_D7 <- subset(resshrink_AC_I_D7, padj < 0.01)
sigtabshrink_AC_I_D7 = cbind(as(sigtabshrink_AC_I_D7, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabshrink_AC_I_D7), ], "matrix"))
resshrink_AC_I_D10 <- lfcShrink(dds, contrast=c("Description", "algae_infected_D10", "control_infected_D10"),type="ashr")
sigtabshrink_AC_I_D10 <- subset(resshrink_AC_I_D10, padj < 0.01)
sigtabshrink_AC_I_D10 = cbind(as(sigtabshrink_AC_I_D10, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabshrink_AC_I_D10), ], "matrix"))
resshrink_AC_U_D7 <- lfcShrink(dds, contrast=c("Description", "algae_uninfected_D7", "control_uninfected_D7"),type="ashr")
sigtabshrink_AC_U_D7 <- subset(resshrink_AC_U_D7, padj < 0.01)
sigtabshrink_AC_U_D7 = cbind(as(sigtabshrink_AC_U_D7, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabshrink_AC_U_D7), ], "matrix"))
resshrink_AC_U_D10 <- lfcShrink(dds, contrast=c("Description", "algae_uninfected_D10", "control_uninfected_D10"),type="ashr")
sigtabshrink_AC_U_D10 <- subset(resshrink_AC_U_D10, padj < 0.01)
sigtabshrink_AC_U_D10 = cbind(as(sigtabshrink_AC_U_D10, "data.frame"), as(tax_table(phyloclass)[rownames(sigtabshrink_AC_U_D10), ], "matrix"))


conditions <- c("A_UI_D10","A_UI_D7", "C_UI_D10", "C_UI_D7", 
                "A_U_D710", "A_I_D710", "C_U_D710", "C_I_D710",
                "AC_I_D7", "AC_I_D10", "AC_U_D7", "AC_U_D10")

reslist <- list(resA_UI_D10, resA_UI_D7, resC_UI_D10, resC_UI_D7, 
                resA_U_D710, resA_I_D710, resC_U_D710, resC_I_D710,
                resAC_I_D7, resAC_I_D10, resAC_U_D7, resAC_U_D10)

sigtabs <- list(sigtabA_UI_D10, sigtabA_UI_D7, sigtabC_UI_D10, sigtabC_UI_D7, 
                sigtabA_U_D710, sigtabA_I_D710, sigtabC_U_D710, sigtabC_I_D710,
                sigtabAC_I_D7, sigtabAC_I_D10, sigtabAC_U_D7, sigtabAC_U_D10)

sigtabshrink <- list(sigtabshrink_A_UI_D10, sigtabshrink_A_UI_D7, sigtabshrink_C_UI_D10, sigtabshrink_C_UI_D7, 
                     sigtabshrink_A_U_D710, sigtabshrink_A_I_D710, sigtabshrink_C_U_D710, sigtabshrink_C_I_D710,
                     sigtabshrink_AC_I_D7, sigtabshrink_AC_I_D10, sigtabshrink_AC_U_D7, sigtabshrink_AC_U_D10)


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
volcano_path <- "plots/all_chicken_plots_micro.png"
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


# Plot shrunk lfc
i = 1
plot3_list <- list()
while (i <= length(sigtabshrink)){
  x = tapply(sigtabshrink[[i]]$log2FoldChange, sigtabshrink[[i]]$Rank3, function(x) max(x))
  x = sort(x, TRUE)
  sigtabshrink[[i]]$Rank3 = factor(as.character(sigtabshrink[[i]]$Rank3), levels=names(x))
  # Genus order
  x = tapply(sigtabshrink[[i]]$log2FoldChange, sigtabshrink[[i]]$Rank4, function(x) max(x))
  x = sort(x, TRUE)
  sigtabshrink[[i]]$Rank4 = factor(as.character(sigtabshrink[[i]]$Rank4), levels=names(x))
  
  p <- ggplot(sigtabshrink[[i]], aes(x=Rank4, y=log2FoldChange, color=Rank3)) + labs(title=conditions[i])+ geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) 
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  plot3_list[[i]] <- p
  i = i + 1
}
plot_path <- "plots/all_plots_shrinklfc.png"
png(plot_path, height = 1200, width = 1800)
multiplot(plotlist = plot3_list, layout = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), nrow = 3, byrow = TRUE))
dev.off()
