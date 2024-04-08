##########################################################################
##########################################################################
# Project: Egr1 regulates regenerative senescence and cardiac repair
# Script purpose: Analyze the processed scRNA_seq data by Cell Ranger of 10x Genomics
# Usage example: 
# Author: Avital Sarusi-Portuguez (avital.sarusi-portuguez@weizmann.ac.il)
####     Lingling Zhang (lingling.zhang@weizmann.ac.il)
##########################################################################
##########################################################################
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(openxlsx)

### Upload saved Seurat object (named seurObj_all_score)
seurObj_all_score <- readRDS("./output/Seurat object/seurObj_all_score.rds")
DimPlot(seurObj_all_score, reduction = "umap", label = T)


### Create Dotplot
genes_of_interest <- c("Postn", "Pdgfra","Trem2","Ddr2","Pecam1","Mki67","Ccr2", "Acta2",
                       "Cd68", "Arg1", "Erbb4","Tnnt2","Timd4","Hbb-bs","S100a9", "Cthrc1",
                       "Krt19", "Satb1","Rgs5", "Lyve1", "Cadm2")
DotPlot(seurObj_all_score, features = genes_of_interest) + RotatedAxis() +
  theme(axis.text.y = element_text(size = 8))  # Adjust the text size for better fit


### Define Sen cells by applying a score threshold
# Set an Opt2 score (Gene list 1) threshold for seurObj_all_score; score expression scale: from -0.4 to 0.8
FetchData(seurObj_all_score_all_score,vars = "Opt2")
seurObj_all_score_all_score$Opt2_pos <- FetchData(seurObj_all_score_all_score,vars = "Opt2")$Opt2 > 0.3
DimPlot(seurObj_all_score_all_score,group.by = "Opt2_pos", cols = c("gray","firebrick1"),
        split.by = "Group", order = T)





### Recluster cardiac fibroblasts (FBs) clusters
Idents(seurObj_all_score)
# Levels: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
# How to subset cluster(s) and re-cluster them
# e.g. Subset FB clusters,seurObj is Seurat object,  0,1,3,5,10,15 are FB clusters
seurObj_all_FB <- subset(seurObj_all_score, idents = c(0,1,3,5,10,15))
DimPlot(seurObj_all_FB, reduction = "umap", label = TRUE)
# Re-cluster FB
seurObj_all_FB <- FindVariableFeatures(object = seurObj_all_FB, selection.method = "vst", nfeatures = 2000)
vars.to.regress <- c("nCount_RNA","percent.mt")
seurObj_all_FB <- ScaleData(object = seurObj_all_FB, vars.to.regress = vars.to.regress)
# Then do the PCA
seurObj_all_FB <- RunPCA(seurObj_all_FB, features = VariableFeatures(object = seurObj_all_FB))
# Elbow Plot
ElbowPlot(seurObj_all_FB)
# Change the number of PC, the UMAP will be changed__New
n_PC <- 25
seurObj_all_FB<- FindNeighbors(object = seurObj_all_FB, dims = 1:n_PC, reduction='pca',verbos=F)
seurObj_all_FB<- FindClusters(object = seurObj_all_FB, resolution = 0.5,verbos=F)
seurObj_all_FB<- RunUMAP(object = seurObj_all_FB, dims = 1:n_PC,verbos=F)
num_clusters <- table(seurObj_all_FB$seurat_clusters)
DimPlot(object = seurObj_all_FB, reduction = "umap", label = TRUE)


### Define Sen cells in WT_PBS or WT_Agrin groups
seurObj_all_FB_WT_Agrin <- subset(seurObj_all_FB, subset = Group %in% c("WT_Agrin"))
# Set an Opt2 score (Gene list 1) threshold for seurObj_all_FB_WT_Agrin; score expression scale: from -0.4 to 0.8
FetchData(seurObj_all_FB_WT_Agrin,vars = "Opt2")
seurObj_all_FB_WT_Agrin$Opt2_pos <- FetchData(seurObj_all_FB_WT_Agrin,vars = "Opt2")$Opt2 > 0.3
DimPlot(seurObj_all_FB_WT_Agrin,group.by = "Opt2_pos", cols = c("gray","firebrick1"),
        split.by = "Group", order = T)


### Violin Plot to show the gene in Sen_Cell Vs Non-Sen cells
#seurObj_all_FB_WT_Agrin
Non_sen_cells <- subset(seurObj_all_FB_WT_Agrin, Opt2_pos == FALSE)
Sen_cells <- subset(seurObj_all_FB_WT_Agrin, Opt2_pos == TRUE)
# Combine TRUE_Sen cells and FALSE_Non-sen cells
#seurObj_all_FB_WT_Agrin
seurObj_all_FB_WT_Agrin_Sen <- merge(y = Non_sen_cells, x = Sen_cells )
# First, ensure that 'Opt2_pos' is a factor and has the levels 'TRUE' and 'FALSE'
seurObj_all_FB_WT_Agrin_Sen$Opt2_pos <- factor(seurObj_all_FB_WT_Agrin_Sen$Opt2_pos, levels = c("FALSE", "TRUE"))
# Now, specify colors for each level
colors <- c("Non_sen_cells" = "gray", "Sen_cells" = "red")
# Specify the genes you want to plot  
genes_of_interest <- c("Il6", "Apoe")  # Replace with your actual genes
# Create a violin plot for the specified genes with the defined colors
VlnPlot(seurObj_all_FB_WT_Agrin_Sen, features = genes_of_interest, group.by = "Opt2_pos",
        combine = TRUE, split.by = "Opt2_pos", pt.size = 0, cols = colors)



####For creating Volcano plot
library(readxl)
library(ggrepel)
library(dplyr)
library(ggplot2)
# Read the Excel file
For_Volcano <- read_excel("./output/Opt2_0.3_WT_Agrin_FB_T vs F_MAST_For GSEA.xlsx")

# Create the "DE" column using mutate
For_Volcano <- For_Volcano %>%
  mutate(DE = ifelse(avg_log2FC > 1 & p_val_adj < 0.05, "UP",
                     ifelse(avg_log2FC < -1 & p_val_adj < 0.05, "DOWN", "NA")))

For_Volcano <- For_Volcano %>% mutate(p_val_adj_mod = -log10(p_val_adj))

For_Volcano$p_val_adj_mod[is.infinite(For_Volcano$p_val_adj_mod)] <- max(For_Volcano$p_val_adj_mod[!is.infinite(For_Volcano$p_val_adj_mod)])+5
p <- For_Volcano %>%
  
  ggplot(aes(x = avg_log2FC, y = p_val_adj_mod, color = DE)) +
  
  geom_point(size = 0.8) + scale_color_manual(values = c("blue","gray","red")) +
  
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "None") +
  
  geom_vline(xintercept = c(-1,1), linetype = "dashed") + geom_hline(yintercept = 1.3, linetype = "dashed") +
  
  ylab("-log10(p_val_adj)")

# Add labels for specific genes
p <- p +
  geom_text(data = For_Volcano %>% filter(Gene %in% c("Il6","Slc10a6","Ccl7")),
            aes(x = avg_log2FC, y = p_val_adj_mod, label = Gene),
            hjust = -0.2, vjust = -0.5, size = 3, color = "black")


#Save as pdf (maybe you need to adjuste the width and hight)

pdf(file = "Volcano.pdf",width = 4,height = 4) 

print(p)

dev.off()




####For creating Scatter plot
library(readxl)
library(ggrepel)
library(dplyr)
library(ggplot2)
# Read the Excel file
Scatter_Pot <- read_excel("./output/Opt2_0.3_WT_Agrin_FB_T vs F_MAST_For GSEA.xlsx")
# Convert 'Pathways' to factor with desired levels (reversed order)
Scatter_Pot$Pathways <- factor(Scatter_Pot$Pathways, levels = rev(unique(Scatter_Pot$Pathways)))

# Order the dataframe by the 'Pathways' column
Scatter_Pot <- Scatter_Pot[order(Scatter_Pot$Pathways),]

# Create the ggplot
IPA_Scatter_Pot <- ggplot(Scatter_Pot, aes(x = `-log10(p-value)`, y = Pathways)) +
  geom_point(aes(color = `z-score`), size = 4) +
  scale_color_continuous(low = "green4", high = "red2", limits = c(-3.5, 3.5),
                         guide = guide_colorbar(reverse = TRUE)) +
  theme(legend.position = "bottom") +
  theme_bw() +
  coord_cartesian(xlim = c(0, 5.6)) +  # Set x-axis limits
  geom_vline(xintercept = 0.6, linetype = "dashed", color = "black")  # Add a dashed line at x = 1.3

# Save as pdf
pdf(file = "IPA_Scatter_Pot.pdf", width = 8, height = 4) 
print(IPA_Scatter_Pot)
dev.off()











