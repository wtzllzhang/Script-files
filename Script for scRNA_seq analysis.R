##########################################################################
##########################################################################
# Project: Egr1 regulates regenerative senescence and cardiac repair
# Script purpose: Analyze the processed scRNA_seq data by Cell Ranger of 10x Genomics
# Usage example: 
# Author: Avital Sarusi-Portuguez (avital.sarusi-portuguez@weizmann.ac.il)
####     Lingling Zhang (lingling.zhang@weizmann.ac.il)
##########################################################################
##########################################################################
#load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(openxlsx)
library(reshape2)
library(future)
library(hdf5r)

### Use parallel option
plan("multisession", workers = 4)
options(future.globals.maxSize = 2000 * 1024^2)

### cutoffs for QC
filter_mt = 25
n_MAD = 4

###################
## Create Seurat ##
###################
# Create Seurat object for each sample, with metadata and filtering for low quality cells
# done iteratively
cellranger_dir <- c(Exp1 = "/gpfs/units/bioinformatics/projects/wis/tzahore/INCPMPM-19250/single_cell_10x/SingleCell_output/cellranger",
                    Exp2 = "/gpfs/units/bioinformatics/projects/wis/tzahore/INCPMPM-20028/single_cell_10x/SingleCell_output/cellranger")

samples_names <- c( "944"="Egr1_WT_Agrin",
                    "957"="Egr1_WT_PBS_1",
                    "960"="Egr1_WT_PBS_2",
                    "965"="Egr1_KO_PBS_1",
                    "966"="Egr1_KO_PBS_2",
                    "968"="Egr1_KO_Agrin",
                    "903"="Egr1_KO_Agrin_1",
                    "910"="Egr1_WT_PBS",
                    "912"="Egr1_KO_Agrin_2",
                    "915"="Egr1_KO_PBS",
                    "917"="Egr1_WT_Agrin_1",
                    "919"="Egr1_WT_Agrin_2")

seurObj_list= list()
meta.data.all <- list()
thresholds.all <- list()
for (Exp in c("Exp1","Exp2")){
  samples <- list.dirs(cellranger_dir[Exp],recursive = F,full.names = F)
  samples <- samples[samples %in% names(samples_names)]
  samples_file <- paste0(cellranger_dir[Exp],"/",samples, "/outs/filtered_feature_bc_matrix.h5")
  for (i in 1:length(samples_file)){
    new_name <- samples_names[samples[i]]
    count.data = Read10X_h5(samples_file[i])
    seurObj = CreateSeuratObject(counts = count.data, min.cells = 0, min.features = 0,
                                 project=new_name)
    seurObj[['percent.mt']] = PercentageFeatureSet(object=seurObj,pattern = "^mt-")
    seurObj <- RenameCells(seurObj,new.names = paste0(colnames(x = seurObj),"_",new_name))
    seurObj$Exp <- Exp
    Meta_to_add <- seurObj@meta.data %>% select(orig.ident) %>% tidyr::separate(orig.ident, into = c("gene","cond","treat"),extra = "drop") %>%
      select(cond,treat)
    seurObj <- AddMetaData(seurObj,Meta_to_add)
    seurObj$Group <- paste0(seurObj$cond,"_", seurObj$treat)
    #Save info before filter for plotting
    meta.data.all[[new_name]] <- seurObj@meta.data
    #find cutoffs and subset
    meta.data.thresh <- seurObj@meta.data %>% group_by(orig.ident) %>% summarise_if(is.numeric,c(mad,median,mean))
    thresholds <- meta.data.thresh %>% transmute(Sample = orig.ident,
                                                 percent.mt_MADS_above = percent.mt_fn2+ percent.mt_fn1*n_MAD,
                                                 nFeature_MADS_above = nFeature_RNA_fn2+ nFeature_RNA_fn1*n_MAD,
                                                 nFeature_MADS_below = nFeature_RNA_fn2- nFeature_RNA_fn1*n_MAD,
                                                 nCount_MADS_below = nCount_RNA_fn2- nCount_RNA_fn1*n_MAD,
                                                 nCount_MADS_above = nCount_RNA_fn2+ nCount_RNA_fn1*n_MAD)
    thresholds[thresholds<0]<-0
    thresholds.all[[new_name]] <- thresholds
    seurObj <- subset(x= seurObj,subset =  nFeature_RNA < thresholds$nFeature_MADS_above & 
                        nFeature_RNA >= thresholds$nFeature_MADS_below &
                        nCount_RNA < thresholds$nCount_MADS_above & 
                        nCount_RNA >= thresholds$nCount_MADS_below &
                        percent.mt <= min(thresholds$percent.mt_MADS_above,filter_mt))
    seurObj_list[[new_name]] <- seurObj
  }
  rm(count.data)
}

###merge
seurObj <- seurObj_list[[1]]

if (length(samples_names)>1){
  seurObj <- merge(x = seurObj_list[[1]], y = seurObj_list[2:length(samples_names)], project = "Merged")
}

#################################
## Data reduction and Clusters ##
#################################
seurObj <- NormalizeData(object = seurObj)
seurObj <- FindVariableFeatures(object = seurObj, selection.method = "vst", nfeatures = 2000)
vars.to.regress <- c("nCount_RNA","percent.mt")
seurObj <- ScaleData(object = seurObj, vars.to.regress = vars.to.regress)
seurObj <- RunPCA(seurObj, features = VariableFeatures(object = seurObj),verbos=F)
#decide on number of PCs and resolution for clusters
n_PC <- 22
res_to_use <- 0.5
seurObj<- FindNeighbors(object = seurObj, dims = 1:n_PC, reduction='pca',verbos=F)
seurObj<- FindClusters(object = seurObj, resolution = res_to_use,verbos=F)
seurObj<- RunUMAP(object = seurObj, dims = 1:n_PC,verbos=F)

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











