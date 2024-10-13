
# https://www.science.org/doi/10.1126/science.aba7721

rm(list = ls())

#setwd()
getwd()

library(tidyverse)
library(Matrix)
library(Seurat)
library(biomaRt)
library(scDblFinder)
library(Nebulosa)

source("chromopack.R")

#############################
###   genes of interest   ###
#############################

# genes of interest
int_genes_df <- read.delim("data/int_genes.csv", sep = ";")
int_genes <- int_genes_df$genesofinterest

# cell type markers
markers_df <- read.delim("data/celltypes_markers.csv", sep = ";")
celltypes_markers <- markers_df$marker


#############################
###    data importation   ###
#############################
# counts and metadata were subset for only cerebrum cells and 
# "H27432" fetus was removed due to 47XY +18. File saved as sparse matrix RDS.
# features with same ensembl id were added together
# Ensembl gene id were also change to gene symbol

meta <- read.delim("data/meta.tsv", sep = "\t")
meta <- meta %>% 
  column_to_rownames("sample")

cts <- readRDS("data/counts.rds")

scbrain <- CreateSeuratObject(counts = cts)
scbrain <- AddMetaData(object = scbrain, metadata = meta)
Idents(scbrain) <- "Fetus_id"

rm(cts, meta) # cleaning


############################
###   filtering cells    ###
############################

scbrain[["percent.mt"]] <- PercentageFeatureSet(scbrain, pattern = "^MT-") # mitochondrial reads
scbrain[["percent.ribo"]] <- PercentageFeatureSet(scbrain, pattern = "^RP[SL]") # Ribosomal percentage

# Plotting: features, reads, and percent mt reads before QC
violin <- VlnPlot(scbrain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 2, pt.size = 0.1, alpha = 0.05)
ggsave("results/exploratory_analysis/before_violin.png", plot = violin, height = 9, width = 9)

# Plotting: nCount vs nFeatures and nCount vs percent.mt before QC
plot1 <- FeatureScatter(scbrain, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scbrain, feature1 = "nCount_RNA", feature2 = "percent.ribo")
plot3 <- FeatureScatter(scbrain, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot4 <- FeatureScatter(scbrain, feature1 = "nFeature_RNA", feature2 = "percent.mt")
combined_plot <- (plot1|plot2)/(plot3|plot4)
ggsave("results/exploratory_analysis/before_scatter.png", plot = combined_plot, height = 8, width = 11)


################ Filtering
scbrain <- subset(scbrain, subset = nCount_RNA > 200 & nCount_RNA < 60000) # reads per cell
scbrain <- subset(scbrain, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000) # features per cell
scbrain <- subset(scbrain, subset = percent.mt < 10) # Percentage of mitochondrial genes 

# doublets detection
sce <- as.SingleCellExperiment(scbrain) # Convert Seurat object to SingleCellExperiment object
sce <- scDblFinder(sce) # Run scDblFinder
scbrain$doublet <- sce$scDblFinder.class # Add doublet information to Seurat object metadata
scbrain <- subset(scbrain, subset = doublet == "singlet") # removing doublets


# Plotting: features, reads, and percent mt reads after QC
violin <- VlnPlot(scbrain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 2, pt.size = 0.1, alpha = 0.05)
ggsave("results/exploratory_analysis/after_violin.png", plot = violin, height = 9, width = 9)

# Plotting: nCount vs nFeatures and nCount vs percent.mt after QC
plot1 <- FeatureScatter(scbrain, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scbrain, feature1 = "nCount_RNA", feature2 = "percent.ribo")
plot3 <- FeatureScatter(scbrain, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot4 <- FeatureScatter(scbrain, feature1 = "nFeature_RNA", feature2 = "percent.mt")
combined_plot <- (plot1|plot2)/(plot3|plot4)
ggsave("results/exploratory_analysis/after_scatter.png", plot = combined_plot, height = 8, width = 11)


###############################
###   filtering features    ###
###############################

scbrain <- scbrain[!grepl("MT-", rownames(scbrain)),] # Removing mitochondrial genes

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # features information
all_features <- getBM(
  attributes = c("external_gene_name","gene_biotype","chromosome_name","start_position","end_position"),
  filters = "external_gene_name",
  values = rownames(scbrain),
  mart = ensembl
)

all_features <- all_features %>% 
  filter(
    complete.cases(.), # dropping features with NAs
    gene_biotype == "protein_coding", # only protein_coding genes
    chromosome_name %in% c(as.character(seq(1, 22))), # only chromosomal genes 
    !duplicated(external_gene_name) | !duplicated(external_gene_name, fromLast = TRUE) # choosing one of duplicated gene names
  ) %>% 
  mutate(
    chromosome_name = factor(chromosome_name, levels = c(as.character(seq(1, 22)))),
    gene_length = end_position - start_position,
    avg_position = (end_position + start_position)/2
  )

# filtering features from seurat object based on aforementioned filters
scbrain <- scbrain[rownames(scbrain) %in% all_features$external_gene_name,]

# Removing genes with zero expression
counts_matrix <- scbrain[["RNA"]]$counts
scbrain[["RNA"]]$counts <- counts_matrix[rowSums(counts_matrix) > 0, ]

# updating all features
all_features <- all_features %>% 
  filter(external_gene_name %in% rownames(scbrain))


########################
###    Processing    ###
########################

scbrain <- NormalizeData(scbrain)

scbrain <- FindVariableFeatures(scbrain) # default 2000 features

top10 <- head(VariableFeatures(scbrain), 10) # Highlight the 10 most highly variable genes

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(scbrain) + ggtitle("All genes")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + ggtitle("Top 10 highly variable genes")
ggsave("results/exploratory_analysis/variableFeatures.png", plot = plot1 + plot2, width = 10, height = 5)

scbrain <- ScaleData(scbrain, features = rownames(scbrain))

scbrain <- RunPCA(scbrain, features = VariableFeatures(object = scbrain))

elbow <- ElbowPlot(scbrain, ndims = 50) + # help define number of PCs
  theme(
    plot.background = element_rect(fill = "white"),  # Set background color to white
    panel.background = element_rect(fill = "white"),  # Set panel background color to white
    panel.border = element_blank(),  # Remove panel borders
    axis.line = element_line(color = "black"),  # Set axis line color
    legend.background = element_rect(fill = "white"),  # Set legend background color to white
  )
ggsave("results/exploratory_analysis/elbow.png", plot = elbow, height = 4, width = 4)

scbrain <- RunUMAP(scbrain, dims = 1:10)

scbrain <- FindNeighbors(scbrain, dims = 1:10, k.param = 30) # 10PCs, 30KNN

scbrain <- FindClusters(scbrain, resolution = 0.1)

Idents(scbrain) <- "seurat_clusters"

#########################################
###      Batch effect verification    ###
#########################################

# UMAP
for (i in c("seurat_clusters","Fetus_id","Sex","Development_day")) {
  clusters <- DimPlot(scbrain, reduction = "umap", group.by = i, pt.size = 0.1, label = TRUE, raster=FALSE) +
    ggtitle("") +
    labs(color = i) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    coord_fixed() +
    theme(
      legend.text = element_text(size = 26),  # Adjust legend text size and style
      legend.title = element_text(size = 28, face = "bold"),  # Adjust legend text size and style
      axis.text = element_text(size = 16),    # Adjust axis text size and style
      axis.title = element_text(size = 22, face = "bold"),   # Adjust axis title size and style
      axis.line = element_line(linewidth = 1.5), # Customize axis lines
      axis.ticks = element_line(linewidth = 1.5),                 # Adjust axis ticks size
    )
  
  ggsave(paste0("results/exploratory_analysis/umap_",i,".png"), plot = clusters, width = 12, height = 12)
}

# Bar plots
for (i in c("Fetus_id","Sex","Development_day")) {
  aux <- scbrain@meta.data %>% 
    dplyr::select(seurat_clusters, !!sym(i)) %>% 
    group_by(seurat_clusters, !!sym(i)) %>% 
    summarise(num_cells = n())%>%
    group_by(seurat_clusters) %>%
    mutate(pct = num_cells / sum(num_cells) * 100) %>% 
    ungroup()
  
  # absolute value bar plot
  bar_plot <- barfunc(aux, x_axis = "seurat_clusters", y_axis = "num_cells", fill_groups = i, rotate_x = 45)
  ggsave(paste0("results/exploratory_analysis/bar_",i,".png"), plot = bar_plot, width = 7, height = 7)
  
  # pct value bar plot
  bar_plot <- barfunc(aux, x_axis = "seurat_clusters", y_axis = "pct", fill_groups = i, annot_loc = "center", annot_text = "num_cells", rotate_x = 45)
  ggsave(paste0("results/exploratory_analysis/pctbar_",i,".png"), plot = bar_plot, width = 7, height = 7)
}


####################################
###       Checking markers       ###
####################################

# Verifying assay type. Should be RNA
DefaultAssay(scbrain)

# Changing identity to cluster column
Idents(scbrain) <- "seurat_clusters"
DefaultAssay(scbrain) <- "RNA"

# hline annotation
auxmark <- markers_df %>% 
  group_by(celltype) %>% 
  summarise(n_markers = n()) %>% 
  mutate(celltype = factor(celltype, levels = unique(markers_df$celltype))) %>% 
  arrange(celltype) %>%  
  mutate(
    y = length(unique(scbrain@meta.data$seurat_clusters))*1.05,
    yend = y,
    x = 1, # (just for the loop to work)
    xend = n_markers, # provisory (just for the loop to work)
    colour = "#000000"
  ) %>% 
  as.data.frame()
for (i in 2:nrow(auxmark)) {
  auxmark[i,"x"] <- auxmark[i-1,"xend"] + 1
  auxmark[i,"xend"] <- auxmark[i,"x"] + auxmark[i,"n_markers"] -1
}

# canonical markers
dot_plot <- DotPlot(
  assay = "RNA",
  object = scbrain, 
  features = celltypes_markers, 
  cols = c("blue","red"),
  cluster.idents = T,
  scale = T 
)+
  theme(
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_blank(),  # Remove panel borders
    axis.line = element_line(color = "black"),  # Set axis line color
    legend.background = element_rect(fill = "white"),  # Set legend background color to white
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
    axis.text.y = element_text(size = 17, face = "bold"),
    #legend.position = "none" # no legend
  ) +
  coord_cartesian(clip = "off")
for (i in 1:nrow(auxmark)) {
  dot_plot <- dot_plot + 
    annotate(
      "segment", # type of annotation to be added
      x = auxmark[i,"x"],
      y = auxmark[i,"y"],
      xend = auxmark[i,"xend"],
      yend = auxmark[i,"yend"],
      #colour = auxmark[i,"colour"],
      
    ) +
    annotate(
      "text", # type of annotation to be added
      x = (auxmark[i,"x"] + auxmark[i,"xend"])/2, 
      y = auxmark[i,"y"] + 0.3, 
      label = auxmark[i,"celltype"], 
      #hjust = 1.1, 
      #vjust = -0.5, 
      color = "black",
      fontface = "bold",
      size = 4
    )
}

ggsave("results/exploratory_analysis/canonical_markers.png", plot = dot_plot, height = 5, width = 22)


###############################
###       Annotating        ###
###############################
Idents(scbrain) <- "seurat_clusters"

scbrain <- RenameIdents(scbrain, c(
  `0` = "Exc",
  `1` = "Inh",
  `2` = "NPC",
  `3` = "MSN",
  `4` = "Exc",
  `5` = "Inh",
  `6` = "RGL",
  `7` = "Ast",
  `8` = "Inh",
  `9` = "End"
  )
)

scbrain@meta.data <- scbrain@meta.data %>%
  mutate(
    seurat_clusters = case_when(
      seurat_clusters == "0" ~ "Exc",
      seurat_clusters == "1" ~ "Inh",
      seurat_clusters == "2" ~ "NPC",
      seurat_clusters == "3" ~ "MSN", # medium spiny neurons
      seurat_clusters == "4" ~ "Exc",
      seurat_clusters == "5" ~ "Inh",
      seurat_clusters == "6" ~ "RGL",
      seurat_clusters == "7" ~ "Ast",
      seurat_clusters == "8" ~ "Inh",
      seurat_clusters == "9" ~ "End"
    )
  )


###################################
###       Annotated UMAP        ###
###################################

# remove legend, increase label size and BOLD
umap_plot <- DimPlot(
  scbrain,
  group.by = "seurat_clusters",
  reduction = "umap", 
  
  label = F, 
  pt.size = 0.1, 
  raster = FALSE 
  #label.size = 12, 
  #cols = c("Oli" = brewer.pal(7, "Set3")[1], "Exc" = brewer.pal(8, "Set3")[8], "OPC" = brewer.pal(7, "Set3")[3], 
  #         "Inh" = brewer.pal(7, "Set3")[4], "Ast" = brewer.pal(7, "Set3")[7], "Mic" = brewer.pal(7, "Set3")[6], 
  #         "End" = brewer.pal(7, "Set3")[5])
) +
  NoAxes() + # removes axis
  theme(
    legend.position = "none" # Removes the legend
  ) + 
  ggtitle(NULL)
ggsave("results/exploratory_analysis/ANNOTATED_umap.png", plot = umap_plot, height = 12, width = 12)


####################################################
###      Repeating Batch effect verification     ###
####################################################

clusters <- DimPlot(scbrain, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = TRUE, raster=FALSE) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  coord_fixed() +
  theme(
    legend.text = element_text(size = 26),  # Adjust legend text size and style
    legend.title = element_text(size = 28, face = "bold"),  # Adjust legend text size and style
    axis.text = element_text(size = 16),    # Adjust axis text size and style
    axis.title = element_text(size = 22, face = "bold"),   # Adjust axis title size and style
    axis.line = element_line(linewidth = 1.5), # Customize axis lines
    axis.ticks = element_line(linewidth = 1.5),                 # Adjust axis ticks size
  )
ggsave("results/exploratory_analysis/ANNOTATED_umap_names.png", plot = clusters, width = 12, height = 12)

# Bar plots
for (i in c("Fetus_id","Sex","Development_day")) {
  aux <- scbrain@meta.data %>% 
    dplyr::select(seurat_clusters, !!sym(i)) %>% 
    group_by(seurat_clusters, !!sym(i)) %>% 
    summarise(num_cells = n())%>%
    group_by(seurat_clusters) %>%
    mutate(pct = num_cells / sum(num_cells) * 100) %>% 
    ungroup()
  
  # absolute value bar plot
  bar_plot <- barfunc(aux, x_axis = "seurat_clusters", y_axis = "num_cells", fill_groups = i, rotate_x = 45)
  ggsave(paste0("results/exploratory_analysis/ANNOTATED_bar_",i,".png"), plot = bar_plot, width = 7, height = 7)
  
  # pct value bar plot
  bar_plot <- barfunc(aux, x_axis = "seurat_clusters", y_axis = "pct", fill_groups = i, annot_loc = "center", annot_text = "num_cells", rotate_x = 45)
  ggsave(paste0("results/exploratory_analysis/ANNOTATED_pctbar_",i,".png"), plot = bar_plot, width = 7, height = 7)
}


##################################################
###          Repeating markers plots           ###
##################################################

Idents(scbrain) <- "seurat_clusters"

# hline annotation
auxmark <- markers_df %>% 
  group_by(celltype) %>% 
  summarise(n_markers = n()) %>% 
  mutate(celltype = factor(celltype, levels = unique(markers_df$celltype))) %>% 
  arrange(celltype) %>%  
  mutate(
    y = length(unique(scbrain@meta.data$seurat_clusters)) + 0.25,
    yend = y,
    x = 1, # (just for the loop to work)
    xend = n_markers, # provisory (just for the loop to work)
    colour = "#000000"
  ) %>% 
  as.data.frame()
for (i in 2:nrow(auxmark)) {
  auxmark[i,"x"] <- auxmark[i-1,"xend"] + 1
  auxmark[i,"xend"] <- auxmark[i,"x"] + auxmark[i,"n_markers"] -1
}


# canonical markers
dot_plot <- DotPlot(
  object = scbrain, 
  features = celltypes_markers, 
  cols = c("blue","red"),
  cluster.idents = T,
  #cluster.features = T, # não existe!!!
  scale = T # zscore 
)+
  theme(
    plot.background = element_rect(fill = "white"),  # Set background color to white
    panel.background = element_rect(fill = "white"),  # Set panel background color to white
    #panel.grid.major = element_blank(),  # Remove major grid lines
    #panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove panel borders
    axis.line = element_line(color = "black"),  # Set axis line color
    legend.background = element_rect(fill = "white"),  # Set legend background color to white
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
    axis.text.y = element_text(size = 17, face = "bold"),
    #legend.position = "none" # no legend
  )
for (i in 1:nrow(auxmark)) {
  dot_plot <- dot_plot + 
    annotate(
      "segment", # type of annotation to be added
      x = auxmark[i,"x"],
      y = auxmark[i,"y"],
      xend = auxmark[i,"xend"],
      yend = auxmark[i,"yend"],
      #colour = auxmark[i,"colour"],
      
    ) +
    annotate(
      "text", # type of annotation to be added
      x = (auxmark[i,"x"] + auxmark[i,"xend"])/2, 
      y = auxmark[i,"y"] + 0.15, 
      label = auxmark[i,"celltype"], 
      #hjust = 1.1, 
      #vjust = -0.5, 
      color = "black",
      fontface = "bold",
      size = 4
    )
}

ggsave("results/exploratory_analysis/ANNOTATED_dotplot.png", plot = dot_plot, height = 8, width = 22)


##############################################
###    Differential expression analysis    ###
##############################################

DElist <- list()

for (i in unique(scbrain$seurat_clusters)){
  DElist[[i]] <- FindMarkers(
    scbrain,
    ident.1 = i,
    logfc.threshold = 0.0,
    min.pct = 0.0,
    min.cells.feature = 0
  )
}

DElist <- imap(DElist, function(i, cell_type){
  i <- i %>% 
    rownames_to_column("gene") %>% 
    mutate(
      celltype = cell_type,
      DEG = case_when(
        avg_log2FC >= 1 & p_val_adj < 0.05 ~ "UP",
        avg_log2FC <= -1 & p_val_adj < 0.05 ~ "DOWN",
        TRUE ~ "NO" # else none of the other conditions are true
      ),
      alterat = case_when(
        avg_log2FC > 0 & p_val_adj < 0.05 ~ "UP",
        avg_log2FC < 0 & p_val_adj < 0.05 ~ "DOWN",
        TRUE ~ "NO" # else none of the other conditions are true
      )
    )
  return(i)
})

DE_df <- do.call(rbind, DElist) %>% 
  filter(gene %in% c("CACNA1H","CACNG4","CREBBP","DAB1","FYN","NCK2","PDK1","RELN","SRC"))



##################################
###       Expression UMAP      ###
##################################

for(i in c("CACNA1H","CACNG4","CREBBP","DAB1","FYN","NCK2","PDK1","RELN","SRC")){
  umap_feature <- plot_density(
    scbrain,
    features = i,
    reduction = "umap",
    size = 0.2
  ) +
    scale_color_viridis_c(option = "inferno") +
    theme(
      legend.position = "none"
    )
  
  ggsave(paste0("results/feature_",i,".png"), plot = umap_feature, height = 5, width = 5, limitsize = FALSE)
}

###############################
###        Dot plot         ###
###############################

dot_plot <- DotPlot(
  scbrain, 
  features = c("CACNA1H","CACNG4","CREBBP","DAB1","FYN","NCK2","PDK1","RELN","SRC"), 
  group.by = 'seurat_clusters',
  cols = c("blue","red"),
  cluster.idents = T,
  scale = T
)+
  labs(
    x = NULL,
    y = NULL
  ) +
  theme(
    plot.background = element_rect(fill = "white"),  # Set background color to white
    panel.background = element_rect(fill = "white"),  # Set panel background color to white
    #panel.grid.major = element_blank(),  # Remove major grid lines
    #panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove panel borders
    axis.line = element_line(color = "black"),  # Set axis line color
    legend.background = element_rect(fill = "white"),  # Set legend background color to white
    axis.text.x = element_text(angle = 45, hjust = 1) # rotate x axis test 45 degrees
    #legend.position = "none" # no legend
  )

ggsave("results/dotplot.png", plot = dot_plot, height = 5, width = 8)


