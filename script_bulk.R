rm(list = ls())     # clears Environment 

library(tidyverse) #do it all
library(DESeq2)  
library(ComplexHeatmap)
library(pheatmap)
library(data.table)
library(org.Hs.eg.db) #Homo sapiens OrgDb
library(biomaRt) #gene annotations
library(RColorBrewer) # for pretty colors
library(ggrepel)
library(umap)
library(mdp)
library(kableExtra)
library(readxl) #excel
library(writexl) #excel
library(VennDiagram) #venn digrams
library(clusterProfiler) #functional analysis
library(pathview)
library(cowplot)
library(GOplot)
library(enrichplot)
library(gage)
library(gageData)
library(ggraph) # network plot
library(igraph)
library(sva) #combatseq
library(webshot2)
library(karyoploteR) #GWAS plots
library(vsn)
library(RColorBrewer)
library(scales)
library(ReactomePA) #REACTOME
library(metap) # for pvalues
library(devtools) # for downloading github repositories
library(WGCNA)
library(CEMiTool) 
library(gt)


#source("~/chromopack.R")

if (!file.exists("results")){   
  dir.create("results")
  dir.create("results/exploratory_analysis")
  dir.create("results/tables")
}

##############################
###     counts matrix      ###
##############################

# counts matrix
normal <- read.delim("data/counts.csv", sep = ",", header = FALSE, row.names = 1)

# features meta
all_features <- read.delim("data/features.csv", sep = ",")

normal <- normal %>% 
  set_names(sapply(1:524, function(i) paste0("S", i), simplify = FALSE)) %>% # changing samples names
  mutate(features = all_features$gene_symbol) %>% # adding features names
  group_by(features) %>%
  summarise(across(everything(), sum)) %>% # sum features with same name
  ungroup() %>% 
  column_to_rownames("features") %>% 
  mutate_all(as.numeric)

# biomart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
all_features <- getBM(
  attributes = c("entrezgene_id", "external_gene_name", "gene_biotype", "chromosome_name", "start_position", "end_position"),
  filters = "external_gene_name",
  values = rownames(normal),
  mart = ensembl
)

all_features <- all_features %>% 
  filter(!duplicated(external_gene_name)) %>% # remove duplicates
  filter(gene_biotype == "protein_coding") %>% #only protein coding genes
  filter(chromosome_name %in% c(as.character(seq(1, 22)), "X", "Y")) %>% # only human chr genes
  mutate(chromosome_name = factor(chromosome_name, levels = c(as.character(seq(1, 22)), "X", "Y"))) %>% # ordering
  drop_na()

# filtering features from normal based on aforementioned filters
normal <- normal %>% 
  filter(rownames(.) %in% all_features$external_gene_name)


##############################
###        metadata        ###
##############################

# importing metadata
meta <- read.delim("data/metadata.csv", sep = ",", row.names = 1)

meta <- meta %>% 
  mutate(samples = colnames(normal)) %>% 
  column_to_rownames("samples") %>% 
  mutate_all(factor) %>% 
  dplyr::rename("region" = "structure_acronym") %>% 
  mutate(
    samples = rownames(.),
    # grouping
    grouped_age = case_when(
      age %in% c("8 pcw","9 pcw","12 pcw") ~ "8_12PCW",
      age %in% c("13 pcw","16 pcw","17 pcw","19 pcw") ~ "13_19PCW",
      age %in% c("21 pcw","24 pcw","25 pcw","26 pcw","35 pcw","37 pcw") ~ "20_40PCW",
      age %in% c("4 mos","10 mos","1 yrs","2 yrs") ~ "0_2YEARS",
      age %in% c("3 yrs","4 yrs","8 yrs") ~ "3_10YEARS",
      age %in% c("11 yrs","13 yrs","15 yrs","18 yrs","19 yrs","21 yrs") ~ "11_22YEARS",
      age %in% c("23 yrs","30 yrs","36 yrs","37 yrs","40 yrs") ~ "23_40YEARS"
    ),
    grouped_region = case_when(
      region %in% c("M1C-S1C", "DFC", "MFC", "OFC", "VFC", "M1C") ~ "Frontal_lobe",
      region %in% c("A1C", "ITC", "TCx", "STC") ~ "Temporal_lobe",
      region %in% c("IPC", "S1C", "PCx") ~ "Parietal_lobe",
      region %in% c("V1C", "Ocx") ~ "Occipital_lobe",
      region %in% c("AMY", "HIP", "STR", "MD", "DTH") ~ "Subcortical_structures",
      region %in% c("CB","CBC","URL") ~ "Cerebellum",
      region %in% c("MGE","LGE","CGE") ~ "Primordial_structures"
    ),
    # ordering
    age = factor(age, levels = c("8 pcw","9 pcw","12 pcw","13 pcw","16 pcw","17 pcw","19 pcw","21 pcw","24 pcw","25 pcw",
                                 "26 pcw","35 pcw","37 pcw","4 mos","10 mos","1 yrs","2 yrs","3 yrs","4 yrs","8 yrs","11 yrs",
                                 "13 yrs","15 yrs","18 yrs","19 yrs","21 yrs","23 yrs","30 yrs","36 yrs","37 yrs","40 yrs")),
    region = factor(region, levels = c("M1C-S1C", "DFC", "MFC", "OFC", "VFC", "M1C","A1C", "ITC", "TCx", "STC","IPC", "S1C", 
                                       "PCx","V1C", "Ocx","AMY", "HIP", "STR", "MD", "DTH","CB","CBC","MGE","URL","LGE","CGE")),
    grouped_age = factor(grouped_age, levels = c("8_12PCW", "13_19PCW", "20_40PCW", "0_2YEARS", "3_10YEARS", "11_22YEARS", "23_40YEARS")),
    grouped_region = factor(grouped_region, levels = c("Frontal_lobe", "Temporal_lobe", "Parietal_lobe", "Occipital_lobe", 
                                                       "Subcortical_structures", "Cerebellum", "Primordial_structures"))
  )


##############################
###    genes of interest   ###
##############################

# importing genes of interest
intGenes_df <- read.delim("data/int_genes.csv", sep = ";")
intGenes <- intGenes_df$genesofinterest

# patient mutations
patients <- read.delim("data/patientsmutations.csv", sep = ";")
patients$patient <- factor(patients$patient, levels = c("F2688-1","F10832-1","F11463-1","1-1098-003","2-1259-004","5-5057-003","7-0276-003","AU2168301","AU3756301","AU4027306"))
patients_genes <- unlist(lapply(patients$mutated_genes, function(x) unlist(strsplit(x, "\\|"))))
patients_genes <- patients_genes[!duplicated(patients_genes)]

# SFARI genes
sfari <- read.delim("data/SFARI_genes.csv", sep = ",")
sfari <- sfari[!isNA(sfari$gene.score),] # filtering only genes with scores
#sfari$gene.score <- factor(sfari$gene.score)
sfari_genes <- sfari$gene.symbol

# Bruno genes
bruno <- read.delim("data/bruno_genes.csv", sep = ";")
bruno_genes <- bruno$gene

##############################
###     quality control    ###
##############################

normal <- normal %>% 
  filter(rowSums(.) > 10) %>% # at least 10 counts in all samples
  filter(rowMeans(. == 0) < 0.3) # equal to zero in less than 30% of samples

# genes must be expressed in all combinations of grouped_region and grouped_age
# this important when calculating log2fc and wilcoxon (Inf, -Inf, Nan if mean == 0)
aux <- normal
aux["grouped_age",] <- meta[, "grouped_age"]
aux["grouped_region",] <- meta[, "grouped_region"]
aux <- aux %>% 
  t() %>% 
  as.data.frame() %>% 
  pivot_longer(cols = -c(grouped_age, grouped_region), names_to = "gene", values_to = "expression") %>% 
  as.data.frame() %>% 
  mutate(expression = as.numeric(expression)) %>% 
  group_by(grouped_age, grouped_region, gene) %>%
  summarise(sum_expression = sum(expression)) %>% 
  ungroup() %>% 
  as.data.frame() %>% 
  filter(sum_expression == 0) %>% 
  filter(!duplicated(gene)) %>% 
  pull(gene)

normal <- normal %>% 
  filter(!(rownames(.) %in% aux)) # expressed in all combinations

# Adjusting all_features
all_features <- all_features %>% 
  filter(external_gene_name == rownames(normal))

# Registering sample genes for downstream functional analysis
sample_genes <- rownames(normal)
all_features$entrezgene_id <- as.character(all_features$entrezgene_id)
sample_genes_entrez <- all_features[all_features$external_gene_name %in% sample_genes,]$entrezgene_id


###################################
###     exploratory analysis    ###
###################################

########### data distribution

ggsave("results/exploratory_analysis/age_region.png", plot = sampDistri_func(meta, "age","region"), width = 30, height = 20)

ggsave("results/exploratory_analysis/groupedAge_groupedRegion.png", plot = sampDistri_func(meta, "grouped_age","grouped_region"), width = 30, height = 20)

ggsave("results/exploratory_analysis/groupedRegion_gender.png", plot = sampDistri_func(meta, "grouped_region", "gender"), width = 30, height = 15)

ggsave("results/exploratory_analysis/groupedAge_gender.png", plot = sampDistri_func(meta, "grouped_age", "gender"), width = 20, height = 10) 

ggsave("results/exploratory_analysis/groupedRegion_donorName.png", plot = sampDistri_func(meta, "grouped_region", "donor_name"), width = 15, height = 30) 

ggsave("results/exploratory_analysis/groupedAge_donorName.png", plot = sampDistri_func(meta, "grouped_age", "donor_name"), width = 15, height = 30) 

ggsave("results/exploratory_analysis/age_groupedRegion.png", plot = sampDistri_func(meta, "age", "grouped_region"), width = 45, height = 15)

ggsave("results/exploratory_analysis/groupedAge_region.png", plot = sampDistri_func(meta, "grouped_age", "region"), width = 20, height = 30)
# NOTE: y axis is order in decreasing number of samples (not factor levels order)


############# sample to sample distance

sampleDists <- dist(t(normal))
sampleDistMatrix <- as.matrix(sampleDists)

sampToSamp <- pheatmap(
  sampleDistMatrix,
  main = "Sample to sample distances",
  annotation_col = meta[,c("grouped_region","grouped_age")],
  annotation_row = meta[,c("grouped_region","grouped_age")],
  annotation_names_col = TRUE,
  annotation_names_row = FALSE,
  show_colnames = FALSE,
  show_rownames = FALSE,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  treeheight_row = 0, # hide dendograms
  treeheight_col = 0, # hide dendograms
  col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
  angle_col = 45,
  annotation_colors = list(
    grouped_age = c("8_12PCW" = brewer.pal(7, "Set3")[1], "13_19PCW" = brewer.pal(7, "Set3")[2], "20_40PCW" = brewer.pal(7, "Set3")[3], "0_2YEARS" = brewer.pal(7, "Set3")[4], "3_10YEARS" = brewer.pal(7, "Set3")[5], "11_22YEARS" = brewer.pal(7, "Set3")[6], "23_40YEARS" = brewer.pal(7, "Set3")[7]), 
    grouped_region = c("Frontal_lobe" = brewer.pal(7, "Dark2")[1], "Temporal_lobe" = brewer.pal(7, "Dark2")[2], "Parietal_lobe" = brewer.pal(7, "Dark2")[3], "Occipital_lobe" = brewer.pal(7, "Dark2")[4], "Subcortical_structures" = brewer.pal(7, "Dark2")[5], "Cerebellum" = brewer.pal(7, "Dark2")[6], "Primordial_structures" = brewer.pal(7, "Dark2")[7])
  )
)
ggsave("results/exploratory_analysis/sampToSamp.png", plot = sampToSamp, height = 8, width = 10)


################ dimensional reduction
##### identify highly variable features 

hvf <- data.frame(features = rownames(normal)) %>% 
  mutate(
    features = rownames(normal),
    mean = rowMeans(normal),
    variance = apply(normal, 1, var),
    dispersion = variance/mean,
    log_mean = log1p(mean),
    log_dispersion = log1p(dispersion)
  )

fit <- loess(hvf$log_dispersion ~ hvf$log_mean)

hvf <- hvf %>% 
  mutate(
    expected_dispersion = predict(fit, log_mean),
    diff_dispersion = log_dispersion - expected_dispersion
  ) %>% 
  arrange(-diff_dispersion)

# regression plot
png("results/exploratory_analysis/variable_features_plot.png", width=800, height=600)
plot(hvf$log_mean, hvf$log_dispersion, pch=20, cex=0.6, col="blue", 
     xlab="Log-transformed Mean Expression", ylab="Log-transformed Dispersion",
     main="Loess Fit of Log-transformed Dispersion vs. Mean Expression")
# Add the fitted trend line
lines(sort(hvf$log_mean), hvf$expected_dispersion[order(hvf$log_mean)], col="red", lwd=2)
dev.off()

hvf_2000 <- hvf %>% 
  slice(1:2000) %>% 
  pull(features)

hvf_5000 <- hvf %>% 
  slice(1:5000) %>% 
  pull(features)

##### UMAP
umap_result <- umap(t(normal %>% filter(rownames(.) %in% hvf_2000)))

# Get the UMAP results
umap_data <- data.table(umap_result$layout)
colnames(umap_data) <- c("UMAP1", "UMAP2")
umap_data <- umap_data %>% 
  mutate(samples = colnames(normal)) %>% 
  left_join(meta, by = "samples")

# Plot UMAP results
for(i in c("donor_name", "gender", "age", "grouped_age", "region", "grouped_region")){
  umap_plot <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = !!sym(i))) +
    geom_point(size = 1.1) +
    theme(
      #legend.position = "none",
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.line = element_line()
    )
  ggsave(paste0("results/exploratory_analysis/umap_", i, ".png"), umap_plot, height = 5, width = 6)
}

###############################################
###     differential expression analysis    ###
###############################################

DEres <- list()

# For each grouped_region
for (gr in setdiff(levels(meta$grouped_region), "Primordial_structures")) {
  local_meta <- meta[meta$grouped_region == gr,]
  local_meta <- mutate_all(local_meta, factor) # VERY IMPORTANT!!!
  local_normal <- normal[,meta$grouped_region == gr]
  res_df <- normal[,c(), drop = FALSE]
  # pvalue and log2fc
  for (sb1 in levels(local_meta$grouped_age)) {
    for (sb2 in setdiff(levels(local_meta$grouped_age), sb1)) {
      # pvalue
      res_df[[paste0("pval_", sb1, "_vs_", sb2)]] <- apply(local_normal, 1, function(row){
        wilcox.test(row[local_meta$grouped_age == sb1], row[local_meta$grouped_age == sb2])$p.value
      })
      # padj
      res_df[[paste0("padj_", sb1, "_vs_", sb2)]] <- p.adjust(res_df[[paste0("pval_", sb1, "_vs_", sb2)]], method = "BH")
      # log2fc
      res_df[[paste0("log2fc_", sb1, "_vs_", sb2)]] <- apply(local_normal, 1, function(row){
        log2(mean(row[local_meta$grouped_age == sb1])/mean(row[local_meta$grouped_age == sb2]))
      }) # what if one of the means is zero??? (filtered genes)
    }
  }
  for (sb in levels(local_meta$grouped_age)) {
    # meta p
    res_df[[paste0("meta_pval_", sb)]] <- apply(res_df[,grepl(paste0("padj_", sb), colnames(res_df))], 1, function(row){sumlog(row)$p})
    
    # median log2fc
    res_df[[paste0("med_log2fc_", sb)]] <- apply(res_df[,grepl(paste0("log2fc_", sb), colnames(res_df))], 1, function(row){median(row)})
  }
  
  DEres[[gr]] <- res_df
}

# For each grouped_age
for (gr in levels(meta$grouped_age)) {
  local_meta <- meta[meta$grouped_age == gr,]
  local_meta <- mutate_all(local_meta, factor) # VERY IMPORTANT!!!
  local_normal <- normal[,meta$grouped_age == gr]
  res_df <- normal[,c(), drop = FALSE]
  # pvalue and log2fc
  for (sb1 in levels(local_meta$grouped_region)) {
    for (sb2 in setdiff(levels(local_meta$grouped_region), sb1)) {
      # pvalue
      res_df[[paste0("pval_", sb1, "_vs_", sb2)]] <- apply(local_normal, 1, function(row){
        wilcox.test(row[local_meta$grouped_region == sb1], row[local_meta$grouped_region == sb2])$p.value
      })
      # padj
      res_df[[paste0("padj_", sb1, "_vs_", sb2)]] <- p.adjust(res_df[[paste0("pval_", sb1, "_vs_", sb2)]], method = "BH")
      # log2fc
      res_df[[paste0("log2fc_", sb1, "_vs_", sb2)]] <- apply(local_normal, 1, function(row){
        log2(mean(row[local_meta$grouped_region == sb1])/mean(row[local_meta$grouped_region == sb2]))
      }) # what if one of the means is zero???
    }
  }
  for (sb in levels(local_meta$grouped_region)) {
    # meta p
    res_df[[paste0("meta_pval_", sb)]] <- apply(res_df[,grepl(paste0("padj_", sb), colnames(res_df))], 1, function(row){sumlog(row)$p})
    
    # median log2fc
    res_df[[paste0("med_log2fc_", sb)]] <- apply(res_df[,grepl(paste0("log2fc_", sb), colnames(res_df))], 1, function(row){median(row)})
  }
  
  DEres[[gr]] <- res_df
}


######################
###     boxplot    ###
######################

features_to_plot <- c("CACNA1H","RELN","FYN","CACNG4","CREBBP","DAB1","NCK2","PDK1","SRC")
dodge_lateral <- 1


##### comparing grouped_region

for(i in levels(meta$grouped_age)){
  # subseting for genes of interest
  local_normal <- normal %>% 
    filter(rownames(.) %in% features_to_plot) %>% 
    dplyr::select(meta[meta$grouped_age == i,"samples"])
  local_meta <- meta %>% 
    filter(grouped_age == i) %>% 
    mutate_all(factor)
  
  # annotation (metap and med_log2fc)
  annots <- DEres[[i]][rownames(DEres[[i]]) %in% rownames(local_normal),grepl("meta|med", colnames(DEres[[i]]))]
  for (j in levels(local_meta$grouped_region)) {
    annots[[j]] <- paste0("p=", format(annots[[paste0("meta_pval_",j)]], scientific = TRUE, digits = 1), "\nfc=", round(annots[[paste0("med_log2fc_",j)]], 1))
  }
  annots <- annots[,colnames(annots) %in% levels(local_meta$grouped_region)]
  annots$gene <- rownames(annots)
  annots <- annots %>% pivot_longer(cols = -gene, names_to = "group", values_to = "annot")
  annots <- as.data.frame(annots)
  
  # adding a condition row
  local_normal["grouped_region",] <- meta[meta$grouped_age == i, "grouped_region"]
  
  # transposing
  local_normal <- as.data.frame(t(local_normal))
  
  # long format
  long_int <- local_normal %>% pivot_longer(cols = -grouped_region, names_to = "gene", values_to = "expression")
  long_int <- as.data.frame(long_int)
  long_int$expression <- as.numeric(long_int$expression)
  
  # re-level grouped_region to improve order
  long_int$grouped_region <- factor(long_int$grouped_region, levels = c("Frontal_lobe", "Temporal_lobe", "Parietal_lobe", "Occipital_lobe", "Subcortical_structures", "Cerebellum", "Primordial_structures"))
  
  # adding annotation label
  long_int <- long_int %>%
    group_by(grouped_region, gene) %>%
    mutate(annot = replace(expression, expression != max(expression), NA))
  long_int <- as.data.frame(long_int)
  long_int$annot[!isNA(long_int$annot)] <- annots[match(paste(long_int[!isNA(long_int$annot),"grouped_region"], long_int[!isNA(long_int$annot),"gene"]), paste(annots$group,annots$gene)),"annot"]
  
  # Plotting
  box_plot <- ggplot(long_int, aes(x = gene, y = expression, fill = grouped_region, label = annot)) +
    geom_point( #jitters scatters randomly everytime
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = dodge_lateral), 
      color = "#550000", 
      size = 1, 
      alpha = 0.5
    )+
    geom_boxplot(
      position = position_dodge(width = dodge_lateral), # spacing between boxes
      alpha = 0.5, # because points are behind decrease opacity
      #color = rep(as.list(brewer.pal(7, "Dark2")), 12)
    )+
    geom_text(
      aes(label = annot, y = expression*1.3), #y = expression + 5 # y = max(expression) + 5
      position = position_dodge(width = dodge_lateral),
      color = "black",
      size = 2.8
    )+
    facet_wrap(~gene, scales = "free")+ # separated box plots
    labs(
      x = NULL,
      y = "RPKM normalized expression"
    )+ 
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid = element_line(color = "gray", linetype = "dotted"),
      strip.background = element_rect(fill= brewer.pal(7, "Set3")[6]),
      axis.text.x = element_blank(), #removes bottom gene name
      panel.grid.major.x = element_blank(), # removes vertical lines
      panel.grid.minor.x = element_blank(),
      strip.text = element_text(size = 12, face = "bold")  # Increase facet title size and make it bold
    )
  
  ggsave(paste0("results/box_", i, ".png"), plot = box_plot, height = 8, width = 13)
}


##### comparing grouped_age

# setdiff é o complementar
for(i in setdiff(levels(meta$grouped_region), "Primordial_structures")){
  # subseting for genes of interest
  #local_normal <- normal[rownames(normal) %in% intGenes, meta$grouped_region == i]
  local_normal <- normal[rownames(normal) %in% features_to_plot, meta$grouped_region == i]
  local_meta <- meta[meta$grouped_region == i,]
  local_meta <- mutate_all(local_meta, factor)
  
  # annotation (metap and med_log2fc)
  annots <- DEres[[i]][rownames(DEres[[i]]) %in% rownames(local_normal),grepl("meta|med", colnames(DEres[[i]]))]
  for (j in levels(local_meta$grouped_age)) {
    annots[[j]] <- paste0("p=", format(annots[[paste0("meta_pval_",j)]], scientific = TRUE, digits = 1), "\nfc=", round(annots[[paste0("med_log2fc_",j)]], 1))
  }
  annots <- annots[,colnames(annots) %in% levels(local_meta$grouped_age)]
  annots$gene <- rownames(annots)
  annots <- annots %>% pivot_longer(cols = -gene, names_to = "group", values_to = "annot")
  annots <- as.data.frame(annots)
  
  # adding a condition row
  local_normal["grouped_age",] <- meta[meta$grouped_region == i, "grouped_age"]
  
  # transposing
  local_normal <- as.data.frame(t(local_normal))
  
  # long format
  long_int <- local_normal %>% pivot_longer(cols = -grouped_age, names_to = "gene", values_to = "expression")
  long_int <- as.data.frame(long_int)
  long_int$expression <- as.numeric(long_int$expression)
  
  # re-leveling grouped_age to order correctly
  long_int$grouped_age <- factor(long_int$grouped_age, levels = c("8_12PCW", "13_19PCW", "20_40PCW", "0_2YEARS", "3_10YEARS", "11_22YEARS", "23_40YEARS"))
  
  # adding annotation label
  long_int <- long_int %>%
    group_by(grouped_age, gene) %>%
    mutate(annot = replace(expression, expression != max(expression), NA))
  long_int <- as.data.frame(long_int)
  long_int$annot[!isNA(long_int$annot)] <- annots[match(paste(long_int[!isNA(long_int$annot),"grouped_age"], long_int[!isNA(long_int$annot),"gene"]), paste(annots$group,annots$gene)),"annot"]
  
  # Plotting
  box_plot <- ggplot(long_int, aes(x = gene, y = expression, fill = grouped_age, label = annot)) +
    geom_point( #jitters scatters randomly everytime
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = dodge_lateral), 
      color = "#550000", 
      size = 1, 
      alpha = 0.5
    )+
    geom_boxplot(
      position = position_dodge(width = dodge_lateral),
      alpha = 0.5, # because points are behind decrease opacity
      #fill = rep(as.list(brewer.pal(7, "Dark2")), 12)
    )+ # position = position_dodge(width = 0.9)
    geom_text( # To show all labels
      aes(label = annot, y = expression + expression/3), #y = expression + 5
      position = position_dodge(width = dodge_lateral),
      #max.overlaps = Inf, 
      color = "black", 
      size = 2.8,
      #angle = 45
    )+
    facet_wrap(~gene, scales = "free")+ # separated box plots
    labs(
      x = NULL,
      y = "RPKM normalized expression"
    )+ 
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid = element_line(color = "gray", linetype = "dotted"),
      strip.background = element_rect(fill= brewer.pal(7, "Set3")[7]),
      axis.text.x = element_blank(), #removes bottom gene name
      panel.grid.major.x = element_blank(), # removes vertical lines
      panel.grid.minor.x = element_blank(),
      strip.text = element_text(size = 12, face = "bold")  # Increase facet title size and make it bold
    )
  
  ggsave(paste0("results/box_", i, ".png"), plot = box_plot, height = 8, width = 13)
}


######################
###     CEMiTool   ###
######################

# Running WGCNA
cem <- cemitool(
  expr = log2(normal %>% filter(rownames(.) %in% hvf_5000) + 1),
  #network_type = "signed", # without the absolute value of the correlation
  #tom_type = "signed",
  filter = FALSE, # all genes
  plot = FALSE, 
  verbose = FALSE
)

aux <- data.frame(sum = rowSums(normal %>% filter(rownames(.) %in% hvf_5000)))


############################
###     modules sizes    ###
############################

# modules
mod <- cem@module %>% 
  filter(modules != "Not.Correlated") %>% 
  mutate(
    modules = gsub("M", "", modules),
    modules = as.numeric(modules),
    modules = as.factor(modules) # do not remove!
  ) %>% 
  left_join(sfari %>% select(gene.symbol, gene.score, syndromic, number.of.reports), by = c("genes" = "gene.symbol")) %>% # adding sfari
  rename("score" = "gene.score", "reports" = "number.of.reports") %>% 
  mutate(
    score = case_when(is.na(score) ~ "NO", TRUE ~ as.character(score)),
    score = factor(score, levels = c("1", "2", "3", "NO")),
    reports = case_when(is.na(reports) ~ 0, TRUE ~ reports),
    syndromic = case_when(is.na(syndromic) ~ "NO", TRUE ~ as.character(syndromic)),
    syndromic = factor(syndromic, levels = c("1", "0", "NO"))
  )

# Amount of genes per type per module
long_mod <- mod %>%
  group_by(modules) %>%
  summarise(num_genes = n()) %>%
  ungroup()

source("C:/Users/ppoli/Desktop/bioinfo/helder/degchromodensity/chromopack.R")
# Plot the bar plot
bar_plot <- barfunc(
  long_mod, "modules", "num_genes", annot_loc = "top", 
  fill_col = "#5577bbff", size_yaxis_text = 0,
  x_title = "Modules", y_title = "Number of Features")
ggsave("results/modules_sizes.png", plot = bar_plot, height = 3, width = 9)



###############################
###     SFARI enrichment    ###
############################### 
### syndromic genes can have score 1,2 or 3

annot_df <- data.frame(modules = levels(mod$modules))
for (st_name in c("SFARI","1", "2","3","Syndromic","Non_syndromic")) {
  if(st_name %in% c("Syndromic","Non_syndromic")){
    cl <- "syndromic"
    if(st_name == "Syndromic"){st <- c("1")}else{st <- c("0")}
  }else{
    cl <- "score"
    if(st_name == "SFARI"){st <- c("1","2","3")}else{st <- c(st_name)}
  }
  annot_df[[st_name]] <- apply(annot_df, 1, function(row){
    M <- row[["modules"]]
    a <- mod %>% filter(modules == M, !!sym(cl) %in% st) %>% nrow # set in module
    b <- mod %>% filter(!!sym(cl) %in% st) %>% nrow() # total set genes
    c <- mod %>% filter(modules == M, !(!!sym(cl) %in% st)) %>% nrow() # not set in module
    d <- mod %>% filter(!(!!sym(cl) %in% st)) %>% nrow() # total non-set genes
    
    # Hypergeometric test parameters
    k <- a # successes in sample (sfari in module)
    m <- b # successes in population (total sfari)
    n <- c + a # total in sample (total in module)
    N <- b + d # total population (total genes)
    
    pval <- phyper(k - 1, m, N - m, n, lower.tail = FALSE)
    return(pval)
  })
}

annot_df <- annot_df %>%
  column_to_rownames("modules") %>% 
  t() %>% 
  as.data.frame()

for (i in colnames(annot_df)) {
  annot_df[[paste0(i, "_padj")]] <- p.adjust(annot_df[[i]], method = "BH")
}

dot <- annot_df %>% 
  t() %>% 
  as.data.frame() %>% 
  filter(!grepl("_padj", rownames(.))) %>% 
  rownames_to_column("modules") %>% 
  pivot_longer(cols = -modules, names_to = "ref", values_to = "pval") %>% 
  mutate(
    #modules = gsub("_padj", "", modules),
    modules = as.numeric(modules),
    modules = as.factor(modules),
    ref = case_when(
      ref == "SFARI" ~ "All SFARI",
      ref == "1" ~ "Score 1",
      ref == "2" ~ "Score 2",
      ref == "3" ~ "Score 3",
      ref == "Syndromic" ~ "Syndromic",
      ref == "Non_syndromic" ~ "Non Syndromic",
    ),
    ref = factor(ref, levels = c("All SFARI", "Score 1", "Score 2", "Score 3", "Syndromic", "Non Syndromic")),
    col = case_when(pval < 0.05 ~ "p<0.05", TRUE ~ "p>=0.05")
  ) %>% 
  filter(pval != 1)


# Plot
dot_plot <- ggplot(dot, aes(x = ref, y = modules, col = col, size = -log10(pval))) +
  geom_point() +
  scale_color_manual(values = c("p>=0.05" = "#444444", "p<0.05" = "#dd4444")) +
  theme_minimal() +
  labs(
    x = NULL,
    y = NULL,
  )+
  theme(
    plot.background = element_rect(fill = "white", color = NA),  # Set background color to white
    panel.background = element_rect(fill = "white", color = NA),  # Set panel background color to white
    panel.border = element_blank(),  # Remove panel borders
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black", linewidth = 0.5),  # Make axes lines bold
    axis.title.x = element_text(face = "bold", size = 15),  # Make axis titles bold
    axis.title.y = element_text(face = "bold", size = 15),  # Make axis titles bold
    axis.text.y = element_text(face = "bold", size = 15),  # Make axis text bold
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1, size = 13), # rotate x axis test 45 degrees
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),  # Make y-axis ticks bold
    legend.text = element_text(size = 11),# Change legend text font size
    legend.title = element_text(size = 13) # Change legend title font size
    #legend.position = "none" # No legend
  ) +
  guides(col = guide_legend(title = NULL)) # remove col legend title
ggsave("results/sfari_on_clusters.png", plot = dot_plot, height = 6, width = 4)




#################### Sfari score bar plot

# Amount of genes per type per module
long_mod <- mod %>%
  group_by(modules, score) %>%
  summarise(num_genes = n()) %>%
  group_by(modules) %>%
  mutate(pct = num_genes / sum(num_genes) * 100) %>% ########### add column just for annotation   ######## título de legenda
  ungroup()

# absolute value bar plot
bar_plot <- barfunc(
  long_mod, "modules", "num_genes", x_title = "Modules", y_title = "Number of Features", fill_groups = "score", 
  annot_loc = "repel", fill_col = list("1" = "#88111188", "2" = "#88881188", "3" = "#11881188", "NO" = "#11118888")
)
ggsave("results/SFARI_abs_score.png", plot = bar_plot, height = 3, width = 9)


# percentage value bar plot
bar_plot <- barfunc(
  long_mod, "modules", "pct", x_title = "Modules", y_title = "Percentage of Features", fill_groups = "score", 
  annot_loc = "F", fill_col = list("1" = "#88111188", "2" = "#88881188", "3" = "#11881188", "NO" = "#11118888")
)
ggsave("results/SFARI_pct_score.png", plot = bar_plot, height = 3, width = 9)

#################### Syndromic bar plot

long_mod <- mod %>%
  group_by(modules, syndromic) %>%
  summarise(num_genes = n()) %>%
  group_by(modules) %>%
  mutate(pct = num_genes / sum(num_genes) * 100) %>% 
  ungroup()

# absolute value bar plot
bar_plot <- barfunc(
  long_mod, "modules", "num_genes", x_title = "Modules", y_title = "Number of Features", fill_groups = "syndromic", 
  annot_loc = "repel", fill_col = list("1" = "#88111188", "0" = "#88881188", "NO" = "#11118888")
)
ggsave("results/SFARI_abs_syndromic.png", plot = bar_plot, height = 3, width = 9)


# percentage value bar plot
bar_plot <- barfunc(
  long_mod, "modules", "pct", x_title = "Modules", y_title = "Percentage of Features", fill_groups = "syndromic", 
  annot_loc = "F", fill_col = list("1" = "#88111188", "0" = "#88881188", "NO" = "#11118888")
)
ggsave("results/SFARI_pct_syndromic.png", plot = bar_plot, height = 3, width = 9)



#######################
###   network plot  ###
#######################

modules_colors <- c( # max 49 colors/modules
  brewer.pal(brewer.pal.info["Set1", "maxcolors"], "Set1"),
  brewer.pal(brewer.pal.info["Dark2", "maxcolors"], "Dark2"),
  brewer.pal(brewer.pal.info["Paired", "maxcolors"], "Paired"),
  brewer.pal(brewer.pal.info["Accent", "maxcolors"], "Accent"),
  brewer.pal(brewer.pal.info["Spectral", "maxcolors"], "Spectral")
)
# subseting
modules_colors <- modules_colors[1:nlevels(mod$modules)]
names(modules_colors) <- levels(mod$modules)


########################## RELN
filt <- rownames(cem@adjacency) %in% (cem@adjacency %>% as.data.frame() %>% dplyr::select("RELN") %>% arrange(-RELN) %>% slice(1:25) %>% rownames(.))
adj <- cem@adjacency[filt, filt]
diag(adj) <- 0 # remove self-loops
adj[adj < 0.027] <- 0 # threshold for better visualization

moduleGraph <- graph_from_adjacency_matrix(as.matrix(adj), mode = "undirected", weighted = TRUE, diag = FALSE) # graph object
gene_to_module <- setNames(mod$modules, mod$genes) # First, create a named vector that maps gene names to their respective modules
V(moduleGraph)$module <- gene_to_module[V(moduleGraph)$name] # Then, assign the module information to the 'module' attribute of the graph's vertices

set.seed(42)
net_plot <- ggraph(moduleGraph, layout = 'fr') + 
  geom_edge_link(aes(width = weight), alpha = 0.2, color = "#999999") +
  geom_node_point(aes(color = module), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, color = "black", size = 3, fontface = "bold") +
  theme_graph() +
  scale_color_manual(values = modules_colors) +
  theme(
    #legend.position = "none"
  )
ggsave("results/network_RELN.png", plot = net_plot, width = 5, height = 5)


########################## CACNA1H
filt <- rownames(cem@adjacency) %in% (cem@adjacency %>% as.data.frame() %>% dplyr::select("CACNA1H") %>% arrange(-CACNA1H) %>% slice(1:25) %>% rownames(.))
adj <- cem@adjacency[filt, filt]
diag(adj) <- 0 # remove self-loops
adj[adj < 0.15] <- 0 # threshold for better visualization

moduleGraph <- graph_from_adjacency_matrix(as.matrix(adj), mode = "undirected", weighted = TRUE, diag = FALSE) # graph object
gene_to_module <- setNames(mod$modules, mod$genes) # First, create a named vector that maps gene names to their respective modules
V(moduleGraph)$module <- gene_to_module[V(moduleGraph)$name] # Then, assign the module information to the 'module' attribute of the graph's vertices

set.seed(42)
net_plot <- ggraph(moduleGraph, layout = 'fr') + 
  geom_edge_link(aes(width = weight), alpha = 0.2, color = "#999999") +
  geom_node_point(aes(color = module), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, color = "black", size = 3, fontface = "bold") +
  theme_graph() +
  scale_color_manual(values = modules_colors) +
  theme(
    #legend.position = "none"
  )
ggsave("results/network_CACNA1H.png", plot = net_plot, width = 5, height = 5)



########################## Module 5
filt <- rownames(cem@adjacency) %in% (mod %>% filter(modules == 5) %>% pull(genes))
adj <- cem@adjacency[filt, filt]
diag(adj) <- 0 # remove self-loops
TOM <- TOMsimilarity(adj)
rownames(TOM) <- rownames(adj)
colnames(TOM) <- colnames(adj)


####################
###   tile plot  ###
####################

# Hub genes
hubs <- get_hubs(cem, 10) #hubs per module

######### M5
tile_data <- mod[mod$modules == 5,"genes", drop = F] %>%
  mutate(
    SFARI = as.factor(case_when(genes %in% sfari_genes ~ 1, TRUE ~ 0)),
    Cohorts = as.factor(case_when(genes %in% patients_genes ~ 1, TRUE ~ 0)),
    IntGenes = as.factor(case_when(genes %in% intGenes ~ 1, TRUE ~ 0)),
    HubGene = as.factor(case_when(genes %in% names(hubs[["M5"]]) ~ 1, TRUE ~ 0))
  ) %>%
  #column_to_rownames("genes") %>% 
  filter(SFARI==1|Cohorts==1|IntGenes==1|HubGene==1) %>% 
  pivot_longer(cols = -genes, names_to = "Category", values_to = "Presence")

tile_plot <- ggplot(tile_data, aes(x = Category, y = genes, fill = Presence)) + 
  geom_tile() +
  scale_fill_manual(values = c("1" = "#bb0000", "0" = "white")) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),  # Set background color to white
    panel.background = element_rect(fill = "white", color = NA),  # Set panel background color to white
    panel.border = element_blank(),  # Remove panel borders
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black", linewidth = 0.5),  # Make axes lines bold
    axis.title.x = element_blank(),  # Make axis titles bold
    axis.title.y = element_blank(),  # Make axis titles bold
    axis.text.y = element_text(face = "bold", size = 10),  # Make axis text bold
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1, size = 10), # rotate x axis test 45 degrees
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),  # Make y-axis ticks bold
    legend.position = "none" # No legend
  )
ggsave("results/tile_M5.png", plot = tile_plot, height = 10, width = 2)

########## M2

tile_data <- mod[mod$modules == 2,"genes", drop = F] %>%
  mutate(
    SFARI = as.factor(case_when(genes %in% sfari_genes ~ 1, TRUE ~ 0)),
    Cohorts = as.factor(case_when(genes %in% patients_genes ~ 1, TRUE ~ 0)),
    IntGenes = as.factor(case_when(genes %in% intGenes ~ 1, TRUE ~ 0)),
    HubGene = as.factor(case_when(genes %in% names(hubs[["M2"]]) ~ 1, TRUE ~ 0))
  ) %>%
  #column_to_rownames("genes") %>% 
  filter(SFARI==1|Cohorts==1|IntGenes==1|HubGene==1) %>% 
  pivot_longer(cols = -genes, names_to = "Category", values_to = "Presence")

tile_plot <- ggplot(tile_data, aes(x = Category, y = genes, fill = Presence)) + 
  geom_tile() +
  scale_fill_manual(values = c("1" = "#bb0000", "0" = "white")) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),  # Set background color to white
    panel.background = element_rect(fill = "white", color = NA),  # Set panel background color to white
    panel.border = element_blank(),  # Remove panel borders
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black", linewidth = 0.5),  # Make axes lines bold
    axis.title.x = element_blank(),  # Make axis titles bold
    axis.title.y = element_blank(),  # Make axis titles bold
    axis.text.y = element_text(face = "bold", size = 10),  # Make axis text bold
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1, size = 10), # rotate x axis test 45 degrees
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),  # Make y-axis ticks bold
    legend.position = "none" # No legend
  )
ggsave("results/tile_M2.png", plot = tile_plot, height = 18, width = 2)



#######################################
###   over representation analysis  ###
#######################################

##########
### M2 ###
##########

# Finding GO enriched sets
overM2 <- enrichGO(
  keyType = "SYMBOL", #there is no preference for using SYMBOL, ENSEMBL or ENTREZ!
  gene = mod %>% filter(modules == 2) %>% pull("genes"),
  universe = rownames(normal), # sample genes. IMPORTANT!!!
  OrgDb = org.Hs.eg.db,
  ont = "all",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

###########
### M5 ###
###########

# Finding GO enriched sets
overM5 <- enrichGO(
  keyType = "SYMBOL", #there is no preference for using SYMBOL, ENSEMBL or ENTREZ!
  gene = mod %>% filter(modules == 5) %>% pull("genes"),
  universe = rownames(normal), # sample genes. IMPORTANT!!!
  OrgDb = org.Hs.eg.db,
  ont = "all",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)





########## Plotting

aux <- rbind(
    overM2@result %>% mutate(Module = "2"),
    overM5@result %>% mutate(Module = "5")
  ) %>% 
  filter(
    ID %in% c('GO:0010975','GO:0030900','GO:0006338','GO:0061564','GO:0007411','GO:0001764','GO:0098978','GO:0000149','GO:0097060','GO:0098984',
            'GO:0022843','GO:1990351','GO:0098982')
  ) %>% 
  mutate(
    Description = case_when(
      duplicated(Description) ~ paste0(Description, " "), # gambiarrinha
      TRUE ~ Description
    )
  )

enrich_plot <- funcenrichplot(aux, highlight = "Module", color_list = modules_colors, x_text_size = 11)
ggsave("results/GO_modules.png", plot = enrich_plot, height = 4, width = 6)

##########
### M2 ###
##########

# Bar plot data
aux_bar <- overM2@result[,c("Description","p.adjust","ONTOLOGY")]
aux_bar$score <- -log10(aux_bar$p.adjust)
aux_bar <- aux_bar[order(aux_bar$score),]
aux_bar$Description <- factor(aux_bar$Description, levels = aux_bar$Description)
aux_bar$ONTOLOGY <- factor(aux_bar$ONTOLOGY)

# Bar plot
bar_plot <- ggplot(aux_bar) +
  geom_col( # fill = ONTOLOGY ######### colorir de acordo com BP, CC, MF
    aes(x = score, y = Description, fill = ONTOLOGY), 
    #fill = "#4422aa77", 
    width = 0.6
  ) +
  labs(
    title = NULL,
    y = NULL,
    x = expression("-log"[10]*"(padjust)")
  ) +
  scale_x_continuous(
    limits = c(0, max(aux_bar$score)),
    breaks = seq(0, max(aux_bar$score), by = 1), 
    expand = c(0, 0), # The horizontal axis does not extend to either side
    position = "top"  # Labels are located on the top
  ) +
  # The vertical axis only extends upwards 
  #scale_y_discrete(expand = expansion(add = c(0, 0.5))) +
  theme(
    # Set background color to white
    panel.background = element_rect(fill = "white"),
    # Set the color and the width of the grid lines for the horizontal axis
    panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.3),
    # Remove tick marks by setting their length to 0
    axis.ticks.length = unit(0, "mm"),
    # Remove the title for both axes
    #axis.title = element_blank(),
    # Only left line of the vertical axis is painted in black
    axis.line.y.left = element_line(color = "black"),
    # Remove labels from the vertical axis
    axis.text.y = element_blank(),
    # But customize labels for the horizontal axis
    axis.text.x = element_text(family = "Econ Sans Cnd", size = 16)
  ) +
  geom_vline(
    xintercept = -log10(0.05),
    linetype = "dotted", 
    color = "#444444"
  ) + 
  scale_fill_manual(
    values = c("BP" = "#4422aa77", "CC" = "#22aa4477", "MF" = "#aa224477")
  ) +
  geom_text(
    aes(0, y = Description, label = Description),
    hjust = 0,
    nudge_x = 0.05,
    colour = "black",
    family = "Econ Sans Cnd",
    size = 2
  )

# Saving
ggsave("results/M2bar.png", plot = bar_plot, height = 4, width = 7)

###########
### M5 ###
###########

# Bar plot data
aux_bar <- overM10@result[,c("Description","p.adjust","ONTOLOGY")]
aux_bar$score <- -log10(aux_bar$p.adjust)
aux_bar <- aux_bar[order(aux_bar$score),]
aux_bar$Description <- factor(aux_bar$Description, levels = aux_bar$Description)
aux_bar$ONTOLOGY <- factor(aux_bar$ONTOLOGY)

# Bar plot
bar_plot <- ggplot(aux_bar) +
  geom_col( # fill = ONTOLOGY ######### colorir de acordo com BP, CC, MF
    aes(x = score, y = Description, fill = ONTOLOGY), 
    #fill = "#4422aa77", 
    width = 0.6
  ) +
  labs(
    title = NULL,
    y = NULL,
    x = expression("-log"[10]*"(padjust)")
  ) +
  scale_x_continuous(
    limits = c(0, max(aux_bar$score)),
    breaks = seq(0, max(aux_bar$score), by = 1), 
    expand = c(0, 0), # The horizontal axis does not extend to either side
    position = "top"  # Labels are located on the top
  ) +
  # The vertical axis only extends upwards 
  #scale_y_discrete(expand = expansion(add = c(0, 0.5))) +
  theme(
    # Set background color to white
    panel.background = element_rect(fill = "white"),
    # Set the color and the width of the grid lines for the horizontal axis
    panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.3),
    # Remove tick marks by setting their length to 0
    axis.ticks.length = unit(0, "mm"),
    # Remove the title for both axes
    #axis.title = element_blank(),
    # Only left line of the vertical axis is painted in black
    axis.line.y.left = element_line(color = "black"),
    # Remove labels from the vertical axis
    axis.text.y = element_blank(),
    # But customize labels for the horizontal axis
    axis.text.x = element_text(family = "Econ Sans Cnd", size = 16)
  ) +
  geom_vline(
    xintercept = -log10(0.05),
    linetype = "dotted", 
    color = "#444444"
  ) + 
  scale_fill_manual(
    values = c("BP" = "#4422aa77", "CC" = "#22aa4477", "MF" = "#aa224477")
  ) +
  geom_text(
    aes(0, y = Description, label = Description),
    hjust = 0,
    nudge_x = 0.05,
    colour = "black",
    family = "Econ Sans Cnd",
    size = 2
  )

# Network plots
net_plot <- cnetplot(overM10) + 
  #labs(title = "GO - Biological Process") +
  theme(plot.title = element_text(hjust = 0.5))

# Circle plots (careful with order)
circ_plot <- cnetplot(overM10, circular = TRUE, colorEdge = TRUE) + 
  #labs(title = "GO - Biological Process") +
  theme(plot.title = element_text(hjust = 0.5))

# Saving
ggsave("results/correlation/functional/M10bar.png", plot = bar_plot, height = 6, width = 7)



