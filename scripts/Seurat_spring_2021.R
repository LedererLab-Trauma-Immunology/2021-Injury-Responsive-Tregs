library(dplyr)
library(Seurat)
library(patchwork)
library(progress)
library(ggplot2)
library(stringr)

treg.data <- Read10X(data.dir = "C:/Users/bh719/Dropbox (Partners HealthCare)/Sally and Jim/FeiData of scRNAseq  Treg project/10xProcessedData of day 7 after burn/GenomeSeq/outs/filtered_feature_bc_matrix")
treg <- CreateSeuratObject(counts = treg.data,project = 'treg7day',min.cells = 3, min.features = 200)

#Pre Process
treg[["percent.mt"]] <- PercentageFeatureSet(treg, pattern = "^mt-") # for some stupiddd reason, I had to change "^MT-" to "^mt-"

#VlnPlot(treg, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(treg, features = c("percent.mt")) + geom_hline(yintercept = 6, linetype="dashed", color = "red")
VlnPlot(treg, features = c("nFeature_RNA")) + geom_hline(yintercept = 201, linetype="dashed", color = "red") + geom_hline(yintercept = 3000, linetype="dashed", color = "red")
abline(h = 5)
#Unique features over 4000, under 500.  

treg <- subset(treg, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 6) #values that make more sense to me
#treg <- subset(treg, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)  #default in Seurat tutorial

treg <- NormalizeData(treg)

#Identify variable features
treg <- FindVariableFeatures(treg, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
#top25 <- head(VariableFeatures(treg), 25)

# plot variable features with and without labels
#plot1 <- VariableFeaturePlot(treg)
#plot1 <- LabelPoints(plot = plot1, points = top25, repel = TRUE)
#plot1

#Scale Data
all.genes <- rownames(treg)
treg <- ScaleData(treg, features = all.genes)

#PCA
treg <- RunPCA(treg, features = VariableFeatures(object = treg))

#print(treg[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(treg, dims = 1:2, reduction = "pca")

#DimPlot(treg, reduction = "pca")

#DimHeatmap(treg, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(treg, dims = 1:15, cells = 500, balanced = TRUE)

#Determine 'dimensionality' of the dataset
#treg <- JackStraw(treg, num.replicate = 100)
#treg <- ScoreJackStraw(treg, dims = 1:20)

#JackStrawPlot(treg, dims = 1:25)

#ElbowPlot(treg)

#Cluster
treg <- FindNeighbors(treg, dims = 1:20, nn.method = 'rann')
treg <- FindClusters(treg, resolution = .5)
#treg <- RunUMAP(treg, dims = 1:20)
treg <- RunTSNE(treg,dims = 1:20)
DimPlot(treg, reduction = "tsne", label = TRUE)

#Show that cluster 6 is B cells
cluster6.markers <- FindMarkers(treg,ident.1 = 6, min.pct = 0.25)
#cluster6.markers <- cluster6.markers[order(abs(cluster6.markers$avg_logFC),decreasing = TRUE),]
cluster6.markers <- cluster6.markers[order(cluster6.markers$avg_log2FC,decreasing = TRUE),]
head(cluster6.markers, n = 25)
# Cd74, Cd79, Ms4a1, Cd19 - B cell markers

#Remove B cells and redo clustering
treg <- subset(treg, subset = seurat_clusters != 6)
treg <- FindVariableFeatures(treg, selection.method = "vst", nfeatures = 2000)

#Metascape
#top2000 <- head(VariableFeatures(treg), 2000)
#sink("sc2000.txt")
#writeLines(unlist(lapply(top2000, paste, collapse=" ")))
#sink()
#

#plot1 <- VariableFeaturePlot(treg)
#plot1 <- LabelPoints(plot = plot1, points = top25, repel = TRUE)
#plot1
treg <- RunPCA(treg, features = VariableFeatures(object = treg))
treg <- RunTSNE(treg,dims = 1:10)
DimPlot(treg, reduction = "tsne", label = TRUE)

#pcas <- c(10,11,12,13,14,15)
#for (i in pcas){
#  treg <- RunTSNE(treg, dims = 1:i)
#  DimPlot(treg, reduction = "tsne", label = TRUE) 
#}


#treg <- JackStraw(treg, num.replicate = 100)
#treg <- ScoreJackStraw(treg, dims = 1:20)
#JackStrawPlot(treg, dims = 1:15) #kinda worthless
#ElbowPlot(treg, ndims = 30)
treg <- FindNeighbors(treg, dims = 1:25, nn.method = 'rann')
treg <- FindClusters(treg, resolution = 2.0)
DimPlot(treg, reduction = "tsne", label = TRUE)

#VDJ
tcr_folder <- "C:/Users/bh719/Dropbox (Partners HealthCare)/Sally and Jim/FeiData of scRNAseq  Treg project/10xProcessedData of day 7 after burn/12-13-multi/12-13_multi/outs/"
treg <- add_clonotype(tcr_folder,treg)
exp_cells <- WhichCells(treg, expression = Cell.State == 'Expanded')
DimPlot(treg, cells.highlight = exp_cells, reduction = 'tsne',sizes.highlight = 1.5,label.size = 12) + scale_color_manual(labels = c("Non Expanded","Expanded"),values = c("grey","red"))

#Sample data
treg <- f_add_sample(treg)
DimPlot(treg, group.by = "sample")
DimPlot(treg, cells.highlight = exp_cells,split.by = "sample", reduction = 'tsne',sizes.highlight = 1.5,label.size = 12) + scale_color_manual(labels = c("Non Expanded","Expanded"),values = c("grey","red"))
table(FetchData(treg, vars = c('sample'), cells = exp_cells))
uninjured <- WhichCells(treg, expression = sample == "Uninjured")
seven_D <- WhichCells(treg, expression = sample == '7D after Injury')
exp_un <- intersect(exp_cells, uninjured)
exp_seven <- intersect(exp_cells, seven_D)
DimPlot(treg, cells.highlight = exp_un, reduction = 'tsne',sizes.highlight = 1.5,label.size = 12) + scale_color_manual(labels = c("Non Expanded","Expanded"),values = c("grey","red"))
DimPlot(treg, cells.highlight = exp_seven, reduction = 'tsne',sizes.highlight = 1.5,label.size = 12) + scale_color_manual(labels = c("Non Expanded","Expanded"),values = c("grey","red"))

Idents(treg) <- treg$sample
DimPlot(treg, cols = c('#00BFC4','#F8766D'))


#Match expanded cells to clusters
exp_clusters <- FetchData(treg, vars = c('seurat_clusters'), cells = exp_cells)
x <- table(exp_clusters)
t <- table(treg$seurat_clusters)
exp_per <- x / t
exp_per
exp_clust <- names(exp_per[exp_per > 0.15])

#Classify clusters as expanded or non / active or not
treg <- f_add_act_pheno(treg, exp_clust)
Idents(treg) <- treg$act_pheno
DimPlot(treg, reduction = 'tsne',cols = rev(ggplotColours(2)))

#Compare expanded clusters vs expanded cells
Idents(treg) <- treg$Cell.State
exp_markers <- FindMarkers(treg, ident.1 = "Expanded", min.pct = 0.25)

Idents(treg) <- treg$act_pheno
pheno_markers <- FindMarkers(treg, ident.1 = "Expanded Phenotype")
pheno_markers <- FindMarkers(treg, ident.1 = "Expanded Phenotype", test.use = "MAST")
pheno_markers_all <- FindMarkers(treg, ident.1 = "Expanded Phenotype", logfc.threshold = 0, min.pct = 0.01)

#Vplots
vplots <- VlnPlot(treg, features = c("Klrg1",'Tigit','Icos','Itgae','S100a6','Cd44'), split.by = "sample", combine = FALSE, pt.size = 0)
VlnPlot(treg, features = c("Klrg1"), split.by = "sample", combine = FALSE, pt.size = 1)
wrap_plots(plots = vplots, ncol = 3)
#Width: 423, Height: 796
VlnPlot(treg, features = c("Klrg1")) + NoLegend()
VlnPlot(treg, features = c("Tigit")) + NoLegend()
VlnPlot(treg, features = c("Icos")) + NoLegend()
VlnPlot(treg, features = c("Itgae")) + NoLegend()
VlnPlot(treg, features = c("S100a6")) + NoLegend()
VlnPlot(treg, features = c("Cd44")) + NoLegend()

#Volcano Plots
library(EnhancedVolcano)
treg_uninjured <- subset(treg, subset = sample == "Uninjured")
treg_7D <- subset(treg, subset = sample == "7D after Injury")
pheno_markers_uninjured <- FindMarkers(treg_uninjured, ident.1 = "Expanded Phenotype", logfc.threshold = 0, min.pct = 0.01)
pheno_markers_7D <- FindMarkers(treg_7D, ident.1 = "Expanded Phenotype", logfc.threshold = 0, min.pct = 0.01)

EnhancedVolcano(pheno_markers_uninjured, lab = row.names(pheno_markers_uninjured), x = 'avg_logFC', y = 'p_val')


saveRDS(treg, file = "C:/Users/bhanc/Dropbox (Partners HealthCare)/Harvard CyTof/for Brandon/Sally/treg.rds")
load("C:/Users/bhanc/Dropbox (Partners HealthCare)/Harvard CyTof/for Brandon/Sally/treg.rds") # failed

cluster1.markers <- FindMarkers(treg, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 10)

VlnPlot(treg, features = c("Izumo1r"))

#Genes
FeaturePlot(treg, features = c("Klrg1",'Tigit','Icos','Itgae','S100a6','Cd44'), reduction = 'tsne',order = TRUE)

#topGo
library(topGO)
library(mygene)
sig_pheno_markers <- pheno_markers[pheno_markers$avg_logFC > 1,]
sig_genes <- row.names(sig_pheno_markers)
all_genes <- row.names(treg)
gene_list <- factor(as.integer(all_genes %in% sig_genes))
names(gene_list) <- all_genes

res <- queryMany(all_genes,scopes = 'symbol', fields=c('entrezgene','ensembl.gene','go','description'),species = 'mouse')
res <- res[!duplicated(res$query),]
ezgene_to_GO <- list()
for (i in 1:length(res$entrezgene)){
  ezgene_to_GO[[res$entrezgene[[i]]]] <- res$go.BP[[i]]$id
}

for (i in 1:length(gene_list)){
  n = names(gene_list)[i]
  index = match(n,res$query)
  names(gene_list)[i] <- res$entrezgene[index]
}

GOdata <- new("topGOdata", ontology = "BP", allGenes = gene_list, annot = annFUN.gene2GO, gene2GO = ezgene_to_GO)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
alllRes <- GenTable(GOdata, classicFisher = resultFisher, ranksOf = "classicFisher", topNodes = 20)
View(alllRes)

GO_fig <- head(allRes, n = 10)
barplot(-log10(as.numeric(GO_fig$classicFisher)),names.arg = GO_fig$Term, horiz = TRUE)

barplot(rev(-log10(as.numeric(GO_fig$classicFisher))),names.arg = rev(GO_fig$Term), horiz = TRUE) #weird labels
barplot(rev(-log10(as.numeric(GO_fig$classicFisher))), horiz = TRUE)

fig4_GO_list = list()

for (i in 1:length(GO_fig$Term)){
  term <- GO_fig$Term[i]
  term_genes <- c()
  for (j in 1:length(sig_genes)){
    index <- match(sig_genes[j],res$query)
    terms <- res$go.BP[[index]]$term
    if(term %in% terms){
      term_genes <- c(term_genes, sig_genes[j])
    }
  }
  fig4_GO_list[[term]] <- term_genes
}

# org.Mm

GOdata <- new("topGOdata", ontology = "BP", allGenes = gene_list, annotationFun = annFUN.org, mapping = "org.Mm.eg")

#Metascape
sink("sig.txt")
writeLines(unlist(lapply(sig_genes, paste, collapse=" ")))
sink()

sink("background.txt")
writeLines(unlist(lapply(all_genes, paste, collapse=" ")))
sink()







# VDJ 

tcr_folder <- "C:/Users/bh719/Dropbox (Partners HealthCare)/Sally and Jim/FeiData of scRNAseq  Treg project/10xProcessedData of day 7 after burn/12-13-multi/12-13_multi/outs/"
treg <- add_clonotype(tcr_folder,treg)
treg <- f_add_tsne(treg)

clone_counts <- f_clone_counts(treg)

exp_tregs <- subset(treg,subset = clonotype_id == "clonotype1")

clonotypeN <- WhichCells(treg, expression = clonotype_id == "clonotype1")

exp_cells_2 <- WhichCells(treg, expression = clonotype_id %in% sprintf("clonotype%s",seq(1:266)))
exp_cells_3 <- WhichCells(treg, expression = clonotype_id %in% sprintf("clonotype%s",seq(1:88)))
exp_cells <- WhichCells(treg, expression = Cell.State == 'Expanded')

DimPlot(treg, cells.highlight = exp_cells, reduction = 'tsne')
DimPlot(treg, cells.highlight = exp_cells, reduction = 'umap',sizes.highlight = 1.5,label.size = 12) + scale_color_manual(labels = c("Non Expanded","Expanded"),values = c("grey","red"))
DimPlot(treg, reduction = "tsne",label = TRUE,label.size = 6,repel = TRUE) + scale_color_manual(values = c("#FF9999","#9B9EFF","#9FF6FF","#BD6E6E","#874E4E","#CCF2BC","#FDCF52","#A969F2","#F29B69","#85EBEB","gray"))

Idents(treg) <- treg$tsne_fig

tsne_fig_list <- list()
for (i in 1:length(levels(treg))){
  level <- levels(treg)[order(levels(treg))][i]
  tsne_fig_list[[level]] <- WhichCells(treg, idents = level)
}

tsne_fig_sizes <- c(0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,2.25)
tsne_fig_colors <- c("#DC143C","#00FFFF","pale green","dark orange","royal blue","burly wood","spring green","gold","light pink","light sky blue","gray","plum")
tsne_fig_colors_dark <- c("red","#00FFFF","pale green","dark orange","royal blue","burly wood","spring green","gold","light pink","steel blue","light sky blue","plum")
DimPlot(treg, reduction = "tsne",label = FALSE, cells.highlight = tsne_fig_list,sizes.highlight = tsne_fig_sizes,shuffle = TRUE) + scale_color_manual(values = tsne_fig_colors)
DimPlot(treg, reduction = "tsne",label = FALSE, cells.highlight = tsne_fig_list,sizes.highlight = tsne_fig_sizes) + scale_color_manual(values = tsne_fig_colors_dark) + DarkTheme()

#Custom color clusters, no expanded TCRs 
tsne_fig_colors_clusters <- c("#00FFFF","pale green","dark orange","royal blue","burly wood","spring green","gold","light pink","light sky blue","gray","plum")
tsne_fig_sizes_clusters <- c(0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55)
DimPlot(treg, reduction = "tsne",label = TRUE, repel = TRUE, sizes.highlight = tsne_fig_sizes_clusters,shuffle = TRUE) + scale_color_manual(values = rev(tsne_fig_colors_clusters))

Idents(treg) <- treg@meta.data$Cell.State
Idents(treg) <- treg$seurat_clusters

DimPlot(treg, cells.highlight = exp_cells, reduction = 'tsne',sizes.highlight = 1.5,label = TRUE) + scale_color_manual(labels = c("Non Expanded","Expanded"),values = c("grey","red"))
TSNEPlot(treg,cells.highlight = exp_cells,sizes.highlight = 1.5,label = TRUE) + scale_color_manual(labels = c("Non Expanded","Expanded"),values = c("grey","red"))

#cluster 9
cluster9.markers <- FindMarkers(treg,ident.1 = 9,ident.2 = c(2),min.pct = 0.25)
cluster9.markers <- cluster9.markers[order(abs(cluster9.markers$avg_logFC),decreasing = TRUE),]

exp_markers <- FindMarkers(treg, ident.1 = clone_counts$clones[1:231])
exp_markers <- FindMarkers(treg, ident.1 = 'Expanded')

nine <- WhichCells(treg, expression = seurat_clusters == 9)
exp_nine_clones <- unique(treg$clonotype_id[which(names(treg$clonotype_id) %in% intersect(nine,exp_cells))])

treg$seurat_clusters[which(treg$clonotype_id %in% exp_nine_clones)]

#cluster 7

cluster7.markers <- FindMarkers(treg,ident.1 = 7,ident.2 = 2,min.pct = 0.25)
cluster7.markers <- cluster7.markers[order(abs(cluster7.markers$avg_logFC),decreasing = TRUE),]

seven <- WhichCells(treg, expression = seurat_clusters == 7)

exp_seven_clones <- treg$clonotype_id[which(names(treg$clonotype_id) %in% intersect(seven,exp_cells))]

exp_seven_clones <- unique(exp_seven_clones)

seven_clones_distribution <- treg$seurat_clusters[which(treg$clonotype_id %in% exp_seven_clones)]

c7_vs_c3 <- FindMarkers(treg,ident.1 = 7,ident.2 = 3,min.pct = 0.25)
c7_vs_c3 <- c7_vs_c3[order(abs(c7_vs_c3$avg_logFC),decreasing = TRUE),]
#Hmgb2 "reported to be involved in the final ligation step in V(D)J recombination."


c = 0
for (i in 1:length(seven_clones_distribution)){if (seven_clones_distribution[i] == 10){c = c + 1}}


#cluster 3
c3 <- FindMarkers(treg, ident.1 = 3, min.pct = 0.25)
c3 <- c3[order(abs(c3$avg_logFC), decreasing = TRUE),]

#cluster 6
cluster6.markers <- FindMarkers(treg,ident.1 = 6, min.pct = 0.25)
cluster6.markers <- cluster6.markers[order(abs(cluster6.markers$avg_logFC),decreasing = TRUE),]


#cluster 8
cluster8.markers <- FindMarkers(treg,ident.1 = 8, min.pct = 0.25)
cluster8.markers <- cluster8.markers[order(abs(cluster8.markers$avg_logFC),decreasing = TRUE),]

#cluster 10
c10 <- FindMarkers(treg,ident.1 = 10, min.pct = 0.25)
c10 <- c10[order(abs(c10$avg_logFC), decreasing = TRUE),]



#sample data
treg <- f_add_sample(treg)

Idents(treg) <- treg$sample

TSNEPlot(treg,pt.size = 1)

FetchData(treg,vars = c('Cell.State','sample'),cells = nine)

          