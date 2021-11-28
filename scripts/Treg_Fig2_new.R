library(readxl)
library(readr)
library(sjmisc)
library(ggplot2)
library(hash)
library(pheatmap)
library(EnhancedVolcano)

STAR_gene_counts <- read_csv("C:/Users/bh719/Dropbox (Partners HealthCare)/Harvard CyTof/for Brandon/Sally/STAR_Gene_Counts.csv")

#2 genes have duplicate entries, 1-Mar and 2-Mar
STAR_gene_counts <- STAR_gene_counts[!duplicated(STAR_gene_counts$Gene_ID),]

#Convert fully capitalized gene Ids such that only the first letter is capitalized 

Cap <- function(g){
  g <- paste(toupper(substring(g,1,1)), tolower(substring(g,2)), sep = '')
  return(g)
}

STAR_gene_counts$Gene_ID <- sapply(STAR_gene_counts$Gene_ID, Cap)

#DESeq2
library(DESeq2)
library(apeglm)
library("writexl")

colmetadat <- data.frame(Injury = c(rep("Uninjured",9),rep("7D after Injury",10)),CD44 = c(rep("CD44 High",4),rep("CD44 Low",5),rep("CD44 High",4),rep("CD44 Low",6)))
row.names(colmetadat) <- colnames(STAR_gene_counts)[2:length(colnames(STAR_gene_counts))]
STAR_gene_counts <- STAR_gene_counts[!duplicated(STAR_gene_counts$Gene_ID),]
gene_row <- STAR_gene_counts$Gene_ID
cmat <- STAR_gene_counts
cmat <- cmat[,!(names(cmat) %in% c('Gene_ID'))]
row.names(cmat) <- gene_row

dds <- DESeqDataSetFromMatrix(countData = cmat,colData = colmetadat,design = ~ CD44 + Injury)
dds <- dds[rowSums(counts(dds)) >= 10,]
dds$group <- factor(paste0(dds$CD44,dds$Injury))
design(dds) <- ~ group

dds$group = factor(dds$group, levels = c('CD44 LowUninjured','CD44 HighUninjured','CD44 Low7D after Injury','CD44 High7D after Injury'))
dds <- DESeq(dds)
high_low <- lfcShrink(dds, coef = "group_CD44.HighUninjured_vs_CD44.LowUninjured", type = "apeglm")
low7D_low <- lfcShrink(dds, coef = "group_CD44.Low7D.after.Injury_vs_CD44.LowUninjured", type = "apeglm")
dds$group <- factor(dds$group, levels = c('CD44 HighUninjured','CD44 High7D after Injury','CD44 LowUninjured','CD44 Low7D after Injury'))
dds <- DESeq(dds)
high7D_high <- lfcShrink(dds,coef = "group_CD44.High7D.after.Injury_vs_CD44.HighUninjured", type = "apeglm")
dds$group <- factor(dds$group, levels = c('CD44 Low7D after Injury','CD44 High7D after Injury','CD44 LowUninjured','CD44 HighUninjured'))
dds <- DESeq(dds)
high7D_low7D <- lfcShrink(dds,coef = "group_CD44.High7D.after.Injury_vs_CD44.Low7D.after.Injury", type = "apeglm")

curr_res <- high_low
curr_df <- data.frame("GeneID" = rownames(curr_res), "AveExpr" = curr_res$baseMean,"Log2FC" = curr_res$log2FoldChange, "P Value" = curr_res$pvalue, "P Adjust" = curr_res$padj)
curr_df <- curr_df[order(abs(curr_df$Log2FC), decreasing = TRUE),]
high_low_df <- data.frame("GeneID" = curr_df$GeneID, signif(curr_df[!names(curr_df) == "GeneID"], digits = 3))

curr_res <- low7D_low
curr_df <- data.frame("GeneID" = rownames(curr_res), "AveExpr" = curr_res$baseMean,"Log2FC" = curr_res$log2FoldChange, "P Value" = curr_res$pvalue, "P Adjust" = curr_res$padj)
curr_df <- curr_df[order(abs(curr_df$Log2FC), decreasing = TRUE),]
low7D_low_df <- data.frame("GeneID" = curr_df$GeneID, signif(curr_df[!names(curr_df) == "GeneID"], digits = 3))

curr_res <- high7D_high
curr_df <- data.frame("GeneID" = rownames(curr_res), "AveExpr" = curr_res$baseMean,"Log2FC" = curr_res$log2FoldChange, "P Value" = curr_res$pvalue, "P Adjust" = curr_res$padj)
curr_df <- curr_df[order(abs(curr_df$Log2FC), decreasing = TRUE),]
high7D_high_df <- data.frame("GeneID" = curr_df$GeneID, signif(curr_df[!names(curr_df) == "GeneID"], digits = 3))

curr_res <- high7D_low7D
curr_df <- data.frame("GeneID" = rownames(curr_res), "AveExpr" = curr_res$baseMean,"Log2FC" = curr_res$log2FoldChange, "P Value" = curr_res$pvalue, "P Adjust" = curr_res$padj)
curr_df <- curr_df[order(abs(curr_df$Log2FC), decreasing = TRUE),]
high7D_low7D_df <- data.frame("GeneID" = curr_df$GeneID, signif(curr_df[!names(curr_df) == "GeneID"], digits = 3))

write_xlsx(high_low_df,"C:/Users/bhanc/Dropbox (Partners HealthCare)/Harvard CyTof/for Brandon/Sally/DESeq2/DESeq2_Uninjured_High_vs_Uninjured_Low.xlsx")
write_xlsx(low7D_low_df,"C:/Users/bhanc/Dropbox (Partners HealthCare)/Harvard CyTof/for Brandon/Sally/DESeq2/DESeq2_7D_Low_vs_Uninjured_Low.xlsx")
write_xlsx(high7D_high_df,"C:/Users/bhanc/Dropbox (Partners HealthCare)/Harvard CyTof/for Brandon/Sally/DESeq2/DESeq2_7D_High_vs_Uninjured_High.xlsx")
write_xlsx(high7D_low7D_df,"C:/Users/bhanc/Dropbox (Partners HealthCare)/Harvard CyTof/for Brandon/Sally/DESeq2/DESeq2_7D_High_vs_7D_Low.xlsx")


#Enhanced Volcano

#low7D vs low

l7dl_genes <- c('HSPA1A','HSPA1B','SUN2','ZFP658','PIK3IP1','S100A6','IL7R','DDIT2','CCNB2','CXCR6','CCR10','KLRG1','HRH4','GAS2L3','LGALS3')
l7dl_genes <- sapply(l7dl_genes, Cap)

keyvals.colour <- ifelse(
  low7D_low$log2FoldChange < -1 & low7D_low$pvalue < 1e-5, '#6666FF',
  ifelse(low7D_low$log2FoldChange > 1 & low7D_low$pvalue < 1e-5, '#FF6666',
         'grey')
)
keyvals.colour[is.na(keyvals.colour)] <- 'grey'
names(keyvals.colour)[keyvals.colour == '#6666FF'] <- 'Low'
names(keyvals.colour)[keyvals.colour == 'grey'] <- 'Mid'
names(keyvals.colour)[keyvals.colour == '#FF6666'] <- 'High'

#H 615, W 763
EnhancedVolcano(low7D_low,lab = rownames(low7D_low), title = "",subtitle = '',
                x = 'log2FoldChange',y = 'pvalue',
                labSize = 6.0, colCustom = keyvals.colour,
                titleLabSize = 1,subtitleLabSize = 1,axisLabSize = 22,
                boxedLabels = TRUE,drawConnectors = TRUE,
                selectLab = l7dl_genes,
                xlab = bquote(~Log[2]~ 'Fold Change'))

#high 7D vs high

h7dh_genes <- c('HSPA1A','DLK1','IRS2','CCL5','CCR10','SATB1','CXCR5','CCL8','CD160','CLU','CD9','ISG15','CD22','TBX21')
h7dh_genes <- sapply(h7dh_genes, Cap)

keyvals.colour <- ifelse(
  high7D_high$log2FoldChange < -1 & high7D_high$pvalue < 1e-5, '#6666FF',
  ifelse(high7D_high$log2FoldChange > 1 & high7D_high$pvalue < 1e-5, '#FF6666',
         'grey')
)
keyvals.colour[is.na(keyvals.colour)] <- 'grey'
names(keyvals.colour)[keyvals.colour == '#6666FF'] <- 'Low'
names(keyvals.colour)[keyvals.colour == 'grey'] <- 'Mid'
names(keyvals.colour)[keyvals.colour == '#FF6666'] <- 'High'

EnhancedVolcano(high7D_high,lab = rownames(high7D_high), title = "",subtitle = '',
                x = 'log2FoldChange',y = 'pvalue',
                labSize = 6.0, colCustom = keyvals.colour,
                titleLabSize = 1,subtitleLabSize = 1,axisLabSize = 22,
                boxedLabels = TRUE,drawConnectors = TRUE,
                selectLab = h7dh_genes,
                xlab = bquote(~Log[2]~ 'Fold Change'))

#high vs low

keyvals.colour <- ifelse(
  high_low$log2FoldChange < -1 & high_low$pvalue < 1e-5, '#6666FF',
  ifelse(high_low$log2FoldChange > 1 & high_low$pvalue < 1e-5, '#FF6666',
         'grey')
)
keyvals.colour[is.na(keyvals.colour)] <- 'grey'
names(keyvals.colour)[keyvals.colour == '#6666FF'] <- 'Low'
names(keyvals.colour)[keyvals.colour == 'grey'] <- 'Mid'
names(keyvals.colour)[keyvals.colour == '#FF6666'] <- 'High'

EnhancedVolcano(high_low,lab = rownames(high_low), title = "",subtitle = '',
                x = 'log2FoldChange',y = 'pvalue',
                labSize = 6.0, colCustom = keyvals.colour,
                titleLabSize = 1,subtitleLabSize = 1,axisLabSize = 22,
                boxedLabels = TRUE,drawConnectors = TRUE,
                selectLab = hl_genes,
                xlab = bquote(~Log[2]~ 'Fold Change'))

EnhancedVolcano(high_low, lab = rownames(high_low), title = "", subtitle = "",
                x = 'log2FoldChange',y = 'pvalue',
                colCustom = keyvals.colour,
                xlab = bquote(~Log[2]~ 'Fold Change'))

#7Dhigh vs 7D low
keyvals.colour <- ifelse(
  high7D_low7D$log2FoldChange < -1 & high7D_low7D$pvalue < 1e-5, '#6666FF',
  ifelse(high7D_low7D$log2FoldChange > 1 & high7D_low7D$pvalue < 1e-5, '#FF6666',
         'grey')
)
keyvals.colour[is.na(keyvals.colour)] <- 'grey'
names(keyvals.colour)[keyvals.colour == '#6666FF'] <- 'Low'
names(keyvals.colour)[keyvals.colour == 'grey'] <- 'Mid'
names(keyvals.colour)[keyvals.colour == '#FF6666'] <- 'High'


EnhancedVolcano(high7D_low7D, lab = rownames(high_low), title = "", subtitle = "",
                x = 'log2FoldChange',y = 'pvalue',
                colCustom = keyvals.colour,
                xlab = bquote(~Log[2]~ 'Fold Change'))


#PCA
vsd <- vst(dds,blind = FALSE)


pcaData <- plotPCA(vsd, intgroup=c("Injury","CD44"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = Injury, shape = CD44)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme(text = element_text(size = 20))

#scale_fill_discrete(name = "CD44", labels = c("High","Loow")) #this did not work :(



#Genes with high variance 

var_list <- f_get_var(assay(vsd))
var_list <- var_list[order(var_list,decreasing = TRUE)]
gene_list <- head(names(var_list),2000)





#Upset
library(mygene)
library(UpSetR)
library(hash)
library(sjmisc)

all_genes <- STAR_gene_counts$Gene_ID[!duplicated(STAR_gene_counts$Gene_ID)]

res <- queryMany(gene_list,scopes = 'symbol', fields=c('entrezgene','ensembl.gene','go','description'),species = 'mouse')
res <- res[!duplicated(res$query),]

genes_term_map <- f_genes_term_map(res,gene_list)
genes_func_map <- f_genes_func_map(genes_term_map,f_upset_dic(f_upset_funcs()))
listInput <- f_listInput(genes_func_map)
upset(fromList(listInput), order.by = "freq",mainbar.y.label = 'Gene Count: Intersection',sets.x.label = 'Gene Count: Gene Ontology',nsets = 9,set_size.scale_max = 1000,set_size.show = TRUE,set_size.angles = 0,text.scale = 1.9)



#Heat Map
#Cytokines and Cell Surface Markers
gene_list <- toupper(gene_list) #5/25/2021
diff_cytokines_cd <- gene_list[grep("^CCL|^CXC|^IFN|^IL|^TNF|CD40LG|FASL|CD70|TGFB",gene_list)] 
treg_genes <- c("TGFB","IL10","ENTPD1","NT5E","LAG3",'TIGIT','CTLA4','ITGAE','KLRG1','ICOS','IL10RA','FGL2','HAVCR2','CD83')
treg_genes <- sapply(treg_genes, Cap)

gene_list_heat_cyto <- c(f_diff_cytokines(gene_list),f_treg_genes())
gene_list_heat_cyto <- gene_list_heat_cyto[gene_list_heat_cyto %in% gene_list]
gene_list_heat_cyto <- gene_list_heat_cyto[!duplicated(gene_list_heat_cyto)]
gene_list_heat_cyto <- sapply(gene_list_heat_cyto, Cap) #5/25/2021

curr_heat <- STAR_gene_counts[STAR_gene_counts$Gene_ID %in% gene_list_heat_cyto,]
heat_names <- curr_heat$Gene_ID
curr_heat <- curr_heat[,!(names(curr_heat) %in% c('Gene_ID'))]
row.names(curr_heat) <- heat_names
#curr_heat <- log(curr_heat + 1, base = 2)
curr_heat <- t(apply(curr_heat, 1, cal_z_score))

Treg <- f_treg_annot(toupper(gene_list_heat_cyto))
names(Treg) <- sapply(names(Treg), Cap)
gene_row_annot <- data.frame(Treg = Treg)
annot_colors <- list(Treg = c('Treg Activity' = '#20B2AA','Other' = '#DCDCDC'),
                     Injury = c('Uninjured' = '#474443','7D after Injury' = '#FFA500'),
                     CD44 = c('CD44 High' = '#9932CC', 'CD44 Low' = '#FFB6C1'))
CD44 <- data.frame(CD44 = c(rep("CD44 High",4),rep("CD44 Low",5),rep("CD44 High",4),rep("CD44 Low",6)))
row.names(CD44) <- colnames(curr_heat)
Injury <- data.frame(Injury = c(rep("Uninjured",9),rep("7D after Injury",10)))
row.names(Injury) <- colnames(curr_heat)
sample_col_annot <- cbind(CD44,Injury)

row.names(curr_heat) <- sapply(row.names(curr_heat), Cap)

pheatmap(curr_heat, annotation_row = gene_row_annot,annotation_colors = annot_colors,cutree_cols = 2,annotation_col = sample_col_annot,annotation_names_col = FALSE,annotation_names_row = FALSE,show_colnames = FALSE,fontsize_row = 7.5,main ='Cytokines & Cell Surface Markers',border_color = NA, treeheight_row = 0)

#Heat Map
#Transcription Factors
diff_transcription_facs <- f_transcription_genes(res,gene_list)
diff_transcription_facs <- diff_transcription_facs[1:66]
diff_transcription_facs <- sapply(diff_transcription_facs, Cap) #5/25/2021

curr_heat <- STAR_gene_counts[STAR_gene_counts$Gene_ID %in% diff_transcription_facs,]
heat_names <- curr_heat$Gene_ID
curr_heat <- curr_heat[,!(names(curr_heat) %in% c('Gene_ID'))]
row.names(curr_heat) <- heat_names
curr_heat <- t(apply(curr_heat, 1, cal_z_score))
annot_colors <- list(Injury = c('Uninjured' = '#474443','7D after Injury' = '#FFA500'),
                     CD44 = c('CD44 High' = '#9932CC', 'CD44 Low' = '#FFB6C1'))
CD44 <- data.frame(CD44 = c(rep("CD44 High",4),rep("CD44 Low",5),rep("CD44 High",4),rep("CD44 Low",6)))
row.names(CD44) <- colnames(curr_heat)
Injury <- data.frame(Injury = c(rep("Uninjured",9),rep("7D after Injury",10)))
row.names(Injury) <- colnames(curr_heat)
sample_col_annot <- cbind(CD44,Injury)
row.names(curr_heat) <- sapply(row.names(curr_heat), Cap)
pheatmap(curr_heat,annotation_colors = annot_colors,cutree_cols = 2,annotation_col = sample_col_annot,annotation_names_col = FALSE,annotation_names_row = FALSE,show_colnames = FALSE,fontsize_row = 7.5,main ='Transcription Factors',border_color = NA, treeheight_row = 0)



