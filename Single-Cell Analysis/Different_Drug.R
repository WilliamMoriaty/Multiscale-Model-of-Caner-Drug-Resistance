library(dplyr)
library(Seurat)
library(patchwork)
library(ggpubr)
library(RColorBrewer)
setwd("~/Desktop/Rcode")
scRNAdptlist<-list()
# 
diffdrugD0.data <- read.table(gzfile("./GSM3972657_D0.dge.txt.gz"), row.names = 1, header = TRUE, sep = "\t")
scRNAdptlist[[1]] <- CreateSeuratObject(counts = diffdrugD0.data, project = "D0", min.cells = 3, min.features = 200)
diffdrug_ERL.data <- read.table(gzfile("./Different_drug/GSM4494350_PC9D3_ERL1_expression_matrix.txt.gz"), row.names = 1, header = TRUE, sep = "\t")
scRNAdptlist[[2]] <- CreateSeuratObject(counts = diffdrug_ERL.data, project = "ERL", min.cells = 3, min.features = 200)
diffdrug_ERL_CRI.data <- read.table(gzfile("./Different_drug/GSM4494352_PC9D3_CRI_ERL1_expression_matrix.txt.gz"), row.names = 1, header = TRUE, sep = "\t")
scRNAdptlist[[3]] <- CreateSeuratObject(counts = diffdrug_ERL_CRI.data, project = "ERL_CRI", min.cells = 3, min.features = 200)
diffdrug_Eto.data <- read.table(gzfile("./Different_drug/GSM4494354_PC9Day3_Eto_expression_matrix.txt.gz"), row.names = 1, header = TRUE, sep = "\t")
scRNAdptlist[[4]] <- CreateSeuratObject(counts = diffdrug_Eto.data, project = 'Eto', min.cells = 3, min.features = 200)

scRNAdpt<-merge(scRNAdptlist[[1]], y=c(scRNAdptlist[[2]], scRNAdptlist[[3]],scRNAdptlist[[4]]))
scRNAdpt                
table(scRNAdpt@meta.data$orig.ident)

#
diffdrug <- subset(scRNAdpt, subset = nFeature_RNA > 800 & nFeature_RNA < 7500)
table(diffdrug @meta.data$orig.ident)
# Normalization
diffdrug <- NormalizeData(diffdrug, normalization.method = "LogNormalize", scale.factor = 10000)

#
diffdrug <- FindVariableFeatures(diffdrug, selection.method = "vst", nfeatures = 2000)
#Scaling the data
all.genes <- rownames(diffdrug)

diffdrug <- ScaleData(diffdrug, features = all.genes)
#Perform linear dimensional reduction
diffdrug <- RunPCA(diffdrug, features = VariableFeatures(object = diffdrug))

diffdrug <- RunUMAP(diffdrug, dims = 1:15)
ElbowPlot(diffdrug)

diffdrug <- FindNeighbors(diffdrug, dims = 1:15)
diffdrug <- FindClusters(diffdrug, resolution = 1.2) 
head(Idents(diffdrug), 6)

plot<-DimPlot(object = diffdrug, reduction = "umap", group.by = "orig.ident")
plot<-plot+labs(x = "", y = "", title = "")+theme(legend.text = element_text(size = 12)) + guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme(axis.text = element_blank()) + ## 删去所有刻度标签
  theme(axis.ticks = element_blank())+ ## 删去所有刻度线
  theme(axis.line = element_line(color = "white"))+labs(x = "", y = "")

plot


# read cellular response to stress response gene set
file2<-"cellular_response_to_stress.xlsx"
cellular_response_to_stress_gene<-readxl::read_xlsx(file2)
View(cellular_response_to_stress_gene)
gene<-as.list(cellular_response_to_stress_gene)
# revis bug Multilayer not use in Genssay
diffdrug <- JoinLayers(diffdrug)

diffdrug<-AddModuleScore(object = diffdrug,features = gene,ctrl = 100,name = 'CD_Features')
colnames(diffdrug@meta.data)
colnames(diffdrug@meta.data)[6]<-'Epigenetic instability'
FeaturePlot(diffdrug,features = "Epigenetic instability")+ 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))+theme(axis.text = element_blank()) + ## 删去所有刻度标签
  theme(axis.ticks = element_blank())+## 删去所有刻度线
  theme(axis.line = element_line(color = "white"))+labs(x = "", y = "")
##
g2m.genes <- Seurat::cc.genes$g2m.genes
g2m.genes <- rownames(diffdrug)[toupper(rownames(diffdrug)) %in% g2m.genes] # genes in dataset
gene<-list(g2m.genes)
diffdrug <- JoinLayers(diffdrug)

diffdrug<-AddModuleScore(object = diffdrug,features = gene,ctrl = 100,name = 'CD_Features')
colnames(diffdrug@meta.data)
colnames(diffdrug@meta.data)[7]<-'Cell_Proliferation'
FeaturePlot(diffdrug,features = "Cell_Proliferation")+ 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))+theme(axis.text = element_blank()) + ## 删去所有刻度标签
  theme(axis.ticks = element_blank())+## 删去所有刻度线
  theme(axis.line = element_line(color = "white"))+labs(x = "", y = "")

# read diffdrug DE gene set
file2<-"DTP_recluster_upgene.csv"
diffdrug_gene<-read.csv(file2)
diffdrug_gene_top20<-as.data.frame(diffdrug_gene[1:20,1])
colnames(diffdrug_gene_top20)[1]<-"gene"
View(diffdrug_gene_top20)
gene<-as.list(diffdrug_gene_top20)
# revis bug Multilayer not use in Genssay
diffdrug <- JoinLayers(diffdrug)

diffdrug<-AddModuleScore(object = diffdrug,features = diffdrug_gene_top20,ctrl = 100,name = 'CD_Features')
colnames(diffdrug@meta.data)
colnames(diffdrug@meta.data)[8]<-'DTP'
FeaturePlot(diffdrug,features = "DTP")+ 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 100, name = "RdBu")))+theme(axis.text = element_blank()) + ## 删去所有刻度标签
  theme(axis.ticks = element_blank())+## 删去所有刻度线
  theme(axis.line = element_line(color = "white"))+labs(x = "", y = "")


violindt=data.frame( diffdrug@meta.data[1], diffdrug@meta.data[6])

mycomparisons <- list(c("D0", "Eto"))
#violindt<-subset(violindt, !(violindt$orig.ident == c("D11_ERL")))
ggviolin(violindt,x='orig.ident',y='Epigenetic_instability',alpha = 1.0, width = 0.3,fill='orig.ident',
         palette = c("#9ECABE","#F6F5B4","#2F528F","#E3AD68"),draw_quantiles = c( .50))+
  stat_compare_means(comparisons = mycomparisons)+
  labs(x="",y="Gene Score")+theme(text=element_text(face = "bold",family = 'Arial'))+
  theme(text = element_text(size = 15))+NoLegend()

feature.use<-cellular_response_to_stress_gene$gene

gene_cell_exp <- AverageExpression(diffdrug,features = feature.use) 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)

#complexheatmap作图
library(ComplexHeatmap)
#顶部细胞类型注释
df <- data.frame(colnames(gene_cell_exp))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,#细胞名/cluster
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('D0'="#9ECABE",
                                                  'ERL'="#F6F5B4",
                                                  'ERL-CRI'="#2F528F",
                                                  "Eto"="#E3AD68")))#颜色设置
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
Heatmap(marker_exp,
        cluster_rows = T,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = T,
        column_title = NULL,
        heatmap_legend_param = list(
          title=' '),
        col = colorRampPalette(c("#0000EF","white","#CD2626"))(100),
        border = 'black',
        rect_gp = gpar(col = "black", lwd = 1),
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 12),
        top_annotation = top_anno)

