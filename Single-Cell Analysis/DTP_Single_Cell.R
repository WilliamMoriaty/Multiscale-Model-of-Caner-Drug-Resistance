library(dplyr)
library(Seurat)
library(patchwork)
library(ggpubr)
#setwd("~/Desktop/Rcode/")
scRNAdptlist<-list()
# 使用read.table()函数从txt.gz格式的文件中读取数据，并将第一列作为行名
dtpD0.data <- read.table(gzfile("./GSM3972657_D0.dge.txt.gz"), row.names = 1, header = TRUE, sep = "\t")
scRNAdptlist[[1]] <- CreateSeuratObject(counts = dtpD0.data, project = "D0", min.cells = 3, min.features = 200)
dtpD1.data <- read.table(gzfile("./GSM3972658_D1.dge.txt.gz"), row.names = 1, header = TRUE, sep = "\t")
scRNAdptlist[[2]] <- CreateSeuratObject(counts = dtpD1.data, project = "D1", min.cells = 3, min.features = 200)
dtpD2.data <- read.table(gzfile("./GSM3972659_D2.dge.txt.gz"), row.names = 1, header = TRUE, sep = "\t")
scRNAdptlist[[3]] <- CreateSeuratObject(counts = dtpD2.data, project = "D2", min.cells = 3, min.features = 200)
dtpD4.data <- read.table(gzfile("./GSM3972660_D4.dge.txt.gz"), row.names = 1, header = TRUE, sep = "\t")
scRNAdptlist[[4]] <- CreateSeuratObject(counts = dtpD4.data, project = "D4", min.cells = 3, min.features = 200)
dtpD9.data <- read.table(gzfile("./GSM3972661_D9.dge.txt.gz"), row.names = 1, header = TRUE, sep = "\t")
scRNAdptlist[[5]] <- CreateSeuratObject(counts = dtpD9.data, project = "D9", min.cells = 3, min.features = 200)
dtpD11.data <- read.table(gzfile("./GSM3972662_D11.dge.txt.gz"), row.names = 1, header = TRUE, sep = "\t")
scRNAdptlist[[6]] <- CreateSeuratObject(counts = dtpD11.data, project = "D11", min.cells = 3, min.features = 200)

scRNAdpt<-merge(scRNAdptlist[[1]], y=c(scRNAdptlist[[2]], scRNAdptlist[[3]],scRNAdptlist[[4]], scRNAdptlist[[5]],scRNAdptlist[[6]]))
scRNAdpt                
table(scRNAdpt@meta.data$orig.ident)

scRNAdpt[["percent.mt"]] <- PercentageFeatureSet(scRNAdpt, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(scRNAdpt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 去除线粒体基因表达比例过高的细胞，和一些极值细胞
dtp <- subset(scRNAdpt, subset = nFeature_RNA > 800 & nFeature_RNA < 7500)
table(dtp@meta.data$orig.ident)
# Normalization
dtp <- NormalizeData(dtp, normalization.method = "LogNormalize", scale.factor = 10000)

#我们使用默认参数，即“vst”方法选取1000个高变基因
dtp <- FindVariableFeatures(dtp, selection.method = "vst", nfeatures = 1000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dtp), 10)
#Scaling the data
all.genes <- rownames(dtp)
dtp <- ScaleData(dtp, features = all.genes)
#Perform linear dimensional reduction
dtp <- RunPCA(dtp, features = VariableFeatures(object = dtp))

dtp <- RunUMAP(dtp, dims = 1:12)
ElbowPlot(dtp)

dtp <- FindNeighbors(dtp, dims = 1:12)
dtp <- FindClusters(dtp, resolution = 1.0) 
head(Idents(dtp), 6)

#UMAP
DimPlot(object = dtp, reduction = "umap", label = T)

plot<-DimPlot(object = dtp, reduction = "umap", group.by = "orig.ident",order = c("D11","D9","D4","D2","D1","D0"))
plot<-plot+labs(x = "", y = "", title = "")+theme(legend.text = element_text(size = 12)) + guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme(axis.text = element_blank()) + ## 删去所有刻度标签
  theme(axis.ticks = element_blank())## 删去所有刻度线
plot
feature.use <- c("TOP2A","ITGA6","MKI67","INHBA","CYP1B1","SPARC")
FeaturePlot(dtp, features = feature.use)

Idents(dtp) <- plyr::mapvalues(Idents(dtp), from = levels(dtp), to = c("DTP","Naive","Naive","DTP","Resistant","Naive","Naive","DTP","Resistant","Naive","DTP"))
#DimPlot(object = dtp, reduction = "umap")
Idents(dtp) <- factor(Idents(dtp), levels = c("Naive","DTP","Resistant"))
plot<-DimPlot(dtp, reduction = "umap")
library(ggplot2)
plot <- plot +labs(x = "", y = "", title = "")+
  theme(legend.text = element_text(size = 12)) + guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme(axis.text = element_blank()) + ## 删去所有刻度标签
  theme(axis.ticks = element_blank())## 删去所有刻度线

plot

DotPlot(dtp,features = feature.use,cols = c("lightgrey","red"))

#计算平均表达量
gene_cell_exp <- AverageExpression(dtp,
                                   features = feature.use) 
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
                             col = list(class = c('Naive'="#9ECABE",
                                                  'DTP'="#F6F5B4",
                                                  'Transitional DTP'="#2F528F",
                                                  "Resistant"="#E3AD68")))#颜色设置
#数据标准化缩放一下
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
Heatmap(marker_exp,
        cluster_rows = T,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = T,
        column_title = NULL,
        heatmap_legend_param = list(
          title=' '),
        col = colorRampPalette(c("#0000EF","black","#FDFE00"))(100),
        border = 'black',
        rect_gp = gpar(col = "black", lwd = 1),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        top_annotation = top_anno)


# find markers
dtp <- JoinLayers(dtp)
DTP.markers <- FindMarkers(dtp, ident.1 = "DTP",only.pos = TRUE)
DTP.markers<- subset(DTP.markers, subset = p_val_adj<0.05 & avg_log2FC > 0.25)
DTP.upgene<- data.frame(rownames(DTP.markers))
write.table(DTP.upgene,"DTP_upgene.csv",row.names=FALSE,col.names=TRUE,sep=",")


Naive.markers <- FindMarkers(dtp, ident.1 = "Naive",only.pos = TRUE)
top20 <- head(Naive.markers, 20)
top20<- data.frame(rownames(top20))
write.table(top20,"Naive_marker.csv",row.names=FALSE,col.names=TRUE,sep=",")

DTP.markers <- FindMarkers(dtp, ident.1 = "DTP",only.pos = TRUE)
top20 <- head(DTP.markers, 20)
top20<- data.frame(rownames(top20))
write.table(top20,"DTP_marker.csv",row.names=FALSE,col.names=TRUE,sep=",")


#data = read.table("gene.txt",header = F,sep=",")#读取文本数据
#write.table(data,"cellular_response_to_stress.csv",row.names=FALSE,col.names=TRUE,sep=",")

# G2 and M phases gene set in cell cycle
g2m.genes <- Seurat::cc.genes$g2m.genes
g2m.genes <- rownames(dtp)[toupper(rownames(dtp)) %in% g2m.genes] # genes in dataset
gene<-list(g2m.genes)
dtp <- JoinLayers(dtp)

dtp<-AddModuleScore(object = dtp,features = gene,ctrl = 100,name = 'CD_Features')
colnames(dtp@meta.data)
colnames(dtp@meta.data)[7]<-'Cell_Proliferation_Score'

mycomparisons <- list(c("Naive", "DTP"), c("DTP","Resistant"),c("Naive", "Resistant"))

violindt = data.frame(Celltype = Idents(dtp), dtp@meta.data[7])
violindt<-subset(violindt, !(violindt$Celltype == "Transitional DTP"))  
ggviolin(violindt,x='Celltype',y='Cell_Proliferation_Score',fill='Celltype',palette = c('lancet'),draw_quantiles = c( .50))+
  stat_compare_means(comparisons = mycomparisons)+labs(x="",y="Gene expression")+NoLegend()+theme(text=element_text(
    face = "bold",family = 'Arial'))+theme(text = element_text(size = 15))

# read cellular response to stress response gene set
file2<-"cellular_response_to_stress.xlsx"
cellular_response_to_stress_gene<-readxl::read_xlsx(file2)
View(cellular_response_to_stress_gene)
gene<-as.list(cellular_response_to_stress_gene)
# revis bug Multilayer not use in Genssay
dtp <- JoinLayers(dtp)

dtp<-AddModuleScore(object = dtp,features = gene,ctrl = 100,name = 'CD_Features')
colnames(dtp@meta.data)
colnames(dtp@meta.data)[8]<-'Cellular_response_to_stress_Score'

mycomparisons <- list(c("Naive", "DTP"), c("DTP","Resistant"))

violindt = data.frame(Celltype = Idents(dtp), dtp@meta.data[8])
violindt<-subset(violindt, !(violindt$Celltype == "Transitional DTP"))  
ggviolin(violindt,x='Celltype',y='Cellular_response_to_stress_Score',fill='Celltype',palette = c('lancet'),draw_quantiles = c( .50))+
  stat_compare_means(comparisons = mycomparisons)+labs(x="",y="Gene Score")+NoLegend()+theme(text=element_text(
    face = "bold",family = 'Arial'))+theme(text = element_text(size = 15))

# S phase gene set in cell cylce
s.genes <- Seurat::cc.genes$s.genes
s.genes <- rownames(dtp)[toupper(rownames(dtp)) %in% s.genes] # genes in dataset


# Adaptation_stress_response
file3<-"Adaptation_stress_response.xlsx"
gene<-readxl::read_xlsx(file3)
View(gene)
gene<-as.list(gene)
# revis bug Multilayer not use in Genssay
dtp <- JoinLayers(dtp)

dtp<-AddModuleScore(object = dtp,features = gene,ctrl = 100,name = 'CD_Features')
colnames(dtp@meta.data)
colnames(dtp@meta.data)[10]<-'Adaptation_stress_response_Score'

mycomparisons <- list(c("Naive", "DTP"), c("DTP","Resistant"))

violindt = data.frame(Celltype = Idents(dtp), dtp@meta.data[10])
violindt<-subset(violindt, !(violindt$Celltype == "Transitional DTP"))  
ggviolin(violindt,x='Celltype',y='Adaptation_stress_response_Score',fill='Celltype',palette = c('lancet'),draw_quantiles = c( .50))+
  stat_compare_means(comparisons = mycomparisons)+labs(x="",y="Gene expression")+NoLegend()+theme(text=element_text(
    face = "bold",family = 'Arial'))+theme(text = element_text(size = 15))

