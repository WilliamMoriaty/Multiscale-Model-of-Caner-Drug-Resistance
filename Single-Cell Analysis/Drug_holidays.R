library(dplyr)
library(Seurat)
library(patchwork)
library(ggpubr)
library(RColorBrewer)
#setwd("Desktop/Rcode")
scRNAdptlist<-list()
# 
dtpD0.data <- read.table(gzfile("./Drug Holidays/GSM3972669_D0.dge.txt.gz"), row.names = 1, header = TRUE, sep = "\t")
scRNAdptlist[[1]] <- CreateSeuratObject(counts = dtpD0.data, project = "D0", min.cells = 3, min.features = 200)

dtpD2.data <- read.table(gzfile("./Drug Holidays/GSM3972670_D2.dge.txt.gz"), row.names = 1, header = TRUE, sep = "\t")
scRNAdptlist[[2]] <- CreateSeuratObject(counts = dtpD2.data, project = "D2", min.cells = 3, min.features = 200)

dtpD11_ERL.data <- read.table(gzfile("./Drug Holidays/GSM3972671_D11_ERL.dge.txt.gz"), row.names = 1, header = TRUE, sep = "\t")
scRNAdptlist[[3]] <- CreateSeuratObject(counts = dtpD11_ERL.data, project = "D11_ERL", min.cells = 3, min.features = 200)

dtpD19_DMSO.data <- read.table(gzfile("./Drug Holidays/GSM3972672_D19DMSO.dge.txt.gz"), row.names = 1, header = TRUE, sep = "\t")
scRNAdptlist[[4]] <- CreateSeuratObject(counts = dtpD19_DMSO.data, project = "D19_DMSO", min.cells = 3, min.features = 200)

dtpD19_ERL.data <- read.table(gzfile("./Drug Holidays/GSM3972673_D19_ERL.dge.txt.gz"), row.names = 1, header = TRUE, sep = "\t")
scRNAdptlist[[5]] <- CreateSeuratObject(counts = dtpD19_ERL.data, project = "D19_ERL", min.cells = 3, min.features = 200)

scRNAdpt<-merge(scRNAdptlist[[1]], y=c(scRNAdptlist[[2]], scRNAdptlist[[3]],scRNAdptlist[[4]], scRNAdptlist[[5]]))
#scRNAdpt<-merge(scRNAdptlist[[1]], y=c(scRNAdptlist[[2]], scRNAdptlist[[3]],scRNAdptlist[[4]]))

scRNAdpt 

table(scRNAdpt@meta.data$orig.ident)

scRNAdpt[["percent.mt"]] <- PercentageFeatureSet(scRNAdpt, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(scRNAdpt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 
drugholiday <- subset(scRNAdpt, subset = nFeature_RNA > 800 & nFeature_RNA < 7500)
table(drugholiday @meta.data$orig.ident)
# Normalization
drugholiday <- NormalizeData(drugholiday, normalization.method = "LogNormalize", scale.factor = 10000)

#
drugholiday <- FindVariableFeatures(drugholiday, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(drugholiday), 10)
#Scaling the data
all.genes <- rownames(drugholiday)

drugholiday <- ScaleData(drugholiday, features = all.genes)
#Perform linear dimensional reduction
drugholiday <- RunPCA(drugholiday, features = VariableFeatures(object = drugholiday))

drugholiday <- RunUMAP(drugholiday, dims = 1:15)
ElbowPlot(drugholiday)

drugholiday <- FindNeighbors(drugholiday, dims = 1:15)
drugholiday <- FindClusters(drugholiday, resolution = 1.2) 
head(Idents(drugholiday), 6)

DimPlot(object = drugholiday, reduction = "umap", label = T)+theme(axis.text = element_blank()) + ## 
  theme(axis.ticks = element_blank())+## 
  theme(axis.line = element_line(color = "white"))+labs(x = "", y = "")+theme(text=element_text(face = "bold",family = 'Arial'))+
  theme(text = element_text(size = 15))


dtp<-drugholiday
# read DTP DE gene set
file2<-"DTP_recluster_upgene.csv"
DTP_gene<-read.csv(file2)
DTP_gene_top20<-as.data.frame(DTP_gene[1:20,1])
colnames(DTP_gene_top20)[1]<-"gene"
View(DTP_gene_top20)
gene<-as.list(DTP_gene_top20)
# revis bug Multilayer not use in Genssay
dtp <- JoinLayers(dtp)

dtp<-AddModuleScore(object = dtp,features = DTP_gene_top20,ctrl = 100,name = 'CD_Features')
colnames(dtp@meta.data)
colnames(dtp@meta.data)[7]<-'DTP'
FeaturePlot(dtp,features = "DTP")+ 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))+theme(axis.text = element_blank()) + ## 
  theme(axis.ticks = element_blank())+## 
  theme(axis.line = element_line(color = "white"))+labs(x = "", y = "")
# naive
file2<-"Naive_recluster_upgene.csv"
Naive_gene<-read.csv(file2)
Naive_gene_top20<-as.data.frame(Naive_gene[1:20,1])
colnames(Naive_gene_top20)[1]<-"gene"
View(Naive_gene_top20)
gene<-as.list(Naive_gene_top20)
# revis bug Multilayer not use in Genssay
dtp <- JoinLayers(dtp)

dtp<-AddModuleScore(object = dtp,features = gene,ctrl = 100,name = 'CD_Features')
colnames(dtp@meta.data)
colnames(dtp@meta.data)[8]<-'Naive'
FeaturePlot(dtp,features = "Naive")+ 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))+theme(axis.text = element_blank()) + ## 
  theme(axis.ticks = element_blank())+## 
  theme(axis.line = element_line(color = "white"))+labs(x = "", y = "")

# Resistant
file2<-"Resistant_recluster_upgene.csv"
Resistant_gene<-read.csv(file2)
Resistant_gene_top20<-as.data.frame(Resistant_gene[1:20,1])
colnames(Resistant_gene_top20)[1]<-"gene"
View(Resistant_gene_top20)
gene<-as.list(Resistant_gene_top20)
# revis bug Multilayer not use in Genssay
dtp <- JoinLayers(dtp)

dtp<-AddModuleScore(object = dtp,features = gene,ctrl = 100,name = 'CD_Features')
colnames(dtp@meta.data)
colnames(dtp@meta.data)[9]<-'Resistant'
FeaturePlot(dtp,features = "Resistant")+ 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))+theme(axis.text = element_blank()) + ## 
  theme(axis.ticks = element_blank())+## 
  theme(axis.line = element_line(color = "white"))+labs(x = "", y = "")

plot<-DimPlot(object = drugholiday, reduction = "umap", group.by = "orig.ident",order = c("D19_ERL","D19_DMSO","D11_ERL","D2","D0"))
plot<-plot+labs(x = "", y = "", title = "")+theme(legend.text = element_text(size = 12)) + guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme(axis.text = element_blank()) + ## 
  theme(axis.ticks = element_blank())+
  theme(axis.line = element_line(color = "white"))+theme(text=element_text(face = "bold",family = 'Arial'))+
  theme(text = element_text(size = 15))
  ## 

plot




Idents(dtp) <- plyr::mapvalues(Idents(dtp), from = levels(dtp), 
                               to = c("Resistant","Resistant","Naive","Naive","Naive","Resistant","Naive","Naive","DTP","Resistant","Naive","DTP","Resistant"))
Idents(dtp) <- factor(Idents(dtp), levels = c("Naive","DTP","Resistant"))
plot<-DimPlot(dtp, reduction = "umap")+labs(x = "", y = "", title = "")+theme(legend.text = element_text(size = 12)) + guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme(axis.text = element_blank()) + ## 
  theme(axis.ticks = element_blank())+## 
  theme(axis.line = element_line(color = "white"))
plot

# read cellular response to stress response gene set
file2<-"cellular_response_to_stress.xlsx"
cellular_response_to_stress_gene<-readxl::read_xlsx(file2)
View(cellular_response_to_stress_gene)
gene<-as.list(cellular_response_to_stress_gene)
# revis bug Multilayer not use in Genssay
dtp <- JoinLayers(dtp)

dtp<-AddModuleScore(object = dtp,features = gene,ctrl = 100,name = 'CD_Features')
colnames(dtp@meta.data)
colnames(dtp@meta.data)[10]<-'Epigenetic_instability'



violindt=data.frame( dtp@meta.data[1], dtp@meta.data[10])

mycomparisons <- list(c("D2", "D19_ERL"))
violindt<-subset(violindt, !(violindt$orig.ident == c("D11_ERL")))
ggviolin(violindt,x='orig.ident',y='Epigenetic_instability',alpha = 1.0, width = 0.3,fill='orig.ident',
         palette = c("#9ECABE","#F6F5B4","#2F528F","#E3AD68"),draw_quantiles = c( .50))+
  stat_compare_means(comparisons = mycomparisons)+
  labs(x="",y="Gene Score")+theme(text=element_text(face = "bold",family = 'Arial'))+
  theme(text = element_text(size = 15))+NoLegend()

##########
scRNAdptlist<-list()

dtpD0.data <- read.table(gzfile("./Drug Holidays/GSM3972669_D0.dge.txt.gz"), row.names = 1, header = TRUE, sep = "\t")
scRNAdptlist[[1]] <- CreateSeuratObject(counts = dtpD0.data, project = "D0", min.cells = 3, min.features = 200)

dtpD2.data <- read.table(gzfile("./Drug Holidays/GSM3972670_D2.dge.txt.gz"), row.names = 1, header = TRUE, sep = "\t")
scRNAdptlist[[2]] <- CreateSeuratObject(counts = dtpD2.data, project = "D2", min.cells = 3, min.features = 200)

dtpD19_DMSO.data <- read.table(gzfile("./Drug Holidays/GSM3972672_D19DMSO.dge.txt.gz"), row.names = 1, header = TRUE, sep = "\t")
scRNAdptlist[[3]] <- CreateSeuratObject(counts = dtpD19_DMSO.data, project = "D19_DMSO", min.cells = 3, min.features = 200)

dtpD19_ERL.data <- read.table(gzfile("./Drug Holidays/GSM3972673_D19_ERL.dge.txt.gz"), row.names = 1, header = TRUE, sep = "\t")
scRNAdptlist[[4]] <- CreateSeuratObject(counts = dtpD19_ERL.data, project = "D19_ERL", min.cells = 3, min.features = 200)

scRNAdpt<-merge(scRNAdptlist[[1]], y=c(scRNAdptlist[[2]], scRNAdptlist[[3]],scRNAdptlist[[4]]))
#scRNAdpt<-merge(scRNAdptlist[[1]], y=c(scRNAdptlist[[2]], scRNAdptlist[[3]],scRNAdptlist[[4]]))

scRNAdpt 

table(scRNAdpt@meta.data$orig.ident)

scRNAdpt[["percent.mt"]] <- PercentageFeatureSet(scRNAdpt, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(scRNAdpt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


drugholiday <- subset(scRNAdpt, subset = nFeature_RNA > 800 & nFeature_RNA < 7500)
table(drugholiday @meta.data$orig.ident)
# Normalization
drugholiday <- NormalizeData(drugholiday, normalization.method = "LogNormalize", scale.factor = 10000)

feature.use<-cellular_response_to_stress_gene$gene

gene_cell_exp <- AverageExpression(drugholiday,features = feature.use) 
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
                                                  'D2'="#F6F5B4",
                                                  'D19-DMSO'="#2F528F",
                                                  "D19-ERL"="#E3AD68")))#颜色设置
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




#### find cell numbers in different cell types
# data=matrix(c(1:20),ncol=4,byrow=TRUE)
# colnames(data) = c('Sample','Naive','DTP','Resistant')
# data[1:5]<-c('D0','D2(+Erl)','D11(+Erl,before holiday)','D19(still on holiday)','D19(+Erl,after holiday)')
# cell_type<- as.data.frame(data)
# holiday=data.frame(CellType=Idents(dtp),dtp@meta.data[1])
# D0_naive=subset(holiday,holiday$orig.ident==c("D0") & holiday$CellType==c("Naive"))
# D0_DTP=subset(holiday,holiday$orig.ident==c("D0") & holiday$CellType==c("DTP"))
# D0_Resist=subset(holiday,holiday$orig.ident==c("D0") & holiday$CellType==c("Resistant"))
# # 
# cell_type$Naive[1]=dim(D0_naive)[1]/762
# cell_type$DTP[1]=dim(D0_DTP)[1]/762
# cell_type$Resistant[1]=dim(D0_Resist)[1]/762
# # 
# D2_number=144;
# D2_naive=subset(holiday,holiday$orig.ident==c("D2") & holiday$CellType==c("Naive"))
# D2_DTP=subset(holiday,holiday$orig.ident==c("D2") & holiday$CellType==c("DTP"))
# D2_Resist=subset(holiday,holiday$orig.ident==c("D2") & holiday$CellType==c("Resistant"))
# # 
# cell_type$Naive[2]=dim(D2_naive)[1]/D2_number
# cell_type$DTP[2]=dim(D2_DTP)[1]/D2_number
# cell_type$Resistant[2]=dim(D2_Resist)[1]/D2_number
# # 
# D11_number=720;
# D11_naive=subset(holiday,holiday$orig.ident==c("D11_ERL") & holiday$CellType==c("Naive"))
# D11_DTP=subset(holiday,holiday$orig.ident==c("D11_ERL") & holiday$CellType==c("DTP"))
# D11_Resist=subset(holiday,holiday$orig.ident==c("D11_ERL") & holiday$CellType==c("Resistant"))
# 
# cell_type$Naive[3]=dim(D11_naive)[1]/D11_number
# cell_type$DTP[3]=dim(D11_DTP)[1]/D11_number
# cell_type$Resistant[3]=dim(D11_Resist)[1]/D11_number
# # 
# D19D_number=706;
# D19D_naive=subset(holiday,holiday$orig.ident==c("D19_DMSO") & holiday$CellType==c("Naive"))
# D19D_DTP=subset(holiday,holiday$orig.ident==c("D19_DMSO") & holiday$CellType==c("DTP"))
# D19D_Resist=subset(holiday,holiday$orig.ident==c("D19_DMSO") & holiday$CellType==c("Resistant"))
# 
# cell_type$Naive[4]=dim(D19D_naive)[1]/D19D_number
# cell_type$DTP[4]=dim(D19D_DTP)[1]/D19D_number
# cell_type$Resistant[4]=dim(D19D_Resist)[1]/D19D_number
# # 
# D19E_number=719;
# D19E_naive=subset(holiday,holiday$orig.ident==c("D19_ERL") & holiday$CellType==c("Naive"))
# D19E_DTP=subset(holiday,holiday$orig.ident==c("D19_ERL") & holiday$CellType==c("DTP"))
# D19E_Resist=subset(holiday,holiday$orig.ident==c("D19_ERL") & holiday$CellType==c("Resistant"))
# 
# cell_type$Naive[5]=dim(D19E_naive)[1]/D19E_number
# cell_type$DTP[5]=dim(D19E_DTP)[1]/D19E_number
# cell_type$Resistant[5]=dim(D19E_Resist)[1]/D19E_number
# 
# ggplot(cell_type, aes(x = Sample, y = Naive, group = 1)) + geom_line()
# 
# write.table(cell_type,"holidays.csv",row.names=FALSE,col.names=TRUE,sep=",")


