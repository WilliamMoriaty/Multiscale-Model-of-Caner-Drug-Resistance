library(dplyr)
library(Seurat)
library(patchwork)
library(ggpubr)
library(RColorBrewer)
setwd("~/Desktop/Rcode")
scRNAdptlist<-list()

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

# Remove cells with a high proportion of mitochondrial gene expression, and some extreme cells
dtp <- subset(scRNAdpt, subset = nFeature_RNA > 800 & nFeature_RNA < 7500)
table(dtp@meta.data$orig.ident)
# Normalization
dtp <- NormalizeData(dtp, normalization.method = "LogNormalize", scale.factor = 10000)

#Default value，“vst“method to choose 1000 high variable genes
dtp <- FindVariableFeatures(dtp, selection.method = "vst", nfeatures = 1000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dtp), 10)
#Scaling the data
all.genes <- rownames(dtp)
dtp <- ScaleData(dtp, features = all.genes)

dtp <- RunPCA(dtp, features = VariableFeatures(object = dtp))

ElbowPlot(dtp)

dtp <- FindNeighbors(dtp, dims = 1:20)
dtp <- FindClusters(dtp, resolution = 3.0) 
head(Idents(dtp), 6)

file2<-"cellular_response_to_stress.xlsx"
cellular_response_to_stress_gene<-readxl::read_xlsx(file2)
View(cellular_response_to_stress_gene)
cellular_response_to_stress_gene <- rownames(dtp)[toupper(rownames(dtp)) %in% cellular_response_to_stress_gene$gene] # genes in dataset
cellular_response_to_stress_gene
gene<-list(cellular_response_to_stress_gene)

# revis bug Multilayer not use in Genssay
dtp <- JoinLayers(dtp)

dtp<-AddModuleScore(object = dtp,features = gene,ctrl = 100,name = 'CD_Features')
colnames(dtp@meta.data)
colnames(dtp@meta.data)[7]<-'Epigenetic_instability'
## pseudobulk gene expression per cell-type

getPseudobulk <- function(mat, celltype) {
  mat.summary <- do.call(cbind, lapply(levels(celltype), function(ct) {
    cells <- names(celltype)[celltype==ct]
    pseudobulk <- rowSums(mat[, cells])
    return(pseudobulk)
  }))
  colnames(mat.summary) <- levels(celltype)
  return(mat.summary)
}
mat.summary <- getPseudobulk(dtp[["RNA"]]$counts,dtp@active.ident)
matr.bulk=data.matrix(mat.summary)


library(GSVA)
file2<-"cellular_response_to_stress.xlsx"
cellular_response_to_stress_gene<-readxl::read_xlsx(file2)
View(cellular_response_to_stress_gene)
cellular_response_to_stress_gene <- rownames(dtp)[toupper(rownames(dtp)) %in% cellular_response_to_stress_gene$gene] # genes in dataset
cellular_response_to_stress_gene
gene1<-list(cellular_response_to_stress_gene)
# gsva<-gsva(expr = matr.bulk,gene)
# ?gsva

##IA
file2<-"TIME_IA.xlsx"
TIME_IA<-readxl::read_xlsx(file2)

IA.genes <- rownames(dtp)[toupper(rownames(dtp)) %in% TIME_IA$gene_IA] # genes in dataset
IA.genes
gene2<-list(IA.genes)
#ISM
file2<-"TIME_ISM.xlsx"
TIME_ISM<-readxl::read_xlsx(file2)
ISM.genes <- rownames(dtp)[toupper(rownames(dtp)) %in% TIME_ISM$gene_ISM] # genes in dataset
ISM.genes
gene3<-list(ISM.genes)
#ISS
file2<-"TIME_ISS.xlsx"
TIME_ISS<-readxl::read_xlsx(file2)
ISS.genes <- rownames(dtp)[toupper(rownames(dtp)) %in% TIME_ISS$gene_ISS] # genes in dataset
ISS.genes
gene4<-list(ISS.genes)

#IE
file2<-"TIME_IE.xlsx"
TIME_IE<-readxl::read_xlsx(file2)
IE.genes <- rownames(dtp)[toupper(rownames(dtp)) %in% TIME_IE$gene_IE] # genes in dataset
IE.genes
gene5<-list(IE.genes)
#IR
file2<-"TIME_IR.xlsx"
TIME_IR<-readxl::read_xlsx(file2)
IR.genes <- rownames(dtp)[toupper(rownames(dtp)) %in% TIME_IR$gene_IR] # genes in dataset
IR.genes
gene6<-list(IR.genes)
gene<-c(gene1,gene2,gene3,gene4,gene5,gene6)

library(psych)
gsva<-gsva(expr = matr.bulk,gene)
rownames(gsva)<-c("Epigen_Instability","IA","ISM","ISS","IE","IR")
gsva<-t(gsva)
cor_pearson <- cor(gsva, method = 'pearson')
cor_pearson<-data.frame(cor_pearson)
cor_spearman<-cor(gsva,method="spearman")

P11<-cor.test(gsva[,1],gsva[,2],method = 'pearson')$p.value
P12<-cor.test(gsva[,1],gsva[,3],method = 'pearson')$p.value
P13<-cor.test(gsva[,1],gsva[,4],method = 'pearson')$p.value
P14<-cor.test(gsva[,1],gsva[,5],method = 'pearson')$p.value
P15<-cor.test(gsva[,1],gsva[,6],method = 'pearson')$p.value


write.table(cor_pearson,"Pearson_TIME.csv",row.names=TRUE,col.names=TRUE,sep=",")
