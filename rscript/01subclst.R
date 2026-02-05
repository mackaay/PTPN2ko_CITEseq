load(file = "./matrix_cellranger.RData")
load(file = "./matrix_cellranger_genename.RData")
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)
library(plyr)
library(RColorBrewer)
pal.group <- c("lightgrey", "mediumpurple")
names(pal.group) <- c("WT", "KO")

pal.type <- brewer.pal(7,"Dark2")
names(pal.type) <- c("Tcm", "Tem","Teff", "Ttdem", "Tscm", "Tpex", "Ttex")


#####0. all #####






#####1. subset samples WT and KO #####
idx <- grep("WT", cmo$Assignment)
keep <- intersect(colnames(mat), cmo$Barcode[idx])
sample.mat <- mat[,keep]
adt <- as.matrix(mat[24047:24054,])
rownames(adt) <- gsub("_TotalSeqB", "", rownames(adt))
hto <- as.matrix(mat[24055:nrow(mat),])

cbmc <- CreateSeuratObject(counts = sample.mat[1:24046,])
cbmc[["percent.mt"]] <- PercentageFeatureSet(cbmc, pattern = "MT-")
cbmc[["percent.rb"]] <- PercentageFeatureSet(cbmc, pattern = "^RP[SL]")

VlnPlot(cbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(cbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
cbmc <- subset(cbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)

cbmc <- FindVariableFeatures(cbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(cbmc), 10)
plot1 <- VariableFeaturePlot(cbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

####ADT
cbmc[["ADT"]] <- CreateAssayObject(counts = adt[,keep])
# Validate that the object now contains multiple assays
Assays(cbmc)
cbmc <- NormalizeData(cbmc, assay = "ADT", normalization.method = "CLR", margin = 2)

###HTOdemux 
cbmc[["HTO"]] <- CreateAssayObject(counts = hto[1:2,keep])
cbmc <- NormalizeData(cbmc, assay = "HTO", normalization.method = "CLR")
cbmc <- HTODemux(cbmc, assay = "HTO", positive.quantile = 0.99)
# Global classification results
table(cbmc$HTO_classification.global)
# Group cells based on the max HTO signal
Idents(cbmc) <- "HTO_maxID"
cbmc[["HTO.ident"]] <- Idents(cbmc)



##Finding cluster by RNA
DefaultAssay(cbmc) <- "RNA"
# standard log-normalization
cbmc <- NormalizeData(cbmc)
# choose ~1k variable features
cbmc <- FindVariableFeatures(cbmc)
# standard scaling (no regression)
cbmc <- ScaleData(cbmc)
# Run PCA, select PCs for tSNE visualization and graph-based clustering
cbmc <- RunPCA(cbmc, verbose = TRUE)
ElbowPlot(cbmc, ndims = 50)
# Cluster the cells using the first 25 principal components.
cbmc <- FindNeighbors(cbmc, dims = 1:20)
cbmc <- FindClusters(cbmc, resolution = 0.2) #0.8 for 14 subclusters

#cbmc <- RunTSNE(cbmc, dims = 1:25, method = "FIt-SNE")
cbmc <- RunUMAP(cbmc, dims = 1:20)
DimPlot(cbmc, label = TRUE, reduction = "umap") + NoLegend()+ ggtitle("WT")

# Find the markers that define each cluster, and use these to annotate the
# clusters, we use max.cells.per.ident to speed up the process
cbmc.rna.markers <- FindAllMarkers(cbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cbmc.rna.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(cbmc, features = top10$gene) + NoLegend()+ggtitle("WT")

#feature.names = read.delim("./dataset/features.tsv.gz", header = FALSE,  stringsAsFactors = FALSE)
#for (i in 0:13) {
  cbmc.rna.markers %>% filter(cluster == i)
genemarker <- cbmc.rna.markers[cbmc.rna.markers$cluster == i,"gene"]
anno.name <- feature.names[feature.names$V1 %in% genemarker,]
anno.name <- anno.name[match(genemarker,  anno.name$V1),]
cbmc.rna.markers[cbmc.rna.markers$cluster == i,"gene"] <- anno.name$V2





# Now, we will visualize CD14 levels for RNA and protein By setting the default assay, we can
# visualize one or the other
DefaultAssay(cbmc) <- "ADT"
p1 <- FeaturePlot(cbmc, "CD45RA-TotalSeqB", cols = c("lightgrey", "darkgreen")) + ggtitle("CD45RA protein")
DefaultAssay(cbmc) <- "RNA"
p2 <- FeaturePlot(cbmc, "PTPRC") + ggtitle("CD45RA RNA")
p1 | p2
DefaultAssay(cbmc) <- "ADT"
p1 <- FeaturePlot(cbmc, "CD279-TotalSeqB", cols = c("lightgrey", "darkgreen")) + ggtitle("CD279 protein")
DefaultAssay(cbmc) <- "RNA"
p2 <- FeaturePlot(cbmc, "PDCD1") + ggtitle("CD279 RNA")
p1 | p2
DefaultAssay(cbmc) <- "ADT"
p1 <- FeaturePlot(cbmc, "CD8-TotalSeqB", cols = c("lightgrey", "darkgreen")) + ggtitle("CD8 protein")
DefaultAssay(cbmc) <- "RNA"
p2 <- FeaturePlot(cbmc, "CD8A") + ggtitle("CD8 RNA")
p1 | p2
DefaultAssay(cbmc) <- "ADT"
p1 <- FeaturePlot(cbmc, "CD62L-TotalSeqB", cols = c("lightgrey", "darkgreen")) + ggtitle("CD62L protein")
DefaultAssay(cbmc) <- "RNA"
p2 <- FeaturePlot(cbmc, "SELL") + ggtitle("CD62L RNA")
p1 | p2


#ident 
cbmc[["cluster.ident"]] <- cbmc[["HTO.ident"]] 
Idents(cbmc) <- cbmc[["HTO.ident"]] 

pdf("./plots/01UMAP_KO_HTOident.pdf", width = 5, height = 5)
DimPlot(cbmc, label = TRUE, reduction = "umap")
RidgePlot(cbmc, assay = "HTO", features = rownames(cbmc[["HTO"]]), ncol = 2)
dev.off()



cbmp.wt <- cbmc
rna.markers.wt <- cbmc.rna.markers
save(cbmp.wt, rna.markers.wt, file = "./seurat_wt.RData")
save(cbmp.wt, file = "./seurat_wt.rds")


cbmp.ko <- cbmc
rna.markers.ko <- cbmc.rna.markers
save(cbmp.ko, rna.markers.ko, file = "./seurat_ko.RData")
save(cbmp.ko, file = "./seurat_ko.rds")


###cell proportion 
load("./seurat_wt.rds")
load("./seurat_ko.rds")

cbmc <- cbmp.wt
cbmc <- cbmp.ko

DefaultAssay(cbmc) <- "ADT"
cart <- subset(cbmc, subset =  CART>0) 
cart <- cart[,which(cart[["HTO.ident"]]== "KO2")]

FeatureScatter(cart, feature1 = "CD4", feature2 = "CD8") 
cd4t <- subset(cart, subset =  CD4 >0.6 & CD8 < 1) # CD4 T
cd8t <- subset(cart, subset =  CD4 <1 & CD8 > 1.2) # CD8 T

DefaultAssay(cart) <- "RNA"
ncol(subset(cart, subset = ENSG00000175354 > 0.1)) / ncol(cart)#PTPN2 KO
FeatureScatter(cart, feature1 = "CD45RA", feature2 = "ENSG00000081237") #CD45RA ADT > 
FeatureScatter(cart, feature1 = "CD62L", feature2 = "ENSG00000188404") # CD62L ADT > 
FeatureScatter(cart, feature1 = "CD25", feature2 = "ENSG00000134460") #CD25 ADT > 
FeatureScatter(cart, feature1 = "CD366", feature2 = "ENSG00000135077") # CD366 ADT > 
FeatureScatter(cart, feature1 = "CD279", feature2 = "ENSG00000188389") # CD279 ADT >

FeatureScatter(cd8t, feature1 = "CD45RA", feature2 = "ENSG00000081237") #CD45RA ADT > 1.25
FeatureScatter(cd8t, feature1 = "CD62L", feature2 = "ENSG00000188404") # CD62L ENSG00000188404 > 0
FeatureScatter(cd8t, feature1 = "CD25", feature2 = "ENSG00000134460") #CD25 ENSG00000134460 > 0
FeatureScatter(cd8t, feature1 = "CD366", feature2 = "ENSG00000135077") # CD366 ENSG00000135077 > 0
FeatureScatter(cd8t, feature1 = "CD279", feature2 = "ENSG00000188389") # CD279 ENSG00000188389 >0

FeatureScatter(cd8t, feature1 = "ENSG00000188404", feature2 = "ENSG00000126353") # #TCF7
FeatureScatter(cd8t, feature1 = "CD8", feature2 = "ENSG00000175354") #PTPN2



DefaultAssay(cd8t) <- "RNA"
ncol(subset(cd8t, subset = adt_CD45RA <1.25 & ENSG00000188404 > 0.1)) / ncol(cd8t) #CD45RA-/lo CD62Lhi/+ central/memory T cells
ncol(subset(cd8t, subset = adt_CD45RA < 1.25 & ENSG00000188404 < 0.1)) / ncol(cd8t) #CD45RA-/lo CD62L-/lo effector/memory T cells
ncol(subset(cd8t, subset = adt_CD45RA < 1.25 & ENSG00000134460 >0.1)) / ncol(cd8t) #CD45RA-/lo CD25+/hi effector T cells
ncol(subset(cd8t, subset = adt_CD45RA >1.25 & ENSG00000188404 <0.1)) / ncol(cd8t) #CD45RA+/hi CD62L-/lo terminally differentiated effector/memory T cells
ncol(subset(cd8t, subset = adt_CD45RA > 1.25 & ENSG00000188404 > 0.1)) / ncol(cd8t) # CD45RA+/hi CD62L+/hi stem cell like memory T cells
ncol(subset(cd8t, subset = ENSG00000188389 > 0.1 & ENSG00000135077 < 0.1)) / ncol(cd8t)#PD-1 intermediate/hi Tim-3lo/- progenitor exhausted T cells
ncol(subset(cd8t, subset = ENSG00000188389 > 0.1 & ENSG00000135077 > 0.1)) / ncol(cd8t) #PD-1 +/hi Tim-3+/hi terminally exhausted T cells

ncol(subset(cd8t, subset = ENSG00000081059 > 0.1 & ENSG00000135077 > 0.1)) / ncol(cd8t)#TCF-1 hi/inter TIM-3+ =T PEX
ncol(subset(cd8t, subset = ENSG00000175354 > 0.1)) / ncol(cd8t)#PTPN2 KO

FeatureScatter(cd4t, feature1 = "CD45RA", feature2 = "ENSG00000081237") #CD45RA ADT > 1.25
DefaultAssay(cd4t) <- "RNA"
ncol(subset(cd4t, subset = adt_CD45RA <1.25 & ENSG00000188404 > 0.1)) / ncol(cd4t) #CD45RA-/lo CD62Lhi/+ central/memory T cells
ncol(subset(cd4t, subset = adt_CD45RA < 1.25 & ENSG00000188404 < 0.1)) / ncol(cd4t) #CD45RA-/lo CD62L-/lo effector/memory T cells
ncol(subset(cd4t, subset = adt_CD45RA < 1.25 & ENSG00000134460 >0.1)) / ncol(cd4t) #CD45RA-/lo CD25+/hi effector T cells
ncol(subset(cd4t, subset = adt_CD45RA >1.25 & ENSG00000188404 <0.1)) / ncol(cd4t) #CD45RA+/hi CD62L-/lo terminally differentiated effector/memory T cells
ncol(subset(cd4t, subset = adt_CD45RA > 1.25 & ENSG00000188404 > 0.1)) / ncol(cd4t) # CD45RA+/hi CD62L+/hi stem cell like memory T cells
ncol(subset(cd4t, subset = ENSG00000188389 > 0.1 & ENSG00000135077 < 0.1)) / ncol(cd4t)#PD-1 intermediate/hi Tim-3lo/- progenitor exhausted T cells
ncol(subset(cd4t, subset = ENSG00000188389 > 0.1 & ENSG00000135077 > 0.1)) / ncol(cd4t) #PD-1 +/hi Tim-3+/hi terminally exhausted T cells



###cell proportion strong positive CAR-T
load("./seurat_wt.rds")
load("./seurat_ko.rds")

cbmc <- cbmp.wt
cbmc <- cbmp.ko

DefaultAssay(cbmc) <- "ADT"
FeatureScatter(cbmc, feature1 = "CD8", feature2 = "CART") 
Idents(cbmc) <- cbmc[["seurat_clusters"]]
RidgePlot(cbmc, features = "CART")
RidgePlot(cbmc, features = "CD8")

cart <- subset(cbmc, subset =  CART>0.2) 
cart <- cart[,which(cart[["HTO.ident"]]== "KO2")]

FeatureScatter(cart, feature1 = "CD4", feature2 = "CD8") 
cd4t <- subset(cart, subset =  CD4 >0.75 & CD8 < 1) # CD4 T
cd8t <- subset(cart, subset =  CD4 <0.5 & CD8 > 1.2) # CD8 T

DefaultAssay(cd8t) <- "RNA"
FeatureScatter(cd8t, feature1 = "adt_CD45RA", feature2 = "SELL")
ncol(subset(cd8t, subset = adt_CD45RA <1.25 & SELL > 0.1)) / ncol(cd8t) #CD45RA-/lo CD62Lhi/+ central/memory T cells
ncol(subset(cd8t, subset = adt_CD45RA < 1.25 & SELL < 0.1)) / ncol(cd8t) #CD45RA-/lo CD62L-/lo effector/memory T cells
ncol(subset(cd8t, subset = adt_CD45RA < 1.25 & IL2RA >0.1)) / ncol(cd8t) #CD45RA-/lo CD25+/hi effector T cells
ncol(subset(cd8t, subset = adt_CD45RA >1.25 & SELL <0.1)) / ncol(cd8t) #CD45RA+/hi CD62L-/lo terminally differentiated effector/memory T cells
ncol(subset(cd8t, subset = adt_CD45RA > 1.25 & SELL > 0.1)) / ncol(cd8t) # CD45RA+/hi CD62L+/hi stem cell like memory T cells
ncol(subset(cd8t, subset = PDCD1 > 0.1 & HAVCR2 < 0.1)) / ncol(cd8t)#PD-1 intermediate/hi Tim-3lo/- progenitor exhausted T cells
ncol(subset(cd8t, subset = PDCD1 > 0.1 & HAVCR2 > 0.1)) / ncol(cd8t) #PD-1 +/hi Tim-3+/hi terminally exhausted T cells


FeatureScatter(cd4t, feature1 = "CD45RA", feature2 = "SELL") #CD45RA ADT > 1.25
DefaultAssay(cd4t) <- "RNA"
ncol(subset(cd4t, subset = adt_CD45RA <1.25 & SELL > 0.1)) / ncol(cd4t) #CD45RA-/lo CD62Lhi/+ central/memory T cells
ncol(subset(cd4t, subset = adt_CD45RA < 1.25 & SELL < 0.1)) / ncol(cd4t) #CD45RA-/lo CD62L-/lo effector/memory T cells
ncol(subset(cd4t, subset = adt_CD45RA < 1.25 & IL2RA >0.1)) / ncol(cd4t) #CD45RA-/lo CD25+/hi effector T cells
ncol(subset(cd4t, subset = adt_CD45RA >1.25 & SELL <0.1)) / ncol(cd4t) #CD45RA+/hi CD62L-/lo terminally differentiated effector/memory T cells
ncol(subset(cd4t, subset = adt_CD45RA > 1.25 & SELL > 0.1)) / ncol(cd4t) # CD45RA+/hi CD62L+/hi stem cell like memory T cells
ncol(subset(cd4t, subset = PDCD1 > 0.1 & HAVCR2 < 0.1)) / ncol(cd4t)#PD-1 intermediate/hi Tim-3lo/- progenitor exhausted T cells
ncol(subset(cd4t, subset = PDCD1 > 0.1 & HAVCR2 > 0.1)) / ncol(cd4t) #PD-1 +/hi Tim-3+/hi terminally exhausted T cells




####2. Rename by ADT and RNA types ####
load(file = "seurat_all.rds")

DefaultAssay(cbmc) <- "ADT"
FeatureScatter(cbmc, feature1 = "CD8", feature2 = "CART") 
Idents(cbmc) <- cbmc[["seurat_clusters"]]

cart <- subset(cbmc, subset =  CART>0) 

FeatureScatter(cart, feature1 = "CD4", feature2 = "CD8") 
cd4t <- subset(cart, subset =  CD4 >0.75 & CD8 < 1) # CD4 T
cd8t <- subset(cart, subset =  CD4 <0.5 & CD8 > 1.2) # CD8 T


cd8t@assays$ADT@counts[c("CD45RA"),]
cd8t@assays$RNA@counts[c("SELL","IL2RA","PDCD1", "HAVCR2"),]
adt_rna <- as.matrix(rbind(cd8t@assays$ADT@counts[c("CD45RA"),], cd8t@assays$RNA@counts[c("SELL","IL2RA","PDCD1", "HAVCR2"),]))
rownames(adt_rna)[1] <- "CD45RA"
cd8t[["Supervised"]] <- CreateAssayObject(counts = adt_rna)

# Validate that the object now contains multiple assays
Assays(cd8t)
cd8t <- NormalizeData(cd8t, assay = "Supervised", normalization.method = "CLR", margin = 2)
DefaultAssay(cd8t) <- "Supervised"
cd8t <- FindVariableFeatures(cd8t)
cd8t <- ScaleData(cd8t)
cd8t <- RunPCA(cd8t, verbose = TRUE)
ElbowPlot(cd8t, ndims = 10)
cd8t <- FindNeighbors(cd8t, dims = 1:4)
cd8t <- FindClusters(cd8t, resolution = 0.01) #0.8 for 14 subclusters
Idents(cd8t) <- "Supervised_snn_res.0.01"
cd8t <- RunUMAP(cd8t, dims = 1:4)
DimPlot(cd8t, label = TRUE, reduction = "umap") + NoLegend()+ ggtitle("")


DefaultAssay(cd8t) <- "RNA"
cd8t@meta.data$rename <- "T"
idx <- which(colnames(cd8t) %in% colnames(subset(cd8t, subset = adt_CD45RA <1.25 & SELL > 0.1)))
cd8t@meta.data$rename[idx] <- "Tcm"
idx <- which(colnames(cd8t) %in% colnames(subset(cd8t, subset = adt_CD45RA < 1.25 & SELL < 0.1)))
cd8t@meta.data$rename[idx] <- "Tem"
idx <- which(colnames(cd8t) %in% colnames(subset(cd8t, subset = adt_CD45RA < 1.25 & IL2RA >0.1)))
cd8t@meta.data$rename[idx] <- "Teff"
idx <- which(colnames(cd8t) %in% colnames(subset(cd8t, subset = adt_CD45RA >1.25 & SELL <0.1)))
cd8t@meta.data$rename[idx] <- "Ttdem"
idx <- which(colnames(cd8t) %in% colnames(subset(cd8t, subset = adt_CD45RA > 1.25 & SELL > 0.1)))
cd8t@meta.data$rename[idx] <- "Tscm"
idx <- which(colnames(cd8t) %in% colnames(subset(cd8t, subset = PDCD1 > 0.1 & HAVCR2 < 0.1)))
cd8t@meta.data$rename[idx] <- "Tpex"
idx <- which(colnames(cd8t) %in% colnames(subset(cd8t, subset = PDCD1 > 0.1 & HAVCR2 > 0.1)))
cd8t@meta.data$rename[idx] <- "Ttex"
cd8t@meta.data$rename
Idents(cd8t) <- "rename"
png("./plots/01DimPlot_rename_CD8_ADTRNAcluster.png", width = 12, height = 10, units = "cm", res = 300)
DimPlot(cd8t, label = F, reduction = "umap", cols = pal.type) + ggtitle("CD8")
dev.off()

DefaultAssay(cd8t) <- "RNA"
######Target Genes######
stat3 <- c("FOSL2", "NTRK3", "IL16", "DUSP10", "CASP1", "FNDC9", "PSD3", "FLT1", "TSHZ2", "RORA", "KCNK1", "NAPSA","FNBP1L", 
           "CALCB", "IL1R1", "COL6A3", "CCR6", "IL24", "HSD11B1", "IFNL1", "IL23R", "MAF", "PALLD", "HOPX", "IFIT2", "GPR87",
           "BCAR3", "IL1R2","DHRS9", "IGFBP6", "PRG4","CXCR5", "GZMB", "IL12RB2", "GBP4","TNFSF13B","TNFSF11","CSF2", "AFF3",
           "TNFSF4", "TMEM200A", "PTGER2", "SLC12A8", "B3GALNT1", "LAIR2", "PHLDA1", "NR4A2",  "QPCT")
stat1 <- c("IDO1", "CXCL10", "IRF1", "GBP5", "ITK", "ICAM1", "APOL6", "TMEM140", "PARP9", "GBP2", "CTSS", "CXCL9", "TRIM21", 
           "GBP1", "CCL2", "DTX3L", "PSMB9", "CCDC68", "APOL1", "PDCD1LG2")
stat5<- c("IL2RA", "CISH", "CDK6", "LRRC32", "MAP3K8", "LTA", "OSM", "LIF", "NABP1", "GZMB")
immune_response <- c("AIM2", "APOL1", "C1S", "C3", "C4A", "CCL2", "CIITA", "CTSS",  "CXCL10", "CXCL9", "GBP1", "GBP2", "GBP5", "GCH1", "HLA-E",
                     "ICAM1", "IFI35", "IL7", "IL4R", "LYN", "ORAI1", "PDCD1LG2",
                     "PSMB8", "PSMB9", "RNF19B", "TAP1", "TAP2", "BCL6", "C1S", "C3", "C4A", "F2RL1", "FYN", "ICAM1", "IDO1",
                     "IL4R", "IL7", "LYN", "PDCD1LG2", "PVR", "TGFB2")
memory_regulation <- c("BCL6", "CD46", "FGL2", "IL23A", "PCK1","TNFSF4")
Idents(cd8t) <- cd8t@meta.data$HTO.ident
cd8t <- RenameIdents(object = cd8t, 'KO1' = "KO")
cd8t <- RenameIdents(object = cd8t, 'KO2' = "KO")
cd8t <- RenameIdents(object = cd8t, 'WT1' = "WT")
cd8t <- RenameIdents(object = cd8t, 'WT2' = "WT")
cd8t@meta.data$group <- Idents(cd8t)


library(ggpubr)
vp_case1 <- function(gene_signature, file_name, test_sign){
  plot_case1 <- function(signature, y_max = NULL){
    VlnPlot(tcell, features = signature,
            pt.size = 0.1, 
            #group.by = "Response", 
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) + stat_compare_means(comparisons = test_sign, label = "p.signif")
  }
  plot_list <- list()
  y_max_list <- list()
  for (gene in gene_signature) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]]) # get the max no. for each gene
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 1) )
  }
  cowplot::plot_grid(plotlist = plot_list)
  file_name <- paste0(file_name, "_r.png")
  ggsave(file_name, width = 14, height = 8)
}
vp_case2 <- function(gene_signature, file_name, test_sign){
  plot_case1 <- function(signature, y_max = NULL){
    VlnPlot(tcell, features = signature,
            pt.size = 0.1, 
            #group.by = "Response", 
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) + stat_compare_means(comparisons = test_sign, label = "p.signif")
  }
  plot_list <- list()
  y_max_list <- list()
  for (gene in gene_signature) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]]) # get the max no. for each gene
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 1) )
  }
  cowplot::plot_grid(plotlist = plot_list)
  file_name <- paste0(file_name, "_r.png")
  ggsave(file_name, width = 16, height = 16)
}
comparisons <- list(c("WT", "KO"))

Idents(cd8t) <- "rename"
tcell <- cd8t[,which(Idents(cd8t) == "Tcm")]
Idents(tcell) <- "group"
gene_sig <- stat1
vp_case1(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_Tcm_stat1", test_sign = comparisons)
gene_sig <- stat3
vp_case2(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_Tcm_stat3", test_sign = comparisons)
gene_sig <- stat5
vp_case1(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_Tcm_stat5", test_sign = comparisons)

Idents(cd8t) <- "rename"
tcell <- cd8t[,which(Idents(cd8t) == "Tem")]
Idents(tcell) <- "group"
gene_sig <- stat1
vp_case1(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_Tem_stat1", test_sign = comparisons)
gene_sig <- stat3
vp_case2(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_Tem_stat3", test_sign = comparisons)
gene_sig <- stat5
vp_case1(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_Tem_stat5", test_sign = comparisons)

Idents(cd8t) <- "rename"
tcell <- cd8t[,which(Idents(cd8t) == "Teff")]
Idents(tcell) <- "group"
gene_sig <- stat1
vp_case1(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_Teff_stat1", test_sign = comparisons)
gene_sig <- stat3
vp_case2(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_Teff_stat3", test_sign = comparisons)
gene_sig <- stat5
vp_case1(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_Teff_stat5", test_sign = comparisons)

Idents(cd8t) <- "rename"
tcell <- cd8t[,which(Idents(cd8t) == "Ttdem")]
Idents(tcell) <- "group"
gene_sig <- stat1
vp_case1(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_Ttdem_stat1", test_sign = comparisons)
gene_sig <- stat3
vp_case2(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_Ttdem_stat3", test_sign = comparisons)
gene_sig <- stat5
vp_case1(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_Ttdem_stat5", test_sign = comparisons)

Idents(cd8t) <- "rename"
tcell <- cd8t[,which(Idents(cd8t) == "Tscm")]
Idents(tcell) <- "group"
gene_sig <- stat1
vp_case1(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_Tscm_stat1", test_sign = comparisons)
gene_sig <- stat3
vp_case2(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_Tscm_stat3", test_sign = comparisons)
gene_sig <- stat5
vp_case1(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_Tscm_stat5", test_sign = comparisons)


Idents(cd8t) <- "rename"
tcell <- cd8t[,which(Idents(cd8t) == "Tcm")]
Idents(tcell) <- "group"
png("./plots/01VlnPlot_Tcm_CD8_STAT3.png", width = 50, height = 70, units = "cm", res = 300)
VlnPlot(tcell, features = stat3, ncol = 5)
dev.off()
png("./plots/01VlnPlot_Tcm_CD8_STAT1.png", width = 50, height = 30, units = "cm", res = 300)
VlnPlot(tcell, features = stat1, ncol = 5)
dev.off()
png("./plots/01VlnPlot_Tcm_CD8_STAT5.png", width = 50, height = 20, units = "cm", res = 300)
VlnPlot(tcell, features = stat5, ncol = 5) + stat_compare_means(comparisons = as.character( tcell@meta.data$group), label = "p.signif")
dev.off()
png("./plots/01VlnPlot_Tcm_CD8_immune_response.png", width = 50, height = 70, units = "cm", res = 300)
VlnPlot(tcell, features = immune_response, ncol = 5)
dev.off()

tcell <- cd8t[,which(Idents(cd8t) == "Tem")]
Idents(tcell) <- "group"
png("./plots/01VlnPlot_Tem_CD8_STAT3.png", width = 50, height = 70, units = "cm", res = 300)
VlnPlot(tcell, features = stat3, ncol = 5)
dev.off()
png("./plots/01VlnPlot_Tem_CD8_STAT1.png", width = 50, height = 30, units = "cm", res = 300)
VlnPlot(tcell, features = stat1, ncol = 5)
dev.off()
png("./plots/01VlnPlot_Tem_CD8_STAT5.png", width = 50, height = 20, units = "cm", res = 300)
VlnPlot(tcell, features = stat5, ncol = 5)
dev.off()
png("./plots/01VlnPlot_Tem_CD8_immune_response.png", width = 50, height = 70, units = "cm", res = 300)
VlnPlot(tcell, features = immune_response, ncol = 5)
dev.off()


tcell <- cd8t[,which(Idents(cd8t) == "Teff")]
Idents(tcell) <- "group"
png("./plots/01VlnPlot_Teff_CD8_STAT3.png", width = 50, height = 70, units = "cm", res = 300)
VlnPlot(tcell, features = stat3, ncol = 5)
dev.off()
png("./plots/01VlnPlot_Teff_CD8_STAT1.png", width = 50, height = 30, units = "cm", res = 300)
VlnPlot(tcell, features = stat1, ncol = 5)
dev.off()
png("./plots/01VlnPlot_Teff_CD8_STAT5.png", width = 50, height = 20, units = "cm", res = 300)
VlnPlot(tcell, features = stat5, ncol = 5)
dev.off()
png("./plots/01VlnPlot_Teff_CD8_immune_response.png", width = 50, height = 70, units = "cm", res = 300)
VlnPlot(tcell, features = immune_response, ncol = 5)
dev.off()

tcell <- cd8t[,which(Idents(cd8t) == "Ttdem")]
Idents(tcell) <- "group"
png("./plots/01VlnPlot_Ttdem_CD8_STAT3.png", width = 50, height = 70, units = "cm", res = 300)
VlnPlot(tcell, features = stat3, ncol = 5)
dev.off()
png("./plots/01VlnPlot_Ttdem_CD8_STAT1.png", width = 50, height = 30, units = "cm", res = 300)
VlnPlot(tcell, features = stat1, ncol = 5)
dev.off()
png("./plots/01VlnPlot_Ttdem_CD8_STAT5.png", width = 50, height = 20, units = "cm", res = 300)
VlnPlot(tcell, features = stat5, ncol = 5)
dev.off()
png("./plots/01VlnPlot_Ttdem_CD8_immune_response.png", width = 50, height = 70, units = "cm", res = 300)
VlnPlot(tcell, features = immune_response, ncol = 5)
dev.off()

tcell <- cd8t[,which(Idents(cd8t) == "Tscm")]
Idents(tcell) <- "group"
png("./plots/01VlnPlot_Tscm_CD8_STAT3.png", width = 50, height = 70, units = "cm", res = 300)
VlnPlot(tcell, features = stat3, ncol = 5)
dev.off()
png("./plots/01VlnPlot_Tscm_CD8_STAT1.png", width = 50, height = 30, units = "cm", res = 300)
VlnPlot(tcell, features = stat1, ncol = 5)
dev.off()
png("./plots/01VlnPlot_Tscm_CD8_STAT5.png", width = 50, height = 20, units = "cm", res = 300)
VlnPlot(tcell, features = stat5, ncol = 5)
dev.off()
png("./plots/01VlnPlot_Tscm_CD8_immune_response.png", width = 50, height = 70, units = "cm", res = 300)
VlnPlot(tcell, features = immune_response, ncol = 5)
dev.off()

tcell <- cd8t[,which(Idents(cd8t) == "Tpex")]
Idents(tcell) <- "group"
png("./plots/01VlnPlot_Tpex_CD8_STAT3.png", width = 50, height = 70, units = "cm", res = 300)
VlnPlot(tcell, features = stat3, ncol = 5)
dev.off()
png("./plots/01VlnPlot_Tpex_CD8_STAT1.png", width = 50, height = 30, units = "cm", res = 300)
VlnPlot(tcell, features = stat1, ncol = 5)
dev.off()
png("./plots/01VlnPlot_Tpex_CD8_STAT5.png", width = 50, height = 20, units = "cm", res = 300)
VlnPlot(tcell, features = stat5, ncol = 5)
dev.off()
png("./plots/01VlnPlot_Tpex_CD8_immune_response.png", width = 50, height = 70, units = "cm", res = 300)
VlnPlot(tcell, features = immune_response, ncol = 5)
dev.off()

tcell <- cd8t[,which(Idents(cd8t) == "Ttex")]
Idents(tcell) <- "group"
png("./plots/01VlnPlot_Ttex_CD8_STAT3.png", width = 50, height = 70, units = "cm", res = 300)
VlnPlot(tcell, features = stat3, ncol = 5)
dev.off()
png("./plots/01VlnPlot_Ttex_CD8_STAT1.png", width = 50, height = 30, units = "cm", res = 300)
VlnPlot(tcell, features = stat1, ncol = 5)
dev.off()
png("./plots/01VlnPlot_Ttex_CD8_STAT5.png", width = 50, height = 20, units = "cm", res = 300)
VlnPlot(tcell, features = stat5, ncol = 5)
dev.off()
png("./plots/01VlnPlot_Ttex_CD8_immune_response.png", width = 50, height = 70, units = "cm", res = 300)
VlnPlot(tcell, features = immune_response, ncol = 5)
dev.off()



DefaultAssay(cd4t) <- "RNA"
DimPlot(cd4t, label = TRUE, reduction = "umap") + NoLegend()
Idents(cd4t)
cd4t@meta.data$rename <- "T"
idx <- which(colnames(cd4t) %in% colnames(subset(cd4t, subset = adt_CD45RA <1.25 & SELL > 0.1)))
cd4t@meta.data$rename[idx] <- "Tcm"
idx <- which(colnames(cd4t) %in% colnames(subset(cd4t, subset = adt_CD45RA < 1.25 & SELL < 0.1)))
cd4t@meta.data$rename[idx] <- "Tem"
idx <- which(colnames(cd4t) %in% colnames(subset(cd4t, subset = adt_CD45RA < 1.25 & IL2RA >0.1)))
cd4t@meta.data$rename[idx] <- "Teff"
idx <- which(colnames(cd4t) %in% colnames(subset(cd4t, subset = adt_CD45RA >1.25 & SELL <0.1)))
cd4t@meta.data$rename[idx] <- "Ttdem"
idx <- which(colnames(cd4t) %in% colnames(subset(cd4t, subset = adt_CD45RA > 1.25 & SELL > 0.1)))
cd4t@meta.data$rename[idx] <- "Tscm"
idx <- which(colnames(cd4t) %in% colnames(subset(cd4t, subset = PDCD1 > 0.1 & HAVCR2 < 0.1 )))
cd4t@meta.data$rename[idx] <- "Tpex"
idx <- which(colnames(cd4t) %in% colnames(subset(cd4t, subset = PDCD1 > 0.1 & HAVCR2 > 0.1)))
cd4t@meta.data$rename[idx] <- "Ttex"
cd4t@meta.data$rename
Idents(cd4t) <- "rename"
png("./plots/01DimPlot_rename_CD4_WT.png", width = 12, height = 10, units = "cm", res = 300)
DimPlot(cd4t, label = F, reduction = "umap", cols = pal.type) + ggtitle("WT CD4")
dev.off()



cbmc <- cbmp.ko

DefaultAssay(cbmc) <- "ADT"
FeatureScatter(cbmc, feature1 = "CD8", feature2 = "CART") 
Idents(cbmc) <- cbmc[["seurat_clusters"]]
RidgePlot(cbmc, features = "CART")
RidgePlot(cbmc, features = "CD8")

cart <- subset(cbmc, subset =  CART>0) 

FeatureScatter(cart, feature1 = "CD4", feature2 = "CD8") 
cd4t <- subset(cart, subset =  CD4 >0.75 & CD8 < 1) # CD4 T
cd8t <- subset(cart, subset =  CD4 <0.5 & CD8 > 1.2) # CD8 T

DefaultAssay(cd8t) <- "RNA"
DimPlot(cd8t, label = TRUE, reduction = "umap") + NoLegend()
Idents(cd8t)

cd8t@meta.data$rename <- "T"
idx <- which(colnames(cd8t) %in% colnames(subset(cd8t, subset = adt_CD45RA <1.25 & SELL > 0.1)))
cd8t@meta.data$rename[idx] <- "Tcm"
idx <- which(colnames(cd8t) %in% colnames(subset(cd8t, subset = adt_CD45RA < 1.25 & SELL < 0.1)))
cd8t@meta.data$rename[idx] <- "Tem"
idx <- which(colnames(cd8t) %in% colnames(subset(cd8t, subset = adt_CD45RA < 1.25 & IL2RA >0.1)))
cd8t@meta.data$rename[idx] <- "Teff"
idx <- which(colnames(cd8t) %in% colnames(subset(cd8t, subset = adt_CD45RA >1.25 & SELL <0.1)))
cd8t@meta.data$rename[idx] <- "Ttdem"
idx <- which(colnames(cd8t) %in% colnames(subset(cd8t, subset = adt_CD45RA > 1.25 & SELL > 0.1)))
cd8t@meta.data$rename[idx] <- "Tscm"
idx <- which(colnames(cd8t) %in% colnames(subset(cd8t, subset = PDCD1 > 0.1 & HAVCR2 < 0.1)))
cd8t@meta.data$rename[idx] <- "Tpex"
idx <- which(colnames(cd8t) %in% colnames(subset(cd8t, subset = PDCD1 > 0.1 & HAVCR2 > 0.1)))
cd8t@meta.data$rename[idx] <- "Ttex"
cd8t@meta.data$rename
Idents(cd8t) <- "rename"
png("./plots/01DimPlot_rename_CD8_KO.png", width = 12, height = 10, units = "cm", res = 300)
DimPlot(cd8t, label = F, reduction = "umap", cols = pal.type) + ggtitle("KO CD8")
dev.off()

DefaultAssay(cd4t) <- "RNA"
DimPlot(cd4t, label = TRUE, reduction = "umap") + NoLegend()
Idents(cd4t)
cd4t@meta.data$rename <- "T"
idx <- which(colnames(cd4t) %in% colnames(subset(cd4t, subset = adt_CD45RA <1.25 & SELL > 0.1)))
cd4t@meta.data$rename[idx] <- "Tcm"
idx <- which(colnames(cd4t) %in% colnames(subset(cd4t, subset = adt_CD45RA < 1.25 & SELL < 0.1)))
cd4t@meta.data$rename[idx] <- "Tem"
idx <- which(colnames(cd4t) %in% colnames(subset(cd4t, subset = adt_CD45RA < 1.25 & IL2RA >0.1)))
cd4t@meta.data$rename[idx] <- "Teff"
idx <- which(colnames(cd4t) %in% colnames(subset(cd4t, subset = adt_CD45RA >1.25 & SELL <0.1)))
cd4t@meta.data$rename[idx] <- "Ttdem"
idx <- which(colnames(cd4t) %in% colnames(subset(cd4t, subset = adt_CD45RA > 1.25 & SELL > 0.1)))
cd4t@meta.data$rename[idx] <- "Tscm"
idx <- which(colnames(cd4t) %in% colnames(subset(cd4t, subset = PDCD1 > 0.1 & HAVCR2 < 0.1)))
cd4t@meta.data$rename[idx] <- "Tpex"
idx <- which(colnames(cd4t) %in% colnames(subset(cd4t, subset = PDCD1 > 0.1 & HAVCR2 > 0.1)))
cd4t@meta.data$rename[idx] <- "Ttex"
cd4t@meta.data$rename
Idents(cd4t) <- "rename"
png("./plots/01DimPlot_rename_CD4_KO.png", width = 12, height = 10, units = "cm", res = 300)
DimPlot(cd4t, label = F, reduction = "umap", cols = pal.type) + ggtitle("KO CD4")
dev.off()








####2. Rename by ADT types ####
load(file = "./seurat_all.rds")
DefaultAssay(cbmc) <- "ADT"
cart <- subset(cbmc, subset =  CART>0) 

FeatureScatter(cart, feature1 = "CD4", feature2 = "CD8") 
cd4t <- subset(cart, subset =  CD4 >0.6 & CD8 < 1) # CD4 T
cd8t <- subset(cart, subset =  CD4 <0.5 & CD8 > 1.2) # CD8 T

#ADT cluster
DefaultAssay(cd8t) <- "ADT"
cd8t <- NormalizeData(cd8t)
cd8t <- FindVariableFeatures(cd8t)
cd8t <- ScaleData(cd8t)
cd8t <- RunPCA(cd8t, verbose = TRUE)
ElbowPlot(cd8t, ndims = 20)
cd8t <- FindNeighbors(cd8t, dims = 1:7)
cd8t <- FindClusters(cd8t, resolution = 0.5) #0.8 for 14 subclusters
cd8t <- RunUMAP(cd8t, dims = 1:7)
Idents(cd8t) <- "ADT_snn_res.0.5"
png("./plots/01Umap_CD8_ADTcluster.png", width = 10, height = 10, res = 300, units = "cm")
DimPlot(cd8t, label = TRUE, reduction = "umap") + NoLegend() + ggtitle("CD8")
dev.off()

png("./plots/01FeaturePlot_ADTcluster_CD8.png", width = 40, height = 10, units = "cm", res = 300)
p1 <- FeaturePlot(cd8t, "CD45RA") + ggtitle("CD45RA")
p2 <- FeaturePlot(cd8t, "CD62L") + ggtitle("CD62L")
p3 <- FeaturePlot(cd8t, "CD25") + ggtitle("CD25")
p4 <- FeaturePlot(cd8t, "CD366") + ggtitle("CD366")
p5 <- FeaturePlot(cd8t, "CD279") + ggtitle("CD279")
p1 | p2 |p3 | p4 |p5
dev.off()


#RNA cluster
DefaultAssay(cd8t) <- "RNA"
cd8t <- NormalizeData(cd8t)
cd8t <- FindVariableFeatures(cd8t)
cd8t <- ScaleData(cd8t)
cd8t <- RunPCA(cd8t, verbose = TRUE)
ElbowPlot(cd8t, ndims = 20)
cd8t <- FindNeighbors(cd8t, dims = 1:20)
cd8t <- FindClusters(cd8t, resolution = 0.35) #0.8 for 14 subclusters
cd8t <- RunUMAP(cd8t, dims = 1:20)
Idents(cd8t) <- "RNA_snn_res.0.35"
png("./plots/01Umap_CD8_ADT.png", width = 10, height = 10, res = 300, units = "cm")
DimPlot(cd8t, label = TRUE, reduction = "umap") + NoLegend() + ggtitle("CD8")
dev.off()

Idents(cd8t) <- "HTO.ident"
cd8t <- RenameIdents(object = cd8t, 'KO1' = "KO")
cd8t <- RenameIdents(object = cd8t, 'KO2' = "KO")
cd8t <- RenameIdents(object = cd8t, 'WT1' = "WT")
cd8t <- RenameIdents(object = cd8t, 'WT2' = "WT")
cd8t@meta.data$group <- Idents(cd8t)
Idents(cd8t) <- "group"
png("./plots/01Umap_CD8_ADT_WTKO.png", width = 10, height = 10, res = 300, units = "cm")
DimPlot(cd8t, label = F, reduction = "umap", cols = pal.group) + ggtitle("CD8")
dev.off()


DefaultAssay(cd4t) <- "RNA"
cd4t <- NormalizeData(cd4t)
cd4t <- FindVariableFeatures(cd4t)
cd4t <- ScaleData(cd4t)
cd4t <- RunPCA(cd4t, verbose = TRUE)
ElbowPlot(cd4t, ndims = 20)
cd4t <- FindNeighbors(cd4t, dims = 1:20)
cd4t <- FindClusters(cd4t, resolution = 0.35) #0.8 for 14 subclusters
cd4t <- RunUMAP(cd4t, dims = 1:20)
Idents(cd4t) <- "RNA_snn_res.0.35"
png("./plots/01Umap_CD4_ADT.png", width = 10, height = 10, res = 300, units = "cm")
DimPlot(cd4t, label = TRUE, reduction = "umap") + NoLegend() + ggtitle("CD4")
dev.off()

Idents(cd4t) <- "HTO.ident"
cd4t <- RenameIdents(object = cd4t, 'KO1' = "KO")
cd4t <- RenameIdents(object = cd4t, 'KO2' = "KO")
cd4t <- RenameIdents(object = cd4t, 'WT1' = "WT")
cd4t <- RenameIdents(object = cd4t, 'WT2' = "WT")
cd4t@meta.data$group <- Idents(cd4t)
Idents(cd4t) <- "group"
png("./plots/01Umap_CD4_ADT_WTKO.png", width = 10, height = 10, res = 300, units = "cm")
DimPlot(cd4t, label = F, reduction = "umap", cols = pal.group) + ggtitle("CD4")
dev.off()



#####3. cluster by RNA #####
#load(file = "./seurat_wt.RData")
#load(file = "./seurat_ko.RData")
load(file = "seurat_all.rds")
library(RColorBrewer)
pal.group <- c("red", "black")
names(pal.group) <- c("WT", "KO")

#CD8
DefaultAssay(cbmc) <- "ADT"
cart <- subset(cbmc, subset =  CART>0) 
FeatureScatter(cart, feature1 = "CD4", feature2 = "CD8") + geom_hline(yintercept=1.2, linetype="dashed", color = "black")+ geom_vline(xintercept=1, linetype="dashed", color = "black")
cd8t <- subset(cart, subset =  CD4 <1 & CD8 > 1.2) # CD8 T
Idents(cd8t) <- cd8t@meta.data$HTO.ident
cd8t <- RenameIdents(object = cd8t, 'KO1' = "KO")
cd8t <- RenameIdents(object = cd8t, 'KO2' = "KO")
cd8t <- RenameIdents(object = cd8t, 'WT1' = "WT")
cd8t <- RenameIdents(object = cd8t, 'WT2' = "WT")
cd8t@meta.data$group <- Idents(cd8t)
DefaultAssay(cd8t) <- "RNA"
cd8t <- FindVariableFeatures(cd8t)
cd8t <- ScaleData(cd8t)
cd8t <- RunPCA(cd8t, verbose = TRUE)
ElbowPlot(cd8t, ndims = 15)
cd8t <- FindNeighbors(cd8t, dims = 1:10)
cd8t <- FindClusters(cd8t, resolution = 0.3) #5 clusters
#cd8t <- FindClusters(cd8t, resolution = 0.25) # 4 clusters

cd8t <- RunUMAP(cd8t, dims = 1:10)
png("./plots/01DimPlot_CD8_RNA.png", width = 10, height = 8, res = 300, units = "cm")
Idents(cd8t) <- "RNA_snn_res.0.3"
DimPlot(cd8t, label = TRUE, reduction = "umap") + NoLegend()+ ggtitle("CD8")
dev.off()
png("./plots/01DimPlot_CD8_RNA_WTKO.png", width = 10, height = 8, res = 300, units = "cm")
Idents(cd8t) <- "group"
DimPlot(cd8t, label = F, reduction = "umap", cols = pal.group) + ggtitle("CD8") 
dev.off()


save(cd8t, file = "./seurat_cd8t.rds")

load("./seurat_cd8t.rds")

#CD4
DefaultAssay(cbmc) <- "ADT"
cart <- subset(cbmc, subset =  CART>0) 
FeatureScatter(cart, feature1 = "CD4", feature2 = "CD8") 
cd4t <- subset(cart, subset =  CD4 >0.61 & CD8 < 1) # CD8 T
Idents(cd4t) <- cd4t@meta.data$HTO.ident
cd4t <- RenameIdents(object = cd4t, 'KO1' = "KO")
cd4t <- RenameIdents(object = cd4t, 'KO2' = "KO")
cd4t <- RenameIdents(object = cd4t, 'WT1' = "WT")
cd4t <- RenameIdents(object = cd4t, 'WT2' = "WT")
cd4t@meta.data$group <- Idents(cd4t)
DefaultAssay(cd4t) <- "RNA"
cd4t <- FindVariableFeatures(cd4t)
cd4t <- ScaleData(cd4t)
cd4t <- RunPCA(cd4t, verbose = TRUE)
ElbowPlot(cd4t, ndims = 20)
cd4t <- FindNeighbors(cd4t, dims = 1:10)
cd4t <- FindClusters(cd4t, resolution = 0.3) #0.3 for 6 subclusters

cd4t <- RunUMAP(cd4t, dims = 1:10)
png("./plots/01DimPlot_CD4_RNA.png", width = 10, height = 8, res = 300, units = "cm")
Idents(cd4t) <- "RNA_snn_res.0.3"
DimPlot(cd4t, label = TRUE, reduction = "umap") + NoLegend()+ ggtitle("CD4")
dev.off()
png("./plots/01DimPlot_CD4_RNA_WTKO.png", width = 10, height = 8, res = 300, units = "cm")
Idents(cd4t) <- "group"
DimPlot(cd4t, label = F, reduction = "umap", cols = pal.group) + ggtitle("CD4") 
dev.off()



######Find markers for subtypes##### 
cbmc.rna.markers <- FindAllMarkers(cd8t, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cbmc.rna.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(cd8t, features = top10$gene) + NoLegend()+ggtitle("CD8")
write.csv(cbmc.rna.markers, file = "./cd8_markers.csv")

cbmc.rna.markers <- FindAllMarkers(cd4t, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cbmc.rna.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(cd4t, features = top10$gene) + NoLegend()+ggtitle("CD4")
write.csv(cbmc.rna.markers, file = "./cd4_markers.csv")


png("./plots/01FeaturePlot_TF_CD8.png", width = 40, height = 10, units = "cm", res = 300)
p1 <- FeaturePlot(cd8t, "TCF7") + ggtitle("TCF7")
p2 <- FeaturePlot(cd8t, "FOXP1") + ggtitle("FOXP1")
p3 <- FeaturePlot(cd8t, "NFATC1") + ggtitle("NFATC1")
p4 <- FeaturePlot(cd8t, "FOXP3") + ggtitle("FOXP3")
p1 | p2 |p3 | p4
dev.off()
png("./plots/01FeaturePlot_TF_CD4.png", width = 40, height = 10, units = "cm", res = 300)
p1 <- FeaturePlot(cd4t, "TCF7") + ggtitle("TCF7")
p2 <- FeaturePlot(cd4t, "FOXP1") + ggtitle("FOXP1")
p3 <- FeaturePlot(cd4t, "NFATC1") + ggtitle("NFATC1")
p4 <- FeaturePlot(cd4t, "FOXP3") + ggtitle("FOXP3")
p1 | p2 |p3 | p4
dev.off()



RidgePlot(cd8t, assay = "HTO", features = rownames(cd8t[["HTO"]]), ncol = 2)
RidgePlot(cd8t, assay = "RNA", features = c("PTPRC", "SELL"), ncol = 2)
pdf("./plots/cd8.pdf", width = 5, height = 5)
DimPlot(cd8t, label = TRUE, reduction = "umap") + NoLegend()+ ggtitle("CD8")
FeaturePlot(cd8t, "PTPRC") + ggtitle("CD45RO")
FeaturePlot(cd8t, "PTPRC") + ggtitle("CD45RA")
FeaturePlot(cd8t, "SELL") + ggtitle("CD62L")
FeaturePlot(cd8t, "CCR7") + ggtitle("CCR7")
FeaturePlot(cd8t, "IL2RA") + ggtitle("CD25")
FeaturePlot(cd8t, "CD27") + ggtitle("CD27")
FeaturePlot(cd8t, "TCF7") + ggtitle("TCF7")
FeaturePlot(cd8t, "PDCD1") + ggtitle("PD-1")
FeaturePlot(cd8t, "HAVCR2") + ggtitle("TIM3")
FeaturePlot(cd8t, "TIGIT") + ggtitle("TIGIT")
FeaturePlot(cd8t, "LAG3") + ggtitle("LAG3")
FeaturePlot(cd8t, "KLRG1") + ggtitle("KLRG1")
FeaturePlot(cd8t, "IL7R") + ggtitle("IL7R")
FeaturePlot(cd8t, "GZMB") + ggtitle("GZMB")
FeaturePlot(cd8t, "IL2") + ggtitle("IL2")
FeaturePlot(cd8t, "ENTPD1") + ggtitle("CD39")
dev.off()
pdf("./plots/cd4.pdf", width = 5, height = 5)
DimPlot(cd4t, label = TRUE, reduction = "umap") + NoLegend()+ ggtitle("CD4")
FeaturePlot(cd4t, "PTPRC") + ggtitle("CD45RO")
FeaturePlot(cd4t, "PTPRC") + ggtitle("CD45RA")
FeaturePlot(cd4t, "SELL") + ggtitle("CD62L")
FeaturePlot(cd4t, "CCR7") + ggtitle("CCR7")
FeaturePlot(cd4t, "IL2RA") + ggtitle("CD25")
FeaturePlot(cd4t, "CD27") + ggtitle("CD27")
FeaturePlot(cd4t, "TCF7") + ggtitle("TCF7")
FeaturePlot(cd4t, "PDCD1") + ggtitle("PD-1")
FeaturePlot(cd4t, "HAVCR2") + ggtitle("TIM3")
FeaturePlot(cd4t, "TIGIT") + ggtitle("TIGIT")
FeaturePlot(cd4t, "LAG3") + ggtitle("LAG3")
FeaturePlot(cd4t, "KLRG1") + ggtitle("KLRG1")
FeaturePlot(cd4t, "IL7R") + ggtitle("IL7R")
FeaturePlot(cd4t, "GZMB") + ggtitle("GZMB")
FeaturePlot(cd4t, "IL2") + ggtitle("IL2")
dev.off()

Idents(cd8t) <- "RNA_snn_res.0.3"
RidgePlot(cd8t, assay = "RNA", features = "PTPN2")


png("./plots/01FeaturePlot_Memory_CD8.png", width = 40, height = 10, units = "cm", res = 300)
p1 <- FeaturePlot(cd8t, "TCF7") + ggtitle("TCF7")
p2 <- FeaturePlot(cd8t, "CCR7") + ggtitle("CCR7")
p3 <- FeaturePlot(cd8t, "SELL") + ggtitle("SELL")
p4 <- FeaturePlot(cd8t, "IL7R") + ggtitle("IL7R")
p1 | p2 |p3 | p4
dev.off()

png("./plots/01FeaturePlot_Effector_CD8.png", width = 40, height = 10, units = "cm", res = 300)
p1 <- FeaturePlot(cd8t, "IL2RA") + ggtitle("IL2RA")
p2 <- FeaturePlot(cd8t, "CD27") + ggtitle("CD27")
p3 <- FeaturePlot(cd8t, "CD28") + ggtitle("CD28")
p4 <- FeaturePlot(cd8t, "CD69") + ggtitle("CD69")
p1 | p2 |p3 | p4
dev.off()

png("./plots/01FeaturePlot_Exhaustion_CD8.png", width = 70, height = 10, units = "cm", res = 300)
p1 <- FeaturePlot(cd8t, "HAVCR2") + ggtitle("TIM3")
p2 <- FeaturePlot(cd8t, "LAG3") + ggtitle("LAG3")
p3 <- FeaturePlot(cd8t, "PDCD1") + ggtitle("PD-1")
p4 <-FeaturePlot(cd8t, "KLRG1") + ggtitle("KLRG1")
p5 <- FeaturePlot(cd8t, "TIGIT") + ggtitle("TIGIT")
p6 <- FeaturePlot(cd8t, "TOX") + ggtitle("TOX")
p7 <- FeaturePlot(cd8t, "CTLA4") + ggtitle("CTLA4")
p1 | p2 |p3 | p4 | p5 | p6 |p7
dev.off()


png("./plots/01FeaturePlot_Metabolism_CD8.png", width = 40, height = 10, units = "cm", res = 300)
p1 <- FeaturePlot(cd8t, "LDHA") + ggtitle("LDHA")
p2 <- FeaturePlot(cd8t, "SLC2A1") + ggtitle("SLC2A1")
p3 <- FeaturePlot(cd8t, "SLC2A3") + ggtitle("SLC2A3")
p4 <-FeaturePlot(cd8t, "PFKFB3") + ggtitle("PFKFB3")
p1 | p2 |p3 | p4 
dev.off()


png("./plots/01FeaturePlot_markers_CD8.png", width = 50, height = 10, units = "cm", res = 300)
p1 <- FeaturePlot(cd8t, "LDHA") + ggtitle("LDHA")
p2 <- FeaturePlot(cd8t, "TBX21") + ggtitle("TBX21")
p3 <- FeaturePlot(cd8t, "EOMES") + ggtitle("EOMES")
p4 <-FeaturePlot(cd8t, "SLAMF6") + ggtitle("SLAMF6")
p5 <- FeaturePlot(cd8t, "CD45RA") + ggtitle("CD45RA")
p1 | p2 |p3 | p4 |p5
dev.off()

png("./plots/01FeaturePlot_GZM_CD8.png", width = 40, height = 10, units = "cm", res = 300)
p1 <- FeaturePlot(cd8t, "GZMH") + ggtitle("GZMH")
p2 <- FeaturePlot(cd8t, "GZMB") + ggtitle("GZMB")
p3 <- FeaturePlot(cd8t, "NKG7") + ggtitle("NKG7")
p4 <-FeaturePlot(cd8t, "CX3CR1") + ggtitle("CX3CR1")
p1 | p2 |p3 | p4 
dev.off()

######Marker Dotplot######
DefaultAssay(cd4t) <- "RNA"
Idents(cd4t) <-"seurat_clusters"
key.genes <- c("S100A4", "FOXP1", "CXCR6", "LAG3", "KLRG1","IFNG", "PFN1", "GZMB" , "IL7R", "SELL", "FOXP3", "TCF7",
               "PDCD1", "HAVCR2")
png("./plots/01Dotplot_CD4.png", width = 18, height = 15, units = "cm", res = 300)
DotPlot(cd4t, features = key.genes) + RotatedAxis() + ggtitle("CD4")
dev.off()

DefaultAssay(cd8t) <- "RNA"
Idents(cd8t) <-"seurat_clusters"
key.genes <- c("S100A4", "FOXP1", "CXCR6", "LAG3", "KLRG1","IFNG", "PFN1", "GZMB" , "IL7R", "SELL", "FOXP3", "TCF7",
               "PDCD1", "HAVCR2")
DotPlot(cd8t, features = key.genes) + RotatedAxis() + ggtitle("CD8")
dev.off()
key.genes <- c("PTPN2", "S100A4", "FOXP1", "CXCR6", "LAG3", "KLRG1","IFNG", "PFN1", "GZMB" , "IL7R", "SELL", "TCF7","PDCD1", "HAVCR2")
png("./plots/01Dotplot_CD8.png", width = 18, height = 15, units = "cm", res = 300)
DotPlot(cd8t, features = key.genes) + RotatedAxis() + ggtitle("CD8")
dev.off()

DefaultAssay(cd8t) <- "ADT"
VlnPlot(cd8t, "CD45RA") + ggtitle("CD45RA")

selected.gene <- c("CCR7", "IL7R", "SELL", "CX3CR1","TCF7", "NKG7", "PRF1", "GNLY", "GZMH", "GZMB", "IFNG", "LTA", "TNF", "IFIH1", "IFIT3", "IFIT2", "IFIT1", "CENPU", "MCM4", "MKI67", "CENPM", "CENPF", 
                   "TK1", "RRM2",  "EOMES", "HLA-DQB1", "HLA-DPB1", "HLA-DRB5", "JUN", "JUND", "REL", "RELB", "NFKBIA", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ATP6", 
                   "PDCD1", "LAG3", "TIGIT", "HAVCR2", "TOX", "KLRG1", "CTLA4")
png("./plots/01DotPlot_CD8_markers.png", width = 12.5, height = 25, res = 300, units = "cm")
DotPlot(cd8t, features = selected.gene) + RotatedAxis() + ggtitle("CD8")+ coord_flip()
dev.off()


selected.gene <- c("ALDOA", "BPGM", "ENO1", "ENO2", "GAPDH", "GPI", "HK1", "HK2","HKDC1", "PFKL", "PFKM", "PGAM1", "PGAM2", "PGAM4", "PGK1", "PKLR", "PKM", "TPI1")
png("./plots/01DotPlot_CD8_glycolysis.png", width = 12.5, height = 15, res = 300, units = "cm")
DotPlot(cd8t, features = selected.gene) + RotatedAxis() + ggtitle("Glycolysis")+ coord_flip()
dev.off()


######Target Genes######
stat3 <- c("FOSL2", "NTRK3", "IL16", "DUSP10", "CASP1", "FNDC9", "PSD3", "FLT1", "TSHZ2", "RORA", "KCNK1", "NAPSA","FNBP1L", 
           "CALCB", "IL1R1", "COL6A3", "CCR6", "IL24", "HSD11B1", "IFNL1", "IL23R", "MAF", "PALLD", "HOPX", "IFIT2", "GPR87",
           "BCAR3", "IL1R2","DHRS9", "IGFBP6", "PRG4","CXCR5", "GZMB", "IL12RB2", "GBP4","TNFSF13B","TNFSF11","CSF2", "AFF3",
           "TNFSF4", "TMEM200A", "PTGER2", "SLC12A8", "B3GALNT1", "LAIR2", "PHLDA1", "NR4A2",  "QPCT")
stat1 <- c("IDO1", "CXCL10", "IRF1", "GBP5", "ITK", "ICAM1", "APOL6", "TMEM140", "PARP9", "GBP2", "CTSS", "CXCL9", "TRIM21", 
           "GBP1", "CCL2", "DTX3L", "PSMB9", "CCDC68", "APOL1", "PDCD1LG2")
stat5<- c("IL2RA", "CISH", "CDK6", "LRRC32", "MAP3K8", "LTA", "OSM", "LIF", "NABP1", "GZMB")
immune_response <- c("AIM2", "APOL1", "C1S", "C3", "C4A", "CCL2", "CIITA", "CTSS",  "CXCL10", "CXCL9", "GBP1", "GBP2", "GBP5", "GCH1", "HLA-E",
                     "ICAM1", "IFI35", "IL7", "IL4R", "LYN", "ORAI1", "PDCD1LG2",
                     "PSMB8", "PSMB9", "RNF19B", "TAP1", "TAP2", "BCL6", "C1S", "C3", "C4A", "F2RL1", "FYN", "ICAM1", "IDO1",
                     "IL4R", "IL7", "LYN", "PDCD1LG2", "PVR", "TGFB2")
memory_regulation <- c("BCL6", "CD46", "FGL2", "IL23A", "PCK1","TNFSF4")

library(ggpubr)
Idents(cd8t) <- "seurat_clusters"
tcell <- cd8t[,which(Idents(cd8t) == "0")]
Idents(tcell) <- "group"
gene_sig <- stat1
vp_case1(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_cluster0_stat1", test_sign = comparisons)
gene_sig <- stat3
vp_case2(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_cluster0_stat3", test_sign = comparisons)
gene_sig <- stat5
vp_case1(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_cluster0_stat5", test_sign = comparisons)

tcell <- cd8t[,which(Idents(cd8t) == "1")]
Idents(tcell) <- "group"
gene_sig <- stat1
vp_case1(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_cluster1_stat1", test_sign = comparisons)
gene_sig <- stat3
vp_case2(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_cluster1_stat3", test_sign = comparisons)
gene_sig <- stat5
vp_case1(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_cluster1_stat5", test_sign = comparisons)

tcell <- cd8t[,which(Idents(cd8t) == "2")]
Idents(tcell) <- "group"
gene_sig <- stat1
vp_case1(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_cluster2_stat1", test_sign = comparisons)
gene_sig <- stat3
vp_case2(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_cluster2_stat3", test_sign = comparisons)
gene_sig <- stat5
vp_case1(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_cluster2_stat5", test_sign = comparisons)

tcell <- cd8t[,which(Idents(cd8t) == "3")]
Idents(tcell) <- "group"
gene_sig <- stat1
vp_case1(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_cluster3_stat1", test_sign = comparisons)
gene_sig <- stat3
vp_case2(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_cluster3_stat3", test_sign = comparisons)
gene_sig <- stat5
vp_case1(gene_signature = gene_sig, file_name = "./plots/01VlnPlot_CD8_cluster3_stat5", test_sign = comparisons)



tcell <- cd8t[,which(Idents(cd8t) == "0")]
Idents(tcell) <- "group"
png("./plots/01VlnPlot_cluster0_CD8_STAT3.png", width = 50, height = 70, units = "cm", res = 300)
VlnPlot(tcell, features = stat3, ncol = 5)+ stat_compare_means(comparisons = group, label = "p.signif")
dev.off()
png("./plots/01VlnPlot_cluster0_CD8_STAT1.png", width = 50, height = 30, units = "cm", res = 300)
VlnPlot(tcell, features = stat1, ncol = 5)
dev.off()
png("./plots/01VlnPlot_cluster0_CD8_STAT5.png", width = 50, height = 20, units = "cm", res = 300)
VlnPlot(tcell, features = stat5, ncol = 5)+ stat_compare_means(comparisons = group, label = "p.signif")
dev.off()
png("./plots/01VlnPlot_cluster0_CD8_immune_response.png", width = 50, height = 70, units = "cm", res = 300)
VlnPlot(tcell, features = immune_response, ncol = 5)
dev.off()

tcell <- cd8t[,which(Idents(cd8t) == "1")]
Idents(tcell) <- "group"
png("./plots/01VlnPlot_cluster1_CD8_STAT3.png", width = 50, height = 70, units = "cm", res = 300)
VlnPlot(tcell, features = stat3, ncol = 5)
dev.off()
png("./plots/01VlnPlot_cluster1_CD8_STAT1.png", width = 50, height = 30, units = "cm", res = 300)
VlnPlot(tcell, features = stat1, ncol = 5)
dev.off()
png("./plots/01VlnPlot_cluster1_CD8_STAT5.png", width = 50, height = 20, units = "cm", res = 300)
VlnPlot(tcell, features = stat5, ncol = 5)
dev.off()
png("./plots/01VlnPlot_cluster1_CD8_immune_response.png", width = 50, height = 70, units = "cm", res = 300)
VlnPlot(tcell, features = immune_response, ncol = 5)
dev.off()


tcell <- cd8t[,which(Idents(cd8t) == "2")]
Idents(tcell) <- "group"
png("./plots/01VlnPlot_cluster2_CD8_STAT3.png", width = 50, height = 70, units = "cm", res = 300)
VlnPlot(tcell, features = stat3, ncol = 5)
dev.off()
png("./plots/01VlnPlot_cluster2_CD8_STAT1.png", width = 50, height = 30, units = "cm", res = 300)
VlnPlot(tcell, features = stat1, ncol = 5)
dev.off()
png("./plots/01VlnPlot_cluster2_CD8_STAT5.png", width = 50, height = 20, units = "cm", res = 300)
VlnPlot(tcell, features = stat5, ncol = 5)
dev.off()
png("./plots/01VlnPlot_cluster2_CD8_immune_response.png", width = 50, height = 70, units = "cm", res = 300)
VlnPlot(tcell, features = immune_response, ncol = 5)
dev.off()

tcell <- cd8t[,which(Idents(cd8t) == "3")]
Idents(tcell) <- "group"
png("./plots/01VlnPlot_cluster3_CD8_STAT3.png", width = 50, height = 70, units = "cm", res = 300)
VlnPlot(tcell, features = stat3, ncol = 5)
dev.off()
png("./plots/01VlnPlot_cluster3_CD8_STAT1.png", width = 50, height = 30, units = "cm", res = 300)
VlnPlot(tcell, features = stat1, ncol = 5)
dev.off()
png("./plots/01VlnPlot_cluster3_CD8_STAT5.png", width = 50, height = 20, units = "cm", res = 300)
VlnPlot(tcell, features = stat5, ncol = 5)
dev.off()
png("./plots/01VlnPlot_cluster3_CD8_immune_response.png", width = 50, height = 70, units = "cm", res = 300)
VlnPlot(tcell, features = immune_response, ncol = 5)
dev.off()



######relabeling##### 
#0 = Teff, 1= Tscm, 2= Tmem, 3 = Tmem_exh, 4= NaiveT
DefaultAssay(cd8t) <- "RNA"
new.cluster.ids <- c("Teff", "Tscm", "Tmem", "Tmem_exh", "NaiveT")
names(new.cluster.ids) <- levels(cd8t)
cd8t <- RenameIdents(cd8t, new.cluster.ids)
cd8t@meta.data$relabel <-  Idents(cd8t)
DotPlot(cd8t, features = key.genes) + RotatedAxis() + ggtitle("CD8")

png("./plots/01DimPlot_CD8_relabel.png", width = 10, height = 10, res = 300, units = "cm")
DimPlot(cd8t, label = T, reduction = "umap",label.size = 6) + ggtitle("CD8") + NoLegend()
dev.off()



######Percentage######
length(which(grepl("KO", cd8t@meta.data$HTO.ident))) / length(cd8t@meta.data$HTO.ident)
length(which(grepl("KO", cd8t@meta.data$HTO.ident)&grepl("0", cd8t@meta.data$RNA_snn_res.0.3))) / length(which(grepl("0", cd8t@meta.data$RNA_snn_res.0.3)))
length(which(grepl("KO", cd8t@meta.data$HTO.ident)&grepl("1", cd8t@meta.data$RNA_snn_res.0.3))) / length(which(grepl("1", cd8t@meta.data$RNA_snn_res.0.3)))
length(which(grepl("KO", cd8t@meta.data$HTO.ident)&grepl("2", cd8t@meta.data$RNA_snn_res.0.3))) / length(which(grepl("2", cd8t@meta.data$RNA_snn_res.0.3)))
length(which(grepl("KO", cd8t@meta.data$HTO.ident)&grepl("3", cd8t@meta.data$RNA_snn_res.0.3))) / length(which(grepl("3", cd8t@meta.data$RNA_snn_res.0.3)))
length(which(grepl("KO", cd8t@meta.data$HTO.ident)&grepl("4", cd8t@meta.data$RNA_snn_res.0.3))) / length(which(grepl("4", cd8t@meta.data$RNA_snn_res.0.3)))

proportion$Group <- gsub("WT1", "WT", proportion$Group)
proportion$Group <- gsub("WT2", "WT", proportion$Group)
proportion$Group <- gsub("KO1", "KO", proportion$Group)
proportion$Group <- gsub("KO2", "KO", proportion$Group)
ggplot(proportion, aes(fill=Group, y=Proportion, x=Cells)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+ggtitle("CD8") +geom_hline(yintercept=0.4855997, linetype="dashed", 
                                             color = "black", size=2)

Idents(cd8t) <- "RNA_snn_res.0.25"
length(which(grepl("KO", cd8t@meta.data$HTO.ident))) / length(cd8t@meta.data$HTO.ident)
length(which(grepl("KO", cd8t@meta.data$HTO.ident)&grepl("0", cd8t@meta.data$RNA_snn_res.0.25))) / length(which(grepl("0", cd8t@meta.data$RNA_snn_res.0.25)))
length(which(grepl("KO", cd8t@meta.data$HTO.ident)&grepl("1", cd8t@meta.data$RNA_snn_res.0.25))) / length(which(grepl("1", cd8t@meta.data$RNA_snn_res.0.25)))
length(which(grepl("KO", cd8t@meta.data$HTO.ident)&grepl("2", cd8t@meta.data$RNA_snn_res.0.25))) / length(which(grepl("2", cd8t@meta.data$RNA_snn_res.0.25)))
length(which(grepl("KO", cd8t@meta.data$HTO.ident)&grepl("3", cd8t@meta.data$RNA_snn_res.0.25))) / length(which(grepl("3", cd8t@meta.data$RNA_snn_res.0.25)))
proportion <- read.csv("./dataset/proportion_4cluster.csv")
png("./plots/01ProportionBarplot_RNAcluster_CD8_rm4.png", width = 10, height = 10, units = "cm", res = 300)
ggplot(proportion, aes(fill=Group, y=Proportion, x=Cell)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+ggtitle("CD8") +geom_hline(yintercept=0.4843581, linetype="dashed", 
                                             color = "black", size=2)
dev.off()



pdf("./plots/01Heatmap_relabel_CD8.pdf", width = 8, height = 5)
jakstat = read.csv("./dataset/ssGSEA_genelist.csv", stringsAsFactors = FALSE, check.names = FALSE)
stat3 <- c("FOSL2", "NTRK3", "IL16", "DUSP10", "CASP1", "FNDC9", "PSD3", "FLT1", "TSHZ2", "RORA", "KCNK1", "NAPSA","FNBP1L", 
           "CALCB", "IL1R1", "COL6A3", "CCR6", "IL24", "HSD11B1", "IFNL1", "IL23R", "MAF", "PALLD", "HOPX", "IFIT2", "GPR87",
           "BCAR3", "IL1R2","DHRS9", "IGFBP6", "PRG4","CXCR5", "GZMB", "IL12RB2", "GBP4","TNFSF13B","TNFSF11","CSF2", "AFF3",
           "TNFSF4", "TTMEM200A", "PTGER2", "SLC12A8", "B3GALNT1", "LAIR2", "PHLDA1", "NR4A2",  "QPCT")

#CD8 Teff 
Tem <- subset(x = cd8t, idents = "Teff")
Idents(Tem) <- "group"
VlnPlot(Tem, features = c("PTPN2")) + ggtitle("Teff PTPN2")

rna.markers <- FindAllMarkers(Tem, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
write.csv(rna.markers, file = "./cd8_cluster0.csv")
rna.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(Tem, features = top10$gene) + NoLegend()+ggtitle("Teff DEG")
DoHeatmap(Tem, features = jakstat$KEGG_JAK_STAT_SIGNALING_PATHWAY) + NoLegend()+ggtitle("Teff JAK-STAT")
DoHeatmap(Tem, features = stat3) + NoLegend()+ggtitle("Teff STAT3 Targets")

#CD8 Tscm
Tem <- subset(x = cd8t, idents = "Tscm")
Idents(Tem) <- "group"
VlnPlot(Tem, features = c("PTPN2")) + ggtitle("Tscm PTPN2")

rna.markers <- FindAllMarkers(Tem, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
write.csv(rna.markers, file = "./cd8_cluster1.csv")
rna.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(Tem, features = top10$gene) + NoLegend()+ggtitle("Tscm DEG")
DoHeatmap(Tem, features = jakstat$KEGG_JAK_STAT_SIGNALING_PATHWAY) + NoLegend()+ggtitle("Tscm JAK-STAT")
DoHeatmap(Tem, features = stat3) + NoLegend()+ggtitle("Tscm STAT3 Targets")


#CD8 Tmem
Tem <- subset(x = cd8t, idents = "Tmem")
Idents(Tem) <- "group"
VlnPlot(Tem, features = c("PTPN2")) + ggtitle("Tmem PTPN2")

rna.markers <- FindAllMarkers(Tem, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
write.csv(rna.markers, file = "./cd8_cluster2.csv")
rna.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(Tem, features = top10$gene) + NoLegend()+ggtitle("Tmem DEG")
DoHeatmap(Tem, features = jakstat$KEGG_JAK_STAT_SIGNALING_PATHWAY) + NoLegend()+ggtitle("Tmem JAK-STAT")
DoHeatmap(Tem, features = stat3) + NoLegend()+ggtitle("Tmem STAT3 Targets")


#CD8 Tmem_exh 
Tem <- subset(x = cd8t, idents = "Tmem_exh")
Idents(Tem) <- "group"
VlnPlot(Tem, features = c("PTPN2")) + ggtitle("Tmem_exh PTPN2")

rna.markers <- FindAllMarkers(Tem, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
write.csv(rna.markers, file = "./cd8_cluster3.csv")
rna.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(Tem, features = top10$gene) + NoLegend()+ggtitle("Tmem_exh DEG")
DoHeatmap(Tem, features = jakstat$KEGG_JAK_STAT_SIGNALING_PATHWAY) + NoLegend()+ggtitle("Tmem_exh JAK-STAT")
DoHeatmap(Tem, features = stat3) + NoLegend()+ggtitle("Tmem_exh STAT3 Targets")


#CD8 NaiveT 
Tem <- subset(x = cd8t, idents = "NaiveT")
Idents(Tem) <- "group"
VlnPlot(Tem, features = c("PTPN2")) + ggtitle("NaiveT PTPN2")

rna.markers <- FindAllMarkers(Tem, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
write.csv(rna.markers, file = "./cd8_cluster4.csv")
rna.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(Tem, features = top10$gene) + NoLegend()+ggtitle("NaiveT DEG")
DoHeatmap(Tem, features = jakstat$KEGG_JAK_STAT_SIGNALING_PATHWAY) + NoLegend()+ggtitle("NaiveT JAK-STAT")
DoHeatmap(Tem, features = stat3) + NoLegend()+ggtitle("NaiveT STAT3 Targets")
dev.off()



#####CD8+ remove cluster4 20230112######
load("./seurat_cd8t.rds")
FeatureScatter(cd8t, feature1 = "CD8", feature2 = "CD4") 
cd8t <- subset(cd8t, idents = "4", invert = T) #remove cluster 4

tiff("./plots/01DimPlot_CD8_RNA_rm4.tif", width = 10, height = 8, res = 300, units = "cm")
Idents(cd8t) <- "RNA_snn_res.0.3"
DimPlot(cd8t, label = F, reduction = "umap") + ggtitle("CD8") #+ NoLegend()
dev.off()
tiff("./plots/01DimPlot_CD8_RNA_WTKO_rm4.tif", width = 10, height = 8, res = 300, units = "cm")
Idents(cd8t) <- "group"
DimPlot(cd8t, label = F, reduction = "umap", cols = c("mediumpurple","red"), pt.size =2) + ggtitle("CD8") 
dev.off()
png("./plots/01DimPlot_CD8_RNAcluster_WTKO_rm4.png", width = 20, height = 8, res = 300, units = "cm")
Idents(cd8t) <- "group"
wt.cd8t <- subset(cd8t, idents = "WT")
Idents(wt.cd8t) <- "RNA_snn_res.0.3"
p1 <- DimPlot(wt.cd8t, label = F, reduction = "umap") + ggtitle("CD8 WT") 
Idents(cd8t) <- "group"
wt.cd8t <- subset(cd8t, idents = "KO")
Idents(wt.cd8t) <- "RNA_snn_res.0.3"
p2 <- DimPlot(wt.cd8t, label = F, reduction = "umap") + ggtitle("CD8 KO") 
p1 | p2
dev.off()


Idents(cd8t) <- "RNA_snn_res.0.3"
cbmc.rna.markers <- FindAllMarkers(cd8t, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cbmc.rna.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
cbmc.rna.markers %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC) -> top30
cbmc.rna.markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top100
DoHeatmap(cd8t, features = top10$gene) + NoLegend()+ggtitle("CD8")
png("./plots/01Heatmap_RNAcluster_CD8_rm4_top100.png", height = 10, width = 10, units = "cm", res = 600)
DoHeatmap(cd8t, features = top100$gene, draw.lines = TRUE, hjust = -10, angle = 90) + 
  NoLegend()+ggtitle("CD8")+ 
  scale_fill_gradientn(colors = c("#4575B4", "#91BFDB","white", "#FC8D59", "#D73027"))+
  theme(axis.text=element_text(size=1), axis.title=element_text(size=14,face="bold"))
dev.off()
write.csv(cbmc.rna.markers, file = "./cd8_markers_RNA_rm4.csv")


png("./plots/01FeaturePlot_TF_CD8_rm4.png", width = 40, height = 10, units = "cm", res = 300)
p1 <- FeaturePlot(cd8t, "TCF7") + ggtitle("TCF7")
p2 <- FeaturePlot(cd8t, "FOXP1") + ggtitle("FOXP1")
p3 <- FeaturePlot(cd8t, "NFATC1") + ggtitle("NFATC1")
p4 <- FeaturePlot(cd8t, "FOXP3") + ggtitle("FOXP3")
p1 | p2 |p3 | p4
dev.off()



RidgePlot(cd8t, assay = "HTO", features = rownames(cd8t[["HTO"]]), ncol = 2)
RidgePlot(cd8t, assay = "RNA", features = c("PTPRC", "SELL"), ncol = 2)
pdf("./plots/cd8_rm4.pdf", width = 5, height = 5)
DimPlot(cd8t, label = TRUE, reduction = "umap") + NoLegend()+ ggtitle("CD8")
FeaturePlot(cd8t, "PTPRC") + ggtitle("CD45RO")
FeaturePlot(cd8t, "PTPRC") + ggtitle("CD45RA")
FeaturePlot(cd8t, "SELL") + ggtitle("CD62L")
FeaturePlot(cd8t, "CCR7") + ggtitle("CCR7")
FeaturePlot(cd8t, "IL2RA") + ggtitle("CD25")
FeaturePlot(cd8t, "CD27") + ggtitle("CD27")
FeaturePlot(cd8t, "TCF7") + ggtitle("TCF7")
FeaturePlot(cd8t, "PDCD1") + ggtitle("PD-1")
FeaturePlot(cd8t, "HAVCR2") + ggtitle("TIM3")
FeaturePlot(cd8t, "TIGIT") + ggtitle("TIGIT")
FeaturePlot(cd8t, "LAG3") + ggtitle("LAG3")
FeaturePlot(cd8t, "KLRG1") + ggtitle("KLRG1")
FeaturePlot(cd8t, "IL7R") + ggtitle("IL7R")
FeaturePlot(cd8t, "GZMB") + ggtitle("GZMB")
FeaturePlot(cd8t, "IL2") + ggtitle("IL2")
FeaturePlot(cd8t, "ENTPD1") + ggtitle("CD39")
dev.off()


Idents(cd8t) <- "RNA_snn_res.0.3"
RidgePlot(cd8t, assay = "RNA", features = "PTPN2")


png("./plots/01FeaturePlot_Memory_CD8_rm4.png", width = 40, height = 10, units = "cm", res = 300)
p1 <- FeaturePlot(cd8t, "TCF7") + ggtitle("TCF7")
p2 <- FeaturePlot(cd8t, "CCR7") + ggtitle("CCR7")
p3 <- FeaturePlot(cd8t, "SELL") + ggtitle("SELL")
p4 <- FeaturePlot(cd8t, "IL7R") + ggtitle("IL7R")
p1 | p2 |p3 | p4
dev.off()

tiff("./plots/01FeaturePlot_ext_CD8_rm4.tif", width = 40, height = 10, units = "cm", res = 300)
p1 <- FeaturePlot(cd8t, "PRDM1") + ggtitle("PRDM1")
p2 <- FeaturePlot(cd8t, "BACH2") + ggtitle("BACH2")
p3 <- FeaturePlot(cd8t, "LEF1") + ggtitle("LEF1")
p4 <- FeaturePlot(cd8t, "PRF1") + ggtitle("PRF1")
p1 | p2 |p3 | p4
dev.off()

png("./plots/01FeaturePlot_Effector_CD8_rm4.png", width = 40, height = 10, units = "cm", res = 300)
p1 <- FeaturePlot(cd8t, "IL2RA") + ggtitle("IL2RA")
p2 <- FeaturePlot(cd8t, "CD27") + ggtitle("CD27")
p3 <- FeaturePlot(cd8t, "CD28") + ggtitle("CD28")
p4 <- FeaturePlot(cd8t, "CD69") + ggtitle("CD69")
p1 | p2 |p3 | p4
dev.off()

png("./plots/01FeaturePlot_Exhaustion_CD8_rm4.png", width = 70, height = 10, units = "cm", res = 300)
p1 <- FeaturePlot(cd8t, "HAVCR2") + ggtitle("TIM3")
p2 <- FeaturePlot(cd8t, "LAG3") + ggtitle("LAG3")
p3 <- FeaturePlot(cd8t, "PDCD1") + ggtitle("PD-1")
p4 <-FeaturePlot(cd8t, "KLRG1") + ggtitle("KLRG1")
p5 <- FeaturePlot(cd8t, "TIGIT") + ggtitle("TIGIT")
p6 <- FeaturePlot(cd8t, "TOX") + ggtitle("TOX")
p7 <- FeaturePlot(cd8t, "CTLA4") + ggtitle("CTLA4")
p1 | p2 |p3 | p4 | p5 | p6 |p7
dev.off()


png("./plots/01FeaturePlot_Metabolism_CD8_rm4.png", width = 40, height = 10, units = "cm", res = 300)
p1 <- FeaturePlot(cd8t, "LDHA") + ggtitle("LDHA")
p2 <- FeaturePlot(cd8t, "SLC2A1") + ggtitle("SLC2A1")
p3 <- FeaturePlot(cd8t, "SLC2A3") + ggtitle("SLC2A3")
p4 <-FeaturePlot(cd8t, "PFKFB3") + ggtitle("PFKFB3")
p1 | p2 |p3 | p4 
dev.off()


png("./plots/01FeaturePlot_markers_CD8_rm4.png", width = 50, height = 10, units = "cm", res = 300)
p1 <- FeaturePlot(cd8t, "LDHA") + ggtitle("LDHA")
p2 <- FeaturePlot(cd8t, "TBX21") + ggtitle("TBX21")
p3 <- FeaturePlot(cd8t, "EOMES") + ggtitle("EOMES")
p4 <-FeaturePlot(cd8t, "SLAMF6") + ggtitle("SLAMF6")
p5 <- FeaturePlot(cd8t, "CD45RA") + ggtitle("CD45RA")
p1 | p2 |p3 | p4 |p5
dev.off()

png("./plots/01FeaturePlot_GZM_CD8_rm4.png", width = 40, height = 10, units = "cm", res = 300)
p1 <- FeaturePlot(cd8t, "GZMH") + ggtitle("GZMH")
p2 <- FeaturePlot(cd8t, "GZMB") + ggtitle("GZMB")
p3 <- FeaturePlot(cd8t, "NKG7") + ggtitle("NKG7")
p4 <-FeaturePlot(cd8t, "CX3CR1") + ggtitle("CX3CR1")
p1 | p2 |p3 | p4 
dev.off()

save(cd8t, file = "./seurat_cd8t_rm4.rds")


Idents(cd8t) <- "RNA_snn_res.0.3"
length(which(grepl("KO", cd8t@meta.data$HTO.ident))) / length(cd8t@meta.data$HTO.ident)
length(which(grepl("KO", cd8t@meta.data$HTO.ident)&grepl("0", cd8t@meta.data$RNA_snn_res.0.3))) / length(which(grepl("0", cd8t@meta.data$RNA_snn_res.0.3)))
length(which(grepl("KO", cd8t@meta.data$HTO.ident)&grepl("1", cd8t@meta.data$RNA_snn_res.0.3))) / length(which(grepl("1", cd8t@meta.data$RNA_snn_res.0.3)))
length(which(grepl("KO", cd8t@meta.data$HTO.ident)&grepl("2", cd8t@meta.data$RNA_snn_res.0.3))) / length(which(grepl("2", cd8t@meta.data$RNA_snn_res.0.3)))
length(which(grepl("KO", cd8t@meta.data$HTO.ident)&grepl("3", cd8t@meta.data$RNA_snn_res.0.3))) / length(which(grepl("3", cd8t@meta.data$RNA_snn_res.0.3)))
proportion <- read.csv("./dataset/proportion_4cluster.csv")
png("./plots/01ProportionBarplot_RNAcluster_CD8_rm4.png", width = 10, height = 10, units = "cm", res = 300)
ggplot(proportion, aes(fill=Group, y=Proportion, x=Cell)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+ggtitle("CD8") +geom_hline(yintercept=0.4843581, linetype="dashed", 
                                             color = "black", size=2)
dev.off()

binom.test(length(which(grepl("KO", cd8t@meta.data$HTO.ident)&grepl("0", cd8t@meta.data$RNA_snn_res.0.3))), length(which(grepl("0", cd8t@meta.data$RNA_snn_res.0.3))), p = 0.4843581, alternative = "two.sided")
binom.test(length(which(grepl("KO", cd8t@meta.data$HTO.ident)&grepl("1", cd8t@meta.data$RNA_snn_res.0.3))), length(which(grepl("1", cd8t@meta.data$RNA_snn_res.0.3))), p = 0.4843581, alternative = "two.sided")
binom.test(length(which(grepl("KO", cd8t@meta.data$HTO.ident)&grepl("2", cd8t@meta.data$RNA_snn_res.0.3))), length(which(grepl("2", cd8t@meta.data$RNA_snn_res.0.3))), p = 0.4843581, alternative = "two.sided")
binom.test(length(which(grepl("KO", cd8t@meta.data$HTO.ident)&grepl("3", cd8t@meta.data$RNA_snn_res.0.3))), length(which(grepl("3", cd8t@meta.data$RNA_snn_res.0.3))), p = 0.4843581, alternative = "two.sided")


selected.gene <- c("IL7R","CCR7",  "SELL", "TCF7", "NKG7", "PRF1", "GNLY", "GZMH", "GZMB", "HLA-DPB1", "HLA-DRB5", "CD69",
                   "IFIH1", "IFIT3", "IFIT2", "IFIT1",
                   "REL", "RELB", "NFKBIA", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ATP6", 
                   "TPI1","PGK1", "LDHA", "ENO1",  "PKM",  "GAPDH", 
                    "LAG3",   "CX3CR1", "KLRG1")
tiff("./plots/01DotPlot_CD8_markers.tif", width = 12, height = 20, res = 300, units = "cm")
DotPlot(cd8t, features = selected.gene) + RotatedAxis() + ggtitle("CD8") + coord_flip()
dev.off()

cd8t@assays$RNA@data["CD25",]

pal.group <- brewer.pal(3, "Dark2")[1:2]
names(pal.group) <- c("WT", "KO")
Idents(cd8t) <- "RNA_snn_res.0.3"
T_RNAcluster <- cd8t[,which(Idents(cd8t) == "0")]
Idents(T_RNAcluster) <- "group"
rna.markers <- FindAllMarkers(T_RNAcluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
rownames(rna.markers[which(rna.markers$p_val_adj < 0.05),])
tiff("./plots/01VlnPlot_WTvsKO_DEG_RNAcluster0_rm4.tif", width = 10, height = 20, units = "cm", res = 300)
VlnPlot(T_RNAcluster, features = rownames(rna.markers[which(rna.markers$p_val_adj < 0.05),]), cols = pal.group, ncol = 1) + 
  stat_compare_means(comparisons = list(c("WT", "KO")), label = "p.signif")
dev.off()

T_RNAcluster <- cd8t[,which(Idents(cd8t) == "1")]
Idents(T_RNAcluster) <- "group"
rna.markers <- FindAllMarkers(T_RNAcluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
rownames(rna.markers[which(rna.markers$p_val_adj < 0.05),])
tiff("./plots/01VlnPlot_WTvsKO_DEG_RNAcluster1_rm4.tif", width = 10, height = 60, units = "cm", res = 300)
VlnPlot(T_RNAcluster, features = rownames(rna.markers[which(rna.markers$p_val_adj < 0.05),]), cols = pal.group, ncol = 1) + 
  stat_compare_means(comparisons = list(c("WT", "KO")), label = "p.signif")
dev.off()

T_RNAcluster <- cd8t[,which(Idents(cd8t) == "2")]
Idents(T_RNAcluster) <- "group"
rna.markers <- FindAllMarkers(T_RNAcluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
rownames(rna.markers[which(rna.markers$p_val_adj < 0.05),])
tiff("./plots/01VlnPlot_WTvsKO_DEG_RNAcluster2_rm4.tif", width = 10, height = 10, units = "cm", res = 300)
VlnPlot(T_RNAcluster, features = rownames(rna.markers[which(rna.markers$p_val_adj < 0.05),]), cols = pal.group, ncol = 1) + 
  stat_compare_means(comparisons = list(c("WT", "KO")), label = "p.signif")
dev.off()

T_RNAcluster <- cd8t[,which(Idents(cd8t) == "3")]
Idents(T_RNAcluster) <- "group"
rna.markers <- FindAllMarkers(T_RNAcluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
rownames(rna.markers[which(rna.markers$p_val_adj < 0.05),])
tiff("./plots/01VlnPlot_WTvsKO_DEG_RNAcluster3_rm4.tif", width = 10, height = 10, units = "cm", res = 300)
VlnPlot(T_RNAcluster, features = rownames(rna.markers[which(rna.markers$p_val_adj < 0.05),]), cols = pal.group, ncol = 1) + 
  stat_compare_means(comparisons = list(c("WT", "KO")), label = "p.signif",  na.rm = TRUE)
dev.off()

#####CD4 T cells ####
load("./seurat_all.rds")
Idents(cart) <- "HTO.ident"
table(Idents(cart))
length(grep("KO", Idents(cart) ) )/ ncol(cart)

DefaultAssay(cart) <- "ADT"
p1 <- FeaturePlot(cart, "CD4", cols = c("lightgrey", "darkgreen")) + ggtitle("CD4 protein")
DefaultAssay(cart) <- "RNA"
p2 <- FeaturePlot(cart, "CD4") + ggtitle("CD4 RNA")
p1 | p2

#CD4
DefaultAssay(cbmc) <- "ADT"
cart <- subset(cbmc, subset =  CART>0) 
FeatureScatter(cart, feature1 = "CD4", feature2 = "CD8") 
cd4t <- subset(cart, subset =  CD4 >0.61 & CD8 < 1) # CD8 T
Idents(cd4t) <- cd4t@meta.data$HTO.ident
cd4t <- RenameIdents(object = cd4t, 'KO1' = "KO")
cd4t <- RenameIdents(object = cd4t, 'KO2' = "KO")
cd4t <- RenameIdents(object = cd4t, 'WT1' = "WT")
cd4t <- RenameIdents(object = cd4t, 'WT2' = "WT")
cd4t@meta.data$group <- Idents(cd4t)
DefaultAssay(cd4t) <- "RNA"
cd4t <- FindVariableFeatures(cd4t)
cd4t <- ScaleData(cd4t)
cd4t <- RunPCA(cd4t, verbose = TRUE)
ElbowPlot(cd4t, ndims = 20)
cd4t <- FindNeighbors(cd4t, dims = 1:20)
cd4t <- FindClusters(cd4t, resolution = c(0.1,0.3,0.5)) #0.3 for 6 subclusters
cd4t <- FindClusters(cd4t, resolution = c(0.1)) #0.3 for 6 subclusters

cd4t <- RunUMAP(cd4t, dims = 1:20)
png("./plots/01DimPlot_CD4_RNA.png", width = 10, height = 8, res = 300, units = "cm")
Idents(cd4t) <- "RNA_snn_res.0.1"
DimPlot(cd4t, label = TRUE, reduction = "umap") + NoLegend()+ ggtitle("RNA_snn_res.0.1")
Idents(cd4t) <- "group"
DimPlot(cd4t, label = TRUE, reduction = "umap") + NoLegend()+ ggtitle("RNA_snn_res.0.1")
Idents(cd4t) <- "RNA_snn_res.0.3"
DimPlot(cd4t, label = TRUE, reduction = "umap") + NoLegend()+ ggtitle("RNA_snn_res.0.3")
Idents(cd4t) <- "RNA_snn_res.0.5"
DimPlot(cd4t, label = TRUE, reduction = "umap") + NoLegend()+ ggtitle("RNA_snn_res.0.5")
dev.off()
VlnPlot(cd4t, features = "PTPN2", pt.size = 0)+ ggtitle("CD4 Tcell PTPN2")
VlnPlot(cd4t, features = "STAT5", pt.size = 0)+ ggtitle("CD4 Tcell PTPN2")
png("./plots/01DimPlot_CD4_RNA_WTKO.png", width = 10, height = 8, res = 300, units = "cm")
Idents(cd4t) <- "group"
DimPlot(cd4t, label = F, reduction = "umap", cols = pal.group) + ggtitle("CD4") 
dev.off()


##find markers for each cluster
Idents(cd4t) <- "RNA_snn_res.0.1"
pbmc.markers <- FindAllMarkers(cd4t[grep("^MT-|^RP[SL]", rownames(cd4t), invert = T),], only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top100
write.csv(top100, file = "./cd4_markers_res.0.1_top100.csv")

DoHeatmap(cd4t, features = top10$gene)


Idents(cd4t) <- "RNA_snn_res.0.3"
pbmc.markers <- FindAllMarkers(cd4t[grep("^MT-|^RP[SL]", rownames(cd4t), invert = T),], only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top100
write.csv(top100, file = "./cd4_markers_res.0.3_top100.csv")

DoHeatmap(cd4t, features = top10$gene)


Idents(cd4t) <- "RNA_snn_res.0.5"
pbmc.markers <- FindAllMarkers(cd4t[grep("^MT-|^RP[SL]", rownames(cd4t), invert = T),], only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top100
write.csv(top100, file = "./cd4_markers_res.0.5_top100.csv")

DoHeatmap(cd4t, features = top10$gene)


##percentage 
group.df <- data.frame(Percentage =rep(0, 2*4), cluster=c(rep(c("0","1","2","3"), 2)))
Idents(cd4t) <- cd4t[["RNA_snn_res.0.1"]] 
tmp <- cd4t[,grep("KO", cd4t@meta.data$group)]
Idents(tmp) <- tmp[["RNA_snn_res.0.1"]] 
table(Idents(tmp))/ncol(tmp)
table(Idents(tmp))
#DimPlot(tmp, cols = pal.rename)+ggtitle("N_Ctrl")
group.df$Percentage[1:4] <- table(Idents(tmp))/ncol(tmp)
tmp <- cd4t[,grep("WT", cd4t@meta.data$group)]
Idents(tmp) <- tmp[["RNA_snn_res.0.1"]] 
table(Idents(tmp))/ncol(tmp)
table(Idents(tmp))
#DimPlot(tmp, cols = pal.rename)+ggtitle("N_Ctrl")
group.df$Percentage[5:8] <- table(Idents(tmp))/ncol(tmp)

group.df$group <- c(rep("KO",4), rep("WT", 4))
group.df$group <- factor(group.df$group, levels = c("WT", "KO" ) )
group.df$Percentage <- round(group.df$Percentage, 4)
group.df$cluster <- factor(group.df$cluster, levels = c("0", "1", "2","3"))

#https://r-charts.com/colors/
png("./plots/31percentage.png", width = 15, height = 15, units = "cm", res = 600)
group.df[,] %>%
  ggplot(data = ., mapping = aes(x = group, y = Percentage, fill = cluster)) +
  geom_col() +
  geom_text(mapping = aes(label = paste0(Percentage*100, "%") ),              # converting the values to percent
            size = 5,                                             # size of the font
            position = position_stack(vjust = 0.5)) +             # positioning in the middle
  #scale_fill_brewer(palette = "Accent") +                           # coloring the plot
  #scale_fill_manual(values = pal.rename)+
  #facet_grid(.~sex) +
  labs(x = "Group",                                              # labelling x axis
       y = "Percentage",                                        # labeling y axis
       title = "Percentage (CD4)",        # title
       fill = "Cluster") +                               # legend
  scale_y_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 90,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  ) +
  geom_hline(yintercept=0.4843581, linetype="dashed", color = "black", size=2)


group.df <- data.frame(Percentage =rep(0, 2*4), group=c(rep(c("WT", "KO"), 4)))
Idents(cd4t) <- cd4t[["RNA_snn_res.0.1"]] 
tmp <- cd4t[,grep("0", Idents(cd4t))]
Idents(tmp) <- tmp[["group"]] 
table(Idents(tmp))/ncol(tmp)
table(Idents(tmp))
group.df$Percentage[1:2] <- table(Idents(tmp))/ncol(tmp)
tmp <- cd4t[,grep("1", Idents(cd4t))]
Idents(tmp) <- tmp[["group"]] 
table(Idents(tmp))/ncol(tmp)
table(Idents(tmp))
group.df$Percentage[3:4] <- table(Idents(tmp))/ncol(tmp)
tmp <- cd4t[,grep("2", Idents(cd4t))]
Idents(tmp) <- tmp[["group"]] 
table(Idents(tmp))/ncol(tmp)
table(Idents(tmp))
group.df$Percentage[5:6] <- table(Idents(tmp))/ncol(tmp)
tmp <- cd4t[,grep("3", Idents(cd4t))]
Idents(tmp) <- tmp[["group"]] 
table(Idents(tmp))/ncol(tmp)
table(Idents(tmp))
group.df$Percentage[7:8] <- table(Idents(tmp))/ncol(tmp)

group.df$cluster <- c("C0", "C0", "C1", "C1", "C2", "C2", "C3", "C3")
group.df$group <- factor(group.df$group, levels = c("WT", "KO" ) )
group.df$Percentage <- round(group.df$Percentage, 4)
group.df$cluster <- factor(group.df$cluster, levels = c("C0","C1","C2","C3"))

Idents(cd4t) <- "group"
tmp <- cd4t[,grep("KO", Idents(cd4t))]
ncol(tmp)/ncol(cd4t)
#https://r-charts.com/colors/
png("./plots/31percentage.png", width = 15, height = 15, units = "cm", res = 600)
group.df[,] %>%
  ggplot(data = ., mapping = aes(x = cluster, y = Percentage, fill = group)) +
  geom_col() +
  geom_text(mapping = aes(label = paste0(Percentage*100, "%") ),              # converting the values to percent
            size = 5,                                             # size of the font
            position = position_stack(vjust = 0.5)) +             # positioning in the middle
  #scale_fill_brewer(palette = "Accent") +                           # coloring the plot
  #scale_fill_manual(values = pal.rename)+
  #facet_grid(.~sex) +
  labs(x = "Cluster",                                              # labelling x axis
       y = "Percentage",                                        # labeling y axis
       title = "Percentage (CD4)",        # title
       fill = "Cluster") +                               # legend
  scale_y_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 90,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  ) +
  geom_hline(yintercept=0.3437713, linetype="dashed", color = "black", size=2)

binom.test(length(which(grepl("KO", cd4t@meta.data$HTO.ident)&grepl("0", cd4t@meta.data$RNA_snn_res.0.1))), length(which(grepl("0", cd4t@meta.data$RNA_snn_res.0.1))), 
           p = 0.3437713, alternative = "two.sided")
binom.test(length(which(grepl("KO", cd4t@meta.data$HTO.ident)&grepl("1", cd4t@meta.data$RNA_snn_res.0.1))), length(which(grepl("1", cd4t@meta.data$RNA_snn_res.0.1))), 
           p = 0.3437713, alternative = "two.sided")



## KO vs WT
Idents(cd4t) <- "group"
pbmc.markers <- FindAllMarkers(cd4t[grep("^MT-|^RP[SL]", rownames(cd4t), invert = T),], only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top100

VlnPlot(cd4t,features = c(pbmc.markers$gene[4:nrow(pbmc.markers)], "PTPN1"), pt.size = 0, ncol = 6)
