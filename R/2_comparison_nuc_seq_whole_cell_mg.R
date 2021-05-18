library(tidyverse)
library(Seurat)
library(data.table)
library(Matrix)
library(clustree)

source(file.path("R", "functions.R"))

if (!file.exists(file.path("data", "seurat_object.RData"))) {
  #load datasets
  #jäkel
  ms_jakel <- read_tsv(file.path("data", "GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt")) %>%
    mutate(cell_type = str_replace_all(Celltypes , c("[1-9]" = "", 
                                                     "_Macrophages" = "", 
                                                     "Endothelial_cells" = "Endothelial-mural", 
                                                     "Pericytes" = "Endothelial-mural",
                                                     "Vasc_smooth_muscle" = "Endothelial-mural",
                                                     "Oligo" = "Oligodendrocytes"), Celltypes),
           Condition = gsub("Ctrl", "Control", .$Condition)) %>%
    filter(!cell_type %in% c("ImOlGs", "COPs", "Immune_cells", "Macrophages"))
  colnames(ms_jakel)[4] <- "diagnosis"
  
  #load counts
  jakel_umis <- fread(file.path("data", "GSE118257_MSCtr_snRNA_ExpressionMatrix_R(1).txt.gz")) %>%
    as.data.frame()
  rownames(jakel_umis) <- jakel_umis$V1
  jakel_umis$V1 <- NULL
  
  ms_jakel$UMIs <- colSums(jakel_umis[,ms_jakel$Detected])
  ms_jakel <- ms_jakel %>%
    mutate(norm_UMIs = 10000*(UMIs-min(UMIs))/(max(UMIs)- min(UMIs))) 
  
  jakel_umis <- jakel_umis[,ms_jakel$Detected[ms_jakel$diagnosis=="Control" & ms_jakel$cell_type =="Microglia"] ] 
  
  #sankowski et al
  load("/Users/romansankowski/Documents/single_cell_analysis/sankowski-et-al-microglia/data/prdata.RData")
  load("/Users/romansankowski/Documents/single_cell_analysis/sankowski-et-al-microglia/data/metadata_ctrl.RData")
  
  whole_micr <- prdata[,df$ID]
  rownames(whole_micr) <- make.unique(gsub("_.*", "", rownames(prdata)))
  boxplot(whole_micr["MID1",])
  
  # schirmer et al data downloaded from https://cells.ucsc.edu/ms/rawMatrix.zip
  schirmer_et_al <- Read10X(file.path("data","schirmer_et_al"))
  
  #schirmer metadata
  schirm <- read_tsv(file.path("data", "Schirmer_MS.txt")) %>%
    mutate(cell_type = str_replace_all(cell_type, c("EN-L2-3-Cntl|EN-L2-3-MS|EN-L4|EN-L5-6|EN-MIX|EN-PYR" = "Neuron", 
                                                    "IN-PV|IN-SST|IN-SV2C|IN-VIP" = "Neuron",
                                                    "OL-Cntl|OL-MS-1|OL-MS-2" = "Oligodendrocytes",
                                                    "Phagocytes" = "Macrophages",
                                                    "Endo cells" = "Endothelial-mural",
                                                    "OPC" = "OPCs"))) %>%
    filter(!cell_type %in% c("Stromal cells", "B cells", "T cells", "Glia-MIX", "Macrophages")) %>%
    mutate(norm_UMIs = 10000*(UMIs-min(UMIs))/(max(UMIs)- min(UMIs))) 
  
  #extract microglia
  schirmer_et_al <- as.data.frame(as.matrix(schirmer_et_al[, schirm$cell[schirm$cell_type == "Microglia"& schirm$diagnosis == "Control"]]))
  boxplot(schirmer_et_al["MID1",])
  
  #velmeshev
  velm <- read_tsv(file.path("data", "Schirmer_Autism.txt")) %>%
    mutate(cell_type = str_replace_all(cluster, c("L2/3|L4|L5/6$|L5/6-CC$|Neu-mat|Neu-NRGN-I$|Neu-NRGN-II" = "Neuron", 
                                                  "IN-PV|IN-SST|IN-SV2C|IN-VIP" = "Neuron",
                                                  "AST-FB|AST-PP" = "Astrocytes",
                                                  "Phagocytes" = "Macrophages",
                                                  "Endothelial" = "Endothelial-mural",
                                                  "OPC" = "OPCs"))) %>%
    filter(!cell_type %in% c("Stromal cells", "B cells", "T cells", "Glia-MIX"))  %>%
    mutate(norm_UMIs = 10000*(UMIs-min(UMIs))/(max(UMIs)- min(UMIs))) 
  
  #raw counts extracted from https://cells.ucsc.edu/?ds=autism
  velm_et_al <- Read10X(file.path("data", "velmeshev_et_al"))
  velm_et_al <- as.data.frame(as.matrix(velm_et_al[,velm$cell[velm$diagnosis=="Control" & velm$cell_type=="Microglia"]]))
  boxplot(velm_et_al["MID1",])
  
#     Gerrits E, Heng Y, Boddeke EWGM, Eggen BJL. Transcriptional profiling of microglia; current state of the art and future perspectives. Glia 2020 Apr;68(4):740-755. PMID: 31846124
  gerrits1 <- ReadMtx(mtx = file.path("data", "GSM4023552_Donor1nuclei_matrix.mtx.gz"), 
                      cells = file.path("data", "GSM4023552_Donor1nuclei_barcodes.tsv.gz"),
                      features = file.path("data", "GSM4023552_Donor1nuclei_genes.tsv.gz"))
  
  gerrits1 <- as.matrix(gerrits1[rowSums(gerrits1) >0, colSums(gerrits1) >200]) 
  
  gerrits2 <- ReadMtx(mtx = file.path("data", "GSM4023555_Donor2nuclei_matrix.mtx.gz"), 
                      cells = file.path("data", "GSM4023555_Donor2nuclei_barcodes.tsv.gz"),
                      features = file.path("data", "GSM4023555_Donor2nuclei_genes.tsv.gz"))
  gerrits2 <- as.matrix(gerrits2[rowSums(gerrits2) >0, colSums(gerrits2) >200]) 
  
  gerrits_et_al <- merge(gerrits1, gerrits2, by=0)
  rownames(gerrits_et_al) <- gerrits_et_al$Row.names
  gerrits_et_al$Row.names <- NULL
  boxplot(gerrits_et_al["MID1",])
  
  #run seurat
  #all
  all <- list(
    "whole_cells" = whole_micr,
    "Jakel_et_al" = jakel_umis,
    "Schirmer_et_al" = schirmer_et_al,
    "Velm_et_al" = velm_et_al,
    "Gerrits_et_al" = gerrits_et_al
  )
  
  all <- all %>% 
    map(function(x) CreateSeuratObject(counts = x, min.cells = 10, min.features = 500))
  
  all <- merge(x=all[[1]], y=c(all[[2]], all[[3]], all[[4]], all[[5]]),add.cell.ids = c("whole_cells", "Jakel_et_al", "Schirmer_et_al", "Velm_et_al", "Gerrits_et_al"))
  
  all$dataset <- case_when(
    grepl("whole_cells", colnames(all))  ~ "whole_cells",
    grepl("Jakel_et_al", colnames(all))  ~ "Jakel_et_al",
    grepl("Schirmer_et_al", colnames(all))  ~ "Schirmer_et_al",
    grepl("Gerrits_et_al", colnames(all))  ~ "Gerrits_et_al",
    T  ~ "Velm_et_al",
    
  )
  
  #normalize dataset
  all <- all %>%
    SCTransform(vars.to.regress = c("dataset"),
                variable.features.n = 10000) %>%
    RunPCA() 
  
  #
  ElbowPlot(all)
  
  all<- all %>% 
    RunUMAP(dims=1:20) %>%
    FindNeighbors(dims=1:20) %>%
    FindClusters(resolution=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1))
  
  #url https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html#seurat-objects
  clustree(all)
  
  ggsave("plots/overview-cluster-resolutions.pdf", useDingbats=F)
  
  all<-FindClusters(all,resolution=.3)
  
  
  DimPlot(all, label = TRUE) + NoLegend()
  DimPlot(all, label = TRUE, group.by = "Region") + NoLegend()
  
  #
  ElbowPlot(all)
  
  all<- all %>% 
    RunUMAP(dims=1:20) %>%
    FindNeighbors(dims=1:20) %>%
    FindClusters(resolution=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1))
  
  #url https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html#seurat-objects
  clustree(all)
  
  ggsave("plots/overview-cluster-resolutions.pdf", useDingbats=F)
  
  all<-FindClusters(all,resolution=.5)
  
  
  DimPlot(all, label = TRUE) + NoLegend()
  DimPlot(all, label = TRUE, group.by = "dataset") + NoLegend()
  
  save(all, file = file.path("data","seurat_object.RData"))
} else  {
  load(file.path("data", "seurat_object.RData"))
}

#plot
#violin plots
all2 <- all
Idents(all2) <- all2$dataset

genes <- toupper(c("Cx3cr1", "Tmem119", "Hexb", "SLC2A5", "OLFML3", "CSF1R", "P2ry12", "Rho", "Mid1", "Scoc","Luc7l3", "Prdx1","FAXC", "CALM2","RBFOX3", "PDGFRA", "MOG", "MAG", "CADPS2", "SYN3", "GAD2", "AQP4", "FGFR3", "MFGE8", "PRDX6", "SOX9", "SLC1A3", "CTPS1", "BMP4", "NEU4", "OPALIN","PLP1"))

violin <- all2[["SCT"]]@data[genes[genes %in% rownames(all2)],] %>%
  t() %>%
  as.data.frame() %>%
  mutate(Cluster=all2$dataset,
         ID=colnames(all2))

violin <- violin %>%
  pivot_longer(CX3CR1:PLP1 ,"Gene", values_to="Expression")

#line plots
gene_line_plot(violin, "MFGE8") +
  scale_color_brewer(palette = "Set1") +
  theme(text = element_text(size=20))

walk(genes, function(x) {
  tryCatch({
   p <- gene_line_plot(violin, x) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    theme(text = element_text(size=10))
  print(p)
  ggsave(file.path("plots", paste0(x,"_gene_line_plot.pdf")), useDingbats=F, width = 20, height = 2)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

walk(genes, function(x) {
  p1 <- ggplot(violin[violin$Gene == x,], aes(Cluster, Expression, fill=Cluster)) +
    geom_violin(scale = "width", draw_quantiles = .5) +
    scale_fill_brewer(palette = "Set1") +
    theme_void() +
    theme(text = element_text(size=10))
  print(p1)
  ggsave(file.path("plots", paste0(x,"_gene_violin_plot.pdf")), useDingbats=F, width = 10, height = 10)
})

all.markers <- FindAllMarkers(all2, only.pos=T)
