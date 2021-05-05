library(tidyverse)
library(Seurat)
library(data.table)

#initialize lists
hmn_umi<- lst()
ms_umi <- lst()

hmn_prop<- lst()
ms_prop <- lst()

#mouse data
mouse <- read_tsv(file.path("data", "GSE60361_Linnarson_Mouse.tsv"))

mouse2 <-  mouse %>%
  mutate(cell_type = str_replace_all(.$level1class, c("interneurons" = "Neuron", "pyramidal SS|pyramidal CA1" = "Neuron", "-ependymal" = "")),
         cell_type = str_to_sentence(cell_type),
         tissue = str_replace_all(tissue, c("ca1hippocampus" = "Hippocampus", "sscortex" = "Cortex"))
  ) %>%
  group_by(tissue, cell_type) %>%
  summarise(median_umi = round(median(`UMI Count`),2),
            sd = sd(`UMI Count`)) %>%
  arrange(median_umi)

mouse2 %>%
  ggplot(aes(reorder(cell_type, median_umi), median_umi, color = cell_type)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=median_umi-sd, ymax=median_umi+sd), width=0,size=1)+
  facet_wrap(~ tissue) +
  coord_flip() +
  scale_color_brewer("Cell type",palette = "Dark2") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(size = 15)) +
  labs(x=element_blank(), y="Median  UMI count") +
  geom_text(aes(label=median_umi), nudge_x=0.2, size=5, color="black")

ggsave(file.path("plots", "mouse_umis.pdf"), useDingbats=F)

#cell proportions
mouse3 <- mouse %>%
  mutate(cell_type = str_replace_all(.$level1class, c("interneurons" = "Neuron", "pyramidal SS|pyramidal CA1" = "Neuron", "-ependymal" = "")),
         cell_type = str_to_sentence(cell_type),
         tissue = str_replace_all(tissue, c("ca1hippocampus" = "Hippocampus", "sscortex" = "Cortex"))
  ) %>%
  group_by(tissue, cell_type) %>%
  summarise(rel_count = n()) %>%
  group_by(tissue) %>%
  mutate(rel_count = rel_count/sum(rel_count)*100)

a <- mouse3 %>%
  ggplot(aes(x=2, y=rel_count, fill=cell_type)) +
  geom_bar(position = 'fill', stat = 'identity', color='black', lwd=0.1) +
  coord_polar(theta='y', start=0) +
  theme_void() +
  facet_wrap(~tissue) +
  scale_fill_brewer(palette = "Dark2") +
  xlim(0.5, 2.5) 
print(a)

ms_umi[["linnarsson"]] <- mouse2
ms_prop[["linnarson"]] <- mouse3

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

ms_jakel2 <- ms_jakel %>%
  group_by(diagnosis, cell_type) %>%
  summarise(median_umi = round(median(UMIs),2),
            sd = sd(UMIs),
            median_norm_umi = round(median(norm_UMIs),2),
            sd_norm_umi = sd(norm_UMIs)) %>%
  arrange(median_umi)

ms_jakel2 %>%
  ggplot(aes(reorder(cell_type, median_umi), median_umi, color = cell_type)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=median_umi-sd, ymax=median_umi+sd), width=0,size=1)+
  facet_wrap(~ diagnosis) +
  coord_flip() +
  scale_color_brewer("Cell type",palette = "Dark2") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(size = 15)) +
  labs(x=element_blank(), y="Median  UMI count") +
  geom_text(aes(label=median_umi), nudge_x=0.2, size=5, color="black")

ms_jakel2 %>%
  ggplot(aes(reorder(cell_type, median_norm_umi), median_norm_umi, color = cell_type)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=median_norm_umi-sd_norm_umi, ymax=median_norm_umi+sd_norm_umi), width=0,size=1)+
  facet_wrap(~ diagnosis) +
  coord_flip() +
  scale_color_brewer("Cell type",palette = "Dark2") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(size = 15)) +
  labs(x=element_blank(), y="Median  UMI count") +
  geom_text(aes(label=median_norm_umi), nudge_x=0.2, size=5, color="black")

#cell proportions
ms_jakel3 <- ms_jakel %>%
  group_by(diagnosis,cell_type) %>%
  summarise(rel_count = n()) %>%
  group_by(diagnosis) %>%
  mutate(rel_count = rel_count/sum(rel_count)*100)

a <- ms_jakel3 %>%
  ggplot(aes(x=2, y=rel_count, fill=cell_type)) +
  geom_bar(position = 'fill', stat = 'identity', color='black', lwd=0.1) +
  coord_polar(theta='y', start=0) +
  theme_void() +
  facet_wrap(~diagnosis) +
  scale_fill_brewer(palette = "Dark2") +
  xlim(0.5, 2.5) 
print(a)

hmn_umi[["Jakel_et_al"]] <- ms_jakel2
hmn_prop[["Jakel_et_al"]] <- ms_jakel3

#schirmer
schirm <- read_tsv(file.path("data", "Schirmer_MS.txt")) %>%
  mutate(cell_type = str_replace_all(cell_type, c("EN-L2-3-Cntl|EN-L2-3-MS|EN-L4|EN-L5-6|EN-MIX|EN-PYR" = "Neuron", 
                                                  "IN-PV|IN-SST|IN-SV2C|IN-VIP" = "Neuron",
                                                  "OL-Cntl|OL-MS-1|OL-MS-2" = "Oligodendrocytes",
                                                  "Phagocytes" = "Macrophages",
                                                  "Endo cells" = "Endothelial-mural",
                                                  "OPC" = "OPCs"))) %>%
  filter(!cell_type %in% c("Stromal cells", "B cells", "T cells", "Glia-MIX", "Macrophages")) %>%
  mutate(norm_UMIs = 10000*(UMIs-min(UMIs))/(max(UMIs)- min(UMIs))) 

schirm2 <- schirm %>%
  group_by(diagnosis,cell_type) %>%
  summarise(median_umi = round(median(UMIs),2),
            sd = sd(UMIs),
            median_norm_umi = round(median(norm_UMIs),2),
            sd_norm_umi = sd(norm_UMIs)) %>%
  arrange(median_umi)

schirm2 %>%
  ggplot(aes(reorder(cell_type, median_umi), median_umi, color = cell_type)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=median_umi-sd, ymax=median_umi+sd), width=0,size=1)+
  facet_wrap(~ diagnosis) +
  coord_flip() +
  scale_color_brewer("Cell type",palette = "Dark2") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(size = 15)) +
  labs(x=element_blank(), y="Median  UMI count") +
  geom_text(aes(label=median_umi), nudge_x=0.2, size=5, color="black")

schirm2 %>%
  ggplot(aes(reorder(cell_type, median_norm_umi), median_norm_umi, color = cell_type)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=median_norm_umi-sd_norm_umi, ymax=median_norm_umi+sd_norm_umi), width=0,size=1)+
  facet_wrap(~ diagnosis) +
  coord_flip() +
  scale_color_brewer("Cell type",palette = "Dark2") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(size = 15)) +
  labs(x=element_blank(), y="Median  UMI count") +
  geom_text(aes(label=median_norm_umi), nudge_x=0.2, size=5, color="black")

#cell proportions
schirm3 <- schirm %>%
  group_by(diagnosis,cell_type) %>%
  summarise(rel_count = n()) %>%
  group_by(diagnosis) %>%
  mutate(rel_count = rel_count/sum(rel_count)*100)

a <- schirm3 %>%
  ggplot(aes(x=2, y=rel_count, fill=cell_type)) +
  geom_bar(position = 'fill', stat = 'identity', color='black', lwd=0.1) +
  coord_polar(theta='y', start=0) +
  theme_void() +
  facet_wrap(~diagnosis) +
  scale_fill_brewer(palette = "Dark2") +
  xlim(0.5, 2.5) 
print(a)

hmn_umi[["Schirmer_et_al."]] <- schirm2
hmn_prop[["Schirmer_et_al."]] <- schirm3

#schirmer
velm <- read_tsv(file.path("data", "Schirmer_Autism.txt")) %>%
  mutate(cell_type = str_replace_all(cluster, c("L2/3|L4|L5/6$|L5/6-CC$|Neu-mat|Neu-NRGN-I$|Neu-NRGN-II" = "Neuron", 
                                                  "IN-PV|IN-SST|IN-SV2C|IN-VIP" = "Neuron",
                                                  "AST-FB|AST-PP" = "Astrocytes",
                                                  "Phagocytes" = "Macrophages",
                                                  "Endothelial" = "Endothelial-mural",
                                                  "OPC" = "OPCs"))) %>%
  filter(!cell_type %in% c("Stromal cells", "B cells", "T cells", "Glia-MIX"))  %>%
  mutate(norm_UMIs = 10000*(UMIs-min(UMIs))/(max(UMIs)- min(UMIs))) 

velm2 <- velm %>%
  group_by(diagnosis,cell_type) %>%
  summarise(median_umi = round(median(UMIs),2),
            sd = sd(UMIs),
            median_norm_umi = round(median(norm_UMIs),2),
            sd_norm_umi = sd(norm_UMIs)) %>%
  arrange(median_umi)

velm2 %>%
  ggplot(aes(reorder(cell_type, median_umi), median_umi, color = cell_type)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=median_umi-sd, ymax=median_umi+sd), width=0,size=1)+
  facet_wrap(~ diagnosis) +
  coord_flip() +
  scale_color_brewer("Cell type",palette = "Dark2") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(size = 15)) +
  labs(x=element_blank(), y="Median  UMI count") +
  geom_text(aes(label=median_umi), nudge_x=0.2, size=5, color="black")

velm2 %>%
  ggplot(aes(reorder(cell_type, median_norm_umi), median_norm_umi, color = cell_type)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=median_norm_umi-sd_norm_umi, ymax=median_norm_umi+sd_norm_umi), width=0,size=1)+
  facet_wrap(~ diagnosis) +
  coord_flip() +
  scale_color_brewer("Cell type",palette = "Dark2") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(size = 15)) +
  labs(x=element_blank(), y="Median  UMI count") +
  geom_text(aes(label=median_norm_umi), nudge_x=0.2, size=5, color="black")


#cell proportions
velm3 <- velm %>%
  group_by(diagnosis,cell_type) %>%
  summarise(rel_count = n()) %>%
  group_by(diagnosis) %>%
  mutate(rel_count = rel_count/sum(rel_count)*100)

a <- velm3 %>%
  ggplot(aes(x=2, y=rel_count, fill=cell_type)) +
  geom_bar(position = 'fill', stat = 'identity', color='black', lwd=0.1) +
  coord_polar(theta='y', start=0) +
  theme_void() +
  facet_wrap(~diagnosis) +
  scale_fill_brewer(palette = "Dark2") +
  xlim(0.5, 2.5) 
print(a)

hmn_umi[["Velmeshev_et_al."]] <- velm2
hmn_prop[["Velmeshev_et_al."]] <- velm3

#bakken et al human
bakken_umi <- readRDS(file.path("data","allen_m1c_2019_ssv4_bakken_et_al_2020.rds"))[["RNA"]]@counts
bakken <- readRDS(file.path("data","allen_m1c_2019_ssv4_bakken_et_al_2020.rds"))@meta.data %>%
  na.omit %>% 
  mutate(
    UMIs = colSums(as.matrix(bakken_umi)[rownames(bakken_umi) %in% rownames(jakel_umis), rownames(bakken)]),
    norm_UMIs = 10000*(UMIs-min(UMIs))/(max(UMIs)- min(UMIs))
  )
bakken$cluster2 <- gsub(" .*", "", bakken$cluster)
bakken$cell_type <- str_replace_all(bakken$cluster2, c("Astro"="Astrocytes",
                                                                            "Endo|Peri|VLMC"="Endothelial-mural",
                                                                            "Exc|Inh"="Neuron",
                                                                            "Oligo"="Oligodendrocytes",
                                                                            "OPC"="OPCs",
                                                                            "Micro"="Microglia"))

bakken2 <- bakken %>%
  group_by(cell_type) %>%
  summarise(median_umi = round(median(UMIs),2),
            sd = sd(UMIs),
            median_norm_umi = round(median(norm_UMIs),2),
            sd_norm_umi = sd(norm_UMIs)) %>%
  arrange(median_umi)

bakken2 %>%
  ggplot(aes(reorder(cell_type, median_umi), median_umi, color = cell_type)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=median_umi-sd, ymax=median_umi+sd), width=0,size=1)+
  coord_flip() +
  scale_color_brewer("Cell type",palette = "Dark2") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(size = 15)) +
  labs(x=element_blank(), y="Median  UMI count") +
  geom_text(aes(label=median_umi), nudge_x=0.2, size=5, color="black")

bakken2 %>%
  ggplot(aes(reorder(cell_type, median_norm_umi), median_norm_umi, color = cell_type)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=median_norm_umi-sd_norm_umi, ymax=median_norm_umi+sd_norm_umi), width=0,size=1)+
  coord_flip() +
  scale_color_brewer("Cell type",palette = "Dark2") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(size = 15)) +
  labs(x=element_blank(), y="Median  UMI count") +
  geom_text(aes(label=median_norm_umi), nudge_x=0.2, size=5, color="black")

#cell proportions
bakken3 <- bakken %>%
  group_by(cell_type) %>%
  summarise(rel_count = n()) %>%
  mutate(rel_count = rel_count/sum(rel_count)*100)

a <- bakken3 %>%
  ggplot(aes(x=2, y=rel_count, fill=cell_type)) +
  geom_bar(position = 'fill', stat = 'identity', color='black', lwd=0.1) +
  coord_polar(theta='y', start=0) +
  theme_void() +
  scale_fill_brewer(palette = "Dark2") +
  xlim(0.5, 2.5) 
print(a)

#hmn_umi[["Bakken_et_al."]] <- bakken2
#hmn_prop[["Bakken_et_al."]] <- bakken3

#yao et al murine
yao <- readRDS(file.path("data","allen_mop_2020_yao_et_al_2020.rds"))@meta.data %>% 
  mutate(tissue = "Cortex")
yao$cluster2 <- gsub(".*_", "", yao$cluster)
yao$cell_type <- case_when(
  yao$cluster2 == "Astro" ~ "Astrocytes",
  yao$cluster2 == "Micro" ~ "Microglia",
  yao$cluster2 %in% c("SMC", "Peri", "VLMC", "Endo") ~ "Endothelial-mural",
  yao$cluster2 == "OPC" ~ "OPCs",
  yao$cluster2 == "Oligo" ~ "Oligodendrocytes",
  yao$cluster2 %in% c("PVM", "CR") ~ "Other",
  TRUE ~ "Neuron"
)
yao <- yao[yao$cell_type != "Other",]

yao2 <- yao %>%
  group_by(cell_type) %>%
  summarise(median_umi = round(median(nCount_RNA),2),
            sd = sd(nCount_RNA)) %>%
  arrange(median_umi)  %>% 
  mutate(tissue = "Cortex")

yao2 %>%
  ggplot(aes(reorder(cell_type, median_umi), median_umi, color = cell_type)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=median_umi-sd, ymax=median_umi+sd), width=0,size=1)+
  coord_flip() +
  scale_color_brewer("Cell type",palette = "Dark2") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(size = 15)) +
  labs(x=element_blank(), y="Median  UMI count") +
  geom_text(aes(label=median_umi), nudge_x=0.2, size=5, color="black")

ggsave(file.path("plots", "yao_umis.pdf"), useDingbats=F)

#cell proportions
yao3 <- yao %>%
  group_by(cell_type) %>%
  summarise(rel_count = n()) %>%
  mutate(rel_count = rel_count/sum(rel_count)*100)  %>% 
  mutate(tissue = "Cortex")

a <- yao3 %>%
  ggplot(aes(x=2, y=rel_count, fill=cell_type)) +
  geom_bar(position = 'fill', stat = 'identity', color='black', lwd=0.1) +
  coord_polar(theta='y', start=0) +
  theme_void() +
  scale_fill_brewer(palette = "Dark2") +
  xlim(0.5, 2.5) 
print(a)

ms_umi[["Yao et al."]] <- yao2
ms_prop[["Yao et al."]] <- yao3

#     Gerrits E, Heng Y, Boddeke EWGM, Eggen BJL. Transcriptional profiling of microglia; current state of the art and future perspectives. Glia 2020 Apr;68(4):740-755. PMID: 31846124
gerrits1 <- ReadMtx(mtx = file.path("data", "GSM4023552_Donor1nuclei_matrix.mtx.gz"), 
                   cells = file.path("data", "GSM4023552_Donor1nuclei_barcodes.tsv.gz"),
                   features = file.path("data", "GSM4023552_Donor1nuclei_genes.tsv.gz"))

gerrits1 <- as.matrix(gerrits1[rowSums(gerrits1) >0, colSums(gerrits1) >200]) %>%
  CreateSeuratObject()

gerrits2 <- ReadMtx(mtx = file.path("data", "GSM4023555_Donor2nuclei_matrix.mtx.gz"), 
                    cells = file.path("data", "GSM4023555_Donor2nuclei_barcodes.tsv.gz"),
                    features = file.path("data", "GSM4023555_Donor2nuclei_genes.tsv.gz"))
gerrits2 <- as.matrix(gerrits2[rowSums(gerrits2) >0, colSums(gerrits2) >200]) %>%
  CreateSeuratObject()

gerrits <- merge(gerrits1, gerrits2)
write_rds(gerrits[["RNA"]]@counts, file = file.path("data", "gerrits_et_al.rds"))

.gerrits <- gerrits %>%
  SCTransform() %>%
  RunPCA() %>%
  RunUMAP(dims=1:15) %>%
  FindNeighbors(dims=1:15) %>%
  FindClusters()

#plot all umis humans
hmn_umi2 <- hmn_umi %>% 
  bind_rows(.id="Dataset") %>%
  filter(diagnosis == "Control")

hmn_umi3 <- hmn_umi2 %>% 
  group_by(diagnosis, cell_type) %>% 
  summarise(mean_umi = round(mean(median_umi),2)) %>%
  filter(diagnosis == "Control")

hmn_umi3  %>%
  ggplot(aes(reorder(cell_type, mean_umi), mean_umi, group="Dataset")) +
  geom_point(data = hmn_umi2, aes(x=reorder(cell_type, median_umi), y=median_umi, fill = Dataset), size=10, pch=21, stroke=.25) +
  stat_summary(fun="mean", geom = "crossbar", width=0.5) +
  coord_flip() +
  scale_fill_brewer("Dataset",palette = "Set2") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(size = 20)) +
  labs(x=element_blank(), y="Mean  UMI count") +
  geom_text(aes(label=mean_umi), nudge_x=0.4, size=5, color="black") 
ggsave(file.path("plots", "umi_dotplot_human.pdf"), useDingbats=F)

#plot all props humans
hmn_prop2 <- hmn_prop %>% 
  bind_rows(.id="Dataset") %>%
  filter(diagnosis == "Control")

hmn_prop3 <- hmn_prop2 %>% 
  group_by(diagnosis, cell_type) %>% 
  summarise(mean_prop = round(mean(rel_count),2)) %>%
  filter(diagnosis == "Control")

hmn_prop3  %>%
  ggplot(aes(reorder(cell_type, mean_prop), mean_prop, group="Dataset")) +
  geom_point(data = hmn_prop2, aes(x=reorder(cell_type, rel_count), y=rel_count, fill = Dataset), size=10, pch=21, stroke=.25) +
  stat_summary(fun="mean", geom = "crossbar", width=0.5) +
  coord_flip() +
  scale_fill_brewer("Dataset",palette = "Set2") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(size = 20)) +
  labs(x=element_blank(), y="Mean % of celltype") +
  geom_text(aes(label=mean_prop), nudge_x=0.4, size=5, color="black") 
ggsave(file.path("plots", "prop_dotplot_human.pdf"), useDingbats=F)


#plot all mice
ms_umi2 <- ms_umi %>% 
  bind_rows(.id="Dataset") %>% 
  filter(tissue == "Cortex")

ms_umi3 <- ms_umi2 %>% 
  group_by(tissue, cell_type) %>% 
  summarise(mean_umi = round(mean(median_umi),2)) %>% 
  filter(tissue == "Cortex")

ms_umi3  %>%
  ggplot(aes(reorder(cell_type, mean_umi), mean_umi, group="Dataset")) +
  geom_point(data = ms_umi2, aes(x=reorder(cell_type, median_umi), y=median_umi, fill = Dataset), size=10, pch=21, stroke=.25) +
  stat_summary(fun="mean", geom = "crossbar", width=0.5) +
  coord_flip() +
  scale_fill_brewer("Dataset",palette = "Set1") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(size = 20)) +
  labs(x=element_blank(), y="Mean  UMI count") +
  geom_text(aes(label=mean_umi), nudge_x=0.4, size=5, color="black") 

ggsave(file.path("plots", "umi_dotplot_mouse.pdf"), useDingbats=F)


Freibur9
#plot all mice
ms_prop2 <- ms_prop %>% 
  bind_rows(.id="Dataset") %>% 
  filter(tissue == "Cortex") 

ms_prop3 <- ms_prop2 %>% 
  dplyr::group_by(Dataset,tissue, cell_type) %>% 
  summarise(mean_prop = round(mean(rel_count),2),
            rel_count = round(rel_count,1)) %>% 
  filter(tissue == "Cortex")

ms_prop3  %>%
  ggplot(aes(reorder(cell_type, mean_prop), mean_prop, group="Dataset")) +
  geom_point(data = ms_prop2, aes(x=reorder(cell_type, rel_count), y=rel_count, fill = Dataset), size=10, pch=21, stroke=.25) +
  geom_point(data = ms_prop2, aes(x=reorder(cell_type, rel_count), y=rel_count, fill = Dataset), size=15, pch="l", stroke=.1) +
  coord_flip() +
  scale_fill_brewer("Dataset",palette = "Set1") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(size = 20)) +
  facet_wrap(~ Dataset) +
  labs(x=element_blank(), y="Mean % of celltype") +
  geom_text(aes(label=rel_count), nudge_x=0.4, size=5, color="black") 

ggsave(file.path("plots", "prop_dotplot_mouse.pdf"), useDingbats=F)

#plot yao et al
ms_prop2 <- ms_prop %>% 
  bind_rows(.id="Dataset") %>% 
  filter(Dataset == "Yao et al.") 

ms_prop3 <- ms_prop2 %>% 
  group_by(cell_type) %>% 
  summarise(mean_prop = round(mean(rel_count),2)) 

ms_prop3  %>%
  ggplot(aes(reorder(cell_type, mean_prop), mean_prop, group="Dataset")) +
  geom_point(data = ms_prop2, aes(x=reorder(cell_type, rel_count), y=rel_count, fill = Dataset), size=10, pch=21, stroke=.25) +
  stat_summary(fun="mean", geom = "crossbar", width=0.5) +
  coord_flip() +
  scale_fill_brewer("Dataset",palette = "Set1") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(size = 20)) +
  labs(x=element_blank(), y="Mean % of celltype") +
  geom_text(aes(label=mean_prop), nudge_x=0.4, size=5, color="black") 

ggsave(file.path("plots", "prop_dotplot_yao_et_al_mouse.pdf"), useDingbats=F)

  
  