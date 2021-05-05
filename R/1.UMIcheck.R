setwd("/home/ngsnpath/Desktop/UMI_microglia/")
setwd("/Users/gianni/Desktop/Temp/2021/reviewRoman/UMI_microglia/")

library(ggplot2)

LinnarsonMouse <- read.delim("GSE60361_Linnarson_Mouse.tsv")
SchirmerAutism <- read.delim("Schirmer_Autism.txt")
JakelMS <- read.delim("GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt")
SchirmerMS <- read.delim("Schirmer_MS.txt")
# Muller_GBM <- read.delim("Muller_GBM.tsv")
# Tanaka <- read.delim("Tanaka_organoids.tsv")


colnames(LinnarsonMouse)[colnames(LinnarsonMouse) == "level1class"] <- "Cell_Type"
ggplot(LinnarsonMouse, aes(x=level1class, y=UMI.Count, fill=level1class)) +
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + facet_wrap(. ~ tissue,scales = "free") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())


colnames(SchirmerAutism)[colnames(SchirmerAutism) == "cluster"] <- "Cell_Type"
ggplot(SchirmerAutism, aes(x=Cell_Type, y=UMIs, fill=Cell_Type)) +
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) +# facet_wrap(. ~ region,scales = "free") # +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

colnames(SchirmerMS)[colnames(SchirmerMS) == "cell_type"] <- "Cell_Type"
ggplot(SchirmerMS, aes(x=Cell_Type, y=UMIs, fill=Cell_Type)) +
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) +# facet_wrap(. ~ region,scales = "free") # +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())


colnames(JakelMS)[colnames(JakelMS) == "Celltypes"] <- "Cell_Type"
ggplot(JakelMS, aes(x=Cell_Type, y=genes, fill=Cell_Type)) +
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) +# facet_wrap(. ~ region,scales = "free") # +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())


 
# colnames(Muller_GBM)[colnames(Muller_GBM) == "Cell.Type.Assignment"] <- "Cell_Type"
# ggplot(Muller_GBM, aes(x=Cell_Type, y=nUMI, fill=Cell_Type)) +
#   geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) +# facet_wrap(. ~ region,scales = "free") # +
#   theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())


# colnames(Tanaka)[colnames(Tanaka) == "Annotation"] <- "Cell_Type"
# Tanaka$UMIsLog <- log2(Tanaka$nCount_RNA +1)
# ggplot(Tanaka, aes(x=Cell_Type, y=UMIsLog, fill=Cell_Type)) +
#   geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) +# facet_wrap(. ~ region,scales = "free") # +
#    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
