#####################################################################
#Titulo: Pipeline para Ilumina HumanHT-12 V4.0 expression beadchip
#GPL10558
#Autor: Jose Antonio Ovando Ricardez
#Fecha: Julio/2021

#####################################################################
# Instalacion de paqueterias
#####################################################################


#####################################################################
# Llamado de librerias
#####################################################################

library(limma)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)

#####################################################################
#Preprocesamiento de los datos (SOLO SI NO TENEMOS LOS ARCHIVOS IDATR)
#####################################################################

# En este caso como no tenemos los archivos IDAT los vamos a
# Extraer directamente del formato .txt

rawData <- read.table(
  paste0('data/GSE35145_non-normalized_expression.txt'),
  header = TRUE, sep = '\t', stringsAsFactors = FALSE, skip = 0)


# Extraemos las columnas con el P-Val del objeto x
detectionpvalues <- rawData[,grep('Detection.Pval', colnames(rawData))]
rawData <- rawData[,-grep('Detection.Pval', colnames(rawData))]

# set rownames and tidy up final expression matrix
probes <- rawData$ID_REF
rawData <- data.matrix(rawData[,2:ncol(rawData)])
rownames(rawData) <- probes
colnames(rawData) <- gsub('^X', '',
                    gsub('\\.AVG_Signal', '', colnames(rawData)))

#Manualmente renombramos cada una de las columnas con el nombre de
#nuestras muestras

colnames(rawData)[1:8] <- c("IPF1","IPF2","IPF3","IPF4","Control1",
                      "Control2","Control3","Control4")

#####################################################################
#Pre-Anotacion
#####################################################################

#Seleccionamos nuestro archivo de anotacion
bgxfile <- 'data/GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt'

#Asignamos a nuestros datos las anotaciones
annot <- illuminaio::readBGX(bgxfile)$probes
annot <- annot[,which(colnames(annot) %in% c('Source','Symbol','Transcript','ILMN_Gene','RefSeq_ID',
                                             'Entrez_Gene_ID','Symbol','Protein_Product','Probe_Id','Probe_Type',
                                             'Probe_Start','Chromosome','Probe_Chr_Orientation','Probe_Coordinates',
                                             'Cytoband', 'Definition', 'Ontology_Component', 'Ontology_Process',
                                             'Ontology_Function', 'Synonyms'))]
annot <- annot[which(annot$Probe_Id %in% rownames(rawData)),]
annot <- annot[match(rownames(rawData), annot$Probe_Id),]


#####################################################################
#Procesamiento de los datos y TargetsFile
#####################################################################

#Seleccionamos nuestro archivo targets en formatro .txt
#Este archivo lo cree usando Excel

targetsfile <- 'Phenodata_IPF_VS_CONTROL.txt'
targetinfo <- readTargets(targetsfile, sep = '\t')
rownames(targetinfo) <- gsub('\\.idat$', '',
                             gsub('^raw/', '', targetinfo$IDATfile))
rawData <- rawData[,match(rownames(targetinfo), colnames(rawData))]
if (!all(colnames(rawData) == rownames(targetinfo)))

  
# Creamos un objeto que almacene las intensidades (E)
project <- new('EListRaw')
project@.Data[[1]] <- 'illumina'
project@.Data[[2]] <- targetinfo
project@.Data[[3]] <- NULL
project@.Data[[4]] <- rawData
project@.Data[[5]] <- NULL
project$E <- rawData
project$targets <- targetinfo
project$genes <- NULL
project$other$Detection <- detectionpvalues  


#####################################################################
#Analisis de calidad
#####################################################################

#Boxplot
par(mar=c(8,8,5,5), cex=0.6, cex.axis=1.0, cex.lab=1.3)

boxplot(project$E,
        main="Boxplot of log2-intensitites for the Raw data",
        xlab="", ylab=bquote(~Log[2]~expression),
        names= targetinfo$Sample,
        col = c("brown1", "brown1","brown1","brown1","darkslategray3"
                ,"darkslategray3","darkslategray3","darkslategray3"),
        las=2,
        boxwex = 0.6,
        staplewex = 0.6,
        outline=F)

#####################################################################
#Normalizacion
#####################################################################

project.bgcorrect.norm <- neqc(project, offset = 16)

#Boxplot
par(mar=c(8,8,5,5), cex=0.6, cex.axis=1.0, cex.lab=1.3)

boxplot(project.bgcorrect.norm$E,
        main="Boxplot of log2-intensitites for the Raw data",
        xlab="", ylab=bquote(~Log[2]~expression),
        names= targetinfo$Sample,
        col = c("brown1", "brown1","brown1","brown1","darkslategray3"
                ,"darkslategray3","darkslategray3","darkslategray3"),
        las=2,
        boxwex = 0.6,
        staplewex = 0.6,
        outline=F)

#####################################################################
# PCA (NormData)
#####################################################################

PCA <- prcomp(t(project.bgcorrect.norm$E), scale = F)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Disease =
                       Biobase::pData(targetinfo)$Group)

ggplot(dataGG, aes(PC1, PC2, label=row.names(dataGG))) +
  geom_point(aes(colour = targetinfo$Group))+geom_text(size=1) +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = c("#fa1d19", "#1937fa"))


#####################################################################
# Anotacion y filtrado de controles
#####################################################################

nnot <- annot[which(annot$Probe_Id %in% rownames(project.bgcorrect.norm)),]
project.bgcorrect.norm <- project.bgcorrect.norm[which(rownames(project.bgcorrect.norm) %in% annot$Probe_Id),]
annot <- annot[match(rownames(project.bgcorrect.norm), annot$Probe_Id),]
project.bgcorrect.norm@.Data[[3]] <- annot
project.bgcorrect.norm$genes <- annot


Control <- project.bgcorrect.norm$genes$Source=="ILMN_Controls"
NoSymbol <- project.bgcorrect.norm$genes$Symbol == ""
isexpr <- rowSums(project.bgcorrect.norm$other$Detection <= 0.05) >= 3

project.bgcorrect.norm.filt <- project.bgcorrect.norm[!Control & !NoSymbol & isexpr, ]
dim(project.bgcorrect.norm)
dim(project.bgcorrect.norm.filt)

#Removemos las columnas que no nos sirven

project.bgcorrect.norm.filt$genes <- project.bgcorrect.norm.filt$genes[,c(
  'Probe_Id',
  'Definition','Ontology_Component','Ontology_Process','Ontology_Function',
  'Chromosome','Probe_Coordinates','Cytoband','Probe_Chr_Orientation',
  'RefSeq_ID','Entrez_Gene_ID','Symbol')]

head(project.bgcorrect.norm.filt$genes)

#Creamos la tabla final con la anotación
AnnotFinal = data.frame(project.bgcorrect.norm.filt$genes$Probe_Id,
                        project.bgcorrect.norm.filt$genes$Symbol,
                        project.bgcorrect.norm.filt$E)

#####################################################################
#Expresion diferencial
#####################################################################

casos <- as.factor(targetinfo$Group)
design = model.matrix(~ 0+casos)
colnames(design) = levels(casos)
design

#Hacemos la matriz de contrastes
cont.matrix<- makeContrasts(IPF-Control, levels=design)
cont.matrix

fit<-lmFit(AnnotFinal[3:10],design)
fit2= contrasts.fit(fit,cont.matrix)
fit3= eBayes(fit2)
fit3
head(fit3)

#Guardar tabla
table_CD = topTable(fit3, coef=1, number=nrow(fit), 
                    sort.by= "none",adjust="fdr")

FinalTable = data.frame(table_CD, AnnotFinal)

write.table(FinalTable, file="DEG_GSE35145.txt", sep="\t", 
            row.names=F, col.names=TRUE, quote=FALSE)

#####################################################################
#Volcano
#####################################################################
EnhancedVolcano(FinalTable,lab = FinalTable$"project.bgcorrect.norm.filt.genes.Symbol",
                x = 'logFC',
                y = 'adj.P.Val',
                ylab = bquote(~adj.P.Val),
                xlab = bquote(~logFC),
                ylim = c(0, 2),
                xlim = c(-5, 5),
                pCutoff = 0.05,
                FCcutoff = 1,
                col=c('#c2bcbc', '#a6a2a2', '#f59289', '#eb4334'),
                cutoffLineType = 'solid',
                cutoffLineCol = 'coral4',
                hlineWidth = 0.2,
                cutoffLineWidth = 0.2,
                axisLabSize = 8,
                pointSize = 1.0,
                colAlpha = 1,
                legendPosition = 'bottom',
                legendLabels = c("NS", expression(LogFC), "adj.P.Val", expression(adj.P.Val ~ and
                                                                                  ~ logFC)),
                legendLabSize = 8,
                title = "          Volcano plot of differential expresion of genes",
                titleLabSize = 13,
                subtitleLabSize = 9,
                captionLabSize = 9,
                border = 'full',
                borderWidth = 0.5,
                borderColour = 'black',
                caption = paste0("Total = ", nrow(FinalTable), " genes"),
                subtitle = "Comparation between IPF vs Control                                        pCutoff = 0.05,     FCcutoff = 1 ",
                labSize = 1.0)

#####################################################################
#Seleccionamos el resultado de la expresion diferencial
#####################################################################

#Seleccionamos genes con valores p-values por debajo 
#del umbral

topgenes = FinalTable[FinalTable[, "adj.P.Val"] < 0.05, ]

dim(topgenes)
head(topgenes)

write.table(topgenes, file="TopGenes_GSE35145.txt", sep="\t", 
            row.names=FALSE, col.names=TRUE, quote=FALSE)

#Distinguimos entre los genes up and down regulated
topups = topgenes[topgenes[, "logFC"] > 1, ]
dim(topups)

write.table(topgenes, file="TopUps_GSE35145.txt", sep="\t", 
            row.names=FALSE, col.names=TRUE, quote=FALSE)

topdowns = topgenes[topgenes[, "logFC"] < -1, ]
dim(topdowns)

write.table(topgenes, file="TopDowns_GSE35145.txt", sep="\t", 
            row.names=FALSE, col.names=TRUE, quote=FALSE)

#####################################################################
#Heatmap
#####################################################################

#Seleccionar los datos
glimpse(FinalTable)
data_filtered <- FinalTable %>% filter(adj.P.Val < 0.05 & (logFC > 1 | logFC < -1 ))
nrow(data_filtered)

#Cálculo de valor Z
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
dataZ <- t(apply(data_filtered[,9:16], 1, cal_z_score))
data_filtered2 <- data.frame(data_filtered[,8], dataZ[,5:8], dataZ[,1:4])
data_filtered3 <- data_filtered2 %>% drop_na

#Seleccionar de color
col_fun <- colorRamp2(seq(min(data_filtered3[,2:5]), max(data_filtered3[,6:9]), length = 3), c("#f1a340", "#f7f7f7", "#998ec3"))
col_fun

#Opciones de guardado
calc_ht_size = function(ht, unit = "inch") {
  pdf()
  ht = draw(ht)
  w = ComplexHeatmap:::width(ht)
  w = convertX(w, unit, valueOnly = TRUE)
  h = ComplexHeatmap:::height(ht)
  h = convertY(h, unit, valueOnly = TRUE)
  dev.on()
  c(w, h)
}

#Legends
lgd <- Legend(col_fun = col_fun, title = "Row Z-Score")

#Heatmap without genes
ht <- Heatmap(as.matrix(data_filtered3[,2:9]),
              name = "Z-Score", column_title = "Differential gene expression heatmap (GSE35145)",  column_title_gp = gpar(fontsize = 17, fontface = "bold"),
              col = col_fun,
              column_order = order(as.numeric(gsub("column", "", colnames(data_filtered3[,2:9])))),
              clustering_distance_rows = "euclidean",
              row_names_gp = gpar(fontsize = 0),
              column_names_gp = gpar(fontsize = 5, fontface = "bold"),
              clustering_method_rows = "ward.D2",
              width = unit(13, "cm"),
              height = unit(15, "cm"))

#Hetmap con genes

ht2 <- Heatmap(as.matrix(data_filtered3[,2:9]),
               name = "Z-Score", column_title = "Differential gene expression heatmap (GSE35145)",  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
               col = col_fun,
               column_order = order(as.numeric(gsub("column", "", colnames(data_filtered3[,2:9])))),
               clustering_distance_rows = "euclidean",
               clustering_method_rows = "ward.D2",
               row_names_max_width = max_text_width(
                 rownames(data_filtered3$data_filtered...8.),
                 gp = gpar(fontsize = 12)),
               row_labels = data_filtered3$data_filtered...8.,
               width = unit(13, "cm"),
               height = unit(16, "cm"),
               row_names_gp = gpar(fontsize = 0.5, fontface = "bold"),
               column_names_gp = gpar(fontsize = 5, fontface = "bold"),
               show_heatmap_legend = TRUE)

#Plot PDF
size <- calc_ht_size(ht2)
size
pdf("test.pdf", width = size[1], height = size[2])
ht2
dev.off()
