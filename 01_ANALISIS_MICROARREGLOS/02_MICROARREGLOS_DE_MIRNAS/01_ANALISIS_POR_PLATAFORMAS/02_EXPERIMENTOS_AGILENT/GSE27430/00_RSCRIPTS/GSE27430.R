#####################################################################
#Titulo: Pipeline para Agilent-019118 Human miRNA Microarray 
#2.0 G4470B (Probe Name version)
#Autor: Jose Antonio Ovando Ricardez
#Fecha: Junio/2021

#####################################################################
# Instalacion de paqueterias
#####################################################################


#####################################################################
# Llamado de librerias
#####################################################################

library(GEOquery)
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
# Abrimos la tabla con el PhenoData (targets)
# Esta tabla la hice fuera de R (En excel)
# Se tomararon los datos directamente del GEO2R
#####################################################################

targetinfo <- readTargets("./02_PHENODATA/Phenodata_IPF_VS_CONTROL.txt",
                          row.names="Sample",
                          sep="")

#####################################################################
# Leemos cada uno de los archivos y especificamos
# en formato .txt
#####################################################################

project <- read.maimages(targetinfo$Sample,
                         path="./01_DATOS/DATOS", 
                         source="agilent",
                         green.only=TRUE)

#####################################################################
# Control de calidad 
#####################################################################
# Boxplot
#####################################################################

pdf("./04_ANALISIS_DE_CALIDAD/BP_RAWDATA_GSE27430.pdf")
boxplot(project$E,
        main="Boxplot of log2-intensitites for the Raw data",
        xlab="", ylab=bquote(~Log[2]~expression),
        names= targetinfo$Filename,
        col = c("brown1", "brown1","brown1","brown1","brown1","brown1",
                "brown1", "brown1","brown1","brown1","brown1","brown1","darkslategray3",
                "darkslategray3", "darkslategray3","darkslategray3","darkslategray3", "darkslategray3",
                "darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3",
                "darkslategray3"),
        las=2,
        boxwex = 0.6,
        staplewex = 0.6,
        outline=F)
dev.off()

#####################################################################
# Coreccion de fondo (Background corecction)
#####################################################################

project.bgc <- backgroundCorrect(project, method="normexp", offset=16)


#####################################################################
# Boxplot (Background corecction)
#####################################################################

pdf("./04_ANALISIS_DE_CALIDAD/BP_BC_GSE27430.pdf")
boxplot(project.bgc$E,
        main="Boxplot of log2-intensitites for the BGC data",
        xlab="", ylab=bquote(~Log[2]~expression),
        names= targetinfo$Filename,
        col = c("brown1", "brown1","brown1","brown1","brown1","brown1",
                "brown1", "brown1","brown1","brown1","brown1","brown1","darkslategray3",
                "darkslategray3", "darkslategray3","darkslategray3","darkslategray3", "darkslategray3",
                "darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3",
                "darkslategray3"),
        las=2,
        boxwex = 0.6,
        staplewex = 0.6,
        outline=F)
dev.off()

#####################################################################
# Normalizacion
#####################################################################

project.NormData <-normalizeBetweenArrays(project.bgc,method=
                                            "cyclicloess")

#####################################################################
# Boxplot (Background corecction)
#####################################################################

pdf("./04_ANALISIS_DE_CALIDAD/BP_NORM_DATA_GSE27430.pdf")
boxplot(project.NormData$E,
        main="Boxplot of log2-intensitites for the BGC data",
        xlab="", ylab=bquote(~Log[2]~expression),
        names= targetinfo$Filename,
        col = c("brown1", "brown1","brown1","brown1","brown1","brown1",
                "brown1", "brown1","brown1","brown1","brown1","brown1","darkslategray3",
                "darkslategray3", "darkslategray3","darkslategray3","darkslategray3", "darkslategray3",
                "darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3",
                "darkslategray3"),
        las=2,
        boxwex = 0.6,
        staplewex = 0.6,
        outline=F)
dev.off()

#####################################################################
# PCA (NormData)
#####################################################################

PCA <- prcomp(t(project.NormData$E), scale = F)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Disease =
                       Biobase::pData(targets)$Group)

pdf("./04_ANALISIS_DE_CALIDAD/PCA_GSE27430.pdf")
ggplot(dataGG, aes(PC1, PC2, label=row.names(dataGG))) +
  geom_text(size=1)+geom_point(aes(colour = targetinfo$Group)) +
  ggtitle("GSE27430 | PCA plot of the log-transformed norm expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.2))+
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = c("#0000ff", "#fb0007")) + theme(plot.title = element_text(face = "bold"))
dev.off()

#####################################################################
# Anotacion
#####################################################################

#Obtenemos los datos de anotacion dek  project.NormData
Annot1 <- project.NormData$genes
Annot2 <- Annot1[1:13737, c(3, 4, 5)]
geneMatrix <- cbind(Annot2, project.NormData$E)

#Filtramos las sondas de control y conservamos solo las que
#tienen un valor de 0
Annot3 <- (geneMatrix[geneMatrix$'ControlType' == "0",])

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

fit<-lmFit(Annot3[4:26],design)
fit2= contrasts.fit(fit,cont.matrix)
fit3= eBayes(fit2)

#Guardar tabla
table_CD = topTable(fit3, coef=1, number=nrow(fit), sort.by= "none",adjust="fdr")

#Creamos la tabla final
FinalTable <- cbind(Annot3[2:3], table_CD, Annot3[4:26])

write.table(FinalTable, file="./05_TABLA_DE_EXPRESION_DIFERENCIAL/DEG_GSE27430.txt", sep="\t",
            row.names=F, col.names=T, quote=F)


#####################################################################
#Volcano
#####################################################################

#Valores de cada color
keyvals <- ifelse(FinalTable$logFC <= -0.5 & FinalTable$adj.P.Val <0.05,
                  '#0000ff',  ifelse(FinalTable$logFC >=0.5 & FinalTable$adj.P.Val <0.05,
                                     '#fb0007', '#e3deeb'))
#Nombre de valores de cada color
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == '#fb0007'] <- 'Up'
names(keyvals)[keyvals == '#e3deeb'] <- 'NS'
names(keyvals)[keyvals == '#0000ff'] <- 'Down'

#Generamos el Volcanoplot
pdf("./06_GRAFICOS_DE_EXPRESION_DIFERENCIAL/Volcano_GSE27430.pdf")
EnhancedVolcano(FinalTable,lab = FinalTable$SystematicName,
                x = 'logFC',
                y = 'adj.P.Val',
                ylab = bquote(~adj.P.Val),
                xlab = bquote(~logFC),
                ylim = c(0, 3),
                xlim = c(-1.5, 1.5),
                pCutoff = 0.05,
                FCcutoff = 0.5,
                colCustom = keyvals,
                cutoffLineType = 'solid',
                cutoffLineCol = 'coral4',
                hlineWidth = 0.2,
                cutoffLineWidth = 0.2,
                axisLabSize = 8,
                pointSize = 1.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabels = c("NS", expression(LogFC), "adj.P.Val", expression(adj.P.Val ~ and
                                                                                  ~ logFC)),
                legendLabSize = 8,
                title = "        GSE27430 | Volcano plot of differential expression",
                titleLabSize = 13,
                subtitleLabSize = 9,
                captionLabSize = 9,
                border = 'full',
                borderWidth = 0.5,
                borderColour = 'black',
                caption = NULL,
                subtitle = "Comparation between IPF vs Control                           pCutoff = 0.05,     FCcutoff = 0.5",
                labSize = 1.0)
dev.off()


##########################################################################################################################################
#Filtrado de datos y resulados
##########################################################################################################################################

#Aqui estableceremos nuestros valores de corte, en este caso usamos 0.05 padj y 0.8 de logFC

##########################################################################################################################################

#Filtrado de datos
data_filtered <- FinalTable %>% filter(adj.P.Val < 0.05 & (logFC > 0.5 | logFC < -0.5 ))
data_filtered <- data_filtered %>% drop_na

#Cantidad de genes expresados diferencialmente
nrow(data_filtered)

#Exportamos la lista filtrada
write.table(data_filtered, file="./07_RESULTADOS/DEG_FILTER_GSE27430.txt", sep="\t", row.names=F, col.names=T, quote=F)

#####################################################################
#Heatmap
#####################################################################

#Calculo de valor Z
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
dataZ <- t(apply(data_filtered[,9:31], 1, cal_z_score))
data_filtered2 <- cbind(data_filtered[,1:8], dataZ[,13:23], dataZ[,1:12])
data_filtered3 <- data_filtered2 %>% drop_na

#Seleccionar color
col_fun <- colorRamp2(seq(min(data_filtered3[,9:31]), max(data_filtered3[,9:31]), length = 3), c("#0000ff", "white", "#fb0007"))
col_fun

rwb <- colorRampPalette(colors = c("blue", "white", "red"))(30)

#Legends
lgd <- Legend(col_fun = col_fun, title = "Row Z-Score")

#Heatmap sin genes
pdf("./06_GRAFICOS_DE_EXPRESION_DIFERENCIAL/Heatmap_GSE27430.pdf")
Heatmap(as.matrix(data_filtered3[,9:31]),
        name = "Z-Score", column_title = "GSE27430 | Differential gene expression heatmap",  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
        col = rwb,
        column_order = order(as.numeric(gsub("column", "", colnames(data_filtered3[,9:31])))),
        clustering_distance_rows = "euclidean",
        row_names_gp = gpar(fontsize = 0),
        column_names_gp = gpar(fontsize = 3, fontface = "bold"),
        clustering_method_rows = "ward.D2",
        width = unit(9, "cm"),
        height = unit(15, "cm"))
dev.off()

#####################################################################
#End
#####################################################################
