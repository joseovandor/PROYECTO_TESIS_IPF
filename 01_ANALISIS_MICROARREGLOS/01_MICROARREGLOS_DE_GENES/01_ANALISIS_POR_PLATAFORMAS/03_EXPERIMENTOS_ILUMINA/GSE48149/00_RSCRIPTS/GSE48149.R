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
library(GEOquery)

#####################################################################
#Preprocesamiento de los datos (SOLO SI NO TENEMOS LOS ARCHIVOS IDATR)
#####################################################################

# En este caso como no tenemos los archivos IDAT los vamos a
# Extraer directamente del formato .xls

rawData <- read.table(
  paste0('./01_DATOS/GSE48149.txt'),
  header = TRUE, sep = '\t', stringsAsFactors = F, skip = 0)


# Extraemos las columnas con el P-Val del objeto
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

rawData1 = data.frame(rawData[,1:22])

colnames(rawData1)[1:22] <- c("IPF1","IPF2","IPF3","IPF4","IPF5",
                              "IPF6","IPF7","IPF8", "IPF9", "IPF10"
                              , "IPF11", "IPF12", "IPF13", "Control1"
                              , "Control2", "Control3", "Control4"
                              , "Control5", "Control6", "Control7"
                              , "Control8", "Control9")

#Estas muestras fueron eliminadas por presentar una mala calidad
#en el PCA
rawData1 <- rawData1[,-(3:4)]

#Volvemos a renombrar
colnames(rawData1)[1:20] <- c("IPF1","IPF2","IPF3","IPF4","IPF5",
                              "IPF6","IPF7","IPF8", "IPF9", "IPF10"
                              , "IPF11", "Control1"
                              , "Control2", "Control3", "Control4"
                              , "Control5", "Control6", "Control7"
                              , "Control8", "Control9")


#####################################################################
#Procesamiento de los datos y TargetsFile
#####################################################################

#Seleccionamos nuestro archivo targets en formatro .txt
#Este archivo lo cree usando Excel

targetsfile <- './02_PHENODATA/Phenodata_IPF_VS_CONTROL.txt'
targetinfo <- readTargets(targetsfile, sep = '\t')
rownames(targetinfo) <- gsub('\\.idat$', '',
                             gsub('^raw/', '', targetinfo$IDATfile))
rawData1 <- rawData1[,match(rownames(targetinfo), colnames(rawData1))]
if (!all(colnames(rawData1) == rownames(targetinfo)))
  

# Creamos un objeto que almacene las intensidades (E)
# Correr linea por linea o da error
  
project <- new('EListRaw')
project@.Data[[1]] <- 'illumina'
project@.Data[[2]] <- targetinfo
project@.Data[[3]] <- NULL
project@.Data[[4]] <- rawData1
project@.Data[[5]] <- NULL
project$E <- rawData1
project$targets <- targetinfo
project$genes <- NULL
project$other$Detection <- detectionpvalues[,1:22]
project$other$Detection <- project$other$Detection[,-(3:4)]

#####################################################################
#Analisis de calidad
#####################################################################

#Boxplot
pdf("./04_ANALISIS_DE_CALIDAD/BP_RAWDATA_GSE48149.pdf")
boxplot(project$E,
        main="Boxplot of log2-intensitites for the Raw data",
        xlab="", ylab=bquote(~Log[2]~expression),
        names= targetinfo$Sample,
        col = c("brown1", "brown1","brown1","brown1","brown1","brown1"
                ,"brown1","brown1","brown1","brown1",
                "brown1","darkslategray3","darkslategray3",
                "darkslategray3","darkslategray3","darkslategray3"
                ,"darkslategray3","darkslategray3","darkslategray3"
                ,"darkslategray3"),
        las=2,
        boxwex = 0.6,
        staplewex = 0.6,
        outline=F)
dev.off()

#####################################################################
#Normalizacion
#####################################################################

project.bgcorrect.norm <- neqc(project, offset = 16)

#Boxplot
pdf("./04_ANALISIS_DE_CALIDAD/BP_NORM_DATA_GSE48149.pdf")
boxplot(project.bgcorrect.norm$E,
        main="Boxplot of log2-intensitites for the Raw data",
        xlab="", ylab=bquote(~Log[2]~expression),
        names= targetinfo$Sample,
        col = c("brown1", "brown1","brown1","brown1","brown1","brown1"
                ,"brown1","brown1","brown1","brown1",
                "brown1","darkslategray3","darkslategray3",
                "darkslategray3","darkslategray3","darkslategray3"
                ,"darkslategray3","darkslategray3","darkslategray3"
                ,"darkslategray3"),
        las=2,
        boxwex = 0.6,
        staplewex = 0.6,
        outline=F)
dev.off()

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

pdf("./04_ANALISIS_DE_CALIDAD/PCA_GSE48149.pdf")
ggplot(dataGG, aes(PC1, PC2, label=row.names(dataGG))) +
  geom_text(size=2)+geom_point(aes(colour = targetinfo$Group)) +
  ggtitle("GSE48149 | PCA plot of the log-transformed norm expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.2))+
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = c("#0000ff", "#fb0007")) + theme(plot.title = element_text(face = "bold"))
dev.off()

#####################################################################
# Anotacion
#####################################################################

gpl <- getGEO('GPL16221')
table <- Table(gpl)[1:22304, c(1, 7, 12)]

normData = data.frame(project.bgcorrect.norm$E)

normData <- cbind(rownames(normData), data.frame(normData, row.names=NULL))
colnames(normData)[1] <- "ID"

AnnotFinal <- merge(table, normData, by.x = "ID", by.y = "ID")


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

fit<-lmFit(AnnotFinal[4:23],design)
fit2= contrasts.fit(fit,cont.matrix)
fit3= eBayes(fit2)
fit3
head(fit3)

#Guardar tabla
table_CD = topTable(fit3, coef=1, number=nrow(fit), sort.by= "none",adjust="fdr")

#Creamos la tabla final
FinalTable = data.frame(AnnotFinal[,1:3], table_CD, AnnotFinal[,4:23])

write.table(FinalTable, file="./05_TABLA_DE_EXPRESION_DIFERENCIAL/DEG_GSE48149.txt", sep="\t",
            row.names=F, col.names=T, quote=F)

#####################################################################
#Volcano
#####################################################################

#Valores de cada color
keyvals <- ifelse(FinalTable$logFC <= -0.8 & FinalTable$adj.P.Val <0.05,
                  '#0000ff',  ifelse(FinalTable$logFC >=0.8 & FinalTable$adj.P.Val <0.05,
                                     '#fb0007', '#e3deeb'))
#Nombre de valores de cada color
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == '#fb0007'] <- 'Up'
names(keyvals)[keyvals == '#e3deeb'] <- 'NS'
names(keyvals)[keyvals == '#0000ff'] <- 'Down'

#Generamos el Volcanoplot
pdf("./06_GRAFICOS_DE_EXPRESION_DIFERENCIAL/Volcano_GSE48149.pdf")
EnhancedVolcano(FinalTable,lab = FinalTable$Symbol,
                x = 'logFC',
                y = 'adj.P.Val',
                ylab = bquote(~adj.P.Val),
                xlab = bquote(~logFC),
                ylim = c(0, 6),
                xlim = c(-4, 4),
                pCutoff = 0.05,
                FCcutoff = 0.8,
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
                title = "     GSE48149 | Volcano plot of differential expression",
                titleLabSize = 13,
                subtitleLabSize = 9,
                captionLabSize = 9,
                border = 'full',
                borderWidth = 0.5,
                borderColour = 'black',
                caption = NULL,
                subtitle = "Comparation between IPF vs Control                           pCutoff = 0.05,     FCcutoff = 0.8",
                labSize = 1.0)
dev.off()

##########################################################################################################################################
#Filtrado de datos y resulados
##########################################################################################################################################

#Aqui estableceremos nuestros valores de corte, en este caso usamos 0.05 padj y 0.8 de logFC

##########################################################################################################################################

#Filtrado de datos
data_filtered <- FinalTable %>% filter(adj.P.Val < 0.05 & (logFC > 0.8 | logFC < -0.8 ))
data_filtered <- data_filtered %>% drop_na

#Cantidad de genes expresados diferencialmente
nrow(data_filtered)

#Exportamos la lista filtrada
write.table(data_filtered, file="./07_RESULTADOS/DEG_FILTER_GSE48149.txt", sep="\t", row.names=F, col.names=T, quote=F)


#####################################################################
#Heatmap
#####################################################################

#Calculo de valor Z
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
dataZ <- t(apply(data_filtered[,10:29], 1, cal_z_score))
data_filtered2 <- cbind(data_filtered[,1:9], dataZ[,12:20], dataZ[,1:11])
data_filtered3 <- data_filtered2 %>% drop_na

#Seleccionar color
col_fun <- colorRamp2(seq(min(data_filtered3[,10:29]), max(data_filtered3[,10:29]), length = 3), c("#0000ff", "white", "#fb0007"))
col_fun

rwb <- colorRampPalette(colors = c("blue", "white", "red"))(30)

#Legends
lgd <- Legend(col_fun = col_fun, title = "Row Z-Score")

#Heatmap sin genes
pdf("./06_GRAFICOS_DE_EXPRESION_DIFERENCIAL/Heatmap_GSE48149.pdf")
Heatmap(as.matrix(data_filtered3[,10:29]),
        name = "Z-Score", column_title = "GSE48149 | Differential gene expression heatmap",  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
        col = rwb,
        column_order = order(as.numeric(gsub("column", "", colnames(data_filtered3[,10:29])))),
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
