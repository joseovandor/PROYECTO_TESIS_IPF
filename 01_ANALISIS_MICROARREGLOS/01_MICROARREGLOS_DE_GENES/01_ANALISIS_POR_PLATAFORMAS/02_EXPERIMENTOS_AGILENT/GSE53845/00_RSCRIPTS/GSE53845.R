#####################################################################
#Titulo: Pipeline para Agilent-014850 Whole Human Genome 
#Whole Human Genome Microarray 4x44K G4112F
#Autor: Jose Antonio Ovando Ricardez
#Fecha: Julio/2021

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
                         source="agilent.median",
                         green.only=F)

#####################################################################
# Control de calidad 
#####################################################################
# Boxplot
#####################################################################

pdf("./04_ANALISIS_DE_CALIDAD/BP_RAWDATA_GSE53845.pdf")
boxplot(project$Rb,
        main="Boxplot of log2-intensitites for the Raw data",
        xlab="", ylab=bquote(~Log[2]~expression),
        names= targetinfo$Filename,
        col = c("brown1", "brown1","brown1","brown1","brown1","brown1","brown1","brown1","brown1","brown1","brown1","brown1","brown1","brown1","darkslategray3"
                ,"brown1","brown1","brown1","brown1","brown1","brown1","darkslategray3","brown1","darkslategray3","darkslategray3","brown1","brown1","brown1","brown1","brown1",
                "brown1","brown1","brown1","brown1","brown1","brown1","brown1","darkslategray3","brown1","darkslategray3","darkslategray3","brown1","darkslategray3","brown1"
                ,"brown1","brown1","brown1",
                "brown1"),
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

pdf("./04_ANALISIS_DE_CALIDAD/BP_BC_GSE53845.pdf")
boxplot(project.bgc$R,
        main="Boxplot of log2-intensitites with BC",
        xlab="", ylab=bquote(~Log[2]~expression),
        names= targetinfo$Filename,
        col = c("brown1", "brown1","brown1","brown1","brown1","brown1","brown1","brown1","brown1","brown1","brown1","brown1","brown1","brown1","darkslategray3"
                ,"brown1","brown1","brown1","brown1","brown1","brown1","darkslategray3","brown1","darkslategray3","darkslategray3","brown1","brown1","brown1","brown1","brown1",
                "brown1","brown1","brown1","brown1","brown1","brown1","brown1","darkslategray3","brown1","darkslategray3","darkslategray3","brown1","darkslategray3","brown1"
                ,"brown1","brown1","brown1",
                "brown1"),
        las=2,
        boxwex = 0.6,
        staplewex = 0.6,
        outline=F)
dev.off()

#####################################################################
# Normalizacion
#####################################################################

project.NormData <-normalizeBetweenArrays(project.bgc,method = 
                                            "quantile")

#####################################################################
# Boxplot (Norm Data)
#####################################################################

pdf("./04_ANALISIS_DE_CALIDAD/BP_NORM_DATA_GSE53845.pdf")
boxplot(project.NormData$A,
        main="Boxplot of log2-intensitites normalized",
        xlab="", ylab=bquote(~Log[2]~expression),
        names= targetinfo$Filename,
        col = c("brown1", "brown1","brown1","brown1","brown1","brown1","brown1","brown1","brown1","brown1","brown1","brown1","brown1","brown1","darkslategray3"
                ,"brown1","brown1","brown1","brown1","brown1","brown1","darkslategray3","brown1","darkslategray3","darkslategray3","brown1","brown1","brown1","brown1","brown1",
                "brown1","brown1","brown1","brown1","brown1","brown1","brown1","darkslategray3","brown1","darkslategray3","darkslategray3","brown1","darkslategray3","brown1"
                ,"brown1","brown1","brown1",
                "brown1"),
        las=2,
        boxwex = 0.6,
        staplewex = 0.6,
        outline=F)
dev.off()

#####################################################################
# PCA (NormData)
#####################################################################

PCA <- prcomp(t(project.NormData$A), scale = F)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Disease =
                       Biobase::pData(targetinfo)$Group)

pdf("./04_ANALISIS_DE_CALIDAD/PCA_GSE53845.pdf")
ggplot(dataGG, aes(PC1, PC2, label=row.names(dataGG))) +
  geom_text(size=0.3)+geom_point(aes(colour = targetinfo$Group)) +
  ggtitle("GSE53845 | PCA plot of the log-transformed norm expression data") +
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
Annot2 <- Annot1[1:45015, c(3, 4, 5)]
geneMatrix <- cbind(Annot2, project.NormData$A)

#Filtramos las sondas de control y conservamos solo las que
#tienen un valor de 0
Annot3 <- (geneMatrix[geneMatrix$'ControlType' == "0",])
Annot4 <- Annot3[-1]
Annot4 <- Annot4[-2]

#Descargamos la informacion extra usando GEOQuery
gpl <- getGEO('GPL6480')
Table <- Table(gpl)[1:45015, c(2, 7, 8)]

AnnotFinal <- merge(Table, Annot4, by.x = "SPOT_ID", by.y = "ProbeName")

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

fit<-lmFit(AnnotFinal[4:51],design)
fit2= contrasts.fit(fit,cont.matrix)
fit3= eBayes(fit2)
fit3
head(fit3)

#Guardar tabla
table_CD = topTable(fit3, coef=1, number=nrow(fit), sort.by= "none",adjust="fdr")

#Creamos la tabla final
FinalTable = data.frame(AnnotFinal[,1:3], table_CD, AnnotFinal[,4:51])

write.table(FinalTable, file="./05_TABLA_DE_EXPRESION_DIFERENCIAL/DEG_GSE53845.txt", sep="\t",
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
pdf("./06_GRAFICOS_DE_EXPRESION_DIFERENCIAL/Volcano_GSE53845.pdf")
EnhancedVolcano(FinalTable,lab = FinalTable$GENE_SYMBOL,
                x = 'logFC',
                y = 'adj.P.Val',
                ylab = bquote(~adj.P.Val),
                xlab = bquote(~logFC),
                ylim = c(0, 15),
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
                title = "     GSE53845 | Volcano plot of differential expression",
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
write.table(data_filtered, file="./07_RESULTADOS/DEG_FILTER_GSE53845.txt", sep="\t", row.names=F, col.names=T, quote=F)

#####################################################################
#Heatmap
#####################################################################

#Esta paso lo hice solo para este experimento para cambiar el nombre de las muestras
nrow(data_filtered)
colnames(data_filtered)[10:57] <- c("IPF1","IPF2","IPF3","IPF4","IPF5","IPF6","IPF7"
                                    ,"IPF8","IPF9","IPF10","IPF11","IPF12","IPF13","IPF14"
                                    ,"Control1","IPF15","IPF16","IPF17","IPF18","IPF19","IPF20"
                                    ,"Control2","IPF21","Control3","Control4","IPF22","IPF23"
                                    ,"IPF24","IPF25","IPF26","IPF27","IPF28","IPF29","IPF30","IPF31"
                                    ,"IPF32","IPF33","Control5","IPF34","Control6","Control7","IPF35"
                                    ,"Control8","IPF36","IPF37","IPF38","IPF39","IPF40")


#Aqui estamos ordenando nuestros datos para crear el heatmap
IPF <- data_filtered %>% dplyr::select(-contains("Control"))
Control <- data_filtered %>% dplyr::select(-contains("IPF"))

data_filtered <- data.frame(Control, IPF[10:49])


#Calculo de valor Z
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
dataZ <- t(apply(data_filtered[,10:57], 1, cal_z_score))
data_filtered2 <- cbind(data_filtered[,1:9], dataZ[,1:8], dataZ[,9:48])
data_filtered3 <- data_filtered2 %>% drop_na

#Seleccionar color
col_fun <- colorRamp2(seq(min(data_filtered3[,10:18]), max(data_filtered3[,10:57]), length = 3), c("#0000ff", "white", "#fb0007"))
col_fun

#Legends
lgd <- Legend(col_fun = col_fun, title = "Row Z-Score")

#Heatmap sin genes
pdf("./06_GRAFICOS_DE_EXPRESION_DIFERENCIAL/Heatmap_GSE53845.pdf")
Heatmap(as.matrix(data_filtered3[,10:57]),
        name = "Z-Score", column_title = "GSE53845 | Differential gene expression heatmap",  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
        col = col_fun,
        column_order = order(as.numeric(gsub("column", "", colnames(data_filtered3[,10:57])))),
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
