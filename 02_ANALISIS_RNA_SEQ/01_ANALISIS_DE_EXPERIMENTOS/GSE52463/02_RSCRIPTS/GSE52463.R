##########################################################################################################################################
#Titulo: Pipeline para expresion diferencial usando conteos de HTSEQ
#Autor: Jose Antonio Ovando Ricardez
#Fecha: Agosto/2021

##########################################################################################################################################
# Instalacion de paqueterias
##########################################################################################################################################

#Paquetes generales de Bioconductor
#BiocManager::install("ComplexHeatmap")
#BiocManager::install("circlize")
#BiocManager::install("EnhancedVolcano")
#BiocManager::install("DESeq2")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("AnnotationDbi")

#Plotting y opciones de coloress
#install.packages("RColorBrewer")
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("tidyr")

##########################################################################################################################################
# Llamado de librerias
##########################################################################################################################################

library(devtools)
library(remotes)

#Paquetes generales de bioconductor
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(EnhancedVolcano)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pheatmap)
library(biomaRt)

#Plotting y opciones de color
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)

##########################################################################################################################################

#Definir la localizacion de la carpeta que contiene los conteos
directory <- "./01_COUNTS"
print(directory)

#Visualizar una lista de los archivos de conteos
sampleFiles <- grep(".txt",list.files(directory),value=TRUE)
print(sampleFiles)

##########################################################################################################################################

#Importar la tabla del Phenodata, esta fue creada manualmente en excel
#Donde se incluye el nombre de las muestras y los grupos a los que corresponden

st <- read.table("03_PHENODATA/Phenodata.txt", header = T)
st

#Indicamos que los grupos que estan almacenados en la columna de Group son un factor
st$Group <- as.factor(st$Group)


#Creamos el objeto de sampleTable en donde vamos a seleccionar como valores
#de condicion usando los nombres de nuestros grupos de nuestro Phenodata.csv

sampleTable <- data.frame(sampleName = st$SampleName,
                          fileName = sampleFiles,
                          condition = st$Group)

print(sampleTable)

##########################################################################################################################################

#Creamos el objeto DEseq con el output de HTSeq

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ condition)

print(dds)


##########################################################################################################################################
#Boxplot raw counts
##########################################################################################################################################

pdf("./04_GRAFICOS/BP_RAW_GSE52463.pdf")
boxplot(assay(dds),
        main="Boxplot of raw counts",
        col = c("darkslategray3","darkslategray3","darkslategray3","darkslategray3"
                ,"darkslategray3","darkslategray3","darkslategray3","brown1"
                ,"brown1","brown1","brown1","brown1","brown1",
                "brown1","brown1"),
        las=2,
        boxwex = 0.6,
        staplewex = 0.6,
        outline=F)
dev.off()

##########################################################################################################################################
#PCA
##########################################################################################################################################

#Normalizamos por rlog para crear el PCA, usamos este metodo porque son pocas muestras
rld <- rlog(dds, blind=FALSE)

#Generamos los datos para el PCA
pcaData <- plotPCA(rld, intgroup = "condition", returnData = TRUE)

#Convertimos las coordenadas a porcentajes
percentVar <- round(100 * attr(pcaData, "percentVar"))

#Generamos el PCA plot

pdf("04_GRAFICOS/PCA_GSE52463.pdf")
ggplot(pcaData, aes(x = PC1, y = PC2, label=pcaData$name, color = condition)) +
  geom_point(size = 3) +
  geom_text(size=1)+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("GSE52463 | PCA plot of the log-transformed norm expression data") +
  theme(plot.title = element_text(hjust = 0.2))+
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = c("#0000ff", "#fb0007")) + theme(plot.title = element_text(face = "bold"))
dev.off()



##########################################################################################################################################
#Boxplot rlog convertion
##########################################################################################################################################

pdf("./04_GRAFICOS/BP_TRANSFORM_GSE52463.pdf")
boxplot(assay(rld),
        main="Boxplot of rlog convertion",
        xlab="", ylab=bquote(~rlog~expression),
        col = c("darkslategray3","darkslategray3","darkslategray3","darkslategray3"
                ,"darkslategray3","darkslategray3","darkslategray3","brown1"
                ,"brown1","brown1","brown1","brown1","brown1",
                "brown1","brown1"),
        las=2,
        boxwex = 0.6,
        staplewex = 0.6,
        outline=F)
dev.off()


##########################################################################################################################################
#Distancia entre muestras
##########################################################################################################################################

sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf("./04_GRAFICOS/DISTANCE_HEATMAP_GSE52463.pdf")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()


##########################################################################################################################################
#Expresiond diferencial
##########################################################################################################################################

#Llevamos acabo el test entre dos condiciones
dds <- DESeq(dds)

res <- results(dds, contrast=c("condition", "IPF", "Control"))


##########################################################################################################################################
#Anotacion
##########################################################################################################################################

#Creamos un dataframe con todos los datos finales

table_final <- data.frame(res)


#Eliminamos todos los valores NA de nuestros datos
table_final <- table_final %>% drop_na

#Agregamos las columnas de anotacion
columns(org.Hs.eg.db)
Annot <- AnnotationDbi::select(
  org.Hs.eg.db, keys=rownames(table_final),
  columns=c("ENSEMBL","SYMBOL","GENENAME","GENETYPE" ), keytype="ENSEMBL")

#Convertimos nuestra matrix en dataframe y convertimos las rownames a una columna
table_final <- cbind(rownames(table_final),
                     data.frame(table_final, row.names=NULL))

#Creamos la tabla final con la anotacion y los valores
table_final <- merge(x=Annot,y=table_final,by.x='ENSEMBL',by.y='rownames(table_final)')

#Eliminamos todos los valores NA de nuestros datos
table_final <- table_final %>% drop_na

#Extraemos las tablas de nuestros datos
write.table(table_final, file="05_TABLAS/DEG_GSE52463.txt", sep="\t", row.names=F, col.names=TRUE, quote=F)

##########################################################################################################################################
#Volcano
##########################################################################################################################################

#Valores de cada color
keyvals <- ifelse(table_final$log2FoldChange <= -1 & table_final$padj <0.05,
                  '#0000ff',  ifelse(table_final$log2FoldChange >=1 & table_final$padj <0.05,
                                     '#fb0007', '#e3deeb'))
#Nombre de valores de cada color
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == '#fb0007'] <- 'Up'
names(keyvals)[keyvals == '#e3deeb'] <- 'NS'
names(keyvals)[keyvals == '#0000ff'] <- 'Down'

#Script del volcano
pdf("04_GRAFICOS/Volcano_GSE52463.pdf")
EnhancedVolcano(table_final,lab = table_final$SYMBOL,
                x = 'log2FoldChange',
                y = 'padj',
                ylab = bquote(~padj),
                xlab = bquote(~log2FoldChange),
                ylim = c(0, 10),
                xlim = c(-5, 5),
                pCutoff = 0.05,
                FCcutoff = 1,
                colCustom = keyvals,
                cutoffLineType = 'solid',
                cutoffLineCol = 'coral4',
                hlineWidth = 0.2,
                cutoffLineWidth = 0.2,
                axisLabSize = 8,
                pointSize = 1.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabels = c("NS", expression(log2FoldChange), "adj.P.Val", expression(padj ~ and
                                                                                  ~ log2FoldChange)),
                legendLabSize = 8,
                title = "     GSE52463 | Volcano plot of differential expression",
                titleLabSize = 13,
                subtitleLabSize = 9,
                captionLabSize = 9,
                border = 'full',
                borderWidth = 0.5,
                borderColour = 'black',
                caption = NULL,
                subtitle = "Comparation between IPF vs Control                            pCutoff = 0.05,     FCcutoff = 1",
                labSize = 1.0)
dev.off()

##########################################################################################################################################
#Filtrado de datos y resulados
##########################################################################################################################################

#Aqui estableceremos nuestros valores de corte, en este caso usamos 0.05 padj y 0.8 de logFC

##########################################################################################################################################

#Filtrado de datos
data_filtered <- table_final %>% filter(padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1 ))

#Cantidad de genes expresados diferencialmente
nrow(data_filtered)

#Exportamos la lista filtrada
write.table(data_filtered, file="./06_RESULTADOS/DEG_FILTER_GSE52463.txt", sep="\t", row.names=F, col.names=T, quote=F)


##########################################################################################################################################
#Heatmap
##########################################################################################################################################

#Obtenemos los datos para hacer el heatmap incluyendo los recuentos normalizados

cnorm <- counts(dds, normalized = T)

table_final <-merge(table_final,cnorm,by='row.names',all=TRUE)
table_final <- table_final %>% drop_na

#Realizamos el filtrado

data_filtered <- table_final %>% filter(padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1 ))
nrow(data_filtered)

#Calculo de valor Z
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
dataZ <- t(apply(data_filtered[,11:25], 1, cal_z_score))
data_filtered2 <- cbind(data_filtered[,1:10], dataZ[,1:7], dataZ[,8:15])
data_filtered3 <- data_filtered2 %>% drop_na

#Seleccionar color
col_fun <- colorRamp2(seq(min(data_filtered3[,11:25]), max(data_filtered3[,11:25]), length = 3), c("#0000ff", "white", "#fb0007"))
col_fun

rwb <- colorRampPalette(colors = c("#0000ff", "white", "#fb0007"))(30)

#Legends
lgd <- Legend(col_fun = col_fun, title = "Row Z-Score")

#Heatmap sin genes
pdf("./04_GRAFICOS/Heatmap_GSE52463.pdf")
Heatmap(as.matrix(data_filtered3[,11:25]),
        name = "Z-Score", column_title = "GSE52463 | Differential gene expression heatmap",  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
        col = rwb,
        column_order = order(as.numeric(gsub("column", "", colnames(data_filtered3[,11:25])))),
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

