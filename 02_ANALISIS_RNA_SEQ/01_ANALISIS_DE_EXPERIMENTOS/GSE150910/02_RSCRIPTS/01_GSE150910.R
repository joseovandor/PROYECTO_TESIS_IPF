##########################################################################################################################################
#Titulo: Pipeline para expresion diferencial usando conteos de GEO
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
library(GEOquery)
library(ComplexHeatmap)
library(circlize)
library(EnhancedVolcano)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)

#Plotting y opciones de color
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)

##########################################################################################################################################
#01.- Cargar phenodata
##########################################################################################################################################

#Cargamos el archivo de Phenodata, el cual contiene los datos sobre las
#muestras, este archivo fue creado en excel, también se puede usar GEOQuery para
#este paso

Phenodata <- read.table("./03_PHENODATA/Phenodata.txt", stringsAsFactors = F, sep = "", header =  T)

##########################################################################################################################################
#02.- Cargar nuestra matriz de conteos
##########################################################################################################################################

#Creamos un objeto que va almacenar nuestra matriz de conteo descargada del GEO,
#aqui los datos venian ya en una sola tabla, en caso de venir las cuentas por separado
#esto puede variar. Es importante solo quedarnos con los grupos que necesitaremos
#En este caso en particular se quitaron todas las columnas que contenian al grupo
#"chp" porque no era parte del estudio.

##########################################################################################################################################

#Objeto que contiene nuestros conteos
raw_counts <- read.delim(("./01_COUNTS/GSE150910_counts.csv"), stringsAsFactors = F, sep = ",", row.names = 1)

#Filtrado con  el paquete dplyr de nuestros grupos que no son de nuestro estudio
raw_counts <- raw_counts %>%
  dplyr::select(-contains("chp"))

#Convertimos el objeto a una matriz
raw_counts <- as.matrix(raw_counts)


##########################################################################################################################################
#03.- Union de nuestro Phenodata con la matriz de conteos
##########################################################################################################################################

#Tenemos que generar un objeto que almacene el Phenodata y conteos, hay varias maneras de hacer
#esto, en este pipeline se muestra una manera sencilla de realizarlo

##########################################################################################################################################

#Indicamos que los nombres de filas del Phenodata van a ser los nombres de la columna SampleName
rownames(Phenodata) <- Phenodata$SampleName

#Aqui verificamos si se han aplicado los cambios y si los nombres de filas del Phenodata y
#el nombre de las columnas de la matriz de conteos sea el mismo

all(rownames(Phenodata) %in% colnames(raw_counts))

all(colnames(raw_counts) %in% rownames(Phenodata))

#Indicamos que los grupos que estan almacenados en la columna de Treatment son un factor
Phenodata$Group <- as.factor(Phenodata$Group)

#Generamos el objeto "dds", el cual contenda la informacion de los conteos, el phenodata y
#los grupos que deseamos evaluar para la posterior expresion diferencial

dds = DESeqDataSetFromMatrix(countData = raw_counts,
                             colData = Phenodata,
                             design = ~Group)

##########################################################################################################################################
#04.- PCA
##########################################################################################################################################

#Parea generar el PCA plot tenemos que normalizar los datos, para esto cuando son menos de 30
#muestras se recomienda transformar por rlog y cuando son mas de 30 es mejor usar vst

##########################################################################################################################################

#Transformamos los datos por vst porque son mas de 30 muestras y es lo recomendado
rld <- vst(dds, blind=FALSE)
print(rld)

#Generamos los datos para el PCA
pcaData <- plotPCA(rld, intgroup = "Group", returnData = TRUE)
print(pcaData)

#Convertimos las coordenadas a porcentajes
percentVar <- round(100 * attr(pcaData, "percentVar"))

#Generamos el PCA plot
pdf("./04_GRAFICOS/PCA_GSE150910.pdf")
ggplot(pcaData, aes(x = PC1, y = PC2, label=pcaData$name, color = Group)) +
  geom_point(size = 3) +
  geom_text(size=1)+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("GSE150910 | PCA plot of the log-transformed norm expression data") +
  theme(plot.title = element_text(hjust = 0.2))+
  scale_shape_manual(values = c(4,15)) +
  scale_color_manual(values = c("#0000ff", "#fb0007")) + theme(plot.title = element_text(face = "bold"))
dev.off()

##########################################################################################################################################
#05.- Expresion diferencial con DEseq2
##########################################################################################################################################

#Aqui utilizamos la funcion de DESeq2 llamada DESeq en la que tenemos que seleccionar nuestro
#objeto "dds" el cual almacena los conteos, el phenodata y los grupos. DESeq2 realiza varias
#tareas en una misma funcion

##########################################################################################################################################

#Llevamos acabo la expresion diferencial con la funcion DESeq
dds <- DESeq(dds)

#Especificamos el contraste entre los grupos que deseemos comparar
res <- results(dds, contrast=c("Group", "IPF", "Control"))

#Generamos un dataframe con los resultados
table_final <- data.frame(res)

#Agregamos las columnas de anotacion
columns(org.Hs.eg.db)

Annot <- AnnotationDbi::select(
  org.Hs.eg.db, keys=rownames(table_final),
  columns=c("ENSEMBL","SYMBOL","GENENAME","GENETYPE" ), keytype="SYMBOL")

#Convertimos nuestra matrix en dataframe y convertimos las rownames a una columna
table_final <- cbind(rownames(table_final),
                     data.frame(table_final, row.names=NULL))

#Creamos la tabla final con la anotacion y los valores
AnnotFinal <- merge(x=Annot,y=table_final,by.x='SYMBOL',by.y='rownames(table_final)')

#Eliminamos los NA
AnnotFinal <- AnnotFinal %>% drop_na


#Exportamos la tabla con los resultados
write.table(table_final, file = "./05_TABLAS/DEG_GSE150910.txt", quote = F, row.names = F, sep="\t", col.names = T)


##########################################################################################################################################
#06.- Filtrado de datos y resulados
##########################################################################################################################################

#Aqui estableceremos nuestros valores de corte, en este caso usamos 0.05 padj y 2 de logFC

##########################################################################################################################################

#Filtrado de datos
data_filtered <- AnnotFinal %>% filter(padj < 0.05 & (log2FoldChange > 1.5 | log2FoldChange < -1.5 ))

#Cantidad de genes expresados diferencialmente
nrow(data_filtered)

#Exportamos la lista filtrada
write.table(data_filtered, file="./06_RESULTADOS/DEG_FILTER_GSE150910.txt", sep="\t", row.names=T, col.names=T, quote=F)

##########################################################################################################################################
#07.- Volcano plot
##########################################################################################################################################

#Aqui haremos el grafico del volcan usando el paquete "Enhanced volcano"

##########################################################################################################################################

#Valores de cada color
keyvals <- ifelse(AnnotFinal$log2FoldChange <= -1.5 & AnnotFinal$padj <0.05,
                  '#0000ff',  ifelse(AnnotFinal$log2FoldChange >=1.5 & AnnotFinal$padj <0.05,
                                     '#fb0007', '#e3deeb'))
#Nombre de valores de cada color
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == '#fb0007'] <- 'Up'
names(keyvals)[keyvals == '#e3deeb'] <- 'NS'
names(keyvals)[keyvals == '#0000ff'] <- 'Down'

#Script del volcano
pdf("04_GRAFICOS/Volcano_GSE150910.pdf")
EnhancedVolcano(AnnotFinal,lab = AnnotFinal$SYMBOL,
                x = 'log2FoldChange',
                y = 'padj',
                ylab = bquote(~padj),
                xlab = bquote(~log2FoldChange),
                ylim = c(0, 75),
                xlim = c(-5, 5),
                pCutoff = 0.05,
                FCcutoff = 1.5,
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
                title = "     GSE150910 | Volcano plot of differential expression",
                titleLabSize = 13,
                subtitleLabSize = 9,
                captionLabSize = 9,
                border = 'full',
                borderWidth = 0.5,
                borderColour = 'black',
                caption = NULL,
                subtitle = "Comparation between IPF vs Control                          pCutoff = 0.05,     FCcutoff = 1.5",
                labSize = 1.0)
dev.off()


##########################################################################################################################################
#08.- Heatmap
##########################################################################################################################################

#Aqui haremos el grafico del volcan usando el paquete "Complexheatmap"

##########################################################################################################################################

#Normalizamos los datos
cnorm <- counts(dds, normalized = T)

#Aqui estamos ordenando nuestros datos para crear el heatmap
IPF <- as.data.frame(cnorm) %>% dplyr::select(-contains("control"))
Control <- as.data.frame(cnorm) %>% dplyr::select(-contains("ipf"))

#Unimos los dos objetos en uno mismo
cnorm_sort <- data.frame(Control,IPF)
colnames(cnorm_sort)

#Unimos los datos de conteos normalizados con la tabla de resultados de la expresion diferencial
table_final_filtered <- data.frame(table_final, cnorm_sort)

#Realizamos un filtrado en base a los puntos de corte establecidos
data_filtered <- table_final_filtered %>% filter(padj < 0.05 & (log2FoldChange > 1.5 | log2FoldChange < -1.5 ))

#Calculo de valor Z
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
dataZ <- t(apply(data_filtered[,8:213], 1, cal_z_score))
data_filtered2 <- cbind(data_filtered[,2:7], dataZ[,1:103], dataZ[,104:206])
data_filtered3 <- data_filtered2 %>% drop_na

#Seleccionar color
col_fun <- colorRamp2(seq(min(data_filtered3[,7:212]), max(data_filtered3[,7:212]), length = 3), c("#0000ff", "white", "#fb0007"))
col_fun

rwb <- colorRampPalette(colors = c("#0000ff", "white", "#fb0007"))(30)

#Legends
lgd <- Legend(col_fun = col_fun, title = "Row Z-Score")

#Heatmap sin genes
pdf("./04_GRAFICOS/Heatmap_GSE150910.pdf")
Heatmap(as.matrix(data_filtered3[,7:212]),
        name = "Z-Score", column_title = "GSE150910 | Differential gene expression heatmap",  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
        col = rwb,
        column_order = order(as.numeric(gsub("column", "", colnames(data_filtered3[,7:212])))),
        clustering_distance_rows = "euclidean",
        row_names_gp = gpar(fontsize = 0),
        column_names_gp = gpar(fontsize = 3, fontface = "bold"),
        clustering_method_rows = "ward.D2",
        width = unit(9, "cm"),
        height = unit(15, "cm"))
dev.off()
