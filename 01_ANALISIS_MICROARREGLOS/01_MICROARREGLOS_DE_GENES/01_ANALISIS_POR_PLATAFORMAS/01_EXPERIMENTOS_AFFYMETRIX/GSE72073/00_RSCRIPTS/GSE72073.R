#####################################################################
#Titulo: Pipeline para Affymetrix Human Transcriptome Array 2.0
#Autor: Jose Antonio Ovando Ricardez
#Fecha: Abril/2021

#####################################################################
# Instalacion de paqueterias
#####################################################################
#install.packages("devtools")
#devtools::install_github("r-lib/remotes")
#install.packages("BiocManager")
#BiocManager::install("oligoClasses")
#BiocManager::install("oligo")
#BiocManager::install("limma")
#BiocManager::install("GEOquery")
#BiocManager::install("ComplexHeatmap")
#BiocManager::install("circlize")
#BiocManager::install("EnhancedVolcano")
#BiocManager::install("arrayQualityMetrics")
#BiocManager::install("ArrayExpress")
#BiocManager::install("pd.hg.u133.plus.2")
#BiocManager::install("hgu133plus2.db")
#Plotting and color options packages
#install.packages("RColorBrewer")
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("tidyr")

#####################################################################
# Llamado de librerias
#####################################################################
library(devtools)
library(remotes)
#General Bioconductor packages
library(oligoClasses)
library(oligo)
library(limma)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(EnhancedVolcano)
library(GEOquery)
library(arrayQualityMetrics)
library(ArrayExpress)
library(hta20transcriptcluster.db)
#Plotting and color options packages)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)


#####################################################################
# Descargamos los datos del GEO
# La muestra GSM595417.CEL fue excluida porque presentaba mala calidad
#####################################################################

untar("./01_DATOS/GSE72073_RAW.tar", exdir = "./01_DATOS/CEL")
list.files("./01_DATOS/CEL")

#####################################################################
# Abrimos la tabla con el PhenoData (pd)
#####################################################################
pd = read.table("./02_PHENODATA/Phenodata_IPF_VS_CONTROL.txt", header=TRUE, as.is=TRUE)
pd[,"Filename"] = paste(pd[,"Group"], pd[,"Replicate"], sep=".")
pd

#####################################################################
# Leyendo archivos CEL
# La muestra GSM595417.CEL fue excluida porque presentaba mala calidad
#####################################################################
celList <- list.celfiles("./01_DATOS/CEL", pattern = '*.CEL.gz', full.names=TRUE, listGzipped=TRUE)
raw_data <- read.celfiles(celList)
pData(raw_data) = pd
sampleNames(raw_data) <- pd$Filename

#####################################################################
#Analisis de calidad (Imagen de cada array)
#####################################################################

for (i in 1:7)
{
  name = paste("./04_ANALISIS_DE_CALIDAD/01_IMAGENES/IMAGEN",i,".jpg",sep="")
  jpeg(name)
  image(raw_data[,i],main=pd$Filename[i])
  dev.off()
}


#####################################################################
#Analisis de calidad con arrayQualityMetrics (raw_data)
#####################################################################

arrayQualityMetrics(expressionset = raw_data[, 1:7],
                    outdir = "./04_ANALISIS_DE_CALIDAD/02_REPORTE_DE_CALIDAD_RAW_DATA",
                    force = TRUE,
                    do.logtransform = TRUE)

#####################################################################
#Normalizacion
#####################################################################

norm_data <- oligo::rma(raw_data)
norm_matrix <- Biobase::exprs(norm_data)

#####################################################################
#Analisis de calidad con arrayQualityMetrics (norm_data)
#####################################################################

arrayQualityMetrics(expressionset = norm_data,
                    outdir = "./04_ANALISIS_DE_CALIDAD/02_REPORTE_DE_CALIDAD_NORM_DATA",
                    force = TRUE)

#####################################################################
#PCA
#####################################################################

PCA <- prcomp(t(norm_matrix), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Group = pData(norm_data)$Group,
                     Sample = pData(norm_data)$Sample,
                     Replicate = pData(norm_data)$Replicate 
)

pdf("./04_ANALISIS_DE_CALIDAD/PCA_GSE72073.pdf")
ggplot(dataGG, aes(PC1, PC2, label=row.names(dataGG))) +
  geom_point(aes( colour = Group))+
  geom_text(size=1)+
  ggtitle("GSE72073 | PCA plot of the log-transformed norm expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.2))+
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = c("#0000ff", "#fb0007")) + theme(plot.title = element_text(face = "bold"))
dev.off()

#####################################################################
#Anotacion
#####################################################################

#Descargamos la anotacion de GEO
#Esto varia entre plataforma y plataforma

Annot <- AnnotationDbi::select(hta20transcriptcluster.db,
                               keys = (featureNames(norm_data)),
                               columns = c("SYMBOL",  'ENSEMBL',"GENENAME"),
                               keytype = "PROBEID")

gene_matrix <- as.data.frame(norm_matrix)
gene_matrix <- tibble::rownames_to_column(gene_matrix, "row_names")

Annot_final <- merge(x= Annot,y= gene_matrix,by.x='PROBEID',by.y='row_names',all=F)
Annot_final <- Annot_final %>% drop_na


#####################################################################
# Expresion diferencial
#####################################################################

casos <- as.factor(pd$Group)
design = model.matrix(~ 0+casos)
colnames(design) = levels(casos)
design

#Hacemos la matriz de contrastes
cont.matrix<- makeContrasts(IPF-Control, levels=design)
cont.matrix

fit<-lmFit(Annot_final[5:11],design)
fit2= contrasts.fit(fit,cont.matrix)
fit3= eBayes(fit2)

#Guardar tabla
table_CD = topTable(fit3, coef=1, number=nrow(fit), sort.by= "none",adjust="fdr")

FinalTable = data.frame(Annot_final[,1:4], table_CD, Annot_final[,5:11])

write.table(FinalTable, file="./05_TABLA_DE_EXPRESION_DIFERENCIAL/DEG_GSE72073.txt", sep="\t",
            row.names=F, col.names=T, quote=F)

#####################################################################
#Volcano
#####################################################################

#Valores de cada color
keyvals <- ifelse(FinalTable$logFC <= -0.8 & FinalTable$P.Value <0.01,
                  '#0000ff',  ifelse(FinalTable$logFC >=0.8 & FinalTable$P.Value <0.01,
                                     '#fb0007', '#e3deeb'))
#Nombre de valores de cada color
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == '#fb0007'] <- 'Up'
names(keyvals)[keyvals == '#e3deeb'] <- 'NS'
names(keyvals)[keyvals == '#0000ff'] <- 'Down'

#Generamos el Volcanoplot
pdf("./06_GRAFICOS_DE_EXPRESION_DIFERENCIAL/Volcano_GSE72073.pdf")
EnhancedVolcano(FinalTable,lab = FinalTable$SYMBOL,
                x = 'logFC',
                y = 'P.Value',
                ylab = bquote(~P.Value),
                xlab = bquote(~logFC),
                ylim = c(0, 5),
                xlim = c(-4, 4),
                pCutoff = 0.01,
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
                legendLabels = c("NS", expression(LogFC), "P.Value", expression(P.Value ~ and
                                                                                  ~ logFC)),
                legendLabSize = 8,
                title = "     GSE72073 | Volcano plot of differential expression",
                titleLabSize = 13,
                subtitleLabSize = 9,
                captionLabSize = 9,
                border = 'full',
                borderWidth = 0.5,
                borderColour = 'black',
                caption = NULL,
                subtitle = "Comparation between IPF vs Control                           pCutoff = 0.01,     FCcutoff = 0.8",
                labSize = 1.0)
dev.off()


##########################################################################################################################################
#Filtrado de datos y resulados
##########################################################################################################################################

#Aqui estableceremos nuestros valores de corte, en este caso usamos 0.05 padj y 0.8 de logFC

##########################################################################################################################################

#Filtrado de datos
data_filtered <- FinalTable %>% filter(P.Value < 0.01 & (logFC > 0.8 | logFC < -0.8 ))
data_filtered <- data_filtered %>% drop_na

#Cantidad de genes expresados diferencialmente
nrow(data_filtered)

#Exportamos la lista filtrada
write.table(data_filtered, file="./07_RESULTADOS/DEG_FILTER_GSE72073.txt", sep="\t", row.names=F, col.names=T, quote=F)

#####################################################################
#Heatmap
#####################################################################

#Calculo de valor Z
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
dataZ <- t(apply(data_filtered[,11:17], 1, cal_z_score))
data_filtered2 <- cbind(data_filtered[,1:10], dataZ[,5:7], dataZ[,1:4])
data_filtered3 <- data_filtered2 %>% drop_na

#Seleccionar color
col_fun <- colorRamp2(seq(min(data_filtered3[,11:17]), max(data_filtered3[,11:17]), length = 3), c("#0000ff", "white", "#fb0007"))
col_fun

rwb <- colorRampPalette(colors = c("#0000ff", "white", "#fb0007"))(30)

#Legends
lgd <- Legend(col_fun = col_fun, title = "Row Z-Score")

#Heatmap sin genes
pdf("./06_GRAFICOS_DE_EXPRESION_DIFERENCIAL/Heatmap_GSE72073.pdf")
Heatmap(as.matrix(data_filtered3[,11:17]),
        name = "Z-Score", column_title = "GSE72073 | Differential gene expression heatmap",  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
        col = rwb,
        column_order = order(as.numeric(gsub("column", "", colnames(data_filtered3[,11:17])))),
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
