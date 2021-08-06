library(stringr)
library(grid)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(tidyverse)
library(pheatmap)
library(corrplot)
library(readxl)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

file=as.matrix(read.table('fastani.txt',fill = T, row.names = 1, col.names = 1:50, skip = 1, na.strings = c("NA")))


file1<-file
one <-c(rep(NA,50))
file1<-cbind(file1,one)

mat<-as.matrix(file1)

mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
file1<-mat
diag(file1)=100


file2<-as.data.frame(file1)
colnames(file2) <-row.names(file2)
columns_names <- rownames(file2)
meta <- readxl::read_excel('Bvulgatus.xlsx',sheet=1,na='NA')



meta<-meta[order(meta$Sample_ID,decreasing = F),]
meta$old_name <- factor(columns_names)  
meta <- meta[order(meta$old_name), ]

heatmap_Ann = HeatmapAnnotation(
  Group=meta$Group, 
  col = list(Group = c("dog"="#ec7014","human"="#01665e","ref"="#80cdc1")))


mycol <- colorRamp2(c(95,99,100),colors = c("#3288bd","white","#d53e4f"))

p2<-Heatmap(as.matrix(file2),
            cluster_rows = T,   
            show_row_names = FALSE, 
            cluster_columns = T,
            show_column_names = FALSE,
            #column_order=anno$name,
            top_annotation = heatmap_Ann,
            col=mycol,
            #column_names_rot = 45,
            name ="ANI"
)
p2


