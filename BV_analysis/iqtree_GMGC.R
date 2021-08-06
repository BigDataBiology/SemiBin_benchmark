library(ggtreeExtra)
library(ggplot2)
library(ggtree)
library(treeio)
library(ggstar)
library(ggnewscale)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

tree <- read.newick('iqtree_gmgc.contree',node.label='support')

library(readxl)
info <- read_excel('gmgc.xlsx',sheet=1,na='NA')


cols<- c("dog_gut"="#ec7014","human_gut"="#01665e","reference"="#80cdc1")

p <- ggtree(tree,branch.length='none') %<+% info +
  geom_tippoint(aes(color=Group),size=1) + 
  scale_color_manual(values=cols)+
  theme(legend.position="none")


p




# Heatmap
d=fortify(tree)
dd = subset(d, isTip)
tree_order<-dd$label[order(dd$y, decreasing=FALSE)]

file <- as.data.frame(read_excel('GMGC_Tonb.xlsx',sheet=1,na='NA',col_names = TRUE))
rownames(file)<-file[,1] 
file<-file[,-1] 

library(tidyverse)
library(pheatmap)
library(corrplot)
library(readxl)
meta <- read_excel('gmgc.xlsx',sheet=1,na='NA')

library('ComplexHeatmap')
library('dplyr')
library(circlize)

rownames(meta)<-meta$Sample_ID

meta<-meta[rev(tree_order),]
file<-file[,rev(tree_order)]

heatmap_Ann = rowAnnotation(
  Group=meta$Group, 
  col = list(Group = c("dog_gut"="#ec7014","human_gut"="#01665e","reference"="#80cdc1")
  ))

mycol <- colorRamp2(c(0,1),colors = c("#3288bd","#d53e4f"))


p2<-Heatmap(as.matrix(t(file)),
            cluster_rows = FALSE,   
            show_row_names = FALSE, 
            cluster_columns = T,
            show_column_names = TRUE,
            #column_order=anno$name,
            left_annotation = heatmap_Ann,
            col=c("0"="grey90","1"="grey50"),
            #column_names_rot = 45,
            name ="genes"
)
p2

