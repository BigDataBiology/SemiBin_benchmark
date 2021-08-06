library(ggtreeExtra)
library(ggplot2)
library(ggtree)
library(treeio)
library(ggstar)
library(ggnewscale)

library(ComplexHeatmap)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

tree <- read.newick('iqtree_crc.contree',node.label='support')

library(readxl)
info <- read_excel('Bvulgatus.xlsx',sheet=1,na='NA')

info <-info[,-1]

cols<-c("dog"="#ec7014","human"="#01665e","ref"="#80cdc1")

node_data<-tree@data$node[tree@data$support>70]
node_data<-node_data[-1]

p <- ggtree(tree) %<+% info +
  geom_nodepoint(aes(subset=node %in%node_data,x=x),color="pink", 
                 alpha=1/2, size=2,na.rm=TRUE,show.legend=TRUE)+
  geom_tippoint(aes(color=Group),size=3) + 
  geom_treescale()+
  scale_color_manual(values=cols)+
  theme(legend.position="none")
   
p



# Heatmap
d=fortify(tree)
dd = subset(d, isTip)
tree_order<-dd$label[order(dd$y, decreasing=FALSE)]
file <- as.data.frame(read_excel('Tonb.xlsx',sheet=1,na='NA',col_names = TRUE))
file<-file[,-1]
rownames(file)<-file[,1] 
file<-file[,-1] 

library(tidyverse)
library(pheatmap)
library(corrplot)
library(readxl)
meta <- read_excel('Bvulgatus.xlsx',sheet=1,na='NA')

library('ComplexHeatmap')
library('dplyr')
library(circlize)

rownames(meta)<-meta$Sample_ID

meta<-meta[rev(tree_order),]
file<-file[,rev(tree_order)]

heatmap_Ann = rowAnnotation(
  Group=meta$Group, 
  col = list(Group = cols))

p2<-Heatmap(as.matrix(t(file)),
            cluster_rows = FALSE,   
            show_row_names = FALSE, 
            cluster_columns = T,
            show_column_names = TRUE,
            #column_order=anno$name,
            left_annotation = heatmap_Ann,
            col=c("1"="grey50","0"="grey90"),
            column_names_rot = 75,
            name ="ANI"
)

##################################################################################3
# split
btuB<-file[c(11,12,13),]
Tonb<-file[c(-11,-12,-13),]

p3<-Heatmap(as.matrix(t(Tonb)),
            cluster_rows = FALSE,    
            show_row_names = FALSE, 
            cluster_columns = T,
            show_column_names = TRUE,
            #column_order=anno$name,
            left_annotation = heatmap_Ann,
            col=c("1"="grey50","0"="grey90"),
            column_names_rot = 75,
            name ="Genes"
)

p4<-Heatmap(as.matrix(t(btuB)),
            cluster_rows = FALSE,  
            show_row_names = FALSE, 
            cluster_columns = T,
            show_column_names = TRUE,
            #column_order=anno$name,
            #left_annotation = heatmap_Ann,
            col=c("1"="grey50","0"="grey90"),
            column_names_rot = 75,
            name ="Genes"
)

p3+p4

