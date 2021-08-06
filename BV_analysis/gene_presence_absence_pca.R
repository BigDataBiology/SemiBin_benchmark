library(vegan)
library(ape)
library(ggplot2)
library(ggrepel)
library(readxl)
library(ggsci)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


gene_presence_absence_tab <- read.delim(file="gene_presence_absence.Rtab",header=T,row.names = 1,  check.names=FALSE, sep="\t")

meta <- read_excel('Bvulgatus.xlsx',sheet=1,na='NA')

gene_presence_absence_tab <- as.data.frame(t(gene_presence_absence_tab))

gene_pca <- prcomp(gene_presence_absence_tab)

pca_sum <- summary(gene_pca)
pc12  <-  as.data.frame(gene_pca$x[, 1:2]) 

pc <- summary(gene_pca)$importance[2,]*100 
palette2 <- pal_npg("nrc")(10)


pc12$group <- meta$Group[order(meta$Sample_ID)]
pc12$id <- row.names(pc12)

p7  <- ggplot(pc12)+
  geom_point(aes(x = PC1, y = PC2, color=group,shape=group, alpha =0.5),size=3.5)+ 
  scale_color_manual(values = palette2) +
  stat_ellipse(data = pc12,aes(x = PC1,y = PC2, fill = group,color= group ), geom = 'polygon',
               level = 0.95, alpha = 0.3, show.legend = T)+
  scale_alpha(guide = 'none') +
  labs(x=paste("PC1 (", format(pc[1], digits=4), "%)", sep=""),
       y=paste("PC2 (", format(pc[2], digits=4), "%)", sep=""))+
  theme(axis.title.x = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5))+ 
  theme(axis.title.y = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5))+ 
  theme(axis.text.x = element_text(size = 10, color = "black", vjust = 0.5, hjust = 0.5))+ 
  theme(axis.text.y = element_text(size = 10, color = "black", vjust = 0.5, hjust = 0.5))+
  theme(legend.title = element_text(color="black", size = 10))+
  theme(legend.text = element_text(color="black", size = 10))+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), plot.margin = unit(c(0.5, 0.25, 0.5, 0.2), "inches"),
        axis.line = element_line(colour = "black"))

p7
