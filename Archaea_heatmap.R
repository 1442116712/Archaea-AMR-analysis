# install ####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
# install.packages("vegan")
library(ComplexHeatmap)
library(circlize)

col_fun = circlize::colorRamp2(c(0,1.5,3), c("white", "#B2182B","#67000D"))
# MAGs Heatmap ####
# input data
data1 <- read.table('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/MAGs.txt',
                    header = TRUE,sep = "\t", row.names = 1)
df1 <- as.matrix(data1) 

## Side annotation
rowsum1 <- read.table('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/MAGs_sum.txt',
                      header = TRUE,sep = "\t", row.names = 1)

bar1 = rowAnnotation(
  ARG_occurrence = anno_barplot(rowsum1,
                                  baseline = 0, bar_width = 0.9, width = unit(1.3, "cm"), gp = gpar(col = "white", fill = "#2171B5"), 
                                  border = F, border_gp = gpar(lwd = 2), 
                                  axis_param = list(side = "bottom", at = c(0,10,20), labels = c("0","10","20"))), # direction = "reverse" - can reverse the direction 
  show_annotation_name = T, annotation_name_gp = gpar(fontsize = 8), annotation_name_side = "bottom", annotation_name_rot = 270)


## Bottom annotation
ARG_type1 <- read.table('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/MAGs_type.txt',header = TRUE,sep = "\t", row.names = 1)

bar2 <- HeatmapAnnotation(
  show_annotation_name = T,annotation_name_gp = gpar(fontsize = 8), annotation_name_side = "left",annotation_name_rot = 0,
  ARG_diversity = anno_barplot(    
    ARG_type1, 
    baseline = 0, bar_width = 0.9, height = unit(1.5, "cm"), gp = gpar(col = "white", fill = "#4292C6"), # bar setting
    border = F, border_gp = gpar(lwd = 2), 
    axis_param = list(direction = "reverse", side = "left", at = c(0,3,6), # axis parameter direction = "reverse", 
                      labels = c("","3","6")),ylim = c(0, 6)), gap = unit(1.5, "mm"))


## Top annotation
prevalence1 <- read.table('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/MAGs_bar.txt',header = TRUE,sep = "\t", row.names = 1)
prevalence1$Percentage <- (prevalence1$ARG / prevalence1$Meta) * 100
prevalence_matrix1 <- as.matrix(prevalence1) 
Metadata1 <- prevalence_matrix1[ ,1]
ARG1 <- prevalence_matrix1[ ,2]

bar3 <- HeatmapAnnotation(
  ARG_prevalence = anno_barplot(
    cbind(prevalence1$Meta, prevalence1$ARG),
    gp = gpar(fill = c("#08519C", "#FF69B4")), 
    beside = TRUE,                        
    bar_width = 0.85, height = unit(1.8, "cm"),
    border = F, 
    axis_param = list(side = "right", at = c(0,1000,2000), labels = c("0","1000","2000"))
  ),
  show_annotation_name = T, 
  annotation_name_gp = gpar(fontsize = 8), 
  annotation_name_side = "right", 
  annotation_name_rot = 0,
  
  Percentage = anno_text(
    paste0(round(prevalence1$Percentage, 0), "%"), 
    location = 0.5,  
    just = "center",            
    rot = 0,                   
    gp = gpar(fontsize = 6, fontface = "bold")     
  )
)

# heatmap setting
Heatmap1 = Heatmap(df1,col = col_fun,
                   name = "log2(MSS-normalized counts + 1)",
                   width = unit(6.8, "cm"), 
                   height = unit(8, "cm"), 
                   cluster_columns = F, 
                   cluster_rows = F, 
                   show_row_dend = F, 
                   show_column_dend = F,
                   row_dend_side = "left",
                   column_dend_height = unit(0.5, "cm"), 
                   row_dend_width = unit(0.5, "cm"), 
                   column_title = "MAGs", 
                   row_title = "",
                   border_gp = gpar(col = "black", lwd = 1),
                   row_title_side = "right", 
                   column_title_side = "bottom", 
                   column_title_rot = F,
                   column_title_gp = gpar(fontsize = 10, fontface = "bold"), 
                   row_title_gp = gpar(fontsize = 10, fontface = "bold"), 
                   right_annotation =bar1, 
                   bottom_annotation = bar2,
                   top_annotation = bar3, 
                   row_names_side = "right", 
                   column_names_side = "bottom", 
                   row_names_gp = gpar(fontsize = 8, fontface = "bold"), 
                   column_names_gp = gpar(fontsize = 8, fontface = "bold"), 
                   show_heatmap_legend = T, 
                   show_column_names = T, 
                   show_row_names = F, 
                   column_gap = unit(0.7, 'mm')) 
Heatmap1



# PCGs Heatmap ####
# input data
data2 <- read.table('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/PCGs.txt',
                    header = TRUE,sep = "\t", row.names = 1)
df2 <- as.matrix(data2) 

## Side annotation
rowsum2 <- read.table('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/PCGs_sum.txt',
                      header = TRUE,sep = "\t", row.names = 1)

bar4 = rowAnnotation(
  ARG_occurrence = anno_barplot(rowsum2,
                                  baseline = 0, bar_width = 0.9, width = unit(1.5, "cm"), gp = gpar(col = "white", fill = "#2171B5"), 
                                  border = F, border_gp = gpar(lwd = 2), 
                                  axis_param = list(at = c(0,50,100), direction = "reverse", side = "bottom",labels = c("0","50","100"))), # direction = "reverse" - can reverse the direction 
  show_annotation_name = T, annotation_name_gp = gpar(fontsize = 8), annotation_name_side = "bottom", annotation_name_rot = 90)

## Bottom annotation
ARG_type2 <- read.table('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/PCGs_type.txt',header = TRUE,sep = "\t", row.names = 1)

bar5 <- HeatmapAnnotation(
  show_annotation_name = F,annotation_name_gp = gpar(fontsize = 8), annotation_name_side = "right",annotation_name_rot = 0,
  ARG_diversity = anno_barplot(    
    ARG_type2, 
    baseline = 0, bar_width = 0.9, height = unit(1.5, "cm"), gp = gpar(col = "white", fill = "#4292C6"), # bar setting
    border = F, border_gp = gpar(lwd = 2), 
    axis_param = list(direction = "reverse", side = "right", at = c(0,3,6), # axis parameter direction = "reverse", 
                      labels = c("0","3","6")),ylim = c(0, 6)), gap = unit(1.5, "mm"))

## Top annotation
prevalence2 <- read.table('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/PCGs_bar.txt',header = TRUE,sep = "\t", row.names = 1)
prevalence2$Percentage <- (prevalence2$ARG / prevalence2$Meta) * 100
prevalence_matrix2 <- as.matrix(prevalence2) 
Metadata2 <- prevalence_matrix2[ ,1]
ARG_detect2 <- prevalence_matrix2[ ,2]

bar6 <- HeatmapAnnotation(
  ARG_prevalence = anno_barplot(
    cbind(prevalence2$Meta, prevalence2$ARG),
    gp = gpar(fill = c("#08519C", "#FF69B4")), 
    beside = TRUE,                        
    bar_width = 0.85, height = unit(1.8, "cm"),
    border = F,
    axis_param = list(side = "left", at = c(0,400,800), labels = c("0","400","800"))
  ),
  show_annotation_name = T, 
  annotation_name_gp = gpar(fontsize = 8), 
  annotation_name_side = "left", 
  annotation_name_rot = 0,
  
  Percentage = anno_text(
    paste0(round(prevalence2$Percentage, 0), "%"), 
    location = 0.5,  
    just = "center",            
    rot = 0,                   
    gp = gpar(fontsize = 6, fontface = "bold")     
  )
)

Heatmap2 = Heatmap(df2,col = col_fun,
                   name = "log2(MSS-normalized counts + 1)",
                   width = unit(5, "cm"), 
                   height = unit(8, "cm"), 
                   cluster_columns = F, 
                   cluster_rows = F, 
                   show_row_dend = F, 
                   show_column_dend = F, 
                   row_dend_side = "left", 
                   column_dend_height = unit(0.5, "cm"), 
                   row_dend_width = unit(0.5, "cm"), 
                   column_title = "PCGs", 
                   border_gp = gpar(col = "black", lwd = 1),
                   row_title = "",
                   row_title_side = "left", 
                   column_title_side = "bottom", 
                   column_title_rot = F, 
                   column_title_gp = gpar(fontsize = 10, fontface = "bold"), 
                   row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                   left_annotation =bar4, 
                   bottom_annotation = bar5, 
                   top_annotation = bar6, 
                   row_names_side = "right", 
                   column_names_side = "bottom", 
                   row_names_gp = gpar(fontsize = 8, fontface = "bold"), 
                   column_names_gp = gpar(fontsize = 8, fontface = "bold"),
                   show_heatmap_legend = F, 
                   show_column_names = T, 
                   show_row_names = F) 

Heatmap2


# Final Heatmap ####
Heatmap0 = Heatmap2 + Heatmap1

Heatmap0


draw(Heatmap0, ht_gap = unit(30, "mm"))

png(
  filename = "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Heatmap_white.png",
  width = 12,       
  height = 6,       
  units = "in",     
  res = 600,        
  bg = "white"      
)

draw(Heatmap0, ht_gap = unit(30, "mm"))

dev.off()


# legend setting  ####
labels <- c("Total genomes in each class","ARG-carrying genomes in each class")
label_colors <- c("#08519C", "#FF69B4")
lgd <- Legend(labels = labels[1:2], labels_gp = gpar(fontsize = 8, fontface = "bold", col = "black"),
              legend_gp = gpar(fill = label_colors[1:2],cex = 0.5, frot = 1.2), 
              title_gp = gpar(fontsize = 10, fontface = "bold", col = "black"),
              title = "Legend",  title_position = "topcenter", border = T,
              gap = unit(1, "cm"), ncol = 1)

draw(lgd)

png(
  filename = "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Legend_HD.png",
  width = 4,     
  height = 2,
  units = "in",
  res = 600,     
  bg = "white",
  type = "cairo"
)

draw(lgd)
dev.off()


lgd = Legend(
  col_fun = col_fun, 
  title = "Log2 (MSS-normalised counts + 1)", 
  at = c(0, 1, 2, 3), 
  labels = c("0", "1", "2", "3"),
  title_position = "topcenter", 
  labels_gp = gpar(fontsize = 10, fontface = "bold"),
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  grid_width = unit(0.6, "cm"),
  legend_height = unit(4, "cm")
)
draw(lgd)
png(
  filename = "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Legend_HD1.png",
  width = 4,     
  height = 2,
  units = "in",
  res = 600,     
  bg = "white",
  type = "cairo"
)

draw(lgd)
dev.off()
