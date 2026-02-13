# install ####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
# install.packages("vegan")
library(ComplexHeatmap)
library(circlize)
library(vegan)

# input data
data1 <- read.table('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/MAGs.txt',
                    header = TRUE,sep = "\t", row.names = 1)
# change the data frame into array matrix
data_matrix1 <- as.matrix(data1) 

# Normalisation
dt1 <- decostand(t(data_matrix1), method = "total")
df1 <- t(dt1)

# Side annotation
rowsum1 <- read.table('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/MAGs_sum.txt',
                      header = TRUE,sep = "\t", row.names = 1)

bar1 = rowAnnotation(
  ARG_distribution = anno_barplot(rowsum1,
                                  baseline = 0, bar_width = 0.9, width = unit(1.3, "cm"), gp = gpar(col = "white", fill = "#20B2AA"), 
                                  border = F, border_gp = gpar(lwd = 2), 
                                  axis_param = list(side = "bottom", at = c(0,10,20), labels = c("0","10","20"))), # direction = "reverse" - can reverse the direction 
  show_annotation_name = T, annotation_name_gp = gpar(fontsize = 8), annotation_name_side = "bottom", annotation_name_rot = 90)


# Bottom annotation
ARG_type1 <- read.table('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/MAGs_type.txt',header = TRUE,sep = "\t", row.names = 1)

bar2 <- HeatmapAnnotation(
  show_annotation_name = F,annotation_name_gp = gpar(fontsize = 8), annotation_name_side = "right",annotation_name_rot = 0,
  ARG_diversity = anno_barplot(    
    ARG_type1, 
    baseline = 0, bar_width = 0.9, height = unit(1.5, "cm"), gp = gpar(col = "white", fill = "#6495ED"), # bar setting
    border = F, border_gp = gpar(lwd = 2), 
    axis_param = list(direction = "reverse", side = "left", at = c(0,3,6), # axis parameter direction = "reverse", 
                      labels = c("")),ylim = c(0, 6)), gap = unit(1.5, "mm"))


# Top annotation
prevalence1 <- read.table('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/MAGs_bar.txt',header = TRUE,sep = "\t", row.names = 1)
prevalence1$Percentage <- (prevalence1$ARG / prevalence1$Meta) * 100
prevalence_matrix1 <- as.matrix(prevalence1) 
Metadata1 <- prevalence_matrix1[ ,1]
ARG1 <- prevalence_matrix1[ ,2]

# 创建行注释
bar3 <- HeatmapAnnotation(
  # 叠加条形图
  ARG_prevalence = anno_barplot(
    cbind(prevalence1$Meta, prevalence1$ARG),
    gp = gpar(fill = c("#3CB371", "#FF69B4")), 
    beside = TRUE,                        # 设置为 FALSE 以实现叠加
    bar_width = 0.85, height = unit(1.8, "cm"),
    border = F,
  ),
  show_annotation_name = F, 
  annotation_name_gp = gpar(fontsize = 8), 
  annotation_name_side = "left", 
  annotation_name_rot = 0,
  
  # 百分比文本
  Percentage = anno_text(
    paste0(round(prevalence1$Percentage, 0), "%"), 
    location = 0.5,  # 将文本放在条形图的上方
    just = "center",            
    rot = 0,                   
    gp = gpar(fontsize = 6, fontface = "bold")     
  )
)


col_fun = circlize::colorRamp2(c(0, 0.5, 1), c("#000080", "white", "red")) # the range should match the number of the color

# heatmap setting
Heatmap1 = Heatmap(df1,col = col_fun,
                   name = "Class in MAGs",
                   width = unit(6.8, "cm"), 
                   height = unit(8, "cm"), # size of the well # heatmap_width (including the text and title)
                   cluster_columns = F, 
                   cluster_rows = F, # data cluster
                   show_row_dend = F, 
                   show_column_dend = F,
                   row_dend_side = "left", # need to cluster first
                   column_dend_height = unit(0.5, "cm"), 
                   row_dend_width = unit(0.5, "cm"), 
                   column_title = "MAGs", 
                   row_title = "AMR Categories",
                   row_title_side = "right", 
                   column_title_side = "bottom", 
                   column_title_rot = F, # title location and rotation
                   column_title_gp = gpar(fontsize = 10, fontface = "bold"), 
                   row_title_gp = gpar(fontsize = 10, fontface = "bold"), # fontface = "bold", fill = "green", col = "black", border = "black"), # title setting
                   right_annotation =bar1, 
                   bottom_annotation = bar2,
                   top_annotation = bar3,  # annotation chart
                   row_names_side = "right", 
                   column_names_side = "bottom", 
                   row_names_gp = gpar(fontsize = 8, fontface = "bold"), 
                   column_names_gp = gpar(fontsize = 8, fontface = "bold"), # location and size of the label name
                   show_heatmap_legend = F, 
                   show_column_names = T, 
                   show_row_names = F, 
                   column_gap = unit(0.7, 'mm')) 
Heatmap1



# Heatmap3
# input data
data2 <- read.table('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/PCGs.txt',
                    header = TRUE,sep = "\t", row.names = 1)
# change the data frame into array matrix
data_matrix2 <- as.matrix(data2) 

# Normalisation
dt2 <- decostand(t(data_matrix2), method = "total")
df2 <- t(dt2)

# Side annotation
rowsum2 <- read.table('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/PCGs_sum.txt',
                      header = TRUE,sep = "\t", row.names = 1)

# bar chart annotation sum by column and row
bar4 = rowAnnotation(
  ARG_distribution = anno_barplot(rowsum2,
                                  baseline = 0, bar_width = 0.9, width = unit(1.5, "cm"), gp = gpar(col = "white", fill = "#20B2AA"), 
                                  border = F, border_gp = gpar(lwd = 2), 
                                  axis_param = list(at = c(0,50,100), side = "bottom",labels = c("0","50","100"))), # direction = "reverse" - can reverse the direction 
  show_annotation_name = F, annotation_name_gp = gpar(fontsize = 8), annotation_name_side = "bottom", annotation_name_rot = 90)


# Bottom annotation
ARG_type2 <- read.table('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/PCGs_type.txt',header = TRUE,sep = "\t", row.names = 1)

bar5 <- HeatmapAnnotation(
  show_annotation_name = T,annotation_name_gp = gpar(fontsize = 8), annotation_name_side = "left",annotation_name_rot = 0,
  ARG_diversity = anno_barplot(    
    ARG_type2, 
    baseline = 0, bar_width = 0.9, height = unit(1.5, "cm"), gp = gpar(col = "white", fill = "#6495ED"), # bar setting
    border = F, border_gp = gpar(lwd = 2), 
    axis_param = list(direction = "reverse", side = "left", at = c(0,3,6), # axis parameter direction = "reverse", 
                      labels = c("0","3","6")),ylim = c(0, 6)), gap = unit(1.5, "mm"))


# Top annotation
prevalence2 <- read.table('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/PCGs_bar.txt',header = TRUE,sep = "\t", row.names = 1)
prevalence2$Percentage <- (prevalence2$ARG / prevalence2$Meta) * 100
prevalence_matrix2 <- as.matrix(prevalence2) 
Metadata2 <- prevalence_matrix2[ ,1]
ARG_detect2 <- prevalence_matrix2[ ,2]

# 创建行注释
bar6 <- HeatmapAnnotation(
  # 叠加条形图
  ARG_prevalence = anno_barplot(
    cbind(prevalence2$Meta, prevalence2$ARG),
    gp = gpar(fill = c("#3CB371", "#FF69B4")), 
    beside = TRUE,                        # 设置为 FALSE 以实现叠加
    bar_width = 0.85, height = unit(1.8, "cm"),
    border = F,
  ),
  show_annotation_name = F, 
  annotation_name_gp = gpar(fontsize = 8), 
  annotation_name_side = "left", 
  annotation_name_rot = 0,
  
  # 百分比文本
  Percentage = anno_text(
    paste0(round(prevalence2$Percentage, 0), "%"), 
    location = 0.5,  # 将文本放在条形图的上方
    just = "center",            
    rot = 0,                   
    gp = gpar(fontsize = 6, fontface = "bold")     
  )
)
# heatmap setting
Heatmap2 = Heatmap(df2,col = col_fun,
                   name = "Pure_class",
                   width = unit(5, "cm"), 
                   height = unit(8, "cm"), # size of the well # heatmap_width (including the text and title)
                   cluster_columns = F, 
                   cluster_rows = F, # data cluster
                   show_row_dend = F, 
                   show_column_dend = F, 
                   row_dend_side = "left", # need to cluster first
                   column_dend_height = unit(0.5, "cm"), 
                   row_dend_width = unit(0.5, "cm"), # dendrogram height and location
                   # rect_gp = gpar(col= "lightgrey",lwd = 0.5), # border color and width of each well, delete if do not want the border
                   #border = T, border_gp = gpar(col= "darkgrey",lwd = 1.2), # color and width of border of the whole map
                   column_title = "PCGs", 
                   row_title = "AMR Categories",
                   row_title_side = "left", 
                   column_title_side = "bottom", 
                   column_title_rot = F, # title location and rotation
                   column_title_gp = gpar(fontsize = 10, fontface = "bold"), 
                   row_title_gp = gpar(fontsize = 10, fontface = "bold"), # fontface = "bold", fill = "green", col = "black", border = "black"), # title setting
                   right_annotation =bar4, 
                   bottom_annotation = bar5, 
                   top_annotation = bar6, # annotation chart
                   row_names_side = "left", 
                   column_names_side = "bottom", 
                   row_names_gp = gpar(fontsize = 8, fontface = "bold"), 
                   column_names_gp = gpar(fontsize = 8, fontface = "bold"), # location and size of the label name
                   show_heatmap_legend = F, 
                   show_column_names = T, 
                   show_row_names = T) 

Heatmap2


Heatmap0 = Heatmap2 + Heatmap1

Heatmap0


draw(Heatmap0, ht_gap = unit(7, "mm"))
draw(lgd, x = unit(20, "cm"), y = unit(26.2, "cm")) # 12 inched x 15 inches
draw(lgd, x = unit(12, "cm"), y = unit(22, "cm")) # 800 Width x 900 height



# legend setting  ####
phylum_labels1 <- c("Total genomes","Resistance genomes")
label_colors1 <- c("#3CB371", "#FF69B4")
lgd3 = Legend(labels = phylum_labels1[1:2], labels_gp = gpar(fontsize = 8, fontface = "bold", col = "black"),
              legend_gp = gpar(fill = label_colors1[1:2],cex = 0.5, frot = 1.2), 
              title_gp = gpar(fontsize = 10, fontface = "bold", col = "black"),
              title = "Legend",  title_position = "topcenter", border = T,
              gap = unit(1, "cm"), ncol = 1)

draw(lgd3)
