#' R script for producing outputs in from **CD4+ cells** from Marcela's project
#' Author: Line Wulff 
#' Date (created): 23-08-18


#' #### ---- Initiate libraries ---- ####
library(flowCore)
library(FlowSOM)
library(Rtsne)
library(ggplot2)
library(stringr)
library(uwot)
library(ggrastr)
library(groupdata2)
library(viridis)

#### ---- variables used throughout script ---- ####
projdir <- getwd()
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)

#' #### ---- Read in data, inspect and transform labels if necessary ---- ####
data_FC <- flowCore::exprs(flowCore::read.FCS(
  filename = paste(projdir,"/inputdata/concat_CD4cells.fcs", sep=""), 
  transformation = FALSE, truncate_max_range = FALSE))

head(data_FC)
dim(data_FC)

# Add sample IDs and groups based on SampleID timestamps from Marcela
# Lowest CI7-1 1-5, middle CI8 1-6, top GF 1-4
data_FC_df <- as.data.frame(data_FC)
data_FC_df$group <- "CI8"
data_FC_df[data_FC_df$SampleID>190000,]$group <- "GF"
data_FC_df[data_FC_df$SampleID<90000,]$group <- "CI7-1"

data_FC_df$sample <- NA
data_FC_df[data_FC_df$SampleID<25000,]$sample <- "1"
data_FC_df[data_FC_df$SampleID<40000 & data_FC_df$SampleID>25000,]$sample <- "2"
data_FC_df[data_FC_df$SampleID<58000 & data_FC_df$SampleID>40000,]$sample <- "3"
data_FC_df[data_FC_df$SampleID<74000 & data_FC_df$SampleID>58000,]$sample <- "4"
data_FC_df[data_FC_df$SampleID<90000 & data_FC_df$SampleID>74000,]$sample <- "5"

data_FC_df[data_FC_df$SampleID<106000 & data_FC_df$SampleID>90000,]$sample <- "1"
data_FC_df[data_FC_df$SampleID<122000 & data_FC_df$SampleID>106000,]$sample <- "2"
data_FC_df[data_FC_df$SampleID<138000 & data_FC_df$SampleID>122000,]$sample <- "3"
data_FC_df[data_FC_df$SampleID<156000 & data_FC_df$SampleID>138000,]$sample <- "4"
data_FC_df[data_FC_df$SampleID<173000 & data_FC_df$SampleID>156000,]$sample <- "5"
data_FC_df[data_FC_df$SampleID<190000 & data_FC_df$SampleID>173000,]$sample <- "6"

data_FC_df[data_FC_df$SampleID<205000 & data_FC_df$SampleID>190000,]$sample <- "1"
data_FC_df[data_FC_df$SampleID<222000 & data_FC_df$SampleID>205000,]$sample <- "2"
data_FC_df[data_FC_df$SampleID<238000 & data_FC_df$SampleID>222000,]$sample <- "3"
data_FC_df[data_FC_df$SampleID>238000,]$sample <- "4"

data_FC_df$sample <- paste(data_FC_df$group,data_FC_df$sample,sep="_")

# plot and check matches
ggplot(data_FC_df,aes(x=SampleID,y=SampleID,colour=sample))+
  geom_point()+#scale_color_manual(values=c("red","blue","cyan","magenta"))+
  #geom_hline(yintercept=c(190000,90000))+
  #geom_vline(xintercept = c(222000,238000))
  theme_classic()

## correct names for protein expression
colnames(data_FC_df) <- sub("^FJComp-","",colnames(data_FC_df))
colnames(data_FC) <- sub("^FJComp-","",colnames(data_FC))
#correct laser channel to protein
channel_prot <- c("GATA-3","viability","CD4","T-bet","CD90.2","TCRb","FOXP3","RORyT","Lin","CD45")
names(channel_prot) <- colnames(data_FC_df)[5:14]

for (i in seq(5,length(channel_prot)+4)){
  colnames(data_FC_df)[i] <- channel_prot[colnames(data_FC_df)[i]]
}


#' select protein marker columns to use for dim. reduc.
# markers should NOT include markers used in gating strategy 
marker_cols <- c("FOXP3","RORyT","GATA-3","T-bet")

# apply arcsinh transformation
# (with standard scale factor of 5 for CyTOF data; alternatively 150 for flow 
# cytometry data; see Bendall et al. 2011, Science, Supplementary Figure S2)
asinh_scale <- 150
data_FC_df[, marker_cols] <- asinh(data_FC_df[, marker_cols] / asinh_scale)

summary(data_FC)

#' #### ---- Running UMAP ---- ####
# possibility of downsampling
# downsample to equal amount of cells in each group
groupn <- c()
for (group in unique(data_FC_df$group)){
  groupn[group] <- nrow(data_FC_df[data_FC_df$group==group,])}
# check group with fewest observations
groupn

down_FC <- downsample(data_FC_df, cat_col = "group")

#n_sub <- 6000
#set.seed(1234)
#ix <- sample(1:length(labels), n_sub)
ix <- rownames(down_FC)

# prepare data for umapr (matrix format required)
data_umap <- down_FC[ix, marker_cols]
data_umap <- as.matrix(data_umap)
dups <- duplicated(data_umap)
data_umap <- data_umap[!dups, ]

umap_emb <- umap(data_umap)

# prepare umap embedding output data for plot
data_plot_umap <- as.data.frame(umap_emb)
colnames(data_plot_umap) <- c("UMAP_1", "UMAP_2")
head(data_plot_umap);dim(data_plot_umap)

# add group labels and merged CI label
data_plot_umap[,"sample"] <- down_FC[ix,]$sample
data_plot_umap[,"group"] <- down_FC[ix,]$group
data_plot_umap$mergedCI <- down_FC$group
data_plot_umap[startsWith(data_plot_umap$mergedCI,"CI"),]$mergedCI <- "CI"

#### Plot each individual group as seperate contour on top of total cells
group_col <- c("CI8"="lightslateblue","CI7-1"="darkslategray4","GF"="gray60")
ggplot(data_plot_umap, aes(x = UMAP_1, y = UMAP_2))+ 
  geom_point_rast(size=0.1, colour="lightgrey")+
  geom_density_2d(data=data_plot_umap[data_plot_umap$group=="GF",],
                  aes(x = UMAP_1, y = UMAP_2, colour=group))+
  scale_color_manual(values = group_col)+
  theme_classic()+
  xlim(c(-13,10))+ylim(c(-10,10))+
  guides(color = guide_legend(override.aes = list(size=3)))#+
  #theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(paste(dato,"CD4cells","UMAP_GF_Contour_sample_plot.pdf",sep="_"), height = 4, width = 4)

ggplot(data_plot_umap, aes(x = UMAP_1, y = UMAP_2))+ 
  geom_point_rast(size=0.1, colour="lightgrey")+
  geom_density_2d(data=data_plot_umap[data_plot_umap$group=="CI7-1",],
                  aes(x = UMAP_1, y = UMAP_2, colour=group))+
  scale_color_manual(values = group_col)+
  theme_classic()+
  xlim(c(-13,10))+ylim(c(-10,10))+
  guides(color = guide_legend(override.aes = list(size=3)))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(paste(dato,"CD4cells","UMAP_CI7-1_Contour_sample_plot.pdf",sep="_"), height = 4, width = 4)

ggplot(data_plot_umap, aes(x = UMAP_1, y = UMAP_2))+ 
  geom_point_rast(size=0.1, colour="lightgrey")+
  geom_density_2d(data=data_plot_umap[data_plot_umap$group=="CI8",],
                  aes(x = UMAP_1, y = UMAP_2, colour=group))+
  scale_color_manual(values = group_col)+
  theme_classic()+
  xlim(c(-13,10))+ylim(c(-10,10))+
  guides(color = guide_legend(override.aes = list(size=3)))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(paste(dato,"CD4cells","UMAP_CI8_Contour_sample_plot.pdf",sep="_"), height = 4, width = 4)


ggplot(data_plot_umap, aes(x = UMAP_1, y = UMAP_2))+ 
  geom_point_rast(colour="lightgrey", size=1)+
  geom_density_2d(data=data_plot_umap[data_plot_umap$mergedCI=="CI",],
                  aes(x = UMAP_1, y = UMAP_2),
                  colour="skyblue3")+
  theme_classic()+
  xlim(c(-13,10))+ylim(c(-10,10))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(paste(dato,"CD4cells","UMAP_CImerged_Contour_sample_plot.pdf",sep="_"), height = 4, width = 4)

#### Colour by individual marker and save pdf versions
for (marker in marker_cols){
  plot_mark <- ggplot(data_plot_umap, aes(x = UMAP_1, y = UMAP_2, colour=data_umap[,marker]))+ 
    geom_point_rast(size=1)+
    scale_colour_viridis_c(option = "plasma")+
    theme_classic()+
    labs(colour=marker)+
    theme(axis.text = element_blank(), axis.ticks = element_blank())
  pdf(paste(dato,"CD4cells","UMAP",marker,"plot.pdf", sep="_"), height = 4, width = 4)
  print(plot_mark)
  dev.off()
}


