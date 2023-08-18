#' R script for producing outputs in from **CD45+ cells** from Marcela's project
#' Author: Line Wulff
#' Date (created): 23-08-18

#### ---- Initiate libraries ---- ####
library(flowCore)
library(FlowSOM)
library(Rtsne)
library(ggplot2)

#### ---- variables used throughout script ---- ####
projdir <- getwd()

#### ---- Read in data, inspect and transform if necessary ---- ####
data_FC <- flowCore::exprs(flowCore::read.FCS(
  filename = paste(projdir,"/inputdata/concat_lymphoidcells.fcs", sep=""), 
  transformation = FALSE, truncate_max_range = FALSE))

head(data_FC)
dim(data_FC)

data_FC_df <- as.data.frame(data_FC)
data_FC_df$group <- "CI8"
data_FC_df[data_FC_df$SampleID>190000,]$group <- "GF"
data_FC_df[data_FC_df$SampleID<90000,]$group <- "CI7-1"

data_FC_df$sample <- NA
data_FC_df[data_FC_df$SampleID<25000,]$sample <- "1"
data_FC_df[data_FC_df$SampleID<40000 & data_FC_df$SampleID>25000,]$group <- "2"
data_FC_df[data_FC_df$SampleID<58000 & data_FC_df$SampleID>40000,]$group <- "3"
data_FC_df[data_FC_df$SampleID<74000 & data_FC_df$SampleID>58000,]$group <- "4"
data_FC_df[data_FC_df$SampleID<90000 & data_FC_df$SampleID>74000,]$group <- "5"

data_FC_df[data_FC_df$SampleID<40000 & data_FC_df$SampleID>90000,]$group <- "2"
data_FC_df[data_FC_df$SampleID<58000 & data_FC_df$SampleID>40000,]$group <- "3"
data_FC_df[data_FC_df$SampleID<74000 & data_FC_df$SampleID>58000,]$group <- "4"
data_FC_df[data_FC_df$SampleID<90000 & data_FC_df$SampleID>74000,]$group <- "5"

ggplot(data_FC_df,aes(x=SampleID,y=SampleID,colour=group))+
  geom_point()+scale_color_manual(values=c("red","blue","cyan","magenta"))+
  geom_hline(yintercept=c(250000,190000,90000))+
  geom_vline(xintercept = c(160000,58000))



