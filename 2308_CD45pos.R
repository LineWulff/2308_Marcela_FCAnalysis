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

ggplot(data_FC_df,aes(x=SampleID,y=SampleID,colour=sample))+
  geom_point()+#scale_color_manual(values=c("red","blue","cyan","magenta"))+
  geom_hline(yintercept=c(190000,90000))+
  geom_vline(xintercept = c(222000,238000))



