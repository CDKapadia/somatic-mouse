library(ggplot2)
library(cowplot)
library(patchwork)
library(ggsci)
library(ggpubr)
library(purrr)
library(data.table)

myColors <- c("#882255","#D62B81",    "#AA4499", "#E259CB",    "#CC6677", "#F7798E",
              "#332288","#6246F1",    "#117733", "#14BF4D",     "#44AA99", "#45E8CD")
names(myColors) <- c("MD7180", "MD7180_MPP", "MD7181","MD7181_MPP", "MD7182","MD7182_MPP", 
                     "MD7183","MD7183_MPP","MD7184","MD7184_MPP", "MD7185","MD7185_MPP")
