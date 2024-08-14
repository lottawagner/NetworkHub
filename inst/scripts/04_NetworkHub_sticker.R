# Creation of the NetworkHub sticker using HexSticker

if (!requireNamespace("hexSticker", quietly = TRUE))
  install.packages("hexSticker")
if (!requireNamespace("magick", quietly = TRUE))
  install.packages("magick")
if (!requireNamespace("sysfonts", quietly = TRUE))
  install.packages("sysfonts")
if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

library(hexSticker)
library(magick)
library(sysfonts)
library(tidyverse)



setwd("/Users/mac/R/2024_Master_IMBEI/NetworkHub/inst/scripts")
nwh_img <- image_read("N_orange.png")


sticker(subplot = nwh_img,
  package = "etwork  ub",
  s_width = 1.6,
  s_height = 1.6,
  s_x = 0.98,
  s_y = 1,
  p_x = 1.2,
  p_y = 1,
  p_color = "orange",
  p_size = 64,
  h_size = 1,
  h_fill = "#d6dffaff",
  h_color = "#91a9f2ff",
  url = "www.bioconductor.org",
  u_size = 18,
  u_color = "#473c8bff",
  dpi = 1200) %>% print(sticker())

##ff5100ff








