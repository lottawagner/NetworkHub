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
nwh_img <- image_read("N_green.png")

sticker(subplot = nwh_img,
  package = "etworkHub",
  s_width = 0.65,
  s_height = 0.65,
  s_x = 0.45,
  s_y = 1,
  p_x = 1.275,
  p_y = 1,
  p_color = "#A2CD5A",
  p_size = 17,
  h_size = 1.6,
  h_fill = "slateblue3",
  h_color = "white") %>% print(sticker())









