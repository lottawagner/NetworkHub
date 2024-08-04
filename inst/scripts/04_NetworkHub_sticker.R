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
nwh_img <- image_read("netzwerk.png")

sticker(subplot = nwh_img,
  package = "NetworkHub",
  s_width = 0.9,
  s_height = 0.9,
  s_x = 1,
  s_y = 0.7,
  p_x = 1,
  p_y = 1.3,
  p_color = "darkslateblue",
  p_size = 22,
  h_size = 1.6,
  h_fill = "darkseagreen",
  h_color = "lemonchiffon") %>% print(sticker())









