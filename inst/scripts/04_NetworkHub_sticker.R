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
nwh_img <- image_read("N_blue.png")

sticker(subplot = nwh_img,
  package = "etworkHub",
  s_width = 0.65,
  s_height = 0.65,
  s_x = 0.45,
  s_y = 1,
  p_x = 1.275,
  p_y = 1,
  p_color = "blue",
  p_size = 17,
  h_size = 1.6,
  h_fill = "tomato1",
  h_color = "#E5E5E5",
  url = "https://github.com/lottawagner/NetworkHub",
  u_size = 3.5,
  u_color = "white") %>% print(sticker())









