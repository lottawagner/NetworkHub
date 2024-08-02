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



nwh_img <- "/Users/mac/R/2024_Master_IMBEI/NetworkHub/inst/scripts/netzwerk.png"


sticker(subplot = nwh_img,
  package = "NetworkHub",
  s_width = 1,
  s_height = 1,
  s_x = 1,
  s_y = 0.75
) %>% print(sticker())






