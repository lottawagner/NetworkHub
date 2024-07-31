# Creation of the NetworkHub sticker using HexSticker

#install.packages("hexSticker")
#install.packages("magick")
#install.packages("sysfonts")
#install.packages("tidyverse")

library(hexSticker)
library(magick)
library(sysfonts)
library(tidyverse)

nwh_img <- image_read("netzwerk.png")

sticker(
  subplot = nwh_img,
  package = "NetworkHub"
) %>%print(sticker())

