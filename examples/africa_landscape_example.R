library(raster)
library(tidyverse)
library(parallel)
library(pracma)
library(contactsimulator)

# Create a landscape
malawi<- raster(xmn=0, xmx=10000, ymn=0, ymx=10000)
res(malawi)<- 70
values(malawi)<- 0
