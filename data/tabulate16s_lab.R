rm(list=ls())
library(stringr)

data_raw <- read.csv("isolates.csv",header=T,sep=",")
data <- data_raw[!is.na(data_raw$Fermenter),]
row.names(data)<- NULL