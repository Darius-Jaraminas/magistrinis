rm(list = ls())
setwd("D:/Dropbox/Dokumentai/Studijos/Magistras/Magistrinis/R kodas")

library("vrtest")
library("ggplot2")
library("foreach")
library("doParallel")

source("Functions.R")

plot(kernel1, 0, 5)