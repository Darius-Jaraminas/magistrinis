library("vrtest")
library("ggplot2")
library("foreach")
library("doParallel")
library("reshape")
library("cubature")

source("Functions.R")


y <- rnorm(500)

ww2 <- calc.wj(n = 500*2, fun = f.theta)

prp <- Sys.time()
wbzint2 <- wb.z2int(y = y, B = 1000, prob = c(0.95), w = ww2)
Sys.time() - prp