rm(list = ls())
setwd("D:/Dropbox/Dokumentai/Studijos/Magistras/Magistrinis/R kodas")

library("vrtest")
library("ggplot2")
library("foreach")
library("doParallel")
library("reshape")
library("cubature")

source("Functions.R")


########################################################
#################### Empirinis reiksmingumas ###########
########################################################

d <- data.frame(test = c("aq5", "aq20", "avr", "zint", "z2int",
                         "rz1", "rz1m2", "rz1m5", "rz1m10"),
                distribution = c("Chi", "Chi", "N", "N", "z.sq",
                                 "N", "Chi.2", "Chi.5", "Chi.10"))
ww <- calc.wj(n = 1000, fun = theta)
ww2 <- calc.wj(n = 1000*2, fun = f.theta)
# Rprof()
# sn <- sim.null(N = 20, t = 200, d = d)
# Rprof(NULL)
# summaryRprof()
# debugonce(t.size)
# debugonce(sim.null)
# debugonce(null.tests)

### z_T(theta)
prp <- Sys.time()
tt <- c(25, 50, 100, 200, 500, 1000)
l <- list()
s <- NULL
#Rprof(filename = "Rprof.out", append = FALSE, interval = 0.01)
for (i in 1:length(tt)){
  sn <- sim.null(N = 2000, t = tt[i], d = d, qzsq = 3.25,
                 z.function = zh,
                 z.int.function = zh.int,
                 z.vec.function = zh.vec,
                 z2int.function = z.sq.int2, w = ww, w2 = ww2)
  l[[paste(tt[i])]] <- sn
  ss <- sn$size
  print(ss)
  ss <- cbind(t = tt[i], ss)
  s <- rbind(s, ss)
  save(l, file = "Rdata/zh.size.detailed2.Rdata")
  save(s, file = "Rdata/zh.size2.Rdata")
}
#Rprof(NULL)
#summaryRprof()
Sys.time() - prp

### tilde{z}_T(theta) - atspari heteroskedastiðkumui statistika
prp <- Sys.time()
tt <- c(25, 50, 100, 200, 500, 1000)
l <- list()
s <- NULL
#Rprof(filename = "Rprof.out", append = FALSE, interval = 0.01)
for (i in 1:length(tt)){
  sn <- sim.null(N = 2000, t = tt[i], d = d, qzsq = 3.25,
                 z.function = z,
                 z.int.function = z.int,
                 z.vec.function = z.vec,
                 z2int.function = z.sq.int, w = ww, w2 = ww2)
  l[[paste(tt[i])]] <- sn
  ss <- sn$size
  print(ss)
  ss <- cbind(t = tt[i], ss)
  s <- rbind(s, ss)
  save(l, file = "Rdata/z.size.detailed2.Rdata")
  save(s, file = "Rdata/z.size2.Rdata")
}
#Rprof(NULL)
#summaryRprof()
Sys.time() - prp

write.csv2(s, "z.size2.csv")

get.empirical.critical.values(object = l, dgp = c("normal", "garch"),
                              alpha = 0.05)


########################################################
#################### Testu galia #######################
########################################################

d <- load(file = "Rdata/dzsq.Rdata")
dzsq <- get(d)
ww <- calc.wj(n = 5000, fun = theta)
ww2 <- calc.wj(n = 5000*2, fun = f.theta)
nsm <- load(file = "Rdata/zh.size.detailed2.Rdata")
l <- get(nsm)
cv <- get.empirical.critical.values(object = l,
                                    dgp = c("normal"), alpha = 0.05)
B <- 500
t <- 200
y <- rnorm(t)
y <- gen.alternative(n = 500, method = "ma12")
# debugonce(alternative.tests2)
# at <- alternative.tests2(y = y, B = B, w = ww, w2 = ww2, dzsq = dzsq)
# at

# cl <- makeCluster(4)
# registerDoParallel(cl)
# prp <- Sys.time()
# (gs <- Gen.Spec(y = y, B = B))
# Sys.time() - prp
# stopCluster(cl)

# debugonce(alternative.tests2)
# prp <- Sys.time()
# (at <- alternative.tests2(y = y, B = B, w = ww, w2 = ww2, dzsq = dzsq,
#                           result = "value", cv = cv))
# Sys.time() - prp


# spow <- sim.a2(N = 500, B = 1000, alt = c("ma12"), tt = c(25),
#                cv = cv)
# write.csv2(spow$power, "power.csv")


smth <- sim.a2(N = 500, B = 1000, cv = cv)
save(smth, file = "Rdata/power.ma2.Rdata")
write.csv2(smth$power, "power.ma2.csv")



  