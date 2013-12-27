rm(list = ls())
setwd("D:/Dropbox/Dokumentai/Studijos/Magistras/Magistrinis/R kodas")

library("vrtest")
library("ggplot2")
library("foreach")
library("doParallel")
library("reshape")
library("cubature")
library("fda")

source("Functions.R")

debugonce(z.matrix)
z(theta = 0.5, rnorm(100))

zm <- z.matrix2(3:1)
zz <- matrix(2:4, nrow = 3, ncol = 3, byrow = TRUE)
zz * zm

debugonce(z.sq.int2)
u <- rnorm(2000)
prp <- Sys.time()
z.sq.int2(u)
Sys.time() - prp

prp <- Sys.time()
z.sq.int(u)
Sys.time() - prpdata(exrates)


y <- generate.y(t = 200, alpha = 0, p = 3)
debugonce(Auto.VR)
Auto.VR(y)
plot.ts(y)

dl1 <- DLtest(exrates[,1], p = 1)
debugonce(p.hat)
p.hat(y = rnorm(200), lags = 10)
debugonce(Auto.DL.test)
Auto.DL.test(y = y, lags = 10)



debugonce(z.matrix)
z.m <- z.matrix(u = 1:1000)

y <- generate.y(t = 200, alpha = 0.01, p = 7)
ggplot(aes(y = y, x = x), data = data.frame(y = y, x = 1:length(y))) +
  geom_line(color = "blue", size = 1)
# debugonce(z)
(zt <- z(u = y, theta = 0.5))

# debugonce(plot.z)
y <- generate.y(t = 400, alpha = -0.6, p = 1)
plot.z(u = y, n = 400)

# Rprof(filename = "Rprof.out", append = FALSE, interval = 0.02)
# zm <- z.matrix(u = 1:2000)
# Rprof(NULL)
# summaryRprof()

debugonce(z.int)
z.int(u = y)

debugonce(monte.zint)
mzint <- monte.zint(t = 200, alpha = 0.5, p = 1, n = 200, B = 300, seed = 777)
ggplot(aes(x = mzint), data = data.frame(mzint = mzint)) + geom_histogram()

bzint <- wild.boot(u = y, B = 500, seed = 7777, test = z.int)
d <- data.frame(boot.dist = bzint$boot.dist, zint = bzint$stat)
p <- ggplot(aes(x = boot.dist), data =d ) + geom_histogram() 
p +  geom_vline(aes(xintercept = zint), color = "green", size = 1.5)



prp <- Sys.time()
bzint <- wild.boot.par(u = y, B = 500, seed = 7777, test = z.int,
                   export = list(z.matrix = z.matrix, wj = wj, theta = theta))
Sys.time() - prp


prp <- Sys.time()
bzint <- wild.boot(u = y, B = 500, seed = 7777, test = z.int)
Sys.time() - prp


prp <- Sys.time()
mzint <- monte.zint.par(t = 200, alpha = 0.5, p = 1, n = 100,
                        B = 200, seed = 777,
                       export = list(generate.y = generate.y,
                                     wild.boot = wild.boot,
                                     z.int = z.int, z.matrix = z.matrix,
                                     wj = wj, theta = theta))
Sys.time() - prp


prp <- Sys.time()
debugonce(monte.zint)
mzint <- monte.zint(t = 200, alpha = 0.1, p = 1, n = 100, B = 100, seed = 777)
Sys.time() - prp


#cvm and ks with projection ####

y <- generate.y(t = 200, alpha = 0.6, p = 4)

debugonce(dl.mod)
dl.mod(y, p = 1, result = "cp")

e <- list(generate.y = generate.y, projection.dl = projection.dl)
wb.dl <- wild.boot.par(u = y, B = 100, test = dl.mod, export = e,
                          make.cluster = TRUE, p = 1, result = "cp")


adm <- auto.dl.mod(y = y, lags = 10)
plot.auto.dl.mod(adm)

adm2 <- auto.dl.mod(y = y, lags = 10, projection = FALSE)
plot.auto.dl.mod(adm2)

#### testing zint under H0
n <- 3000
y <- rnorm(n)

# debugonce(z.int)
# z.int(u = y, interval = c(0,1 - 0.001))

m.zint <- foreach(i = 1:2000, .combine = c) %do% {
  y <- rnorm(n)
  z.int(u = y, interval = c(0,1 - 0.0001))$stat.std
}

var(m.zint)
hist(m.zint)
z.int(u = y, interval = c(0,1 - 0.0001))$sigma

d <- data.frame(m.zint = m.zint)
p <- ggplot(aes(x = m.zint), data = d ) + 
  geom_histogram(binwidth = 0.1, color = "Blue") 
p + geom_density(aes(y = ..scaled..), adjust=1/5, color = "Blue") 


f <- function(x, j){
  sqrt(1-x^2) * x^j
}

te <- seq(from = 0, to = 1, length.out= 1000)

dt1 <- foreach(j = 1:20, .combine = cbind) %do% {
  f(x = te, j = j)
}

dt1 <- data.frame(cbind(te, dt1))
dt1 <- melt(dt1, id.vars = "te")

pdf("functions.pdf")
ggplot(aes(x = te, y = value, group = variable), data = dt1) +
  geom_line(color = "Blue") + facet_wrap(~variable)
dev.off()



w <- NULL
for (i in 0:9){
  w <- c(w, wj(j = i, interval = c(0,1-1e-4)))
}
w <- data.frame(w = w, t = 0:9)
ggplot(aes(x = t, y = w), data = w) + geom_point(size = 4) +
  theme(axis.text = element_text(size =rel(1)),
        axis.title = element_text(size =rel(1)),
        title = element_text(size =rel(2))) +
  scale_x_continuous(breaks = c(0:9)) + xlab("Indeksas j") + ylab("Svoriai") +
  ggtitle("Svoriai w")



### multivariate

debugonce(phi.j)
f <- phi.j(j = 5)

adaptIntegrate(f, lowerLimit = c(0,0), upperLimit = c(1,1))
adaptIntegrate(f, lowerLimit = c(0), upperLimit = c(1))


debugonce(z.multi)
zm <- z.multi(theta = c(0.1, -0.2), u = rnorm(500), rtheta = TRUE)

n <- 4000
pchi <- numeric(200)
for (i in 1:200){
  zm <- z.multi(u = rnorm(n), rtheta = TRUE)
  pchi[i] <- zm$chi.pvalue
}

mean(pchi < 0.05, na.rm = TRUE)

#### square
y <- generate.y(t = 500, alpha = 0.2, p = 1)
y <- gen.null(100, "garch")
wbzsq <- wb.z.sq.int(u = y, B = 100)

(trz <- rz(u = y, m = 50, alpha = 0.05))

(trzsq <- rz(u = y, m = 5, alpha = 0.05, test.f = z.sq))
debugonce(all.tests)
at <- all.tests(y, B = 100)

sn <- sim.null(N = 1000, B = 200, t = 250)

zs <- NULL
for (i in 1:1000){
  zs <- c(zs, z.sq.int2(rnorm(1000)))
}

zs1 <- zs
zs <- zs/0.9999

plot(density(zs, bw = 0.01))
x <- seq(0, 5, length=200)
hx <- dchisq(x, df = 1)
lines(x, hx, type="l", lty=2, xlab="x value",
     ylab="Density", main="Comparison of t Distributions")
abline(v = quantile(zs, 0.95), col = "green")
abline(v = qchisq(0.95, df = 1), col = "red")



### fun
N <- 1000
t <- NULL
for (i in 1:10000){
  e <- rnorm(N)
  e <- e^2
  cs <- cumsum(e)
  cs <- cs/(1:N)
  m <- max(cs)
  t <- c(t, m)
}

plot(density(t))

quantile(t, 0.95)
quantile(t, 0.9)

N <- 200
zsq <- NULL
for (i in 1:200){
  u <- rnorm(N)
  zt <- z.sq.int(u = u)
  zsq <- c(zsq, zt)
}

plot(density(zsq))

quantile(zsq, 0.95)
quantile(zsq, 0.9)
qchisq(0.95, df = 1)
qchisq(0.9, df = 1)




# debugonce(z.cor.sim)
th <- seq(0.1, 0.9, by = 0.1)

(co <- z.cov.sim(theta = th, N = 100, t = 500, dgp = "normal"))
dif <- (z.cov(th) - co)/z.cov(th)

th <- c(0.3, 0.7)
tt <- seq(50, 2000, by = 25)
covs <- NULL
for (i in 1:length(tt)){
  co <- z.cov.sim(theta = th, N = 200, t = tt[i], dgp = "garch", z.function = zh.vec)
  covs <- c(covs, co[1,2])
}
plot(tt, covs - z.cov(th)[1,2], type = "l")
plot(tt, sqrt(tt) * (covs - z.cov(th)[1,2]), type = "l")
plot(tt,(covs - z.cov(th)[1,2]) / sqrt(tt), type = "l")

plot(tt, abs(covs - z.cov(th)[1,2]), type = "l", ylim = c(0, 0.15))
lines(tt, 1/sqrt(tt))

cov.dt <- data.frame(t = tt, dif = abs(covs - z.cov(th)[1,2]), t.sq = 1/sqrt(tt))
names(cov.dt) <- c("Laikas", "Absoliutus skirtumas", "1/(t^(1/2))")
cov.dt <- melt(cov.dt, id = "Laikas")
names(cov.dt) <- c("Laikas", "Kintamasis", "v")


g <- ggplot(aes(x = Laikas, group = Kintamasis), data = cov.dt) +
  geom_line(aes(y = v, linetype = Kintamasis), size = 1) +
  xlab("Laikas") +
  ylab(expression(abs(gamma - hat(gamma)))) +
  ggtitle("Kovariacijos konvergavimas") +
  theme(axis.text = element_text(size = rel(2)),
        axis.title = element_text(size = rel(2.5)),
        plot.title = element_text(size = rel(2.5)),
        legend.text = element_text(size = rel(1.3)),
        legend.title = element_text(size = rel(1.8)))
g <- g + guides(linetype = guide_legend())
g


### pvz zint trûkumai

y <- generate.y(t = 200, alpha = -0.7, p = 2)

ggplot(aes(y = y, x = x), data = data.frame(y = y, x = 1:length(y))) +
  geom_line(color = "blue", size = 0.8)
plot.ts(y, main = "Leiko eilutë", ylab = "y", xlab = "Laikas")
acf(y, main = "Autokoreliacijos", ylab = "Autokoreliacija", xlab = "Vëlavimas",
    ci.col = "white")
p <- 0
while (p < 0.05){
  y <- generate.y(t = 200, alpha = -0.5, p = 2)
  zint <- z.int(y)
  p <- zint$pval
}

p <- numeric(1000)
for (i in 1:1000){
  y <- generate.y(t = 200, alpha = -0.5, p = 2)
  zint <- z.int(y)
  p[i] <- zint$stat
}

pp <- numeric(1000)
for (i in 1:1000){
  y <- rnorm(200)
  zint <- z.int(y)
  pp[i] <- zint$stat
}


# just a hunch
n <- 3
w <- c(0.69, 0.22, 0.09)
pp <- numeric(10000)
for (i in 1:10000){
  pp[i] <-  sum(w * rchisq(n, df = 1))
}

plot(density(zs))
lines(density(pp)$x, density(pp)$y, type="l", lty=2, xlab="x value",
      ylab="Density", main="Comparison of t Distributions")
abline(v = quantile(zs, 0.95), col = "green")
abline(v = quantile(pp, 0.95), col = "red")



##### random theta
n <- 300
b <- 2000
K <- 1000
r <- list()
r[["z"]] <- numeric(K)
r[["zp"]] <- numeric(K)
r[["zbp"]] <- numeric(K)
r[["zbq1"]] <- numeric(K)
r[["zbq2"]] <- numeric(K)
for (i in 1:K){
  theta <- runif(1, min = -1, max = 1)
  # u <- generate.y(t = n, alpha = -0.4, p = 1)
  u <- rnorm(n)
  r[["z"]][i] <- z(theta = theta, u = u)
  r[["zp"]][i] <- 2 * pnorm(abs( r[["z"]][i]), lower.tail=FALSE)
  zb <- numeric(b)
  for (j in 1:b){
    u1 <- rnorm(n)
    zb[j] <- z(theta = theta, u = u1)
  }
  r[["zbq1"]][i] <- quantile(zb, 0.025)
  r[["zbq2"]][i] <- quantile(zb, 0.975)
  r[["zbp"]][i] <- mean(zb > abs(r[["z"]][i])) + mean(zb < -abs(r[["z"]][i]))
}

plot(density(r[["z"]]), xlim = c(-max(abs(r[["z"]])), max(abs(r[["z"]]))), 
     ylim = c(0, 0.41))
x <- seq(-5, 5, length=500)
hx <- dnorm(x)
lines(x, hx, type="l", lty=2, xlab="x value",
      ylab="Density", main="Comparison of t Distributions")
abline(v = quantile(r[["z"]], 0.975), col = "green")
abline(v = quantile(r[["z"]], 0.025), col = "green")
abline(v = qnorm(0.975), col = "red")
abline(v = qnorm(0.025), col = "red")

(apower <- mean(r[["z"]] < qnorm(0.025)) + mean(r[["z"]] > qnorm(0.975)))


#### 
debugonce(gen.alternative)
y <- gen.alternative(n = 200, method = "ar12")
plot.ts(y)
acf(y)
pacf(y)

debugonce(rz)
debugonce(rz.boot)
rz.boot(u = y, m = 1, alpha = 0.05, B = 1000)
wb.z.sq.int(u = y, B = 500)
debugonce(rzsq.boot)
rzsq.boot(u = y, m = 1, alpha = 0.05, B = 1000)


prp <- Sys.time()
p <- alternative.tests(y = y, B = 1000)
Sys.time() - prp

debugonce(z)
zzh <- NULL
for (i in 1:1000){
  y <- gen.null(n = 500, "garch")
  zzh <- c(zzh, z(0.5, y))
}

plot(density(zzh), xlim = c(-max(abs(zzh)), max(abs(zzh))), 
     ylim = c(0, 0.41))
x <- seq(-5, 5, length=500)
hx <- dnorm(x)
lines(x, hx, type="l", lty=2, xlab="x value",
      ylab="Density", main="Comparison of t Distributions")

#### fda

splinebasis = create.bspline.basis(c(0,10), 13)
splinebasis <- create.monomial.basis(c(0,1), 4)
plot(splinebasis)


day.5 <- seq(0,365,5)
logprecav <- CanadianWeather$dailyAv[, , "log10precip"]
dayrange = c(0,365)
daybasis = create.fourier.basis(dayrange, 365)
logprecres.fd <- smooth.basis(day.5, logprecav, daybasis)





##
######## Simulated data example 1: a simple regression smooth  ########
##
#  Warning:  In this and all simulated data examples, your results
#  probably won't be the same as we saw when we ran the example because
#  random numbers depend on the seed value in effect at the time of the
#  analysis.
#
#  Set up 51 observation points equally spaced between 0 and 1
n = 51
argvals = seq(0,1,len=n)
#  The true curve values are sine function values with period 1/2
x = sin(4*pi*argvals)
#  Add independent Gaussian errors with std. dev. 0.2 to the true values
sigerr = 0.2
y = x + rnorm(x)*sigerr
#  When we ran this code, we got these values of y (rounded to two
#  decimals):
y = c(0.27,  0.05,  0.58,  0.91,  1.07,  0.98,  0.54,  0.94,  1.13,  0.64,
      0.64,  0.60,  0.24,  0.15, -0.20, -0.63, -0.40, -1.22, -1.11, -0.76,
      -1.11, -0.69, -0.54, -0.50, -0.35, -0.15,  0.27,  0.35,  0.65,  0.75,
      0.75,  0.91,  1.04,  1.04,  1.04,  0.46,  0.30, -0.01, -0.19, -0.42,
      -0.63, -0.78, -1.01, -1.08, -0.91, -0.92, -0.72, -0.84, -0.38, -0.23,
      0.02)
#  Set up a B-spline basis system of order 4 (piecewise cubic) and with
#  knots at 0, 0.1, ..., 0.9 and 1.0, and plot the basis functions
nbasis = 13
basisobj = create.bspline.basis(c(0,1),nbasis)
plot(basisobj)
#  Smooth the data, outputting only the functional data object for the
#  fitted curve.  Note that in this simple case we can supply the basis
#  object as the "fdParobj" parameter
ys = smooth.basis(argvals=argvals, y=y, fdParobj=basisobj)
Ys = smooth.basis(argvals=argvals, y=y, fdParobj=basisobj,
                  returnMatrix=TRUE)
# Ys[[7]] = Ys$y2cMap is sparse;  everything else is the same


all.equal(ys[-7], Ys[-7])


xfd = ys$fd
Xfd = Ys$fd

#  Plot the curve along with the data
plotfit.fd(y, argvals, xfd)
#  Compute the root-mean-squared-error (RMSE) of the fit relative to the
#  truth
RMSE = sqrt(mean((eval.fd(argvals, xfd) - x)^2))
print(RMSE)  #  We obtained 0.069
#  RMSE = 0.069 seems good relative to the standard error of 0.2.
#  Range through numbers of basis functions from 4 to 12 to see if we
#  can do better.  We want the best RMSE, but we also want the smallest
#  number of basis functions, which in this case is the degrees of
#  freedom for error (df).  Small df implies a stable estimate.
#  Note: 4 basis functions is as small as we can use without changing the
#  order of the spline.  Also display the gcv statistic to see what it
#  likes.
for (nbasis in 4:12) {
  basisobj = create.bspline.basis(c(0,1),nbasis)
  ys = smooth.basis(argvals, y, basisobj)
  xfd = ys$fd
  gcv = ys$gcv
  RMSE = sqrt(mean((eval.fd(argvals, xfd) - x)^2))
  # progress report:
  cat(paste(nbasis,round(RMSE,3),round(gcv,3),"\n"))
}
#  We got RMSE = 0.062 for 10 basis functions as optimal, but gcv liked
#  almost the same thing, namely 9 basis functions.  Both RMSE and gcv
#  agreed emphatically that 7 or fewer basis functions was not enough.
#  Unlike RMSE, however, gcv does not depend on knowing the truth.
#  Plot the result for 10 basis functions along with "*" at the true
#  values
nbasis = 10
basisobj = create.bspline.basis(c(0,1),10)
xfd = smooth.basis(argvals, y, basisobj)$fd
plotfit.fd(y, argvals, xfd)
points(argvals,x, pch="*")
#  Homework:
#  Repeat all this with various values of sigerr and various values of n

n <- 101
t <- 15000
N <- 1
x <- seq(0, 1, length = n)
y <- matrix(0, nrow = n, ncol = N)
prp <- Sys.time()
for (j in 1:N){
  u <- rnorm(t)
  y[, j] <- zh.vec(theta = x, u = u)
}
Sys.time() - prp

prp <- Sys.time()
yy <- zh.vec(theta = x, u = u)
Sys.time() - prp

(round(yy, 2) == round(yyy, 2))



plot.ts(y)

nbasis <- 100
basisobj <- create.bspline.basis(c(0, 1), nbasis)
plot(basisobj)
#  Smooth the data, outputting only the functional data object for the
#  fitted curve.  Note that in this simple case we can supply the basis
#  object as the "fdParobj" parameter
ys <- smooth.basis(argvals = x, y = y, fdParobj = basisobj)
xfd <- ys$fd

#  Plot the curve along with the data
plotfit.fd(y, x, xfd)
#  Compute the root-mean-squared-error (RMSE) of the fit relative to the
#  truth
RMSE <- sqrt(mean((eval.fd(x, xfd) - x)^2))
print(RMSE)  #  We obtained 0.069

for (nbasis in 4:12) {
  basisobj <- create.bspline.basis(c(0, 1), nbasis)
  ys <- smooth.basis(x, y, basisobj)
  xfd <- ys$fd
  gcv <- ys$gcv
  RMSE <- sqrt(mean((eval.fd(x, xfd) - x)^2))
  # progress report:
  cat(paste(nbasis,round(RMSE,3),round(gcv,3),"\n"))
}

z.var <- var.fd(xfd)
z.var.mat <- eval.bifd(x, x, z.var)

kovariacija <- z.var.mat
theta1 <- x
theta2 <- x

def.par <- par(no.readonly = TRUE)
par(mar = c(1, 0.5, 1, 0.5))
layout(matrix(c(1,2,1,2), 2, 2, byrow = TRUE))
persp(theta1, theta2, kovariacija, theta = -60, phi = 50, axes = TRUE,
      zlab = NULL, main = "Kovariacijø funkcija")
# persp(x, x, z.var.mat, theta = 20, phi = 0, axes = TRUE,
#       ticktype = "detailed")
contour(x, x, z.var.mat, nlevels = 14, lwd = 2, labcex = 1,
        main = "Kovariacijø funkcijos konturas", xlab = "theta1",
        ylab = "theta2")
par(def.par)

debugonce(pca.fd)
z.pc <- pca.fd(xfd, 5)
ev <- z.pc$values

plot.ts(ev[1:5])

df.ev <- data.frame(ev = ev[1:5], t = 1:5)
ggplot(aes(y = ev, x = t), data = df.ev) +
  geom_line(size = 0.8) +
  geom_point(size = 4) +
  geom_text(aes(label = round(ev, 3), x = t + 0.25, y = ev + 0.02)) +
  ggtitle("Kovariacijos funkcijos tikrinës reikðmës") +
  theme(axis.text = element_text(size =rel(2), color = "black"),
        axis.title = element_text(size =rel(1.5)),
        title = element_text(size =rel(2))) +
  xlab("i") + ylab(expression(lambda[i]))



debugonce(z.cov)
co <- z.cov(x)

dif <- abs(co - z.var.mat)
max(dif, na.rm = TRUE)


N <- 1000000
nn <- length(ev)
pp <- numeric(N)
for (i in 1:N){
  pp[i] <-  sum(ev * rchisq(nn, df = 1))
}


plot(density(pp, from = 0.001), xlim = c(0, 4))
(q1 <- quantile(pp, 0.95))
abline(v = q1, col = "green")


xx <- density(pp, from = 0.07)$x
xx <- xx[1:150]
df.ev <- data.frame(x = xx,
                    dzsq = density(pp, from = 0.07)$y[1:150],
                    chi = dchisq(xx, df = 1))

df.ev <- melt(df.ev, id = x)
names(df.ev)[2] <- "Kintamasis"


ggplot(aes(y = value, x = x, group = Kintamasis), data = df.ev) +
  geom_line(size = 1.5, aes(linetype = Kintamasis)) +
  geom_vline(xintercept = q1, size = 1.5) +
  geom_vline(xintercept = qchisq(0.95, df = 1), linetype = 2, size = 1.5) +
  ggtitle(bquote(.("Asimptotinio") ~ z[sq] ~ .("pasiskirstymo aproksimacija"))) +
  theme(axis.text = element_text(size =rel(2), color = "black"),
        axis.title = element_text(size = rel(1.5)),
        title = element_text(size =rel(2)),
        legend.text = element_text(size =rel(2.5)),
        legend.key.size = unit(2, "cm")) +
  xlab(expression(x)) + ylab(expression(f[x](x))) + 
  scale_linetype_discrete(labels = c(bquote(z[sq] ~ .("tankis")),
                                     bquote(chi[1]^2 ~ .("tankis"))),
                          name = "") + 
  geom_segment(aes(x=0,y=0,
                   xend=0.07,
                   yend=density(pp, from = 0.07)$y[1]), size = 1.5)

N <- 1000
zsq <- NULL
for (i in 1:100){
  u <- rnorm(N)
  zt <- z.sq.int2(u = u)
  zsq <- c(zsq, zt)
}

plot(density(zsq, adjust =0.4), xlim = c(-0.5, 5), ylim = c(0, 1.1))
lines(density(pp)$x, density(pp)$y, type="l", lty=2, xlab="x value",
      ylab="Density", main="Comparison of t Distributions")
abline(v = quantile(zsq, 0.95), col = "green")
abline(v = quantile(pp, 0.95), col = "red")


### example plots

xx <- seq(0.01, 18, by = 0.05)
df.ev1 <- data.frame(x = xx, Kintamasis = "Asimptotinis",
                    n1 = dchisq(xx, df = 4))
df.ev2 <- data.frame(x = xx + 2, Kintamasis = "Empirinis",
                    n1 = dchisq(xx, df = 4))
df.ev <- rbind(df.ev1, df.ev2)

names(df.ev)[3] <- "value"

shade <- subset(df.ev, x >= qchisq(0.95, df = 4))
sh1 <- shade[shade$Kintamasis == "Asimptotinis", - 2]
sh2 <- shade[shade$Kintamasis != "Asimptotinis", - 2]
sh1 <- rbind(c(qchisq(0.95, df = 4), 0), sh1, c(sh1[nrow(sh1), 1], 0))
sh2 <- rbind(c(qchisq(0.95, df = 4), 0), sh2, c(sh2[nrow(sh2), 1], 0))


g <- ggplot() +
  geom_line(size = 1.5, 
            aes(y = value, x = x, group = Kintamasis,
                linetype = Kintamasis),
            data = df.ev) +
  geom_segment(aes(x=qchisq(0.95, df = 4),y=0,
                   xend=qchisq(0.95, df = 4),
                   yend=dchisq(qchisq(0.95, df = 4) - 2, df = 4)))
g <- g  + geom_polygon(data = sh2, aes(x = x, y = value), alpha = 0.6) +
  geom_text(aes(x = qchisq(0.95, df = 4), y = -0.002,
                label = "Q95%"), size = 6)+
  ggtitle("Empirinis reikðmingumas") +
  theme(axis.text = element_text(size =rel(1.5), color = "black"),
        axis.title = element_text(size =rel(1)),
        title = element_text(size =rel(2)),
        legend.text = element_text(size =rel(1.7)),
        legend.key.size = unit(2, "cm")) +
  xlab("x") + ylab(expression(f[x](x))) + 
  scale_linetype_discrete(labels = c("Asimptotinis tankis", "Empirinis tankis"),
                          name = "")
#########################

xx <- seq(0.01, 15, by = 0.05)
df.ev1 <- data.frame(x = xx, Kintamasis = "Asimptotinis H0",
                     n1 = dchisq(xx, df = 4))
df.ev2 <- data.frame(x = xx + 8, Kintamasis = "Empirinis H1",
                     n1 = dchisq(xx, df = 4))
df.ev <- rbind(df.ev1, df.ev2)

names(df.ev)[3] <- "value"

shade <- subset(df.ev, x >= qchisq(0.95, df = 4))
sh1 <- shade[shade$Kintamasis == "Asimptotinis H0", - 2]
sh2 <- shade[shade$Kintamasis != "Asimptotinis H0", - 2]
sh1 <- rbind(c(qchisq(0.95, df = 4), 0), sh1, c(sh1[nrow(sh1), 1], 0))
sh2 <- rbind(c(qchisq(0.95, df = 4), 0), sh2, c(sh2[nrow(sh2), 1], 0))


g <- ggplot() +
  geom_line(size = 1.5, 
            aes(y = value, x = x, group = Kintamasis,
                linetype = Kintamasis),
            data = df.ev) +
  geom_segment(aes(x=qchisq(0.95, df = 4),y=0,
                   xend=qchisq(0.95, df = 4),yend=dchisq(qchisq(0.95, df = 4) - 7, df = 4)))
g <- g  + geom_polygon(data = sh2, aes(x = x, y = value), alpha = 0.6) +
  geom_text(aes(label = "Q95%", x = qchisq(0.95, df = 4), y = -0.002),
            size = 6)+
  ggtitle("Testo galia") +
  theme(axis.text = element_text(size =rel(1.5), color = "black"),
        axis.title = element_text(size =rel(1)),
        title = element_text(size =rel(2)),
        legend.text = element_text(size =rel(1.7)),
        legend.key.size = unit(2, "cm")) +
  xlab(expression(x)) + ylab(expression(f[x](x))) + 
  scale_linetype_discrete(labels = c(bquote("Tankis prie" ~ H[0]),
                                     bquote("Tankis prie" ~ H[1])),
                          name = "") 



#### something from "simulation studies.R" after null simulation

plot(density(sn$student[, "zint"]),
     xlim = c(-max(abs(sn$student[, "zint"])), max(abs(sn$student[, "zint"]))), 
     ylim = c(0, 0.5), lwd = 4, col = "green")
x <- seq(-5, 5, length=500)
hx <- dnorm(x)
lines(x, hx, type="l", lty=2, xlab="x value",
      ylab="Density", main="Comparison of t Distributions", lwd = 3)
d1 <- density(sn$garch[, "zint"])
lines(d1$x, d1$y, type="l", lty=2, xlab="x value",
      ylab="Density", main="Comparison of t Distributions", lwd = 3,
      col = "blue")
d2 <- density(sn$sv[, "zint"])
lines(d2$x, d2$y, type="l", lty=2, xlab="x value",
      ylab="Density", main="Comparison of t Distributions", lwd = 3,
      col = "red")

sm <- cbind(student = sn$student[, "zint"],
            garch = sn$garch[, "zint"],
            sv = sn$sv[, "zint"])
sm <- melt(sm)
sm <- sm[, -1]
names(sm)[1] <- "dgp"
sm <- as.data.frame(sm)

x <- seq(-5, 5, length = 1000) 
nm <- data.frame(x = x, y = dnorm(x))

v <- ggplot() + geom_line(aes(x = x, y = y), data = nm, size = 2.3) +
  geom_density(aes(x = value, group = dgp, color = dgp), data = sm, size = 1.2) +
  scale_x_continuous(limits = c(-5, 5))
v

y <- gen.null(n = 3000, method = "sv")
acf(y)


plot(density(sn$garch[, "rz1m"]),
     xlim = c(0, max(abs(sn$normal[, "rz1m"]))), 
     ylim = c(0, 0.2), lwd = 4, col = "green")
x <- seq(0, 20, length=500)
hx <- dchisq(x, df = 5)
lines(x, hx, type="l", lty=2, xlab="x value",
      ylab="Density", main="Comparison of t Distributions", lwd = 3)

y <- rnorm(500)
aq <- Auto.Q(y, lags = length(y)/2)

prp <- Sys.time()
wbaq <- wb.auto.q(y = y, B = 1000, prob = 0.95)
Sys.time() - prp

prp <- Sys.time()
avr <- AutoBoot.test(y, nboot = 1000, wild = "Normal", prob = c(0.025, 0.975))
Sys.time() - prp

prp <- Sys.time()
wbzint <- wb.zint(y, B = 1000, prob = c(0.025, 0.975), w = ww)
Sys.time() - prp

ww <- calc.wj(n = 500, fun = theta)
ww2 <- calc.wj(n = 500*2, fun = f.theta)


Rprof(filename = "Rprof.out", append = FALSE, interval = 0.01)
wbzint <- wb.zint(y, B = 1000, prob = c(0.025, 0.975), w = ww)
Rprof(NULL)
summaryRprof()

prp <- Sys.time()
zm <- z.matrix(1:2000)
Sys.time() - prp

prp <- Sys.time()
zm <- zm * zm
Sys.time() - prp


# debugonce(wb.z2int)
prp <- Sys.time()
wbzint2 <- wb.z2int(y = y, B = 1000, prob = c(0.95), w = ww2)
Sys.time() - prp

Rprof(filename = "Rprof.out", append = FALSE, interval = 0.01)
wbzint2 <- wb.z2int(y = y, B = 100, prob = c(0.95), w = ww2)
Rprof(NULL)
summaryRprof()


prp <- Sys.time()
for(i in 1:100){
  zm <- z.matrix2(1:500)
}
Sys.time() - prp

prp <- Sys.time()
for(i in 1:100){
  zm <- z.matrix2.old(1:500)
}
Sys.time() - prp
