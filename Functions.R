DLtest <- function (y, p = 1) {
  ym <- as.matrix(y - mean(y))
  n <- nrow(ym)
  s2 <- sum(ym^2)/(n - p)
  sum3 <- numeric(n - p)
  sum2 <- 0
  for (j in (p + 1):n) {
    sum1 <- 0
    for (i in (p + 1):n) {
      indicate <- 0
      zi <- ym[(i - 1):(i - p), 1]
      zj <- ym[(j - 1):(j - p), 1]
      tem1 <- as.numeric(zi <= zj)
      if (prod(tem1) == 1) 
        indicate <- 1
      sum1 <- sum1 + ym[i, 1] * indicate
    }
    sum2 <- sum2 + sum1^2
    sum3[j - p] <- abs(sum1/sqrt(n - p))
  }
  Cp <- sum2/(s2 * (n - p)^2)
  Kp <- max(sum3)/sqrt(s2)
  return(list(Cpstat = Cp, Kpstat = Kp))
}

Auto.DL.test <- function (y, lags = 10) {
  cpstat <- numeric(lags)
  kpstat <- numeric(lags)
  for (i in 1:lags) {
   dl <- DLtest(y, p = i)
   cpstat[i] <- dl$Cpstat
   kpstat[i] <- dl$Kpstat
  }
  aux <- matrix(1:lags)
  pin <- 2
#   q <- 2.4
#   cmax <- sqrt(T) * max(sqrt(ac3))
#   
#   if (maxro <= sqrt(q * log(T))) 
#     pin <- log(T)
  Lp1 <- cpstat - aux * pin
  Lp2 <- kpstat - aux * pin
  phat1 <- which.max(Lp1)
  phat2 <- which.max(Lp2)
  Tn1 <- cpstat[phat1]
  Tn2 <- kpstat[phat2]
#   pvalue <- 1 - pchisq(Tn, 1)
  return(list(Stat = Tn, Pvalue = pvalue))
}

generate.y <- function(t = 500, alpha = 0.7, p = 2){
  y <- c(rnorm(n = p, sd = 5), numeric(t))
  e <- rnorm(n = t)
  for (i in (p+1):(t+p)){
    y[i] <- alpha * y[i-p] + e[i-p]
  }
  y <- y[-(1:p)]
  return(y)
}

# original z
z <- function(theta, u = y){
  t <- length(u)
  sj <- t(u[-t]) %*% z.matrix(u = u[-1])
  th <- c(1, cumprod(rep(theta, t-2)))
  zt <- sj %*% th
  zt <- (sqrt(1-theta^2)/(var(u) * sqrt(t))) * zt
  return(c(zt))
}

z.vec <- function(theta, u = y){
  n <- length(theta)
  t <- length(u)
  sj <- t(u[-t]) %*% z.matrix(u = u[-1])
  th <- matrix(0, nrow = t-1, ncol = n)
  th[1, ] <- 1
  for (i in 2:nrow(th)){
    th[i, ] <- th[i-1, ] * theta
  }
  zt <- colSums(th * c(sj), na.rm = TRUE)
  zt <- (sqrt(1-theta^2)/(var(u) * sqrt(t))) * zt
  return(zt)
}


zh <- function(theta, u = y){
  t <- length(u)
  zm <- z.matrix(u = u[-1])
  sj <- t(u[-t]) %*% zm
  tau <- t(u[-t])^2 %*% zm^2
#   sj <- sj / (t-1):1
#   tau <- tau / (t-1):1
  sj <- sj / t
  tau <- tau / t
  sj <- sj[-length(sj)]
  tau <- tau[-length(tau)]
  rj <- sj / sqrt(tau) 
  th <- c(1, cumprod(rep(theta, t-2)))
  th <- th[-length(th)]
  zt <- rj %*% th
  zt <- (sqrt(1-theta^2) * sqrt(t)) * zt
  return(zt)
}

zh.vec <- function(theta, u = y){
  n <- length(theta)
  t <- length(u)
  zm <- z.matrix(u = u[-1])
  sj <- t(u[-t]) %*% zm
  tau <- t(u[-t])^2 %*% zm^2
  sj <- sj / t
  tau <- tau / t
  sj <- sj[-length(sj)]
  tau <- tau[-length(tau)]
  rj <- sj / sqrt(tau) 
  th <- matrix(0, nrow = t-1, ncol = n)
#   for (i in 1:n){
#     th[, i] <- c(1, cumprod(rep(theta[i], t-2)))
#   }
  th[1, ] <- 1
  for (i in 2:nrow(th)){
    th[i, ] <- th[i-1, ] * theta
  }
  th <- th[-nrow(th), ]
  zt <- colSums(th * rj, na.rm = TRUE)
  zt <- (sqrt(1-theta^2) * sqrt(t)) * zt
  return(zt)
}

# z <- function(theta, u = y){
#   t <- length(u)
#   sj <- t(u[-t]) %*% z.matrix(u = u[-1])
#   sj <- sj[-length(sj)]
#   sj <- sj/(t-2):1
#   th <- c(1, cumprod(rep(theta, t-2)))
#   th <- th[-length(th)]
#   zt <- sj %*% th
#   zt <- (sqrt(1-theta^2) * sqrt(t)/(var(u))) * zt
#   return(zt)
# }

# z.matrix.old <- function(u){
#   m <- u
#   for(i in 2:length(u)){
#     m <- c(m, u[-c(1:(i-1))], rep(0, i-1))
#   }
#   m <- matrix(m, ncol = length(u), nrow = length(u))
#   return(m)
# }

z.matrix <- function(u){
  n <- length(u)
  m <- matrix(0, ncol = n, nrow = n)
  for(i in 1:n){
    s <- seq(i, by = n-1, length.out = i)
    m[s] <- u[i]
  }
  return(m)
}

plot.z <- function(u, n = 100){
  theta <- seq(from = -0.99, to = 0.99, length = n)
  zt <-  Vectorize(z, vectorize.args = "theta")(theta = theta, u = u)
  dt <- data.frame(theta = theta, z = zt, q1 = qnorm(0.025), q2 = qnorm(0.975))
  p <- ggplot(aes(x = theta, y = z), data = dt)
  p <- p + geom_line(color = "blue", size = 1) 
  p <- p +geom_line(aes(y = q1), color = "green", size = 0.2)
  p <- p +geom_line(aes(y = q2), color = "green", size = 0.2)
  print(p)
}


theta <- function(th, j){
  (1-th^2)^(0.5)*th^j
}

f.theta <- function(th, j){
  (1-th^2)*th^j
}

wj <- function(j, interval = c(0,1), f.t = theta){
  integrate(f.t, lower = interval[1], upper = interval[2], j = j,
            subdivisions = 1000)$value
}

calc.wj <- function(n, fun){
  w <- sapply(0:(n), wj, interval = c(0, 0.9999), f.t = fun)
  return(w)
}


zh.int <- function(u = y, interval = c(0, 0.9999), w = NULL){
  t <- length(u)
  zm <- z.matrix(u = u[-1])
  sj <- t(u[-t]) %*% zm
  tau <- t(u[-t] * u[-t]) %*% (zm * zm)
#   sj <- sj / (t-1):1
#   tau <- tau / (t-1):1
  sj <- sj / t
  tau <- tau / t
  sj <- sj[-length(sj)]
  tau <- tau[-length(tau)]
  rj <- sj / sqrt(tau) 
  if (is.null(w)){
    w <- sapply(0:(t-3), wj, interval = interval)
  } else{
    w <- w[1:(t-2)]
  }
  zint <- rj %*% w
  zint <- sqrt(t) * zint
  sigma.zint <- sum(w * w)
  std <- zint/sqrt(sigma.zint)
  pval <- 2 * pnorm(abs(std), lower.tail=FALSE)
  return(list(stat = c(zint), stat.std = std, sigma = sigma.zint, pval = pval))
}

z.int <- function(u = y, interval = c(0, 0.9999), w = NULL){
  t <- length(u)
  sj <- t(u[-t]) %*% z.matrix(u = u[-1])
  if (is.null(w)){
    w <- sapply(0:(t-2), wj, interval = interval)
  } else{
    w <- w[1:(t-1)]
  }
  zint <- sj %*% w
  zint <- (1/(var(u)*sqrt(t)))*zint
  
  sigma.zint <- sum(w^2)
  std <- zint/sqrt(sigma.zint)
  pval <- 2 * pnorm(abs(std), lower.tail=FALSE)
  return(list(stat = c(zint), stat.std = std, sigma = sigma.zint, pval = pval))
}

monte.zint <- function(t = 200, alpha = 0.5, p = 1, n = 500, B = 500,
                       seed = 777){
  set.seed(seed)
  r <- NULL
  for (i in 1:n){
    y = generate.y(t = t, alpha = alpha, p = p)
    z <- wild.boot(u = y, B = B, test = z.int)
    zint <- c(stat = z$stat, pvalue = z$pvalue)
    r <- rbind(r, zint)
  }
  return(r)
}

monte.zint.par <- function(t = 200, alpha = 0.5, p = 1, n = 500, B = 500,
                       seed = 777, export){
  set.seed(seed)
  cl <- makeCluster(4)
  registerDoParallel(cl)
  r <- foreach (i = 1:n, .combine = rbind, .export = names(export),
                .packages = "foreach") %dopar% {
    y = generate.y(t = t, alpha = alpha, p = p)
    z <- wild.boot(u = y, B = B, test = z.int)
    zint <- c(stat = z$stat, pvalue = z$pvalue)
  }
  stopCluster(cl)
  return(r)
}


wild.boot <- function(u, B = 500, test, ...){
  t <- length(u)
  boot <- foreach(1:B, .combine = c) %do% test(rnorm(t) * u, ...)
  stat <- test(u, ...)
  pvalue = sum(boot > stat)/B
  return(list(boot.dist = boot, pvalue = pvalue, stat = stat))
}

wild.boot.par <- function(u, B = 500, test, export,
                          make.cluster = TRUE, ...){
  t <- length(u)
  if (make.cluster){
    cl <- makeCluster(4)
    registerDoParallel(cl)
    boot <- foreach(1:B, .combine = c, .export = names(export)) %dopar% test(rnorm(t) * u, ...)
    stopCluster(cl)
  }
  else{
    boot <- foreach(1:B, .combine = c, .export = names(export)) %dopar% test(rnorm(t) * u, ...)
  }
  stat <- test(u, ...)
  pvalue = sum(boot > stat)/B
  return(list(boot.dist = boot, pvalue = pvalue, stat = stat))
}


dl.mod <- function(y, p = 1, projection = TRUE, result = c("cp", "kp")){
  ym <- as.matrix(y - mean(y))
  n <- nrow(ym)
  s2 <- sum(ym^2)/(n - p)
  sum3 <- numeric(n - p)
  sum2 <- 0
  for (j in (p + 1):n) {
    if (projection){
      sum1 <- 0
      zj <- ym[(j - 1):(j - p), 1]
      for (i in (p + 1):n) {
        indicate <- 0
        zi <- ym[(i - 1):(i - p), 1]
        pr <- projection.dl(rbind(zi, zj))
        tem1 <- as.numeric(pr[1] <= pr[2])
        if (prod(tem1) == 1) 
          indicate <- 1
        sum1 <- sum1 + ym[i, 1] * indicate
      }
    }
    else{
          y1 <- matrix(0, ncol = p, nrow = n-p)
          for (k in 1:ncol(y1)){
            y1[, k] <- ym[(p-k+1):(n-k)] # why ym, and not just y?????
          }
          zj <- ym[(j - 1):(j - p), 1] # why ym, and not just y?????
          y2 <- matrix(rep(zj, n-p), byrow = TRUE, ncol = p, nrow = n-p)
          indicate <- as.numeric(apply(y1 <= y2, 1, all))
          sum1 <- sum(ym[(p+1):n] * indicate)
    }
    sum2 <- sum2 + sum1^2
    sum3[j - p] <- abs(sum1/sqrt(n - p))
  }
  Cp <- sum2/(s2 * (n - p)^2)
  Kp <- max(sum3)/sqrt(s2)
  if (all(result == "cp")){
    res <- c(Cpstat = Cp)
  }
  else if (all(result == "kp")){
    res <- c(Kpstat = Kp)
  }
  else{
    res <- list(Cpstat = Cp, Kpstat = Kp)
  }
  return(res)
}

projection.dl <- function(x){
  if (nrow(x) != 2){
    stop("NUmber of rows in matrix x is not equal to 2")
  }
  B <- x[1, ]  %*% t(x[2, ])
  e1 <- as.numeric(eigen(B)$vectors[, 1])
  result <- c(x %*% e1)
  return(result)
} 

auto.dl.mod <- function (y, lags = 10, ...) {
  cpstat <- numeric(lags)
  kpstat <- numeric(lags)
  for (i in 1:lags) {
    dl <- dl.mod(y, p = i, ...)
    cpstat[i] <- dl$Cpstat
    kpstat[i] <- dl$Kpstat
  }
#   aux <- matrix(1:lags)
#   pin <- 2
#   #   q <- 2.4
#   #   cmax <- sqrt(T) * max(sqrt(ac3))
#   #   
#   #   if (maxro <= sqrt(q * log(T))) 
#   #     pin <- log(T)
#   Lp1 <- cpstat - aux * pin
#   Lp2 <- kpstat - aux * pin
#   phat1 <- which.max(Lp1)
#   phat2 <- which.max(Lp2)
#   Tn1 <- cpstat[phat1]
#   Tn2 <- kpstat[phat2]
#   #   pvalue <- 1 - pchisq(Tn, 1)
#   return(list(Stat = Tn, Pvalue = pvalue))
  return(list( cpstat =  cpstat,  kpstat =  kpstat))
}

plot.auto.dl.mod <- function(x){
  # needs legend and x axis to be fixed
  d <- data.frame(cp = x$cpstat, kp = x$kpstat, p = 1:length(x$kpstat))
  g <- ggplot(aes(x = p), data = d)
  g <- g + geom_line(aes(y = cp), color = "blue", size = 1) 
  g <- g + geom_line(aes(y = kp), color = "red", size = 1)
  print(g)
}

kernel1 <- function(x){
  y <- (25/(12*(pi^2)*(x^2)))*((sin((6*pi*x)/5))/((6*pi*x)/5) - cos((6*pi*x)/5))
  return(y)
}

sigma.zint.f <- function(x, t){
  ((1-x^2)^(1/2)) * ((1-x^(t-1))/(1-x))
}

sigma.zint <- function(t, interval = c(0,1)){
  integrate(sigma.zint.f, lower = interval[1], upper = interval[2], t = t,
            subdivisions = 100L)$value
}

# 
# phi.j <- function(p, j){
#   fl <- list()
#   for(k in 1:(p-1)){
#     fl[[k]] <- function(x){0}
#   }
#   fl[[p]] <- function(x){1}
#   
#   for (i in  (p+1):(j+p-1)){
#     
#     fl[[i]] <- function(x){
#       l1 <- length(x)
#       ph <- 0
#       for (h in 1:l1){
#         ph <- ph + x[h] * fl[[i - h]](x)
#       }
#       return(ph)
#     }
#     
#   }
#   return(fl)
# }
# 
# phi.rec <- function(x, fx){
#   l1 <- length(x)
#   l2 <- length(fx)
#   if (l1 != l2){
#     stop("Unequal vector lenghts")
#   } else {
#     ph <- 0
#     for (i in 1:l1){
#       ph <- ph + x[i] * fx[[i]](x)
#     }
#   }
#   return(ph)
# }
# 
# phi.i <- function(theta, phi){
#   p1 <- length(theta)
#   p2 <- length(phi)
#   if (p1 != p2){
#     stop("Unequal vector lenghts")
#   } else {
#     phi <- sum(theta*phi)
#   }
#   return(phi)
# }
# 
# 
# 
# f.list <- function(j){
#   fl <- list()
#   for (i in 1:j){
#     fl[[i]] <- phi.j(i-1)
#   }
#   fl <- fl[length(fl):1]
#   return(fl)
# }


phi.j <- function(j){
  if (j < 1){
    phi <- function(x){0}
  } else if (j == 1){
    phi <- function(x){
      #sqrt(1 - sum(x^2))
      1
    }
  } else {
    phi <- function(x){
      l1 <- length(x)
      ph <- 0
      for (h in 1:l1){
        phj <- phi.j(j-h)
        ph <- ph + x[h] * phj(x)
      }
      #ph <- ph * sqrt(1 - sum(x^2))
      return(ph)
    }
  }
  return(phi)
}

z.multi <- function(theta, u, rtheta = FALSE){
  if (rtheta){
    if (!missing(theta)){
      warning("Argument theta ignored.")
    }
    theta <- c(runif(2, -1, 1))
    roots <- polyroot(c(1,theta))
    while (any(Mod(roots) <= 1)){
      theta <- c(runif(2, -1, 1))
      roots <- polyroot(c(1,theta))
    }
  }
  roots <- polyroot(c(1,theta))
  if (any(Mod(roots) <= 1)){
    stop("Process is not stacionary.")
  }
  t <- length(u)
  d <- length(theta)
  sj <- t(u[-t]) %*% z.matrix(u = u[-1])
  phi <- numeric(t-1 + d)
  phi[1:(d-1)] <- 0
  phi[d] <- 1
  for(i in (d+1):(length(phi))){
    phi[i] <- theta %*% phi[(i-1):(i-d)]
  }
  phi <- phi[-c(1:(d-1))]
  sums <- numeric(d)
  for (i in 1:d){
    p <- phi[i:(t-2+i)]
    sums[i] <- sj %*% p
  }
  z <- sqrt(1-sum(theta^2))/(var(u) * sqrt(t)) * sums
  chi <- z %*% z
  chi.p <- 1 - pchisq(chi, 1)
  result <- list(z = z, chi = chi, chi.pvalue = chi.p)
  return(result)
}

gen.null <- function(n, method){
  if (method == "normal"){
    y = rnorm(n = n)
  }
  if (method == "student"){
    y = rt(n = n, df = 5)
  }
  if (method == "garch"){
    y <- numeric(n)
    s <- numeric(n)
    s[1] <- 0.01 + 0.8 * 1
    y[1] <- rnorm(1) * sqrt(0.5 * s[1])
    for (i in 2:n){
      s[i] <- 0.01 + 0.8 * s[i-1] + 0.1 * y[i-1]^2
      y[i] <- rnorm(1) * sqrt(0.5 * s[i])
    }
  }
  if (method == "sv"){
    y <- numeric(n)
    s <- numeric(n)
    s[1] <- 0.9 * 1
    y[1] <- rnorm(1) * exp(s[1])
    for (i in 2:n){
      s[i] <- 0.9 * s[i-1] + 0.3 * rnorm(1)
      y[i] <- rnorm(1) * exp(s[i])
    }
  }
  return(y)
}

sim.null <- function(N, t, d,
                     z.function = zh,
                     z.int.function = zh.int,
                     z2int.function = z.sq.int2,
                     z.vec.function = zh.vec,
                     qzsq = NULL, w = NULL, w2 = NULL){
  smth <- null.tests(rnorm(30))
  l <- length(smth)
  s1 <- NULL
  s2 <- NULL
  s3 <- NULL
  s4 <- NULL
  prp <- Sys.time()
  for (i in 1:N){
    y1 <-  gen.null(n = t, method = "student")
    y2 <-  gen.null(n = t, method = "garch")
    y3 <-  gen.null(n = t, method = "sv")
    y4 <-  gen.null(n = t, method = "normal")
    
    s1 <- c(s1, null.tests(y = y1, z.function = z.function,
                           z.int.function = z.int.function,
                           z.vec.function = z.vec.function,
                           w = w, w2 = w2))
    s2 <- c(s2, null.tests(y = y2, z.function = z.function,
                           z.int.function = z.int.function,
                           z.vec.function = z.vec.function,
                           w = w, w2 = w2))
    s3 <- c(s3, null.tests(y = y3, z.function = z.function,
                           z.int.function = z.int.function,
                           z.vec.function = z.vec.function,
                           w = w, w2 = w2))
    s4 <- c(s4, null.tests(y = y4, z.function = z.function,
                           z.int.function = z.int.function,
                           z.vec.function = z.vec.function,
                           w = w, w2 = w2))
  }
  print(paste0("t: ", t))
  print(Sys.time() - prp)
  s1 <- matrix(s1, ncol = l, byrow = TRUE)
  s2 <- matrix(s2, ncol = l, byrow = TRUE)
  s3 <- matrix(s3, ncol = l, byrow = TRUE)
  s4 <- matrix(s4, ncol = l, byrow = TRUE)
  colnames(s1) <- names(smth)
  colnames(s2) <- names(smth)
  colnames(s3) <- names(smth)
  colnames(s4) <- names(smth)
  size <- matrix(NA, ncol = l, nrow = 4)
  size[1, ] <- t.size(m = s1, d = d, qzsq = qzsq)
  size[2, ] <- t.size(m = s2, d = d, qzsq = qzsq)
  size[3, ] <- t.size(m = s3, d = d, qzsq = qzsq)
  size[4, ] <- t.size(m = s4, d = d, qzsq = qzsq)
  rownames(size) <- c("student", "garch", "sv", "normal")
  colnames(size) <- names(smth)
  result <- list(student = s1, garch = s2, sv = s3, normal = s4, size = size)
  return(result)
}

null.tests <- function(y,
                       z.function = zh,
                       z.int.function = zh.int,
                       z2int.function = z.sq.int2,
                       z.vec.function = zh.vec,
                       w = NULL,
                       w2 = NULL){
  aq20 <- Auto.Q(y, lags = 20)
  aq5 <- Auto.Q(y, lags = 5)
  avr <- Auto.VR(y)
  zint <- z.int.function(u = y, w = w)
  z2int <- z2int.function(u = y, interval = c(0, 0.9999), dens = NULL,
                          w = w2)
  rz1 <- rz(u = y, m = 1, alpha = 0.05, z.function = z.function,
            z.vec.function = z.vec.function)
  rz1m2 <- rz(u = y, m = 2, alpha = 0.05, z.function = z.function,
              z.vec.function = z.vec.function)
  rz1m5 <- rz(u = y, m = 5, alpha = 0.05, z.function = z.function,
              z.vec.function = z.vec.function)
  rz1m10 <- rz(u = y, m = 10, alpha = 0.05, z.function = z.function,
               z.vec.function = z.vec.function)
  r <- c(aq5$Stat, aq20$Stat, avr, zint$stat.std, z2int$stat,
         rz1$stat, rz1m2$stat, rz1m5$stat, rz1m10$stat)
  names(r) <- c("aq5", "aq20", "avr", "zint", "z2int",
                "rz1", "rz1m2", "rz1m5", "rz1m10")
  return(r)
}

t.size <- function(m, d, qzsq = NULL){
  s <- numeric(nrow(d))
  names(s) <- colnames(m)
  for (i in 1:ncol(m)){
    nm <- colnames(m)[i]
    test <- m[, i]
    di <- as.character(d$distribution[d$test == nm])
    if (di == "N"){
      q1 <- qnorm(0.975)
      q2 <- qnorm(0.025)
      s[i] <- mean(test < q2) + mean(test >= q1)
    }
    if (di == "Chi") {
      q1 <- qchisq(0.95, df = 1)
      s[i] <- mean(test >= q1)
    }
    if (length(grep("Chi.", di, fixed = TRUE)) > 0) {
      df <- substring(di, 5, nchar(di))
      df <- as.numeric(df)
      q1 <- qchisq(0.95, df = df)
      s[i] <- mean(test >= q1)
    }
    if (di == "z.sq") {
      if (is.null(qzsq)){
        s[i] <- NA
      } else{
        s[i] <- mean(test >= qzsq)
      }
    }
  }
  return(s)
}



 
z.sq <- function(theta, u, z.function = z){
  zt <- z.function(theta = theta, u = u)^2
  return(zt)
}

# z.sq.int <- function(u, interval = c(0, 0.9999)){
#   require("cubature")
#   zint <- adaptIntegrate(z.sq, lowerLimit = interval[1],
#                          upperLimit = interval[2], u = u)
#   return(zint$integral)
# }

wb.z.sq.int <- function(u, B){
  stat <- z.sq.int(u = u)
  boot <- numeric(B)
  for (i in 1:B){
    u.star <- u * rnorm(length(u))
    boot[i] <- z.sq.int2(u = u.star)
  }
  pval <- mean(boot>=stat)
  r <- list(stat = stat, pval = pval)
  return(r)
}


rz <- function(u, m = 1, alpha = 0.05, theta = NULL, z.function,
               z.vec.function){
  if (m == 1){
    q1 <- qnorm(alpha/2)
    q2 <- qnorm(1-alpha/2)
    if (is.null(theta)){
      theta <- runif(1, min = -1, max = 1)
    }
    zs <- as.numeric(z.function(theta = theta, u = u))
    t <- (zs > q1) & (zs <= q2)
    pval <-  2*pnorm(abs(zs), lower.tail = FALSE)
    re <- list(result = t, stat = zs, theta = theta, pval = pval)
  }
  else{
    scov <- try(solve(z.cov(c(-5,-5))), silent = TRUE)
    while (class(scov) == "try-error"){
      th <- NULL
      th <- runif(m, min = -1, max = 1) 
      zcov <- z.cov(th)
      scov <- try(solve(zcov), silent = TRUE)
    }
    zs <- z.vec.function(theta = th, u = u)
    zs <- t(zs) %*% scov %*% zs
    re <- list(stat = zs, pval = 1 - pchisq(zs, df = m), theta = th)
  }
  return(re)
}

rz.boot <- function(u, m = 1, alpha = 0.05, B, z.function = zh){
  stat <- rz(u = u, m = m, alpha = alpha, z.function = z.function)
  ex <- c("rz", "zh", "z", "z.matrix", "z.cov")
  boot.stat <- foreach (i = 1:B, .combine = c, .packages = "foreach",
                        .export = ex) %dopar% {
                          u.star <- rnorm(length(u))
                          brz <- rz(u = u.star, m = m,
                                    alpha = alpha,
                                    theta = stat$theta,
                                    z.function = z.function)$stat
                        }
  
#   boot.stat <- numeric(B)
#   for (i in 1:B){
#     u.star <- rnorm(length(u))
#     boot.stat[i] <- rz(u = u.star, m = m, alpha = alpha,
#                        theta = stat$theta, z.function = z.function)$stat
#   }
  ### calculation of pvalue !!!!!!!!!!!!! must be checked
  pval <- mean(abs(boot.stat) < c(stat$stat))
  r <- list(stat = stat$stat, pval = pval)
  return(r)
}

rzsq <- function(u, m = 1, alpha = 0.05, theta = NULL, z.function = z){
  q <- qchisq(1-alpha/m, df = 1)
  if (is.null(theta)){
    theta <- runif(1)
  }
  zs <- as.numeric(z.sq(theta = theta, u = u, z.function = z.function))
  t <- zs <= q
  t <- all(t)
  pval <- 1 - pchisq(zs, df = m)
  r <- list(result = t, stat = zs, theta = theta, pval = pval)
  return(r)
}

rzsq.boot <- function(u, m = 1, alpha = 0.05, B){
  stat <- rzsq(u = u, m = m, alpha = alpha)
  boot.stat <- numeric(B)
  for (i in 1:B){
    e <- rnorm(length(u))
    u.star <- u * e
    boot.stat[i] <- rzsq(u = u.star, m = m, alpha = alpha,
                         theta = stat$theta)$stat
  }
  q <- quantile(boot.stat, 1 - alpha)
  ### calculation of pvalue !!!!!!!!!!!!! must be checked
  pval <- mean(stat$stat < boot.stat)
  t <- (stat$stat < q)
  names(t) <- NULL
  r <- list(stat = stat$stat, pval = pval, t = t)
  return(r)
}

z.matrix2.old <- function(u){
  n <- length(u)
  m <- matrix(u, ncol = n, nrow = n)
  zeros <- NULL
  for(i in 1:(n-1)){
    s <- seq(i, by = n-1, length.out = i)
    zeros <- c(zeros, s)
  }
  m[zeros] <- 0
  return(m)
}

z.matrix2 <- function(u){
  n <- length(u)
  m <- matrix(rep(u, each = n), nrow = n, byrow = TRUE)
  zm <- z.matrix(rep(1, n))
  s <- seq(n, by = n-1, length.out = n)
  zm[s] <- 0
  zm <- 1 * (zm == 0)
  m <- m * zm
  return(m)
}

z.sq.int <- function(u = y, interval = c(0, 0.9999), dens = NULL, w = NULL){
  t <- length(u)
  zm <- z.matrix(u = u[-1])
  sj <- t(u[-t]) %*% zm
  sj <- sj / t
  # rhoj <- sj/(t*var(u))
  rhoj <- sj/var(u) 
  
  zm1 <- z.m(rhoj)
  x <- 0:(length(rhoj)-1)
  zm <- z.matrix2(x[(length(x)-1):1])
  z1 <- z.matrix2(rep(1, length(x)-1))
  z1[z1==0] <- -1
  zz <- matrix(x[2:length(x)], nrow = nrow(zm), ncol = ncol(zm),
               byrow = TRUE)
  zz <- zz * z1
  dm <- as.numeric(zz + zm)
  zm2 <- dm[dm > 0] + 1
  if (is.null(w)){
    w <- sapply(0:((t-2)*2), wj, interval = interval, f.t = f.theta)
  } 
  w1 <- w[((0:(t-3))*2 +1)]
  w2 <- w[zm2]
  
  z2int <- t * (rhoj^2 %*% w1 + 2 * zm1 %*% w2)
  pval <- mean(dens >= c(z2int))
  result <- list(stat = c(z2int), pval = pval)
  return(result)
}

z.sq.int2 <- function(u = y, interval = c(0, 0.9999),
                      dens = NULL, w = NULL){
  t <- length(u)
  zm <- z.matrix(u = u[-1])
  sj <- t(u[-t]) %*% zm
  tau <- t(u[-t])^2 %*% zm^2
#   sj <- sj / (t-1):1
#   tau <- tau / (t-1):1
  sj <- sj / t
  tau <- tau / t
  # rhoj <- sj/(t*var(u))
  rhoj <- sj/sqrt(tau)
  rhoj <- rhoj[-length(rhoj)]

  zm1 <- z.m(rhoj)
  x <- 0:(length(rhoj)-1)
  zm <- z.matrix2(x[(length(x)-1):1])
  z1 <- z.matrix2(rep(1, length(x)-1))
  z1[z1==0] <- -1
  zz <- matrix(x[2:length(x)], nrow = nrow(zm), ncol = ncol(zm),
               byrow = TRUE)
  zz <- zz * z1
  dm <- as.numeric(zz + zm)
  zm2 <- dm[dm > 0] + 1
  if (is.null(w)){
    w <- sapply(0:((t-2)*2), wj, interval = interval, f.t = f.theta)
  } 
  w1 <- w[((0:(t-3))*2 +1)]
  w2 <- w[zm2]
  
  z2int <- t * (rhoj^2 %*% w1 + 2 * zm1 %*% w2)
  pval <- mean(dens >= c(z2int))
  result <- list(stat = c(z2int), pval = pval)
  return(result)
}


z.m <- function(x){
  zm <- z.matrix2(x[(length(x)-1):1])
  zz <- matrix(x[2:length(x)], nrow = nrow(zm), ncol = ncol(zm),
               byrow = TRUE)
  dm <- as.numeric(zz * zm)
  dm <- dm[dm != 0]
  return(dm)
}

z.cov <- function(x){
  l <- length(x)
  x1 <- matrix(x, nrow = l, ncol = l)
  x2 <- t(x1)
  co <- sqrt((1 - x1^2) * (1 - x2^2)) / (1 - x1 * x2)
  return(co)
}

z.cov.sim <- function(theta, N, t, dgp = "normal", z.function = z){
  zt <- matrix(NA, nrow = N, ncol = length(theta))
  for (j in 1:N){
    u <- gen.null(n = t, method = dgp)
    zt[j, ] <- z.function(theta = theta, u = u)
  }
  co <- cor(zt)
  return(co)
}

gen.alternative <- function(n, method){
  # set.seed(round(runif(1) * n^3))
  if (method == "ar1"){
    s <- numeric(n)
    e <- rnorm(n)
#     s[1] <- 1
#     for (i in 2:n){
#       s[i] <- 0.01 + 0.8 * s[i-1] + 0.1 * e[i-1]^2
#     }
#     z <- sqrt(s) * e
    
    y <- numeric(n+1)
    y[1] <- 0
    for(i in 2:(n+1)){
      #y[i] <- 0.4 * y[i-1] + z[i-1]
      y[i] <- 0.4 * y[i-1] + e[i-1]
    }
    y <- y[-1]
  }
  if (method == "ar12"){
    e <- rnorm(n)
    y <- numeric(n+12)
    y[1:12] <- 0
    for(i in 13:(n+12)){
      y[i] <- 0.3 * y[i-12] + e[i-12]
    }
    y <- y[-c(1:12)]
  }
  if (method == "ma2"){
    e <- rnorm(n+2)
    y <- numeric(n)
    for(i in 1:n){
      y[i] <- -0.3 * e[i+1] + 0.2 * e[i] + e[i+2]
    }
  }
  if (method == "ma12"){
    e <- rnorm(n+12)
    y <- numeric(n)
    for(i in 1:n){
      y[i] <- 0.5 * e[i] + e[i+12]
    }
  }
  if (method == "bl"){
    e <- rnorm(n+2)
    y <- numeric(n+2)
    y[1:2] <- 1
    for(i in 3:(n+2)){
      y[i] <- 0.2 * y[i-1] * e[i-1] + 0.1 * y[i-2] * e[i-2] + e[i]
    }
    y <- y[-c(1:2)]
  }
  return(y)
}

wb.auto.q <- function(y, B, prob = 0.95, lags = length(y)/2){
  set.seed(12345)
  y <- as.matrix(y)
  LC <- Auto.Q(y, lags = lags)$Stat
  statmat <- matrix(NA, nrow = B, ncol = 1)
  for (i in 1:B) {
    ys <- y * rnorm(nrow(y))
    statmat[i, ] <- Auto.Q(ys, lags = lags)$Stat
  }
  tem <- abs(statmat) > abs(LC)
  tem[tem == "TRUE"] <- 1
  p <- mean(tem)
  CI <- quantile(statmat, prob)
  return(list(test.stat = LC, pval = p, CI = CI))
}

wb.zint <- function(y, B, prob = c(0.025, 0.975), w = NULL){
  set.seed(12345)
  y <- as.matrix(y)
  LC <- zh.int(u = y, w = w)$stat
  statmat <- matrix(NA, nrow = B, ncol = 1)
  for (i in 1:B) {
    ys <- y * rnorm(nrow(y))
    statmat[i, ] <- zh.int(u = ys, w = w)$stat
  }
  tem <- abs(statmat) > abs(LC)
  tem[tem == "TRUE"] <- 1
  p <- mean(tem)
  CI <- quantile(statmat, prob)
  return(list(test.stat = LC, pval = p, CI = CI))
}


wb.z2int <- function(y, B, prob = c(0.95), w = NULL){
  set.seed(12345)
  y <- as.matrix(y)
  LC <- z.sq.int2(u = y, w = w)$stat
  #statmat <- matrix(NA, nrow = B, ncol = 1)
  cl <- makeCluster(4)
  registerDoParallel(cl)
  ex <- c("z.sq.int2", "z.matrix", "z.matrix2", "z.m", "pvalue.zsqint")
  statmat <- foreach (i = 1:B, .combine = c, .packages = "foreach",
                .export = ex) %dopar% {
                  ys <- y * rnorm(nrow(y))
                  stat <- z.sq.int2(u = ys, w = w)$stat
                }
  stopCluster(cl)
#   for (i in 1:B) {
#     ys <- y * rnorm(nrow(y))
#     statmat[i, ] <- z.sq.int2(u = ys, w = w)$stat
#   }
  tem <- abs(statmat) > abs(LC)
  tem[tem == "TRUE"] <- 1
  p <- mean(tem)
  CI <- quantile(statmat, prob)
  return(list(test.stat = LC, pval = p, CI = CI))
}

sim.alternative <- function(N, t, B, h1, w, w2, dzsq, cv = NULL){
  at1 <- alternative.tests2(y = rnorm(t), B = 5, w = w, w2 = w2,
                            dzsq = dzsq, only.names = TRUE)
  p <- matrix(0, nrow = N, ncol = length(at1))
  colnames(p) <- names(at1)
  if (is.null(cv)){
    re <- "pvalue"
  } else{
    pv <- c("avr", "gs", "zint", "z2int", "rz1", "rz1m2", "rz1m5")
    act <- c("aq5", "aq20", "rz1m10")
    re <- "value"
    cv1 <- cv[cv$t == t, act]
  }
  # yy <- gen.alternative(n = t, method = h1)
  # gsb <- Gen.Spec.boot(y = yy, B = B)
  for(i in 1:N){
    yy <- gen.alternative(n = t, method = h1)
    p[i, ] <- alternative.tests2(y = yy, B = B, w = w, w2 = w2,
                                 dzsq = dzsq, gsb = NULL,
                                 result = re)
  }
  if (is.null(cv)){
    power <- apply(p, 2, function(x){mean(x < 0.05)})
  } else{
    p1 <- apply(p[, pv], 2, function(x){mean(x < 0.05)})
    for (i in act){
      p1[[i]] <- mean(p[, i] >= c(cv1[i]))
    }
    power <- p1
  }
  r <- list(all = p, power = power)
  return(r)
}

alternative.tests <- function(y, B){
  wbaq <- wb.auto.q(y = y, B = 1000, prob = 0.95)
  avr <- AutoBoot.test(y, nboot = B, wild = "Normal", prob = c(0.025, 0.975))
  # gs <- Gen.Spec(y = y, B = B)
  gs <- NULL
  zint <- z.int(u = y)
  rz1 <- rz.boot(u = y, m = 1, alpha = 0.05, B = B)
  z2int <- wb.z.sq.int(u = y, B = B)
  rz2 <- rzsq.boot(u = y, m = 1, alpha = 0.05, B = B)
  r <- c(aq$Pvalue, avr$pval, gs, zint$pval, rz1$pval, z2int$pval, rz2$pval)
  names(r) <- c("aq", "avr", "gs", "zint", "rz1", "z2int", "rz2")
  return(r)
}



alternative.tests2 <- function(y, B, w = NULL, w2 = NULL,
                               dzsq = NULL, gsb = NULL,
                               only.names = FALSE,
                               result = c("value")){
  if (!only.names){
    aq5 <- Auto.Q(y = y, lags = 5)
    aq20 <- Auto.Q(y = y, lags = 20)
    avr <- AutoBoot.test(y, nboot = B, wild = "Normal", prob = c(0.025, 0.975))
    if (length(y) <= 201){
      gs <- Gen.Spec(y = y, B = B, CvMexpb = gsb)
    }
    zint <- zh.int(u = y, w = w)
    z2int <- z.sq.int2(u = y, w = w2, dens = dzsq)
    rz1 <- rz(u = y, m = 1, alpha = 0.05, z.function = zh)
    rz1m2 <- rz(u = y, m = 2, alpha = 0.05, z.function = zh,
                z.vec.function = zh.vec)
    rz1m5 <- rz(u = y, m = 5, alpha = 0.05, z.function = zh,
                z.vec.function = zh.vec)
    rz1m10 <- rz(u = y, m = 10, alpha = 0.05, z.function = zh,
                 z.vec.function = zh.vec)
    if (result == "pvalue"){
      if (length(y) <= 201){
        r <- c(aq5$Pvalue, aq20$Pvalue, avr$pval, gs,
               zint$pval, z2int$pval, rz1$pval,
               rz1m2$pval, rz1m5$pval, rz1m10$pval)
      }else{
        r <- c(aq5$Pvalue, aq20$Pvalue, avr$pval,
               zint$pval, z2int$pval, rz1$pval,
               rz1m2$pval, rz1m5$pval, rz1m10$pval)
      }
    }
    else if (result == "value"){
      if (length(y) <= 201){
        r <- c(aq5$Stat, aq20$Stat, avr$pval, gs,
               zint$pval, z2int$pval, rz1$pval,
               rz1m2$pval, rz1m5$pval, rz1m10$stat)
      }else{
        r <- c(aq5$Stat, aq20$Stat, avr$pval,
               zint$pval, z2int$pval, rz1$pval,
               rz1m2$pval, rz1m5$pval, rz1m10$stat)
      }
    }
  }
  else{
    if (length(y) <= 201){
      r <- numeric(10)
    }else{
      r <- numeric(9)
    }
  }
  if (length(y) <= 201){
    names(r) <- c("aq5", "aq20", "avr", "gs", "zint", "z2int",
                  "rz1", "rz1m2", "rz1m5", "rz1m10")
  }else{
    names(r) <- c("aq5", "aq20", "avr", "zint", "z2int",
                  "rz1", "rz1m2", "rz1m5", "rz1m10")
  }
  return(r)
}


compweexp <- function(inf) 
{
  n <- length(inf)
  weiexp <- matrix(1, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in (i + 1):n) {
      if (j > n) 
        break
      aux1 <- (inf[i] - inf[j]) %*% t(inf[i] - inf[j])
      weiexp[i, j] <- exp(-0.5 * aux1)
      weiexp[j, i] <- weiexp[i, j]
    }
  }
  return(weiexp)
}


Mammen <- function(n) 
{
  p <- (sqrt(5) + 1)/(2 * sqrt(5))
  zmat <- rep(1, n) * (-(sqrt(5) - 1)/2)
  u <- runif(n, 0, 1)
  zmat[u > p] <- (sqrt(5) + 1)/2
  return(zmat)
}




Gen.Spec <- function (y, B = 300, CvMexpb) 
{
  #set.seed(12345)
  n <- length(y)
  e <- y - mean(y)
  v <- var(e)
  y1 <- y[1:(n - 1)]
  weiexp <- compweexp(y1)
  CvMexp <- 0
  au <- 1:(n - 1)
  aux <- 1/((au * pi)^2)/(n - au + 1)
  for (j in 1:(n - 1)) {
    CvMexp <- CvMexp + aux[j] *
      t(e[(1 + j):n]) %*% weiexp[1:(n - j), 1:(n - j)] %*% e[(1 + j):n]  
  }
  CvMexp <- CvMexp/v

  

  
#   cl <- makeCluster(4)
#   registerDoParallel(cl)
  if (is.null(CvMexpb)){
    ex <- c("Mammen")
    CvMexpb <- foreach (k = 1:B, .combine = rbind, .packages = "foreach",
                        .export = ex) %dopar% {
                          eb <- e * Mammen(n)
                          eb <- eb - mean(eb)
                          tem <- 0
                          for (j in 1:(n - 1)) {
                            tem <- tem + aux[j] *
                              t(eb[(1 + j):n]) %*% weiexp[1:(n - j), 1:(n - j)] %*% eb[(1 + j):n]
                          }
                          r <- cbind(tem/v > CvMexp, tem/v)
                        }
  }
#   stopCluster(cl)
  
  
#   for (k in 1:B) {
#     eb <- e * Mammen(n)
#     eb <- eb - mean(eb)
#     tem <- 0
# 
#     for (j in 1:(n - 1)) {
#       tem <- tem + aux[j] *
#         t(eb[(1 + j):n]) %*% weiexp[1:(n - j), 1:(n - j)] %*% eb[(1 + j):n]
#     }
#     CvMexpb[k, ] <- cbind(tem/v > CvMexp, tem/v)
#   }
  pboot <- mean(CvMexpb[, 2] > c(CvMexp))
  Critboot <- quantile(CvMexpb[, 2], c(0.9, 0.95, 0.99))
  r <- list(stat = CvMexp, pval = pboot, crit = Critboot)
  return(pboot)
}

Gen.Spec.boot <- function(y, B){
  n <- length(y)
  e <- y - mean(y)
  v <- var(e)
  y1 <- y[1:(n - 1)]
  weiexp <- compweexp(y1)
  CvMexp <- 0
  au <- 1:(n - 1)
  aux <- 1/((au * pi)^2)/(n - au + 1)
  # CvMexpb <- matrix(0, nrow = B, ncol = 2)
  ex <- c("Mammen")
  CvMexpb <- foreach (k = 1:B, .combine = rbind, .packages = "foreach",
                      .export = ex) %dopar% {
                        eb <- e * Mammen(n)
                        eb <- eb - mean(eb)
                        tem <- 0
                        for (j in 1:(n - 1)) {
                          tem <- tem + aux[j] *
                            t(eb[(1 + j):n]) %*% weiexp[1:(n - j), 1:(n - j)] %*% eb[(1 + j):n]
                        }
                        r <- cbind(tem/v > CvMexp, tem/v)
                      }
  return(CvMexpb)
}


AutoBoot.test <- function (y, nboot, wild, prob = c(0.025, 0.975)) 
{
 # set.seed(12345)
  y <- as.matrix(y)
  LC <- Auto.VR(y)
  statmat <- matrix(NA, nrow = nboot, ncol = 1)
  if (wild == "Normal") {
    for (i in 1:nboot) {
      ys <- y * rnorm(nrow(y))
      statmat[i, ] <- Auto.VR(ys)
    }
  }
  if (wild == "Mammen") {
    for (i in 1:nboot) {
      ys <- y * Mammen(nrow(y))
      statmat[i, ] <- Auto.VR(ys)
    }
  }
  if (wild == "Rademacher") {
    for (i in 1:nboot) {
      ys <- y * Rademacher(nrow(y))
      statmat[i, ] <- Auto.VR(ys)
    }
  }
  tem <- abs(statmat) > abs(LC)
  tem[tem == "TRUE"] <- 1
  p <- mean(tem)
  CI <- quantile(statmat, prob)
  return(list(test.stat = LC, pval = p, CI = CI))
}

sim.a1 <- function(ww, ww2, dzsq, alt, N, B,
                   tt = c(25, 50, 100, 200),
                   cv = NULL){
  cl <- makeCluster(4)
  registerDoParallel(cl)
 
  po <- list()
  pow <- NULL
  time <- NULL
  prp <- Sys.time()
  for (i in 1:length(tt)){
    prpi <- Sys.time()
    sa <- sim.alternative(N = N, t = tt[i], B = B, h1 = alt, w =ww,
                          w2 = ww2, dzsq = dzsq, cv = cv)
    print(Sys.time() - prpi)
    time <- c(time, Sys.time() - prpi)
    pow <- rbind(pow, c(tt[i], sa$power))
    print(sa$power)
    po[[i]] <- sa$all
  }
  print(Sys.time() - prp)
  colnames(pow) <- c("t", names(sa$power))
  stopCluster(cl)
  return(list(details = po, power = pow))
}

sim.a2 <- function(N = 500, B = 1000,
                   alt = c("ar1", "ar12", "ma2", "ma12", "bl"),
                   tt = c(25, 50, 100, 200), cv = NULL){
  if (is.null(cv)){
    warning("Calculating power using asimptotic critical values.")
  }
  d <- load(file = "Rdata/dzsq.Rdata")
  dzsq <- get(d)
  ww <- calc.wj(n = 5000, fun = theta)
  ww2 <- calc.wj(n = 5000*2, fun = f.theta)
  det <- list()
  power <- NULL
  for (i in alt){
    s <- sim.a1(ww = ww, ww2 = ww2, dzsq = dzsq, alt = i, N = N, B = B,
                tt = tt, cv = cv)
    det[[i]] <- s$details
    save(det, file = "Rdata/power.detailed2.Rdata")
    power <- rbind(power, data.frame(dgp = i, s$power))
    save(power, file = "Rdata/power2.Rdata")
  }
  return(list(detailed = det, power = power))
}

get.empirical.critical.values <- function(object,
                                          dgp = c("normal", "garch"),
                                          alpha = 0.05){
  r <- NULL
  for (i in names(object)){
    for (j in dgp){
      q <- apply(l[[i]][[j]], 2, quantile, probs = 1 - alpha)
      q <- data.frame(t = as.numeric(i), dgp = j, t(q))
      r <- rbind(r, q)
    }
  }
  return(r)
}