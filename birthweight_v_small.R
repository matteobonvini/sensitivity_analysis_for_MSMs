rm(list = ls())
set.seed(2022)

library(quantreg)
library(qgam)
library(ggplot2)
library(RColorBrewer)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./ctseff.R")
dat <- read.csv("./data/birthweight_data_sensitivity_msm.csv")

y <- dat$dbirwt
a <- dat$T
table(a)
a[a >= 3] <- 3
mean(y[a == 1])
mean(y[a == 2])
mean(y[a == 3])
xnames <- c("dmar", "mwhite", "mblack", "mhispan", "fwhite", "fblack", "fhispan",
            "foreignb", "alcohol", "disllbu", "tripre1", "tripre2", "tripre3",
            "tripre0", "adequac1", "adequac2", "adequac3", "nprevisu", "ddeadkids",
            "dbirmon_2", "dbirmon_3", "dbirmon_4", "dbirmon_5", "dbirmon_6",
            "dbirmon_7", "dbirmon_8", "dbirmon_9", "dbirmon_10", "dbirmon_11",
            "dbirmon_12", "dcntyfipb_3", "dcntyfipb_11", "dcntyfipb_17",
            "dcntyfipb_29", "dcntyfipb_45", "dcntyfipb_49", "dcntyfipb_71",
            "dcntyfipb_77", "dcntyfipb_79", "dcntyfipb_91", "dcntyfipb_101",
            "dcntyfipb_129", "dcntyfipb_133", "dmage", "dfage", "dmeduc", "dfeduc",
            "nprevist", "disllb", "dlivord")
x <- dat[, colnames(dat) %in% xnames]
ncol(x)
x$dfage_zero <- 1 * I(x$dfage == 0)
nunique.x <- apply(x, 2, function(u) length(unique(u)))
cont.x <- colnames(x)[which(nunique.x > 2)]
length(cont.x)
disc.x <- colnames(x)[which(nunique.x == 2)]
length(disc.x)
a0 <- 0:3
table(a)

gam.form <- paste0(paste0("y ~ ", paste0(disc.x, collapse = " + ")),
                   " + dmeduc + dfeduc + dmage + dfage",
                   " + s(nprevist) + s(disllb) + dlivord")
##############################
## Nuisance functions specs ##
##############################

q <- function(y, a, x, new.x, new.a, tau) {
  
  a <- factor(a, ordered = TRUE)
  new.a <- factor(new.a, ordered = TRUE)
  fit <- qgam::mqgam(as.formula(gam.form), data = cbind(y, a, x), qu = tau)
  new.dat <- cbind(a = new.a, new.x)
  out <- do.call("cbind", lapply(sort(tau), 
                                 function(u) qdo(fit, u, predict, 
                                                 newdata = new.dat)))
  return(out)
}
mu <- function(y, a, x) {
  fit <- mgcv::gam(as.formula(gam.form), data = cbind(y = y, a = a, x))
  return(fit)
}
v <- mu

drl.lin <- function(y, a, new.a) {
  a.num <- as.numeric(as.character(a))
  new.a.num <- as.numeric(as.character(new.a))
  dat <- data.frame(y = y, one = 1, a1 = a.num)
  fit <- lm(y ~ -1 + ., data = dat)
  new.a.mat <- cbind(one = 1, a1 = new.a.num)
  out.est <- predict.lm(fit, newdata = as.data.frame(new.a.mat))
  out <- cbind(est = out.est, ll = NA, up = NA)
  return(out)
}
drl.quad <- function(y, a, new.a) {
  a.num <- as.numeric(as.character(a))
  new.a.num <- as.numeric(as.character(new.a))
  dat <- data.frame(y = y, one = 1, a1 = a.num, a2 = a.num^2)
  fit <- lm(y ~ -1 + ., data = dat)
  new.a.mat <- cbind(one = 1, a1 = new.a.num, a2 = new.a.num^2)
  out.est <- predict.lm(fit, newdata = as.data.frame(new.a.mat))
  out <- cbind(est = out.est, ll = NA, up = NA)
  return(out)
}
drl.np <- function(y, a, new.a) {
  dat <- data.frame(y = y, a = factor(a, ordered = TRUE))
  fit <- lm(y ~ -1 + a, data = dat)
  out.est <- predict.lm(fit, 
                        newdata = data.frame(a = factor(new.a, ordered = TRUE)))
  out <- cbind(est = out.est, ll = NA, up = NA)
  return(out)
}
gamma.vals <- seq(1, 1.25, length.out = 5)
# do the HULC
n <- length(y)
s <- sample(rep(1:6, ceiling(n / 6))[1:n])
ests.l.lin <- ests.u.lin <- ests.l.quad <- ests.u.quad <- 
  ests.l.np <- ests.u.np <- array(NA, dim = c(length(a0), length(gamma.vals), 6))
for(i in 1:6) {
  res <- ctseff(a0 = a0, learner = c("dr", "dr", "dr"), gamma = gamma.vals, 
                y = y[s == i], a = a[s == i], x = x[s == i, , drop = FALSE],
                drl = list("dr" = drl.lin, "dr" = drl.quad, "dr" = drl.np), 
                q = q, ps = NULL, mu = mu, v = v, nsplits = 1, 
                double_split = FALSE)
  ests.l.lin[, , i] <- res$est.l[[1]][, 1, ]
  ests.u.lin[, , i] <- res$est.u[[1]][, 1, ]
  ests.l.quad[, , i] <- res$est.l[[2]][, 1, ]
  ests.u.quad[, , i] <- res$est.u[[2]][, 1, ]
  ests.l.np[, , i] <- res$est.l[[3]][, 1, ]
  ests.u.np[, , i] <- res$est.u[[3]][, 1, ]
  print(i)
}

# choose whetner linear, quadratic or nonparam should be plotted
model_type <- "np"
if(model_type == "np") {
  ests.l <- ests.l.np
  ests.u <- ests.u.np
}
if(model_type == "lin") {
  ests.l <- ests.l.lin
  ests.u <- ests.u.lin
}
if(model_type == "quad") {
  ests.l <- ests.l.quad
  ests.u <- ests.u.quad
}
pt.est.l <- apply(ests.l, c(1, 2), mean)
pt.est.u <- apply(ests.u, c(1, 2), mean)
lower <- apply(ests.l, c(1, 2), min)
upper <- apply(ests.u, c(1, 2), max)
colnames(upper) <- colnames(lower) <- 
  colnames(pt.est.u) <- colnames(pt.est.l) <- gamma.vals
lower.ll <- reshape2::melt(lower)
upper.uu <- reshape2::melt(upper)
lower.est <- reshape2::melt(pt.est.l)
upper.est <- reshape2::melt(pt.est.u)
################
# Plot results #
################
gammastar <- gamma.vals[which(apply(lower, 2, max) < apply(upper, 2, min))][1]
plot.est <- data.frame(a = upper.est$Var1, gamma = as.factor(upper.est$Var2), 
                       ub = upper.est$value, lb = lower.est$value)
plot.ci <- data.frame(a = upper.uu$Var1, gamma = as.factor(upper.uu$Var2), 
                      ub = upper.uu$value, lb = lower.ll$value)

colors.vals <- c("black", brewer.pal(7, "Blues"))
names(colors.vals) <- levels(plot.est$gamma)
colors <- c("black", colors.vals[5], "red", colors.vals[7:8])
names(colors) <- levels(plot.est$gamma)
png(paste0("./results/smoking_", model_type, ".png"), 
    units = "cm", height = 10, width = 10, res = 313)
p <- ggplot() +
  labs(x = "Cigarettes / day", y = "Expected birth weight") +
  geom_point(data = plot.est[plot.est$gamma == 1, ], aes(x = a, y = lb),
             col = colors.vals[1]) +
  geom_line(data = plot.est[plot.est$gamma == 1, ],
            aes(x = a, y = ub), col = colors.vals[1]) +
  geom_line(data = plot.ci,
              aes(x = a, y = ub, group = gamma, col = gamma),
              alpha = 1) +
  geom_line(data = plot.ci,
            aes(x = a, y = lb, group = gamma, col = gamma),
            alpha = 1) +
  geom_line(data = plot.ci,
            aes(x = a, y = min(ub[plot.est$gamma == gammastar])),
            col = "pink", alpha = 1, linetype = "dashed") +
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks = 1:length(a0),
                     labels = c("0", "1-5", "6-10", "10+")) +
  theme_bw()
p
dev.off()
