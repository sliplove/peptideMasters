setwd("/home/sliplove/Documents/0Masters/gitcop/")
library(Rcpp)
library(coda)
library(lattice)
library(gtools)
sourceCpp("scorecpp/scoreR.cpp")
set.seed(42)

df <- read.csv("../data/Selected.csv", head = TRUE)[1:21,]
id <- 3
idx <- which(df$ID == id)

mat <- as.matrix(read.table(paste0("../data/mgf_input/matrix/matrix", df$LENGTH[idx], ".txt")))
rule <- as.matrix(read.table(paste0("../data/mgf_input/rule/rule", df$LENGTH[idx], ".txt")))
exp_spectrum <- as.numeric(read.table(paste0("../data/mgf_input/spectrum/spectrum", id, ".txt")))

MASS_PROTON = 1.00728
N_MASS = ncol(mat)
TOTAL_MASS = df$MASS[idx]
MAX_SCORE = df$SCORE.MDDPR.[idx]

exp_spectrum <- sort(exp_spectrum)

modify.mass <- function(mass)
  update_mass(mass, rule, FALSE)

get.score <- function(mass)
  score_peak(exp_spectrum, mat %*% mass, TOTAL_MASS, MASS_PROTON, FALSE)

source("wl.R")
source("mh.R")
source("se.R")
#----------------------------------------------
weights <- wl(0, MAX_SCORE)

s.min <- 0

start.mass <- as.numeric(rdirichlet(1, rep(1, N_MASS)))*TOTAL_MASS
start.score <- get.score(start.mass)
start.score 
start.mass

hit.n.run <- function(weights, start.mass, start.score,
                      step = 50000, min.n = 50000, eps = 0.02,
                      level = 0.95, tracelevel = 1, sss = 0) {
  z <- qnorm((1 + level) / 2)
  mh.res <- mh.weighted(min.n,
                        start.mass = start.mass, start.score = start.score,
                        s.min = sss, w = weights, trace = tracelevel > 1)
  one.traj <- mh.res$traj

  sigma <- list()
  w <- 0
  repeat {
    prob.const <- sum(exp(-weights[one.traj - s.min + 1])) / length(one.traj)
    g <- function(x) (x >= MAX_SCORE) * exp(-weights[x - s.min + 1]) / prob.const
    sigma <- se.obm(one.traj, g)
    w <- 2*z*sigma$se.mean
# 
#     if (tracelevel > 0)
#       cat(sprintf("prob: %e\n", prob.const))

      cat(sprintf("length %d, mean: %e, sdev: %e, w: %e, eps: %g, lambda: %e, prob: %e\n",
                  length(one.traj), sigma$mu, sigma$se.mean,
                  w, w/sqrt(sigma$lambda), sigma$lambda, prob.const))

    if (w / sqrt(sigma$lambda) < eps)
      break

    # FIXME: Inefficient
    mh.res <- mh.weighted(min.n,
                          start.mass = mh.res$mass, start.score = mh.res$score,
                          s.min = sss, w = weights, trace = tracelevel > 1)
    one.traj <- c(one.traj, mh.res$traj)
  }

  list(traj = one.traj, mu = sigma$mu, se = sigma$se.mean, lower = sigma$mu - w/2, upper = sigma$mu + w/2)
}


res.est.unif <- hit.n.run(weights, start.mass = start.mass, start.score = start.score,
                          min.n = 50000, sss = 0)
length(res.est.unif$traj)
tr <- res.est.unif$traj
res.est.unif$mu

table(tr)
plot(hist(tr))
plot(weights, type = 'l')


#------------------Standard MC----------------
pval.est <- function(N, trace = TRUE) {
  res <- rdirichlet(N, rep(1, N_MASS))
  v <- apply(res, 1, function(x) get.score(x*TOTAL_MASS))
  v
}

N <- length(res.est.unif$traj) 
v <- pval.est(N*10)
                                                                                                                    table(v)
table(v)

est <- length(v[v >= MAX_SCORE])/N
est + 1.96*sqrt(est*(1 - est)/N)
est - 1.96*sqrt(est*(1 - est)/N)
