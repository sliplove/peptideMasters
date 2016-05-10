setwd("/home/sliplove/Documents/Masters/github_copy/")
library(Rcpp)
library(coda)
library(pracma)
sourceCpp("scorecpp/scoreR.cpp")
set.seed(42)

MAX_SCORE = 15
TOTAL_MASS =  150
MASS_PROTON = 0
max.W = TOTAL_MASS
sz = 4
exp.spectrum <- seq(10, 150, 10)
rule <- as.matrix(read.table("../tables/rg_tst2"))

modify.mass <- function(mass)
  update_mass(mass, rule, TRUE)

modify.mass(c(10, 20, 40, 80))

mass.to.spectrum <- function(mass)
  sort(unique(c(mass[1], mass[2], mass[3],
                mass[1] + mass[2],
                mass[1] + mass[3],
                mass[1] + mass[4],
                mass[2] + mass[3],
                mass[2] + mass[4],
                mass[3] + mass[4],
                mass[1] + mass[2] + mass[3],
                mass[1] + mass[2] + mass[4],
                mass[1] + mass[3] + mass[4],
                mass[2] + mass[3] + mass[4],
                sum(mass))))

get.score <- function(mass)
  score_peak(exp.spectrum, mass.to.spectrum(mass), TOTAL_MASS, MASS_PROTON, FALSE)

source("wl.R")
source("mh.R")
source("se.R")

#test if everything is correct
weights <- wl(0, MAX_SCORE)

SCORE_ <- 14
s.min <- 0
hit.n.run <- function(weights, start.mass, start.score,
                      step = 50000, min.n = 50000, eps = 0.02,
                      level = 0.95, tracelevel = 1) {
  start.mass <-  diff(c(0, sort(randperm(TOTAL_MASS - 1, sz - 1)), TOTAL_MASS))
  start.score <- get.score(start.mass)
  step = 50000
  min.n = 50000
  eps = 0.02
  level = 0.95
  tracelevel = 1
  
  z <- qnorm((1 + level) / 2)
  mh.res <- mh.weighted(min.n,
                        start.mass = start.mass, start.score = start.score,
                        s.min = 0, w = weights, trace = tracelevel > 1)
  one.traj <- mh.res$traj
  
  sigma <- list()
  w <- 0
  prob.const <- sum(exp(-weights[one.traj - s.min + 1])) / length(one.traj)
  repeat {
    g <- function(x) (x >= SCORE_) * exp(-weights[x - s.min + 1]) / prob.const
    
    sigma <- se.obm(one.traj, g)
    w <- 2*z*sigma$se.mean
    
    if (tracelevel > 0)
      cat(sprintf("length %d, mean: %e, sdev: %e, w: %e, eps: %g, lambda: %e\n",
                  length(one.traj), sigma$mu, sigma$se.mean, w, w/sqrt(sigma$lambda), sigma$lambda))
    
    if (w / sqrt(sigma$lambda) < eps)
      break
    
    # FIXME: Inefficient
    mh.res <- mh.weighted(min.n,
                          start.mass = mh.res$mass, start.score = mh.res$score,
                          s.min = 0, w = weights, trace = tracelevel > 1)
    one.traj <- c(one.traj, mh.res$traj)
  }
  
  list(traj = one.traj, mu = sigma$mu, se = sigma$se.mean, lower = sigma$mu - w/2, upper = sigma$mu + w/2)
}


s.min <- 0
res <- hit.n.run(weights = weights, start.mass = start.mass, start.score = start.score)
return(res$mu)

#=======================Standard MC =====================================
pval.est <- function(N, score.1 = 10, trace = TRUE) {
  v <- numeric(N)
  for (i in 1:N) {
    start.mass <- diff(c(0, sort(randperm(TOTAL_MASS - 1, sz - 1)), TOTAL_MASS))
    v[i] <- get.score(start.mass)
    if (i %% 1000 == 0) {
      if (trace) {
        cat(sprintf("iteration %d\n", i))
      }
    }
  }
  v
}

N <- 1000000
v <- pval.est(N)
est <- (length(v[v >= 15]))/N
est - 1.96*sqrt(est*(1 - est)/N)
est + 1.96*sqrt(est*(1 - est)/N)
  
