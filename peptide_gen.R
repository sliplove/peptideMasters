setwd("/home/sliplove/Documents/Masters/gitcop/")
library(Rcpp)
library(coda)
library(lattice)
library(gtools)
sourceCpp("scorecpp/scoreR.cpp")
set.seed(42)

all.matches = read.table("all_matches.tsv", skip = 1)
all.matches.mass = subset(all_matches, select = c(V2, V5,  V8))
names(all.matches.mass) = c('id', 'maxscore', 'spectrummass')

id <- all.matches.mass$id[5]

mat <- as.matrix(read.table(paste0("./tables/matrix_", id, ".txt")))
rule <- as.matrix(read.table(paste0("./tables/rule_", id, ".txt")))
exp_spectrum <- as.numeric(read.table(paste0("./tables/spectrum_", id, ".txt")))

MASS_PROTON = 1.00728
N_MASS = ncol(mat)
TOTAL_MASS = all.matches.mass$spectrummass[all.matches.mass$id == id]
MAX_SCORE = all.matches.mass$maxscore[all.matches.mass$id == id]

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

hit.n.run <- function(weights, start.mass, start.score,
                      step = 50000, min.n = 50000, eps = 0.02,
                      level = 0.95, tracelevel = 1) {
  z <- qnorm((1 + level) / 2)
  mh.res <- mh.weighted(min.n,
                        start.mass = start.mass, start.score = start.score,
                        s.min = 0, w = weights, trace = tracelevel > 1)
  one.traj <- mh.res$traj

  sigma <- list()
  w <- 0
  prob.const <- sum(exp(-weights[one.traj - s.min + 1])) / length(one.traj)
  repeat {
    g <- function(x) (x >= MAX_SCORE) * exp(-weights[x - s.min + 1]) / prob.const
    
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

start.mass <- as.numeric(rdirichlet(1, rep(1, N_MASS)))*TOTAL_MASS
start.score <- get.score(start.mass)
start.score


res.est.unif <- hit.n.run(weights, start.mass = start.mass, start.score = start.score)
length(res.est.unif$traj)
tr <- res.est.unif$traj
res.est.unif$traj

#------------------Standard MC----------------
pval.est <- function(N, score.1 = 14, trace = TRUE) {
  res <- rdirichlet(N, rep(1, N_MASS))
  v <- apply(res, 1, function(x) get.score(x*TOTAL_MASS))
  v
}

N <- length(res.est.unif$traj) 
v <- pval.est(N)
est <- length(v[v >= SCORE_])/N
est + 1.96*sqrt(est*(1 - est)/N)
est - 1.96*sqrt(est*(1 - est)/N)

