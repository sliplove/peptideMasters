setwd("/home/sliplove/Documents/Masters/github_copy/")
library(Rcpp)
library(coda)
library(lattice)
sourceCpp("scorecpp/scoreR.cpp")
set.seed(42)

MAX_SCORE = 14
MASS_PROTON = 1.00728
TOTAL_MASS =  1007.65

mat <- as.matrix(read.table("../tables/matrix"))
rule <- as.matrix(read.table("../tables/rule_graph"))
sf <-read.table("../tables/exp_spector")
max.W = TOTAL_MASS

surfactin <- c(296.089,  324.153,  327.999,  338.181,  341.846,  359.415,  366.881,
               372.44,   386.08,   395.07,   423.098,  433.021,  437.202,  441.282,
               455.124,  462.093,  464.852,  481.362,  483.64,   508.513,  533.283,
               536.238,  540.697,  549.837,  554.36,   568.234,  569.96,   581.238,
               596.306,  616.903,  631.247,  637.379,  649.104,  667.163,  681.25,
               685.378,  686.424,  698.233,  722.808,  735.46,  746.496,  751.526,
               752.401,  764.41,   780.184,  782.373,  816.009,  820.065,  827.28,
               842.327,  847.555,  848.771,  849.451,  877.359,  893.387,  895.455,
               909.388,  913.777,  933.373,  939.125,  951.488,  959.293,  988.316,
               989.367)
surfactin <- sort(surfactin)

modify.mass <- function(mass)
  update_mass(mass, rule, FALSE)

get.score <- function(mass)
  score_peak(surfactin, mat %*% mass, TOTAL_MASS, MASS_PROTON, FALSE)

source("wl.R")
source("mh.R")
source("se.R")

#-----------------------------------------------------------------------------------------------------------------------
#test if everything is correct
weights <- wl(0, MAX_SCORE)

SCORE_ <- 14
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

generate.est <- function() {
  start.mass <- sample.int(max.W, 8, replace = TRUE)
  start.mass <- start.mass*TOTAL_MASS/sum(start.mass)
  start.score <- get.score(start.mass)
  start.score
  s.min <- 0
  res <- hit.n.run(weights = weights, start.mass = start.mass, start.score = start.score)
  return(res)
}
vec <- replicate(100, generate.est(), simplify = FALSE)

#==========Standard MC=================
pval.est <- function(N, score.1 = 14, trace = TRUE) {
  v <- numeric(N)
  for (i in 1:N) {
    tmp <- sample.int(max.W, 8, replace = TRUE)
    tmp.mass <- (tmp/sum(tmp)*TOTAL_MASS)
    v[i] <- get.score(tmp.mass)
    if (i %% 1000 == 0) {
      if (trace) {
        cat(sprintf("iteration %d\n", i))
      }
    }
  }
  v
  #return((length(v[v >= score.1]))/N)
}

N <- 1000000
v <- pval.est(N)
