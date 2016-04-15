setwd("/home/sliplove/Documents/Masters")
library(Rcpp)
library(coda)
library(lattice)
sourceCpp("scorecpp/scoreR.cpp")
set.seed(42)

MAX_SCORE = 14
MASS_PROTON = 1.00728
TOTAL_MASS =  1007.65

mat <- as.matrix(read.table("matrix"))
rule <- as.matrix(read.table("rule_graph"))

max.W = 100

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

modify.mass <- function(mass) {
  id <- sample.int(nrow(rule), 1)
  beg <-  rule[id, 1] + 1
  end <-  rule[id, 2] + 1
  beg.mass <- mass[beg]
  end.mass <- mass[end]
  delta <- runif(1, min = -beg.mass, max = end.mass)

  mass[beg] <- beg.mass + delta
  mass[end] <- end.mass - delta

  mass
}

get.score <- function(mass)
  score_peak(surfactin, mat %*% mass)

source("wl.R")
source("mh.R")
source("se.R")

#-----------------------------------------------------------------------------------------------------------------------
#test if everything is correct
weights <- wl(0, MAX_SCORE)

start.mass <- sample.int(max.W, 8, replace = TRUE)
start.mass <- start.mass*TOTAL_MASS/sum(start.mass)
start.score <- get.score(start.mass)

s.min <- 0
ww <- exp(weights)
L <- 1000000
one.traj <- mh.weighted(L,
                        start.mass = start.mass, start.score = start.score,
                        s.min = 0, w = weights)

one.traj
plot(cumsum(one.traj < 10)/1:L, lwd = 3, ty = "l", xlab = "number of simulations", ylab = "Cumulative mean")

burn.in <- 2e+05
eff.traj <- one.traj[-(1:burn.in)]
# autocorr.diag(mcmc(one.traj), lags = c(100, 500, 1000))
eff.traj <- eff.traj[seq(1, length(eff.traj), 500)]

# eff.traj <- one.traj
prob.const <- sum(exp(-weights[eff.traj - s.min + 1])) / length(eff.traj)
g <- function(x) (x >= 14) * exp(-weights[x - s.min + 1]) / prob.const
my.est <- mean(g(eff.traj))

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
  repeat {
    prob.const <- sum(exp(-weights[one.traj - s.min + 1])) / length(one.traj)
    g <- function(x) (x >= 14) * exp(-weights[x - s.min + 1]) / prob.const

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

oo <- hit.n.run(weights = weights, start.mass = start.mass, start.score = start.score)

print(my.est)
s.min
#==========Standard MC=================
pval.est <- function(N, score.1 = 10, trace = TRUE) {
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
  return((length(v[v >= score.1]))/N)
}

N <- 100000
pval.est(N)
