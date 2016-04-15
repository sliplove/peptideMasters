setwd("/home/sliplove/Documents/Masters")
library(Biostrings)
library(Rcpp)
library(coda)
library(lattice)
library(BoSSA)
sourceCpp("scorecpp/scoreR.cpp")
set.seed(42)

MAX_SCORE = 14
MASS_PROTON = 1.00728
TOTAL_MASS =  1007.65

mat <- as.matrix(read.table("matrix"))
rule <- read.table("rule_graph")

max.W = 100
sz <- 8

ilya <- c(296.089,  324.153,  327.999,  338.181,  341.846,  359.415,  366.881,
             372.44,  386.08,   395.07,   423.098,  433.021,  437.202,  441.282, 
             455.124,  462.093,  464.852,  481.362,  483.64,   508.513,  533.283, 
             536.238,  540.697,  549.837, 554.36,   568.234,  569.96,   581.238,
             596.306,  616.903,  631.247,  637.379,  649.104,  667.163,  681.25,
             685.378,  686.424,  698.233,  722.808,  735.46,  746.496,  751.526,
             752.401,  764.41,   780.184,  782.373,  816.009,  820.065,  827.28,  
             842.327,  847.555,  848.771,  849.451,  877.359,  893.387,  895.455, 
             909.388,  913.777,  933.373,  939.125,  951.488,  959.293,  988.316, 
             989.367) 

modify.mass <- function(mass) {
  id <- sample(1:sz, 1)
  beg <-  rule[id, 1] + 1
  end <-  rule[id, 2] + 1
  delta <- runif(1, min = -mass[beg], max = mass[end])
  new.mass <- mass
  new.mass[beg] = mass[beg] + delta
  new.mass[end] = mass[end] - delta
  return(new.mass)
}

get.score <- function(mass) 
  get_score(sort(ilya), sort(mat%*%mass))

weight.spector <- function(score, s.min = 0) 
  log(score - s.min + 1)

mh.update <- function(mass, score, lweight.spector) {
  new.mass <- modify.mass(mass)
  new.score <- get.score(new.mass)
  if(new.score > MAX_SCORE) {
    new.score <- MAX_SCORE
  }
  
  lold.weight <- lweight.spector(score)
  lnew.weight <- lweight.spector(new.score)
  
  alpha <- min(1, exp(lnew.weight - lold.weight))
  if (runif(1) < alpha) {
    score <- new.score
    mass <- new.mass
  }
  list(mass = mass, score = score)
}

start.mass <- sample(1:max.W, 8, replace = TRUE)
start.mass <- start.mass*TOTAL_MASS/sum(start.mass)
score <- get.score(start.mass)

mh.usual <- function(N) {
  i <- 1
  v <- numeric(N)
  l <- list(mass = start.mass, score = score)
  while (i <= N) {
    l <- do.call(mh.update, c(l, weight.spector))
    if(l$score > MAX_SCORE) {
      l$score <- MAX_SCORE
    }
    
    v[i] <- l$score
    i <- i + 1
  }
  
  v
}

mh.weighted <- function(N, w, s.min, trace = TRUE) {
  weight.w <- function(score) { 
    w[score - s.min + 1]}
  
  i <- 1
  v <- numeric(N)
  l <- list(mass = start.mass, score = score)
  while (i <= N) {
    l <- do.call(mh.update, c(l, weight.w))
    if(l$score > MAX_SCORE) {
      l$score <- MAX_SCORE
    }
    v[i] <- l$score
    i <- i + 1
    if (i %% 1000 == 0) {
      if (trace) {
        cat(sprintf("iteration %d\n", i))
      }
    }
  }
  v
}

wl.step <- function(s.min, s.max, phi,
                    thr = 1000, w, trace = TRUE) {
  h <- numeric(s.max - s.min + 1)
  if (missing(w))
    w <- log(rep(1, length.out = s.max - s.min + 1))
  lphi <- log(phi)
  weight.wl <- function(score) {
    w[score - s.min + 1]
    }
  
  tmp <- sample(1:max.W, 8, replace = TRUE)
  start.mass <- (tmp/sum(tmp)*TOTAL_MASS)
  score <- get.score(start.mass)
  
  #temp!
  if (score > MAX_SCORE) {
    score <- MAX_SCORE 
  }
  l <- list(mass = start.mass, score = score)
  i <- 1
  repeat {
    l <- do.call(mh.update, c(l, weight.wl))
    #temp!!!!
    if(l$score > MAX_SCORE) {
      l$score <- MAX_SCORE
    }
    
    idx <- l$score - s.min + 1
    h[idx] <- h[idx] + 1
    w[idx] <- w[idx] - lphi
    if (i > thr)
      break
    if (i %% 100 == 0) {
      if (trace) {
        cat(sprintf("iteration %d\n", i))
#         print(w)
#         print(h)
      }
      if (all(h > 0.6 * mean(h)) && all(h > 20)) {
        break
      }
    }
    i <- i + 1
  }
  list(w = w - max(w), h = h, i = i)
}

wl <- function(s.min, s.max,
               phi.start = exp(0.3), phi.end = exp(0.0001),
               thr = 100500,
               trace = TRUE) {
  phi <- phi.start
  
  if (trace) cat(sprintf("wl step: phi=%f\n", phi))
  ll <- wl.step(s.min, s.max, phi, thr = thr, trace = trace)  
  while (phi > phi.end) {
    phi <- sqrt(phi)
    if (trace) cat(sprintf("wl step: phi=%f\n", phi))
    ll <- wl.step(s.min, s.max, phi, thr = thr, ll$w, trace = trace)
  }
  ll$w
}
#-----------------------------------------------------------------------------------------------------------------------
#test if everything is correct 
weights <- wl(0, MAX_SCORE)
s.min <- 0
ww <- exp(weights)
L <- 1000000      
one.traj <- mh.weighted(L, s.min = 0, w = weights)
one.traj
plot(cumsum(one.traj < 10)/1:L, lwd = 3, ty = "l", xlab = "number of simulations", ylab = "Cumulative mean")

burn.in <- 2e+05
eff.traj <- one.traj[-(1:burn.in)]
# autocorr.diag(mcmc(one.traj), lags = c(100, 500, 1000))
eff.traj <- eff.traj[seq(1, length(eff.traj), 500)]

eff.traj <- one.traj
prob.const <- sum(exp(-weights[eff.traj - s.min + 1])) / length(eff.traj)
g <- function(x) ifelse(x >= 10, exp(-weights[x - s.min + 1]), 0) / prob.const
my.est <- mean(g(eff.traj))
my.est
s.min
#==========Standard MC=================
pval.est <- function(N, score.1 = 10, trace = TRUE) {
  v <- numeric(N)
  for (i in 1:N) {
    tmp <- sample(1:max.W, 8, replace = TRUE)
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
