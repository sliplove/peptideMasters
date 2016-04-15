setwd("/home/sliplove/Documents/Masters/")
library(Biostrings)
library(Rcpp)
library(coda)
library(BoSSA)
sourceCpp("scorecpp/scoreR.cpp")
set.seed(42)

MAX_SCORE = 7
TOTAL_MASS =  70
sz = 3
tbl <- read.table("alltstpep.txt")
ilya <- c(10,20,30,40,50,60,70)
rule <- read.table("rg_tst")

modify.mass <- function(mass) {
  id <- sample(1:sz, 1)
  beg <-  rule[id, 1] + 1
  end <-  rule[id, 2] + 1
  delta <- sample((-mass[beg]):mass[end], size = 1, replace = TRUE)
  new.mass <- mass
  new.mass[beg] = mass[beg] + delta
  new.mass[end] = mass[end] - delta
  return(new.mass)
}

mass.to.spectrum <- function(mass)
  sort(unique(c(mass[1], mass[2], mass[3], mass[1] + mass[2], mass[2] + mass[3], mass[1] + mass[3], sum(mass))))

get.score <- function(mass) 
  get_score(ilya, sort(mass.to.spectrum(mass)))


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

x1 <- sample(TOTAL_MASS - 2, 1)
x2 <- sample(TOTAL_MASS - 1 - x1, 1)
x3 <- 70 - x1 - x2
start.mass <- sample(c(x1, x2, x3), 3)
start.mass
mass.to.spectrum(start.mass)
score <- get.score(start.mass)
score
ilya

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
    if (i %% 10000 == 0) {
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
  
  x1 <- sample(TOTAL_MASS - 2, 1)
  x2 <- sample(TOTAL_MASS - 1 - x1, 1)
  x3 <- 70 - x1 - x2
  start.mass <- sample(c(x1, x2, x3), 3)
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
        print(w)
        print(h)
      }
      if (all(h[c(FALSE, TRUE)] > 0.6 * mean(h[c(FALSE, TRUE)])) && all(h[c(FALSE, TRUE)] > 30)) {
        break
      }
    }
    i <- i + 1
  }
  list(w = w - max(w), h = h, i = i)
}

wl <- function(s.min, s.max,
               phi.start = exp(0.6), phi.end = exp(0.00003),
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
# 
weights <- wl(0, MAX_SCORE)
MM <- max(weights[c(FALSE, TRUE)])
weights[c(FALSE, TRUE)] <- weights[c(FALSE, TRUE)] - MM 
# s.min <- 0
ww <- exp(weights)
# 
L <- 1000000      
one.traj <- mh.weighted(L, s.min = 0, w = weights)
plot(cumsum(one.traj < 3)/1:L, lwd = 3, ty = "l", xlab = "number of simulations", ylab = "Cumulative mean")

burn.in <- 1e+05
eff.traj <- one.traj[-(1:burn.in)]
# autocorr.diag(mcmc(one.traj), lags = c(100, 500, 1000))
eff.traj <- eff.traj[seq(1, length(eff.traj), 500)]
length(eff.traj)
# eff.traj <- one.traj
s.min <- 0

prob.const <- sum(exp(-weights[eff.traj - s.min + 1])) / length(eff.traj)
g <- function(x) (x >= 7) * exp(-weights[x - s.min + 1]) / prob.const
my.est <- mean(g(eff.traj))
my.est

res.wl <- se.bm(eff.traj, g)
res.wl$mu + qnorm((1 - .95)/2)*res.wl$se.mean
res.wl$mu + qnorm((1 + .95)/2)*res.wl$se.mean
  
#=======================Standard MC ============================================
pval.est <- function(N, score.1 = 10, trace = TRUE) {
  v <- numeric(N)
  for (i in 1:N) {
    id <- sample(2346, 1 )
    start.mass <- as.numeric(tbl[id,]) 
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
est <- (length(v[v >= 7]))/N
est - 1.96*sqrt(est*(1 - est)/N)
est + 1.96*sqrt(est*(1 - est)/N)
tbl <- data.frame(matrix(NA, 2346, 3))
