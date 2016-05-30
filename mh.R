weight.spector <- function(score, s.min = 0) 
  log(score - s.min + 1)

mh.update <- function(mass, score, lweight.spector) {
  new.mass <- modify.mass(mass)
  new.score <- min(get.score(new.mass), MAX_SCORE)
  
  lold.weight <- lweight.spector(score)
  lnew.weight <- lweight.spector(new.score)
  
  alpha <- min(1, exp(lnew.weight - lold.weight))
  if (runif(1) < alpha) {
    score <- new.score
    mass <- new.mass
  }
  list(mass = mass, score = score)
}

mh.weighted <- function(N, w, start.mass, start.score, s.min, trace = TRUE) {
  weight.w <- function(score) { 
    w[score - s.min + 1]}
  
  i <- 1
  v <- numeric(N)
  l <- list(mass = start.mass, score = start.score)
  while (i <= N) {
    l <- do.call(mh.update, c(l, weight.w))
    l$score <- min(l$score, MAX_SCORE)

    v[i] <- l$score
    # v[i] <- l$mass[1]
    i <- i + 1
    if (i %% 10000 == 0) {
      if (trace) {
        cat(sprintf("iteration %d\n", i))
      }
    }
  }
  c(list(traj = v), l)
}
