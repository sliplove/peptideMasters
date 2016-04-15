wl.step <- function(s.min, s.max, phi,
                    thr = 1000, w, trace = TRUE) {
  h <- numeric(s.max - s.min + 1)
  if (missing(w))
    w <- log(rep(1, length.out = s.max - s.min + 1))
  lphi <- log(phi)

  weight.wl <- function(score) {
    w[score - s.min + 1]
  }

  # FIXME
  tmp <- sample(1:max.W, 8, replace = TRUE)
  start.mass <- (tmp/sum(tmp)*TOTAL_MASS)
  start.score <- get.score(start.mass)

  score <- min(score, MAX_SCORE)
  l <- list(mass = start.mass, score = start.score)
  i <- 1
  repeat {
    l <- do.call(mh.update, c(l, weight.wl))
    score <- min(score, MAX_SCORE)
    idx <- l$score - s.min + 1
    h[idx] <- h[idx] + 1
    w[idx] <- w[idx] - lphi
    if (i > thr)
      break
    if (i %% 100 == 0) {
      if (trace) {
        cat(sprintf("iteration %d\n", i))
        # print(w)
        # print(h)
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
