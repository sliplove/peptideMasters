wl.step <- function(s.min, s.max, phi,
                    thr = 100000, w, trace = TRUE) {
  h <- numeric(s.max - s.min + 1)
  if (missing(w))
    w <- log(rep(1, length.out = s.max - s.min + 1))
  lphi <- log(phi)

  weight.wl <- function(score) {
    w[score - s.min + 1]
  }

  # FIXME
  start.mass <- as.numeric(rdirichlet(1, rep(1, N_MASS))*TOTAL_MASS)
  start.score <- get.score(start.mass)

  start.score <- min(start.score, MAX_SCORE)
  l <- list(mass = start.mass, score = start.score)
  i <- 1
  repeat {
    l <- do.call(mh.update, c(l, weight.wl))
    score <- min(l$score, MAX_SCORE)
    idx <- l$score - s.min + 1
    h[idx] <- h[idx] + 1
    w[idx] <- w[idx] - lphi
    if (i > thr)
      break
    if (i %% 1000 == 0) {
      if (trace) {
        cat(sprintf("iteration %d\n", i))
        print(w)
        print(h)
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
               phi.start = exp(0.6), phi.end = exp(0.0003),
               thr = 100500,
               trace = TRUE) {
  phi <- phi.start
  if (trace) cat(sprintf("wl step: phi=%f\n", phi))
  ll <- wl.step(s.min, s.max, phi, thr = thr, trace = trace)  
  ll$w <- ll$w - max(ll$w)   
  while (phi > phi.end) {
    phi <- sqrt(phi)
    if (trace) cat(sprintf("wl step: phi=%f\n", phi))
    ll <- wl.step(s.min, s.max, phi, thr = thr, ll$w, trace = trace)
    ll$w <- ll$w - max(ll$w)   
    print(ll$w)
  }
  ll$w
  
}
