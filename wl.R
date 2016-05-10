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
  tmp <- sample.int(TOTAL_MASS, 8, replace = TRUE)
  start.mass <- (tmp/sum(tmp)*TOTAL_MASS)
  # start.mass <- diff(c(0, sort(randperm(TOTAL_MASS - 1, sz - 1)), TOTAL_MASS))
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
#       if (all(h[c(FALSE, TRUE)] > 0.6 * mean(h[c(FALSE, TRUE)])) && all(h[c(FALSE, TRUE)] > 20)) {
#         break
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
  # ll$w[c(FALSE, TRUE)] <- ll$w[c(FALSE, TRUE)] - max(ll$w[c(FALSE, TRUE)])   
  while (phi > phi.end) {
    phi <- sqrt(phi)
    if (trace) cat(sprintf("wl step: phi=%f\n", phi))
    ll <- wl.step(s.min, s.max, phi, thr = thr, ll$w, trace = trace)
    # ll$w[c(FALSE, TRUE)] <- ll$w[c(FALSE, TRUE)] - max(ll$w[c(FALSE, TRUE)])   
    print(ll$w)
  }
  ll$w
  
}
