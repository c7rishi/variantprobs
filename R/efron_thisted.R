ET <- function(counts = NULL, r = NULL, N_r = NULL,
               m, t, adj = TRUE) {
  # browser()

  # t <- ceiling(t)

  if(any(is.null(r), is.null(N_r))) {
    tmp <- rle(sort(unname(counts)))
    r <- tmp$values
    N_r <- tmp$lengths
  }


  # N <- sum(as.numeric(r * N_r))
  N_r_obs <- N_r
  names(N_r_obs) <- r


  if (any(length(c(r, N_r)) == 0)) {
    return(NA)
  } else {

    one_to_maxr <- seq_len(max(r))
    N_r <- rep(0, length(one_to_maxr))
    names(N_r) <- one_to_maxr
    N_r[names(N_r_obs)] <- N_r_obs

    # efron thisted estimator

    if (t <= 1) {
      (h_r <- (-1)^(one_to_maxr+1) * t)
    } else if (adj) {
      theta <- 2/(2+t)
      k <- floor(0.5 * logb(m*t^2/(t-1), 3))
      (h_r <- (-1)^(one_to_maxr+1) *
          exp(one_to_maxr * log(t) +
                pbinom(q = one_to_maxr - 1,
                       size = k,
                       prob = theta,
                       lower.tail = FALSE,
                       log.p = TRUE)))
    } else {
      theta <- 1/(1+t)
      k <- floor(0.5 * logb(m*t^2/(t-1), 2))
      (h_r <- (-1)^(one_to_maxr+1) *
          exp(one_to_maxr * log(t) +
                pbinom(q = one_to_maxr - 1,
                       size = k,
                       prob = theta,
                       lower.tail = FALSE,
                       log.p = TRUE)))
    }

    floor(sum(N_r * h_r))
  }
}

chao_est <- function(counts = NULL, r = NULL, N_r = NULL,
                     m) {

  if(any(is.null(r), is.null(N_r))) {
    tmp <- rle(sort(unname(counts)))
    r <- tmp$values
    N_r <- tmp$lengths
  }

  if (any(length(c(r, N_r)) == 0)) {
    return(NA)
  } else {
    # if (missing(names(N_r)))
    names(N_r) <- r
    A <- (m-1)/m

    Q1 <- N_r["1"]
    Q2 <- ifelse(all(N_r["2"] != 0, !is.na(N_r["2"])),
                  N_r["2"],
                  0)

    Sobs <- sum(N_r)
    N0chao <- A * ifelse(Q2 > 0,
                         0.5 * A * Q1^2/Q2,
                         0.5 * A * Q1*(Q1 - 1))

    var_N0chao <- ifelse(Q2 > 0,
                         Q2 * (0.5*A*(Q1/Q2)^2 + A^2*(Q1/Q2)^3 +
                                 0.25*A^2*(Q1/Q2)^4),
                         (0.5*A*Q1*(Q1-1)/2 + 0.25*A^2*Q1*(Q1-1)^2 -
                            0.25*A^2*Q1^4/(N0chao+Sobs))
                         )
    if (any(is.na(var_N0chao), var_N0chao < 0))
      var_N0chao <- NA

    attr(N0chao, "se") <- sqrt(var_N0chao)

    N0chao
  }
}


chao_ice_est <- function(counts = NULL, r = NULL,
                         N_r = NULL, m, m_infreq = NULL,
                         k = 10)
  {

  if (is.null(m_infreq)) m_infreq <- m

  if (any(is.null(r), is.null(N_r))) {
    tmp <- rle(sort(unname(counts)))
    r <- tmp$values
    N_r <- tmp$lengths
  }


  # browser()

  # N <- sum(as.numeric(r * N_r))
  N_r_obs <- N_r
  names(N_r_obs) <- r


  if (any(length(c(r, N_r)) == 0)) {
    return(NA)
  } else {

    one_to_maxr <- seq_len(max(r, 2*k))
    N_r <- rep(0, length(one_to_maxr))
    names(N_r) <- one_to_maxr
    N_r[names(N_r_obs)] <- N_r_obs

    S_freq <- sum(N_r[-(1:k)])

    S_infreq <- sum(N_r[1:k])
    C_infreq <- 1 - N_r["1"]/sum((1:k) * N_r[1:k])

    tmp <- sum(as.numeric((1:k)*N_r[1:k]))

    gamma2_infreq <- max( S_infreq/C_infreq * m_infreq/(m_infreq - 1) *
                            sum(as.numeric((1:k)*(0:(k-1))*N_r[1:k])) / (tmp*(tmp-1)) - 1, 0)


    unname(ceiling(S_freq + S_infreq/C_infreq + N_r[1]/C_infreq*gamma2_infreq -
                     sum(as.numeric(N_r))))

  }
}

find_N_r <- function(counts) {
  tmp <- rle(sort(unname(counts)))
  r <- tmp$values
  N_r <- tmp$lengths
  names(N_r) <- r
  N_r
}



fit_trunc_betabin <- function(freq = NULL, r = NULL, N_r = NULL, m) {

  if(any(is.null(r), is.null(N_r))) {
    tmp <- rle(sort(unname(freq)))
    r <- tmp$values
    N_r <- tmp$lengths
  }

  neg_llik_trunc_beta_binom <- function(alpha, beta) {
    - sum(N_r * emdbook::dbetabinom(x = r, size = m,
                                    shape1 = alpha,
                                    shape2 = beta,
                                    log=TRUE)) +
      sum(N_r) * log(1 - emdbook::dbetabinom(0, size = m,
                                             shape1 = alpha,
                                             shape2 = beta,
                                             log = FALSE))

  }


  fit_full <- optim(c(log(0.1), log(0.1)),
                    fn = function(x)
                      neg_llik_trunc_beta_binom(exp(x[1]), exp(x[2])))



  est <- fit_full$par
  (a_hat <- exp(est[1]))
  (b_hat <- exp(est[2]))

  list(a = a_hat, b = b_hat)
}


parametric_N0_delta <- function(counts = NULL, r = NULL, N_r = NULL,
                                m, t = 1) {
  if(any(is.null(r), is.null(N_r))) {
    tmp <- rle(sort(unname(counts)))
    r <- tmp$values
    N_r <- tmp$lengths
  }

  fit <- tryCatch(fit_trunc_betabin(r = r, N_r = N_r, m = m),
                  error = function(e) list(a = NA, b = NA))

  N0 <- N_r[1]/fit$a * (m + fit$b - 1)/m
  Delta <- N0 * (1 - exp(lbeta(fit$a, (t+1)*m + fit$b) - lbeta(fit$a, fit$b)))
  names(Delta) <- paste("parametric_Delta_", t)

  ceiling(c(N0 = N0, Delta))
}
