# source("efron_thisted.R")
GoodTuring_orig <- function(r = NULL, N_r = NULL, conf = 1.96)
{
  n <- N_r[r > 0]
  n0 <- N_r[r == 0]
  r <- r[r > 0]
  out <- .Call(edgeR:::.cxx_simple_good_turing, r, n, conf)
  names(out) <- c("P0", "proportion")
  out$count <- r
  out$n <- n
  out$n0 <- n0
  out
}


GoodTuring <- function(r = NULL, N_r = NULL,
                       m = NULL, conf = 1.96) {
  tmp <- GoodTuring_orig(r = r, N_r = N_r, conf = conf)
  tmp$proportion <- tmp$proportion * sum(as.numeric(r*N_r))/m
  tmp
}

#' computing good turing probabilities
#'
#' @export
GoodTuringProbs <- function(counts = NULL,
                            r = NULL, N_r = NULL,
                            m = NULL, conf = 1.96,
                            N0min = 0) {

  if (all(is.null(counts), is.null(r), is.null(N_r))) {
    stop("Either provide (a) \'counts\', or (b) \'r\' and \'N_r\'")
  }

  if(any(is.null(r), is.null(N_r))) {
    tmp <- rle(sort(unname(counts)))
    r <- tmp$values
    N_r <- tmp$lengths
  }

  if (length(N_r[r==1]) == 0) {
    p_GT <- NA_real_
  }
  else {
    GT <- GoodTuring(r = r, N_r = N_r, m = m, conf = conf)

    N0est <- max(chao_est(r = r, N_r = N_r, m = m), N0min)

    p0 <- N_r[r == 1]/N0est/m

    # browser()
    # p_atleast_1new <- 1 - exp(N0est*log(1-p0))
    p_atleast_1new <- 1 - exp(-N_r[r == 1]/m)


    if (is.null(counts)) {

    p_GT <- c(p_atleast_1new, p0, GT$proportion)
    names(p_GT) <- c("atleast_1new",
                     "0",
                     GT$count)
    } else {
      tmp_probs <- GT$proportion
      names(tmp_probs) <- GT$count

    }
  }
  p_GT
}
