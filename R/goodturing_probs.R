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
  tmp$proportion <- tmp$proportion * sum(as.numeric(r*N_r))/(m+1)
  tmp
}

#' Computing Good Turing probabilities of encountering gene variants
#' (including hitherto unobserved variants) based on training gene
#' mutation frquencies
#'
#' @inheritParams sgt_Delta
#' @param N0 the total number of unobserved variants. If \code{NULL},
#' \code{N0} is estimated using Chao's formula.
#' @param N0min the minimum value of N0, if known, to be used while
#' estiamting \code{N0}. Ignored if N0 is not \code{NULL}.
#' @examples
#' \dontrun{
#' # load tcga data
#' data("tcga")
#' tcga <- data.table::setDT(tcga)
#'
#' # calculate variant frequencies for KRAS
#' var_freq <- tcga[Hugo_Symbol == "KRAS",
#'                  .(v_f = length(unique(patient_id))),
#'                  by = .(Hugo_Symbol, Variant)
#'                  ]
#' v_f <- var_freq$v_f
#' names(v_f) <- var_freq$HGVSp_Short
#'
#' # calculate cohort size
#' m <- length(unique(tcga$patient_id))
#'
#'
#' # Good Turing estimates
#' goodturing_probs(counts = v_f, m = m)
#' }
#' @export
goodturing_probs <- function(counts = NULL,
                            r = NULL, N_r = NULL,
                            m = NULL, conf = 1.96,
                            N0min = 0,
                            N0 = NULL) {

  estN0 <- is.null(N0)

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

    N0est <- ifelse(estN0,
                    max(chao_est(r = r, N_r = N_r, m = m), N0min),
                    N0)

    p0 <- N_r[r == 1]/N0est/(m+1)

    # browser()
    # p_atleast_1new <- 1 - exp(N0est*log(1-p0))
    p_atleast_1new <- 1 - exp(-N_r[r == 1]/m)


    if (is.null(counts) | is.null(names(counts))) {

      p_GT <- c(p_atleast_1new, p0, GT$proportion)
      names(p_GT) <- c("atleast_1new",
                       "0",
                       GT$count)
    } else {
      tmp_probs <- GT$proportion
      names(tmp_probs) <- GT$count



      p_GT <- unname(c(tmp_probs[as.character(counts)],
                       p0,
                       p_atleast_1new))
      names(p_GT) <- c(names(counts),
                       "each_unseen",
                       "atleast_1new")

    }
  }
  p_GT
}
