#' Smooth Good-Toulmin estimate of \eqn{\Delta(t)}, the (expected)
#' number of new variants in
#' a future (test) cohort that is \eqn{t} times as large as the
#' training cohort
#' @param counts vector of counts or frequencies of the observed variants.
#' @param r unique frequencies.
#' @param N_r frequency of frequency r.
#' @param  m training cohort size.
#' @param t positive scalar. The proportion of the future (test) cohort size to the
#' training cohort size.
#' @param adj logical. Should the Orlitsky et al. adjustment be used?
#' Defaults to \code{TRUE}. Ignored if \code{t < 1}.
#' @note  Either (a) \code{counts}, or (b) \code{r} and
#' \code{N_r} must be provided.
#' @details Computes the original Good Toulmin (1956) estimate of \eqn{\Delta(t)} if \code{t <= 1}. If
#' \code{t > 1}, the Efron-Thisted estimate (if \code{adj = FALSE}) or the
#' Efron-Thisted estimate with Orlitsky et al. (2016) adjustment (if \code{adj = TRUE}) is computed. Also
#' returns an approximate standard error ("se") of the estimate as an attribute, computed using the formula
#' provided in Efron-Thisted (1976, equation 5.2).
#'
#' @references
#' Good, I. J., & Toulmin, G. H. (1956). The number of
#' new species, and the increase in population coverage,
#' when a sample is increased. Biometrika, 43(1–2), 45–63.
#' https://doi.org/10.1093/biomet/43.1-2.45.
#'
#' Efron, B., & Thisted, R. (1976). Estimating the Number
#' of Unseen Species: How Many Words Did Shakespeare
#' Know? Biometrika, 63(3), 435–447.
#' Retrieved from http://www.jstor.org/stable/2335721.
#'
#' Orlitsky, A., Suresh, A. T., & Wu, Y. (2016). Optimal
#' prediction of the number of unseen species. Proceedings
#' of the National Academy of Sciences, 113(47), 13283–13288.
#' https://doi.org/10.1073/pnas.1607774113
#'
#' @examples
#' \dontrun{
#' # load tcga data
#' data("tcga")
#' tcga <- data.table::setDT(tcga)
#'
#' # calculate variant frequencies
#' var_freq <- tcga[,
#'                  .(v_f = length(unique(patient_id))),
#'                  by = .(Hugo_Symbol, Variant)
#'                  ]
#'
#' # calculate cohort size
#' m <- length(unique(tcga$patient_id))
#'
#'
#' # SGT Delta(t) estimate for t = 0.5, 1, 10
#' sgt_Delta(counts = var_freq$v_f, m = m, t = 0.5)
#' sgt_Delta(counts = var_freq$v_f, m = m, t = 1)
#' sgt_Delta(counts = var_freq$v_f, m = m, t = 10)
#' }
#'
#' @export
sgt_Delta <- function(counts = NULL,
                r = NULL, N_r = NULL,
                m, t = 1, adj = TRUE) {

  if (all(is.null(counts), is.null(r), is.null(N_r))) {
    stop("Either provide (a) \'counts\', or (b) \'r\' and \'N_r\'")
  }

  if(any(is.null(r), is.null(N_r))) {
    tmp <- rle(sort(unname(counts)))
    r <- tmp$values
    N_r <- tmp$lengths
  }


  N_r_obs <- N_r
  names(N_r_obs) <- r


  if (any(length(c(r, N_r)) == 0)) {
    return(NA)
  } else {
    one_to_maxr <- seq_len(max(r))
    N_r <- rep(0, length(one_to_maxr))
    names(N_r) <- one_to_maxr
    N_r[names(N_r_obs)] <- N_r_obs

    if (t <= 1) {

      h_r <- - (-t)^one_to_maxr

    } else {

      if (adj) {
        theta <- 2/(2 + t)
        k <- floor(0.5 * logb(m * t^2/(t - 1), 3))
      } else {
        theta <- 1/(1 + t)
        k <- floor(0.5 * logb(m * t^2/(t - 1), 2))
      }

      h_r <- -(-1)^one_to_maxr *
        exp(
          one_to_maxr * log(t) +
            pbinom(
              q = one_to_maxr - 1,
              size = k,
              prob = theta,
              lower.tail = FALSE,
              log = TRUE
            )
        )
    }

    SGT <- sum(h_r * N_r)

    attr(SGT, "se") <- sqrt(sum(h_r^2 * N_r))

    SGT
  }
}
