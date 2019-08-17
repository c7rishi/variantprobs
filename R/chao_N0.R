#' The Chao estimate of the total number of unseen variants
#' based on training gene mutation frequencies
#'
#' @inheritParams sgt_Delta
#'
#' @details Calculates Chao (1987) estimate of the total number of
#' unseen variants. Also provides an approximate standard error ("se") of the
#' estimate as an attribute, computed using the formula provided in Chao (1987).
#'
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
#' # SGT estimate for t = 0.5, 1, 10
#' chao_N0(counts = var_freq$v_f, m = m, t = 0.5)
#' }
#'
#' @export
chao_N0 <- function(counts = NULL, r = NULL, N_r = NULL,
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
    N0chao <- ifelse(Q2 > 0,
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
