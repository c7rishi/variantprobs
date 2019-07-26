# # need an unexported function (.cxx_simple_good_turing)
# # from edgeR.
# # Since CRAN does not allow `:::`, this is a
# # workaround:
#
#
getfrompackage <- function(pkg, name)  {
  # get(fun,
  #     envir = asNamespace(pkg),
  #     inherits = FALSE)
  pkg <- as.character(substitute(pkg))
  name <- as.character(substitute(name))
  get(name, envir = asNamespace(pkg), inherits = FALSE)
}



GoodTuring_orig <- function(r = NULL, N_r = NULL, conf = 1.96)
{
  n <- N_r[r > 0]
  n0 <- N_r[r == 0]
  r <- r[r > 0]
  out <- .Call(
    # edgeR:::.cxx_simple_good_turing,
    getfrompackage("edgeR", ".cxx_simple_good_turing"),
    r, n, conf
  )
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
#' @param N12_imp imputed value of N1 and N2 if either of them is 0.
#' Defaults to 1.
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
#' names(v_f) <- var_freq$Variant
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
                             N0 = NULL,
                             N12_imp = 1)  {


  estN0 <- is.null(N0)

  if (all(is.null(counts), is.null(r), is.null(N_r))) {
    stop("Either provide (a) \'counts\', or (b) \'r\' and \'N_r\'")
  }

  if(any(is.null(r), is.null(N_r))) {
    tmp <- rle(sort(unname(counts)))
    r <- tmp$values
    N_r <- tmp$lengths
  }

  r_1_impute <- r_2_impute <- FALSE

  if (any(length(N_r[r == 1]) == 0, N_r[r == 1] == 0)) {
    N_r <- c(N_r, pmax(N12_imp, 1))
    r <- c(r, 1)
    r_1_impute <- TRUE
  }

  if (any(length(N_r[r == 2]) == 0, N_r[r == 2] == 0)) {
    N_r <- c(N_r, pmax(N12_imp, 1))
    r <- c(r, 2)
    r_2_impute <- TRUE
  }

  ord <- order(r)
  r <- r[ord]
  N_r <- N_r[ord]
  names(N_r) <- r


  GT <- GoodTuring(r = r, N_r = N_r, m = m, conf = conf)

  N0est <- ifelse(estN0,
                  max(chao_N0(r = r, N_r = N_r, m = m), N0min),
                  N0)

  N1_adj <- ifelse(r_1_impute, N12_imp, N_r[r == 1])
  N2_adj <- ifelse(r_2_impute, N12_imp, N_r[r == 2])

  p0 <- N1_adj/N0est/(m+1)

  p_atleast_1new <- 1 - exp(-N1_adj/(m+1))

  # adjust the proportions for 1 & 2
  GT$proportion[1] <- GT$proportion[1] * N_r[r == 1]/N1_adj
  GT$proportion[2] <- GT$proportion[2] * N_r[r == 2]/N2_adj


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


  attr(p_GT, "N0") <- N0est
  p_GT
}
