find_n_r <- function(counts) {
  tmp <- rle(sort(unname(counts)))
  r <- tmp$values
  n_r <- tmp$lengths
  names(n_r) <- r
  n_r
}
