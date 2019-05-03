#' MSK-IMPACT data on somatic non-synonymous single
#' neucleotide variants (SNV)
#'
#' A subset of the publicly available MSK-IMPACT (a targeted clinical gene panel)
#' sequencing data containing non-synonymous SNV mutations. More specifically,
#' the subset of the MSK-IMPACT data with
#' \code{Variant_Type == "SNP"} and \code{ Variant_Classification %in%
#' c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation").}
#' The genes MLL3, MLL2, MLL, and MYCL1 are recoded as
#' KMT2C, KMT2D, KMT2A and MYCL respectively.
#'
#' @format A dataset with 68,919 rows and 5 columns:
#' \describe{
#' \item{Hugo_Symbol}{the gene label}
#' \item{Variant}{the variant label. Obtained by conctatenating
#' the columns labeled
#' 'Hugo_Symbol', 'Chromosome', 'Start_Position', 'Tumor_Seq_Allele1',
#' and 'Tumor_Seq_Allele2' in the original MSK-IMPACT data}
#' \item{patient_id}{the patient (tumor) label. Obtained by extracting
#' the first 9 characters of the column'Tumor_Sample_Barcode' in
#' the original MSK-IMPACT data}
#' \item{Cancer_Type}{the cancer category associated with the tumor.}
#'
#' }
#' @source
#' \describe{
#' \item{original MSK-IMPACT data}{...}
#'
#' }

"impact"
