#' TCGA data with somatic non-synonymous single
#' nucleotide variants (SNV)
#'
#' A subset of the publicly available TCGA whole exome sequencing data
#' containing non-synonymous SNV mutations. More specifically,
#' the subset of the TCGA data with
#' \code{Variant_Type == "SNP"} and \code{ Variant_Classification %in%
#' c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation").}
#'
#' @format A dataset with 1,991,488 rows and 4 columns:
#' \describe{
#' \item{Hugo_Symbol}{the gene label}
#' \item{Variant}{the variant label. Obtained by concatenating
#' the columns labeled
#' 'Hugo_Symbol', 'Chromosome', 'Start_Position', 'Tumor_Seq_Allele1',
#' and 'Tumor_Seq_Allele2' in the original TCGA data}
#' \item{patient_id}{the patient (tumor) label. Obtained by extracting
#' the first 16 characters of the column'Tumor_Sample_Barcode' in
#' the original TCGA data}
#' \item{Cancer_Code}{the cancer category associated with the tumor. Obtained by
#' matching the TSS code (the 6th and 7th characters of 'Tumor_Sample_Barcode'
#' column of the original TCGA data) with publicly available TSS study names.}
#' \item{MS}{dominant Mutation Signature of the tumor.
#' Obtained using Alexandrov 2013 mutation signature
#' calling algorithm. Contains 7 categories Non-hypermutated, Smoking (4),
#' Other, APOBEC (2, 13), UV (7), POLE (10), MMR (6, 15, 20, 26), with the numbers
#' within parenthesis indicating the Signature number in COSMIC-30 list. The category
#' Non-hypermutated corresponds to tumors with total number of mutations <= 500.}
#' }
#' @references
#' Alexandrov, Ludmil B., et al. "Signatures of mutational processes in human cancer."
#' Nature 500.7463 (2013): 415-421.
#' @source
#' https://gdc.cancer.gov/about-data/publications/mc3-2017

"tcga"
