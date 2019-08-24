library(variantprobs)
library(data.table)

data("tcga")

# find frequencies of BRCA and TTN
# variants separately for each cancer type

dt.vf <- tcga[
  Hugo_Symbol %in% c("BRCA", "TTN") &
    !is.na(Cancer_Code)
  ][,
    .(
      Hugo_Symbol = Hugo_Symbol[1],
      v_f = length(unique(patient_id))
    ),
    keyby = .(Variant, Cancer_Code)
    ]


# dcast to get a n_variant x n_cancer
# variant frequency matrix
dt.vf.dcast <- dcast(
  dt.vf[,
        .(Cancer_Code, Hugo_Symbol, Variant, v_f)
        ],
  Hugo_Symbol + Variant ~ Cancer_Code,
  value.var = "v_f",
  fill = 0
)

# calculate no. of patients for each cancer type
canc_npatient <- tcga[
  !is.na(Cancer_Code),
  .(cancer_npatient = length(unique(patient_id))),
  by = Cancer_Code
  ]
npatient <- canc_npatient$cancer_npatient
names(npatient) <- canc_npatient$Cancer_Code


# replace the observed cancer specific
# variant frequencies by the corresponding
# Good Turing probabilities
dt.vf.probs <- copy(dt.vf.dcast)[
  ,
  (names(npatient)) := mapply(
    function(counts, m) {
      tmp <- goodturing_probs(counts = counts, m = m)
      tmp[as.character(counts)]
    },
    .SD,
    npatient,
    SIMPLIFY = FALSE
  ),
  .SDcols = names(npatient)
  ]
