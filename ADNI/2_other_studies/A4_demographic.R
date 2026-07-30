library(readr)
library(dplyr)

prediction_file <- paste0(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/",
  "results_a4_longitudinal_ad_epoch/",
  "a4_adni_brain_mri_ad_epoch_scan_level_predictions.tsv"
)

subjinfo_file <- paste0(
  "/Users/hao/cubic-projects/MULTI/download/A4/Clinical/",
  "Derived_Data/SUBJINFO.csv"
)

pred <- read_tsv(
  prediction_file,
  show_col_types = FALSE
)

subj <- read_csv(
  subjinfo_file,
  show_col_types = FALSE
)

id_col <- intersect(
  c("BID", "participant_id", "PTID"),
  names(pred)
)[1]

a4_demo <- pred |>
  transmute(
    BID = as.character(.data[[id_col]]),
    Age = as.numeric(Age),
    followup = as.numeric(years_since_external_baseline)
  ) |>
  group_by(BID) |>
  arrange(followup, Age, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  left_join(
    subj |>
      transmute(
        BID = as.character(BID),
        SEX = as.numeric(SEX)
      ),
    by = "BID"
  )

a4_followup <- pred |>
  group_by(BID = as.character(.data[[id_col]])) |>
  summarise(
    maximum_followup = max(
      as.numeric(years_since_external_baseline),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

a4_summary <- a4_demo |>
  summarise(
    N = n_distinct(BID),
    mean_age = mean(Age, na.rm = TRUE),
    sd_age = sd(Age, na.rm = TRUE),
    female_n = sum(SEX == 1, na.rm = TRUE),
    female_pct = 100 * mean(SEX == 1, na.rm = TRUE)
  ) |>
  mutate(
    maximum_followup_years = max(
      a4_followup$maximum_followup,
      na.rm = TRUE
    )
  )

print(a4_summary)