#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(scales)
  library(grid)
})

# ============================================================
# 1. Paths
# ============================================================

INPUT_FILE <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Result/LDSC_gc_47_disease_epoch_all_targets.tsv"

OUTDIR <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Figure/ldsc_gc_disease_epoch_DE_organ_network_inline_wordcloud_boxes_spread_deduplicated"

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. Bonferroni threshold
# ============================================================

N_DISEASE_ENDPOINTS <- 527
N_DISEASE_EPOCH_CLOCKS <- 47

BONF_P <- 0.05 / N_DISEASE_ENDPOINTS / N_DISEASE_EPOCH_CLOCKS

message("Bonferroni threshold: P < ", signif(BONF_P, 4))

# ============================================================
# 3. Figure settings
# ============================================================

DISEASE_ORDER <- c("Asthma", "COPD", "Dementia", "MI", "Stroke")
MODALITY_ORDER <- c("MRI", "Proteomics", "Metabolomics")

X_CLOCK <- 0.00
X_DOMAIN <- 1.55

WORD_BOX_XMIN <- 2.25
WORD_BOX_XMAX <- 5.30

DOMAIN_BOX_XMIN <- 2.08
DOMAIN_BOX_XMAX <- 5.48

WORDS_PER_ROW <- 5
DOMAIN_BOX_HEIGHT_FRACTION <- 0.78

organ_colors <- c(
  "Brain"                 = "#0072B2",
  "Eye"                   = "#56B4E9",
  "Heart"                 = "#D55E00",
  "Hepatic"               = "#009E73",
  "Renal"                 = "#8B5A2B",
  "Pulmonary"             = "#7B61A8",
  "Endocrine"             = "#E69F00",
  "Immune"                = "#CC79A7",
  "Skin"                  = "#F0E442",
  "Digestive"             = "#1B9E77",
  "Metabolic"             = "#999933",
  "Adipose"               = "#A6761D",
  "Pancreas"              = "#66A61E",
  "Spleen"                = "#666666",
  "Reproductive female"   = "#E78AC3",
  "Reproductive male"     = "#4D4D4D",
  "Other"                 = "#999999"
)

disease_band_colors <- c(
  "Asthma"   = "#F3F7FB",
  "COPD"     = "#FFF7E6",
  "Dementia" = "#F3FAF3",
  "MI"       = "#FFF0EA",
  "Stroke"   = "#F5F2FA"
)

endpoint_domain_order <- c(
  "Cardiovascular",
  "Respiratory",
  "Endocrine/metabolic",
  "Immune/blood",
  "Neurologic",
  "Psychiatric/neurodevelopmental",
  "Digestive/hepatic",
  "Musculoskeletal/pain",
  "Genitourinary/reproductive",
  "Dermatologic",
  "Cancer/neoplasm",
  "Infectious",
  "Eye/ENT",
  "Pregnancy/perinatal",
  "Other"
)

# ============================================================
# 4. Helper functions
# ============================================================

clean_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

safe_rescale <- function(x, to = c(1, 3)) {
  x <- as.numeric(x)
  
  if (length(x) == 0 || all(!is.finite(x))) {
    return(rep(mean(to), length(x)))
  }
  
  finite_x <- x[is.finite(x)]
  
  if (length(unique(finite_x)) <= 1) {
    out <- rep(mean(to), length(x))
    out[!is.finite(x)] <- mean(to)
    return(out)
  }
  
  out <- scales::rescale(x, to = to, from = range(finite_x, na.rm = TRUE))
  out[!is.finite(out)] <- mean(to)
  out
}

short_modality <- function(x) {
  x <- clean_chr(x)
  
  case_when(
    x == "MRI" ~ "MRI",
    x == "Proteomics" ~ "Prot.",
    x == "Metabolomics" ~ "Metab.",
    TRUE ~ x
  )
}

shorten_organ_label <- function(x) {
  x <- clean_chr(x)
  
  case_when(
    str_detect(x, regex("Reproductive female", ignore_case = TRUE)) ~ "Repro. female",
    str_detect(x, regex("Reproductive male", ignore_case = TRUE)) ~ "Repro. male",
    TRUE ~ x
  )
}

standardize_organ_label <- function(organ_label, organ_key = NA_character_) {
  organ_label <- clean_chr(organ_label)
  organ_key <- clean_chr(organ_key)
  
  case_when(
    str_detect(organ_key, regex("^brain$", ignore_case = TRUE)) ~ "Brain",
    str_detect(organ_key, regex("^eye$", ignore_case = TRUE)) ~ "Eye",
    str_detect(organ_key, regex("^heart$", ignore_case = TRUE)) ~ "Heart",
    str_detect(organ_key, regex("^hepatic$|^liver$", ignore_case = TRUE)) ~ "Hepatic",
    str_detect(organ_key, regex("^renal$|^kidney$", ignore_case = TRUE)) ~ "Renal",
    str_detect(organ_key, regex("^pulmonary$|^lung$", ignore_case = TRUE)) ~ "Pulmonary",
    str_detect(organ_key, regex("^endocrine$", ignore_case = TRUE)) ~ "Endocrine",
    str_detect(organ_key, regex("^immune$", ignore_case = TRUE)) ~ "Immune",
    str_detect(organ_key, regex("^skin$", ignore_case = TRUE)) ~ "Skin",
    str_detect(organ_key, regex("^digestive$", ignore_case = TRUE)) ~ "Digestive",
    str_detect(organ_key, regex("^metabolic$", ignore_case = TRUE)) ~ "Metabolic",
    str_detect(organ_key, regex("^adipose$", ignore_case = TRUE)) ~ "Adipose",
    str_detect(organ_key, regex("^pancreas$", ignore_case = TRUE)) ~ "Pancreas",
    str_detect(organ_key, regex("^spleen$", ignore_case = TRUE)) ~ "Spleen",
    str_detect(organ_key, regex("^reproductive_female$", ignore_case = TRUE)) ~ "Reproductive female",
    str_detect(organ_key, regex("^reproductive_male$", ignore_case = TRUE)) ~ "Reproductive male",
    
    str_detect(organ_label, regex("Brain", ignore_case = TRUE)) ~ "Brain",
    str_detect(organ_label, regex("Eye", ignore_case = TRUE)) ~ "Eye",
    str_detect(organ_label, regex("Heart", ignore_case = TRUE)) ~ "Heart",
    str_detect(organ_label, regex("Hepatic|Liver", ignore_case = TRUE)) ~ "Hepatic",
    str_detect(organ_label, regex("Renal|Kidney", ignore_case = TRUE)) ~ "Renal",
    str_detect(organ_label, regex("Pulmonary|Lung", ignore_case = TRUE)) ~ "Pulmonary",
    str_detect(organ_label, regex("Endocrine", ignore_case = TRUE)) ~ "Endocrine",
    str_detect(organ_label, regex("Immune", ignore_case = TRUE)) ~ "Immune",
    str_detect(organ_label, regex("Skin", ignore_case = TRUE)) ~ "Skin",
    str_detect(organ_label, regex("Digestive", ignore_case = TRUE)) ~ "Digestive",
    str_detect(organ_label, regex("Metabolic", ignore_case = TRUE)) ~ "Metabolic",
    str_detect(organ_label, regex("Adipose", ignore_case = TRUE)) ~ "Adipose",
    str_detect(organ_label, regex("Pancreas", ignore_case = TRUE)) ~ "Pancreas",
    str_detect(organ_label, regex("Spleen", ignore_case = TRUE)) ~ "Spleen",
    str_detect(organ_label, regex("Reproductive female", ignore_case = TRUE)) ~ "Reproductive female",
    str_detect(organ_label, regex("Reproductive male", ignore_case = TRUE)) ~ "Reproductive male",
    TRUE ~ "Other"
  )
}

classify_endpoint_domain <- function(target_id, target_source = NA_character_) {
  x <- toupper(clean_chr(target_id))
  source <- clean_chr(target_source)
  
  case_when(
    source == "PGC" & str_detect(x, "ADHD|ASD|BIP|SCZ|OCD|AUD|AN|MDD|PTSD|TS") ~ "Psychiatric/neurodevelopmental",
    source == "PGC" & str_detect(x, "^AD$|ALZ") ~ "Neurologic",
    source == "PGC" ~ "Psychiatric/neurodevelopmental",
    
    str_detect(x, "^I9_|CARD|HEART|HYPTENS|HYPERT|ANGINA|ARRHY|AF|AORTA|IHD|MI|VTE|TIA|STROKE|CVD|CABG|REVASC|ANTIHYPERT|FG_CVD") ~ "Cardiovascular",
    str_detect(x, "^J10_|ASTHMA|COPD|PULM|BRONCH|INFLUENZA|PNEUMON|RESPIR") ~ "Respiratory",
    str_detect(x, "^E4_|T2D|DM_|DIAB|OBES|ENDO|THYROID|HYPERCHOL|HYPERLIP|LIPOPROT|METABOL|PCOS") ~ "Endocrine/metabolic",
    str_detect(x, "AUTOIMMUNE|^D3_|BLOOD|ANAEMIA|ANEMIA|IMMUNE|COAG|PURPUR|HAEMORRHAGIC|HEMORRHAGIC") ~ "Immune/blood",
    str_detect(x, "^G6_|NEURO|DEMENT|ALZHEIMER|MIGRAINE|EPILEP|SLEEPAPNO|HEADACHE|NERPLEX|ROOTPLEX|XTRAPYR|C_STROKE") ~ "Neurologic",
    str_detect(x, "^F5_|KRA_PSY|DEPRESSION|DEPRESSIO|ANXIETY|SCHIZ|BIPO|MOOD|MENTAL|PERSON|SUBSTANCE|ALCOHOL|STRESS|SOMATOFORM") ~ "Psychiatric/neurodevelopmental",
    str_detect(x, "^K11_|LIVER|IBD|IBS|GASTRO|HERNIA|DIGEST|OES|ILEUS|ORAL|PARODON|GULC|HAEMORR|REFLUX") ~ "Digestive/hepatic",
    str_detect(x, "^M13_|PAIN|DORS|ARTH|OSTEO|GOUT|SCIATICA|KNEE|HIP|LOWBACK|MYALGIA|CERVIC|SHOULDER|ROTATOR|FIBRO|RX_CODEINE|TRAMADOL") ~ "Musculoskeletal/pain",
    str_detect(x, "^N14_|RENAL|KIDNEY|URIN|URETHRA|PROST|OVARY|ENDOMETRIOSIS|FEMALE|MALEGEN|MENORRH|CYSTITIS|PYELONEPHR|PHIMOSIS|BREAST") ~ "Genitourinary/reproductive",
    str_detect(x, "^L12_|SKIN|DERMAT|PSORI|URTIC|CELLULITIS|LICHEN|SEBORR|FOLLICULAR|ABSCESS_CUT") ~ "Dermatologic",
    str_detect(x, "^C3_|^CD2_|CANCER|CARCINOMA|BREAST|PROSTATE|COLORECTAL|BRONCHUS_LUNG|MELANOCYTIC|NEOPLASM|LYMPHOID") ~ "Cancer/neoplasm",
    str_detect(x, "^AB1_|BACT|VIRAL|MYCOSES|SEPSIS|ERYSIPELAS|INFECTION|INFECTIONS|GASTROENTERITIS") ~ "Infectious",
    str_detect(x, "^H7_|^H8_|EYE|GLAUCOMA|CATARACT|RETINA|AMD|VISU|HEARING|EAR|TINNITUS|VERTIGO|CONJUNCT") ~ "Eye/ENT",
    str_detect(x, "^O15_|PREG|PUERP|DELIV|LABOUR|PRETERM|PREECLAMPS|GESTAT") ~ "Pregnancy/perinatal",
    TRUE ~ "Other"
  )
}

pretty_endpoint_label <- function(target_id, target_display = NA_character_) {
  id <- clean_chr(target_id)
  raw_upper <- toupper(id)
  
  mapped <- case_when(
    raw_upper %in% c("AD") ~ "Alzheimer disease",
    raw_upper %in% c("ADHD") ~ "ADHD",
    raw_upper %in% c("ASD") ~ "Autism",
    raw_upper %in% c("BIP") ~ "Bipolar disorder",
    raw_upper %in% c("SCZ") ~ "Schizophrenia",
    raw_upper %in% c("OCD") ~ "OCD",
    raw_upper %in% c("AUD") ~ "Alcohol use disorder",
    raw_upper %in% c("AN") ~ "Anorexia",
    raw_upper %in% c("MDD") ~ "Major depression",
    raw_upper %in% c("PTSD") ~ "PTSD",
    raw_upper %in% c("TS") ~ "Tourette syndrome",
    
    raw_upper %in% c("T2D", "T2D_WIDE", "E4_DM2NASCOMP", "DIABETES", "DM2") ~ "Type 2 diabetes",
    str_detect(raw_upper, "DIAB.*INSUL|INSUL.*DIAB") ~ "Insulin-treated diabetes",
    str_detect(raw_upper, "DM2|T2D|DIAB") ~ "Type 2 diabetes",
    str_detect(raw_upper, "OBES") ~ "Obesity",
    str_detect(raw_upper, "HYPERCHOL|HYPERLIP|LIPOPROT") ~ "Hyperlipidemia",
    str_detect(raw_upper, "THYROID") ~ "Thyroid disease",
    
    str_detect(raw_upper, "HYPTENS|HYPERTENS|HYPERT") ~ "Hypertension",
    str_detect(raw_upper, "HEARTFAIL|HEART_FAIL") ~ "Heart failure",
    str_detect(raw_upper, "MI_STRICT|MYOCARDIAL|\\bMI\\b") ~ "Myocardial infarction",
    str_detect(raw_upper, "IHD|ISCH_HEART") ~ "Ischemic heart disease",
    str_detect(raw_upper, "ANGINA") ~ "Angina",
    str_detect(raw_upper, "ARRHY|AFIB|ATRIAL") ~ "Arrhythmia",
    str_detect(raw_upper, "STROKE") ~ "Stroke",
    str_detect(raw_upper, "VTE|PULMEMB|EMBOL") ~ "Thromboembolism",
    str_detect(raw_upper, "AORTA|ANEURYSM") ~ "Aortic disease",
    
    str_detect(raw_upper, "ASTHMA") ~ "Asthma",
    str_detect(raw_upper, "COPD") ~ "COPD",
    str_detect(raw_upper, "BRONCH") ~ "Bronchitis",
    str_detect(raw_upper, "PNEUMON") ~ "Pneumonia",
    str_detect(raw_upper, "SLEEPAPNO|SLEEP_APNO") ~ "Sleep apnea",
    
    str_detect(raw_upper, "AUTOIMMUNE") ~ "Autoimmune disease",
    str_detect(raw_upper, "ANAEMIA|ANEMIA") ~ "Anemia",
    str_detect(raw_upper, "COAG") ~ "Coagulation disorder",
    
    str_detect(raw_upper, "DEMENT|ALZ") ~ "Dementia",
    str_detect(raw_upper, "MIGRAINE") ~ "Migraine",
    str_detect(raw_upper, "EPILEP") ~ "Epilepsy",
    str_detect(raw_upper, "HEADACHE") ~ "Headache",
    str_detect(raw_upper, "NEURO") ~ "Neurologic disease",
    
    str_detect(raw_upper, "DEPRESSION|DEPRESSIO") ~ "Depression",
    str_detect(raw_upper, "ANXIETY") ~ "Anxiety",
    str_detect(raw_upper, "BIPO") ~ "Bipolar disorder",
    str_detect(raw_upper, "SCHIZ") ~ "Schizophrenia",
    str_detect(raw_upper, "MOOD") ~ "Mood disorder",
    str_detect(raw_upper, "ALCOHOL") ~ "Alcohol-related disorder",
    str_detect(raw_upper, "SUBSTANCE") ~ "Substance use disorder",
    
    str_detect(raw_upper, "LIVER") ~ "Liver disease",
    str_detect(raw_upper, "IBD") ~ "Inflammatory bowel disease",
    str_detect(raw_upper, "IBS") ~ "Irritable bowel syndrome",
    str_detect(raw_upper, "GASTRO") ~ "Gastrointestinal disease",
    str_detect(raw_upper, "REFLUX|OES") ~ "Reflux disease",
    str_detect(raw_upper, "HERNIA") ~ "Hernia",
    
    str_detect(raw_upper, "PAIN") ~ "Pain",
    str_detect(raw_upper, "RX_CODEINE_TRAMADOL|CODEINE|TRAMADOL") ~ "Codeine/tramadol use",
    str_detect(raw_upper, "DORSALGIA") ~ "Dorsalgia",
    str_detect(raw_upper, "LOWBACK") ~ "Low back pain",
    str_detect(raw_upper, "ARTH|OSTEO") ~ "Arthritis",
    str_detect(raw_upper, "GOUT") ~ "Gout",
    str_detect(raw_upper, "SCIATICA") ~ "Sciatica",
    str_detect(raw_upper, "FIBRO") ~ "Fibromyalgia",
    
    str_detect(raw_upper, "RENAL|KIDNEY") ~ "Kidney disease",
    str_detect(raw_upper, "URIN|CYSTITIS|PYELONEPHR") ~ "Urinary tract disease",
    str_detect(raw_upper, "PROST") ~ "Prostate disease",
    str_detect(raw_upper, "OVARY|ENDOMETRIOSIS|MENORRH|FEMALE") ~ "Female reproductive disease",
    str_detect(raw_upper, "MALEGEN|PHIMOSIS") ~ "Male reproductive disease",
    
    str_detect(raw_upper, "PSORI") ~ "Psoriasis",
    str_detect(raw_upper, "DERMAT") ~ "Dermatitis",
    str_detect(raw_upper, "URTIC") ~ "Urticaria",
    str_detect(raw_upper, "CELLULITIS") ~ "Cellulitis",
    str_detect(raw_upper, "SKIN") ~ "Skin disease",
    
    str_detect(raw_upper, "BREAST") ~ "Breast disease",
    str_detect(raw_upper, "PROSTATE.*CANCER|CANCER_PROSTATE") ~ "Prostate cancer",
    str_detect(raw_upper, "COLORECTAL|COLON") ~ "Colorectal cancer",
    str_detect(raw_upper, "BRONCHUS_LUNG|LUNG_CANCER") ~ "Lung cancer",
    str_detect(raw_upper, "MELAN") ~ "Melanoma",
    str_detect(raw_upper, "CANCER|CARCINOMA|NEOPLASM|LYMPHOID") ~ "Cancer",
    
    str_detect(raw_upper, "SEPSIS") ~ "Sepsis",
    str_detect(raw_upper, "BACT") ~ "Bacterial infection",
    str_detect(raw_upper, "VIRAL") ~ "Viral infection",
    str_detect(raw_upper, "MYCOSES") ~ "Fungal infection",
    str_detect(raw_upper, "INFECTION|INFECTIONS") ~ "Infection",
    
    str_detect(raw_upper, "GLAUCOMA") ~ "Glaucoma",
    str_detect(raw_upper, "CATARACT") ~ "Cataract",
    str_detect(raw_upper, "RETINA|AMD") ~ "Retinal disease",
    str_detect(raw_upper, "HEARING|EAR|TINNITUS") ~ "Ear/hearing disease",
    str_detect(raw_upper, "EYE|VISU") ~ "Eye disease",
    
    str_detect(raw_upper, "PREG|GESTAT") ~ "Pregnancy-related disease",
    str_detect(raw_upper, "PREECLAMPS") ~ "Preeclampsia",
    str_detect(raw_upper, "PRETERM") ~ "Preterm birth",
    
    TRUE ~ NA_character_
  )
  
  if (!is.na(mapped) && mapped != "") {
    return(mapped)
  }
  
  label <- raw_upper
  label <- str_replace_all(label, "^FG_", "")
  label <- str_replace_all(label, "^[A-Z]+[0-9]+_", "")
  label <- str_replace_all(label, "^KELA_", "")
  label <- str_replace_all(label, "^RX_", "Medication ")
  label <- str_replace_all(label, "_EXMORE|_STRICT|_INCLAVO|_MORE|_PURCH|_REIMB", "")
  label <- str_replace_all(label, "HYPTENS", "HYPERTENSION")
  label <- str_replace_all(label, "DM2NASCOMP", "TYPE 2 DIABETES")
  label <- str_replace_all(label, "DORSALGIA", "BACK PAIN")
  label <- str_replace_all(label, "_", " ")
  label <- str_replace_all(label, "\\s+", " ")
  label <- str_trim(label)
  label <- str_to_sentence(str_to_lower(label))
  
  if (label == "" || is.na(label)) {
    label <- id
  }
  
  label
}

# ============================================================
# 5. Read and filter data
# ============================================================

if (!file.exists(INPUT_FILE)) {
  stop("Input file does not exist: ", INPUT_FILE)
}

df <- fread(INPUT_FILE) %>% as.data.frame()

required_cols <- c(
  "clock_folder",
  "target_category",
  "target_source",
  "disease_label",
  "organ_key",
  "organ_label",
  "modality",
  "target_id",
  "target_display",
  "gc_mean",
  "gc_std",
  "P"
)

missing_cols <- setdiff(required_cols, names(df))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df_base <- df %>%
  mutate(
    P = as.numeric(P),
    gc_mean = as.numeric(gc_mean),
    gc_std = as.numeric(gc_std),
    disease_label = factor(disease_label, levels = DISEASE_ORDER),
    modality = factor(modality, levels = MODALITY_ORDER),
    epoch_organ = standardize_organ_label(organ_label, organ_key),
    epoch_organ = factor(epoch_organ, levels = names(organ_colors)),
    endpoint_domain = classify_endpoint_domain(target_id, target_source),
    endpoint_domain = factor(endpoint_domain, levels = endpoint_domain_order),
    endpoint_label = as.character(mapply(pretty_endpoint_label, target_id, target_display)),
    source_clean = case_when(
      target_source == "PGC" ~ "PGC",
      TRUE ~ "FinnGen"
    ),
    rg_direction = case_when(
      gc_mean > 0 ~ "Positive",
      gc_mean < 0 ~ "Negative",
      TRUE ~ "Zero"
    ),
    neg_log10_p = -log10(P)
  )

de_df <- df_base %>%
  filter(
    target_category == "Disease_endpoint",
    target_source %in% c("FinnGen", "PGC", "FinnGen_or_other"),
    is.finite(P),
    is.finite(gc_mean),
    P < BONF_P
  )

if (nrow(de_df) == 0) {
  stop("No Bonferroni-significant disease endpoint correlations found at P < 0.05/527/47.")
}

missing_organs <- setdiff(unique(as.character(df_base$epoch_organ)), names(organ_colors))
if (length(missing_organs) > 0) {
  extra_cols <- rep("#999999", length(missing_organs))
  names(extra_cols) <- missing_organs
  organ_colors <- c(organ_colors, extra_cols)
}

sig_out <- file.path(
  OUTDIR,
  "LDSC_gc_47_disease_epoch_DE_FinnGen_PGC_Bonf_0.05_527_47_significant.tsv"
)

fwrite(de_df, sig_out, sep = "\t")

name_map_out <- file.path(
  OUTDIR,
  "LDSC_gc_47_disease_epoch_DE_endpoint_short_name_mapping.tsv"
)

de_df %>%
  distinct(target_id, target_display, endpoint_label, endpoint_domain, target_source) %>%
  arrange(endpoint_domain, endpoint_label, target_id) %>%
  fwrite(name_map_out, sep = "\t")

# ============================================================
# 6. Panel G1: organ-specific endpoint-domain association map
# ============================================================

assoc_summary <- de_df %>%
  mutate(
    epoch_clock_label = paste0(
      shorten_organ_label(as.character(epoch_organ)),
      "\n",
      short_modality(as.character(modality))
    )
  ) %>%
  group_by(
    disease_label,
    epoch_clock_label,
    epoch_organ,
    modality,
    endpoint_domain
  ) %>%
  summarise(
    n_assoc = n(),
    n_unique_targets = n_distinct(target_id),
    mean_rg = mean(gc_mean, na.rm = TRUE),
    median_rg = median(gc_mean, na.rm = TRUE),
    max_abs_rg = max(abs(gc_mean), na.rm = TRUE),
    min_p = min(P, na.rm = TRUE),
    max_neg_log10_p = max(neg_log10_p, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    endpoint_domain = factor(endpoint_domain, levels = rev(endpoint_domain_order)),
    modality_order = match(as.character(modality), MODALITY_ORDER),
    organ_order = match(as.character(epoch_organ), names(organ_colors)),
    clock_x_id = paste(disease_label, epoch_clock_label, sep = "__")
  )

clock_levels <- assoc_summary %>%
  distinct(disease_label, clock_x_id, modality_order, organ_order, epoch_clock_label) %>%
  arrange(disease_label, modality_order, organ_order, epoch_clock_label) %>%
  pull(clock_x_id)

assoc_summary <- assoc_summary %>%
  mutate(clock_x_id = factor(clock_x_id, levels = clock_levels))

assoc_out <- file.path(
  OUTDIR,
  "LDSC_gc_47_disease_epoch_DE_organ_endpoint_domain_summary.tsv"
)

fwrite(assoc_summary, assoc_out, sep = "\t")

p_assoc <- ggplot(
  assoc_summary,
  aes(
    x = clock_x_id,
    y = endpoint_domain
  )
) +
  geom_tile(
    aes(fill = mean_rg),
    color = "white",
    linewidth = 0.35,
    alpha = 0.92
  ) +
  geom_point(
    aes(size = n_assoc, color = epoch_organ),
    shape = 21,
    fill = "white",
    stroke = 0.9,
    alpha = 0.95
  ) +
  geom_text(
    aes(label = n_assoc),
    size = 2.35,
    color = "black"
  ) +
  facet_wrap(
    ~ disease_label,
    scales = "free_x",
    nrow = 1
  ) +
  scale_x_discrete(
    labels = function(x) sub("^.*__", "", x)
  ) +
  scale_y_discrete(drop = TRUE) +
  scale_fill_gradient2(
    low = "#2C5AA0",
    mid = "white",
    high = "#B4442A",
    midpoint = 0,
    name = "Mean LDSC rg"
  ) +
  scale_color_manual(
    values = organ_colors,
    name = "EPOCH organ/system",
    drop = FALSE
  ) +
  scale_size_continuous(
    range = c(2.6, 8.0),
    breaks = pretty_breaks(n = 4),
    name = "No. significant\nassociations"
  ) +
  labs(
    x = NULL,
    y = "Disease endpoint domain"
  ) +
  guides(
    fill = guide_colorbar(order = 1, barwidth = 8, barheight = 0.45),
    color = guide_legend(order = 2, override.aes = list(size = 4, fill = "white"), nrow = 2),
    size = guide_legend(order = 3)
  ) +
  theme_bw(base_size = 10.5) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "#F7F3E8", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 10.5),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 7.4),
    axis.text.y = element_text(size = 8.3),
    axis.title.y = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(4, 6, 4, 6)
  )

# ============================================================
# 7. Panel G2: spread network + deduplicated endpoint word boxes
# ============================================================

all_clock_nodes <- df_base %>%
  filter(!is.na(disease_label), !is.na(modality)) %>%
  distinct(
    clock_folder,
    disease_label,
    modality,
    organ_key,
    organ_label,
    epoch_organ
  ) %>%
  mutate(
    disease_label = factor(disease_label, levels = DISEASE_ORDER),
    modality = factor(modality, levels = MODALITY_ORDER),
    epoch_organ = factor(epoch_organ, levels = names(organ_colors)),
    clock_label = paste0(
      shorten_organ_label(as.character(epoch_organ)),
      " ",
      short_modality(as.character(modality))
    ),
    clock_full_label = paste0(
      as.character(disease_label),
      " | ",
      clock_label
    ),
    modality_order = match(as.character(modality), MODALITY_ORDER),
    organ_order = match(as.character(epoch_organ), names(organ_colors))
  ) %>%
  arrange(disease_label, modality_order, organ_order, clock_label, clock_folder)

disease_gap <- 2.2

clock_counts <- all_clock_nodes %>%
  count(disease_label, name = "n_clocks") %>%
  arrange(disease_label) %>%
  mutate(
    group_index = row_number(),
    group_start_raw = cumsum(lag(n_clocks + disease_gap, default = 0))
  )

clock_nodes <- all_clock_nodes %>%
  left_join(clock_counts, by = "disease_label") %>%
  group_by(disease_label) %>%
  arrange(modality_order, organ_order, clock_label, .by_group = TRUE) %>%
  mutate(
    within_disease_index = row_number(),
    y_raw = group_start_raw + within_disease_index
  ) %>%
  ungroup()

max_y_raw <- max(clock_nodes$y_raw, na.rm = TRUE)

clock_nodes <- clock_nodes %>%
  mutate(
    x_clock = X_CLOCK,
    y_clock = max_y_raw - y_raw + 1
  )

disease_bands <- clock_nodes %>%
  group_by(disease_label) %>%
  summarise(
    y_min = min(y_clock) - 0.48,
    y_max = max(y_clock) + 0.48,
    y_mid = mean(range(y_clock)),
    y_header = max(y_clock) + 0.82,
    .groups = "drop"
  ) %>%
  mutate(
    band_fill = disease_band_colors[as.character(disease_label)]
  )

domain_y_top <- max(clock_nodes$y_clock, na.rm = TRUE)
domain_y_bottom <- min(clock_nodes$y_clock, na.rm = TRUE)

domain_nodes <- data.frame(
  endpoint_domain = factor(endpoint_domain_order, levels = endpoint_domain_order),
  domain_order = seq_along(endpoint_domain_order)
) %>%
  semi_join(
    de_df %>% distinct(endpoint_domain),
    by = "endpoint_domain"
  ) %>%
  arrange(domain_order)

n_domain_present <- nrow(domain_nodes)

if (n_domain_present == 1) {
  domain_nodes$y_domain <- mean(c(domain_y_top, domain_y_bottom))
} else {
  domain_nodes$y_domain <- seq(
    from = domain_y_top,
    to = domain_y_bottom,
    length.out = n_domain_present
  )
}

domain_slot_height <- if (n_domain_present <= 1) {
  domain_y_top - domain_y_bottom
} else {
  abs(diff(sort(domain_nodes$y_domain))[1])
}

domain_box_height <- domain_slot_height * DOMAIN_BOX_HEIGHT_FRACTION

domain_nodes <- domain_nodes %>%
  mutate(
    x_domain = X_DOMAIN,
    y_box_top = y_domain + domain_box_height / 2,
    y_box_bottom = y_domain - domain_box_height / 2
  )

# ------------------------------------------------------------
# Deduplicated word terms
# Collapse by endpoint domain + short disease name + organ color.
# This removes same-color duplicate labels like repeated "Pain".
# ------------------------------------------------------------

word_terms <- de_df %>%
  group_by(endpoint_domain, endpoint_label, epoch_organ) %>%
  summarise(
    n_assoc = n(),
    n_unique_target_ids = n_distinct(target_id),
    target_ids = paste(sort(unique(target_id)), collapse = ";"),
    n_disease_clocks = n_distinct(clock_folder),
    min_p = min(P, na.rm = TRUE),
    max_neg_log10_p = max(neg_log10_p, na.rm = TRUE),
    mean_rg = mean(gc_mean, na.rm = TRUE),
    mean_abs_rg = mean(abs(gc_mean), na.rm = TRUE),
    max_abs_rg = max(abs(gc_mean), na.rm = TRUE),
    dominant_organ = as.character(first(epoch_organ)),
    .groups = "drop"
  ) %>%
  mutate(
    endpoint_domain = factor(endpoint_domain, levels = endpoint_domain_order),
    word = endpoint_label,
    word_weight = n_assoc,
    word_size_raw = log1p(n_assoc) + 0.12 * max_neg_log10_p
  ) %>%
  group_by(endpoint_domain) %>%
  arrange(desc(word_weight), min_p, word, dominant_organ, .by_group = TRUE) %>%
  mutate(
    word_size = safe_rescale(word_size_raw, to = c(1.45, 3.10))
  ) %>%
  ungroup()

word_terms_positioned <- word_terms %>%
  left_join(
    domain_nodes %>%
      select(endpoint_domain, y_box_top, y_box_bottom, y_domain),
    by = "endpoint_domain"
  ) %>%
  group_by(endpoint_domain) %>%
  arrange(desc(word_weight), min_p, word, dominant_organ, .by_group = TRUE) %>%
  mutate(
    word_rank = row_number(),
    n_words_in_domain = n(),
    n_word_rows = ceiling(n_words_in_domain / WORDS_PER_ROW),
    word_row = floor((word_rank - 1) / WORDS_PER_ROW) + 1,
    word_col = ((word_rank - 1) %% WORDS_PER_ROW) + 1,
    row_height = (y_box_top - y_box_bottom) / pmax(n_word_rows, 1),
    col_width = (WORD_BOX_XMAX - WORD_BOX_XMIN) / WORDS_PER_ROW,
    x_word = WORD_BOX_XMIN + (word_col - 0.5) * col_width,
    y_word = y_box_top - (word_row - 0.5) * row_height,
    angle_word = case_when(
      word_rank %% 13 == 0 ~ 10,
      word_rank %% 17 == 0 ~ -10,
      TRUE ~ 0
    )
  ) %>%
  ungroup()

dedup_word_out <- file.path(
  OUTDIR,
  "LDSC_gc_47_disease_epoch_DE_deduplicated_endpoint_words_by_domain_and_organ.tsv"
)

fwrite(word_terms_positioned, dedup_word_out, sep = "\t")

edge_df <- de_df %>%
  group_by(
    clock_folder,
    disease_label,
    epoch_organ,
    modality,
    endpoint_domain
  ) %>%
  arrange(P, .by_group = TRUE) %>%
  summarise(
    n_assoc = n(),
    n_unique_targets = n_distinct(target_id),
    mean_rg = mean(gc_mean, na.rm = TRUE),
    median_rg = median(gc_mean, na.rm = TRUE),
    mean_abs_rg = mean(abs(gc_mean), na.rm = TRUE),
    max_abs_rg = max(abs(gc_mean), na.rm = TRUE),
    min_p = min(P, na.rm = TRUE),
    max_neg_log10_p = max(neg_log10_p, na.rm = TRUE),
    .groups = "drop"
  )

edge_plot_df <- edge_df %>%
  left_join(
    clock_nodes %>%
      select(
        clock_folder,
        x_clock,
        y_clock,
        clock_label,
        clock_full_label
      ),
    by = "clock_folder"
  ) %>%
  left_join(
    domain_nodes %>%
      select(endpoint_domain, x_domain, y_domain),
    by = "endpoint_domain"
  ) %>%
  filter(
    is.finite(y_clock),
    is.finite(y_domain)
  )

edge_out <- file.path(
  OUTDIR,
  "LDSC_gc_47_disease_epoch_DE_clock_to_endpoint_domain_network_edges_all_solid_spread.tsv"
)

fwrite(edge_plot_df, edge_out, sep = "\t")

word_background_df <- domain_nodes %>%
  mutate(
    xmin = DOMAIN_BOX_XMIN,
    xmax = DOMAIN_BOX_XMAX,
    ymin = y_box_bottom,
    ymax = y_box_top
  )

plot_y_min <- min(clock_nodes$y_clock, na.rm = TRUE) - 1.2
plot_y_max <- max(clock_nodes$y_clock, na.rm = TRUE) + 2.8

p_network <- ggplot() +
  geom_rect(
    data = disease_bands,
    aes(
      xmin = -0.60,
      xmax = -0.015,
      ymin = y_min,
      ymax = y_max,
      fill = disease_label
    ),
    color = NA,
    alpha = 0.88
  ) +
  geom_rect(
    data = word_background_df,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    fill = "grey98",
    color = "grey82",
    linewidth = 0.28,
    alpha = 0.96
  ) +
  geom_curve(
    data = edge_plot_df,
    aes(
      x = x_clock,
      y = y_clock,
      xend = x_domain,
      yend = y_domain,
      color = epoch_organ,
      linewidth = n_assoc,
      alpha = abs(mean_rg)
    ),
    curvature = 0.18,
    lineend = "round"
  ) +
  geom_point(
    data = clock_nodes,
    aes(
      x = x_clock,
      y = y_clock,
      color = epoch_organ
    ),
    shape = 21,
    fill = "white",
    stroke = 1.0,
    size = 2.3,
    show.legend = FALSE
  ) +
  geom_text(
    data = clock_nodes,
    aes(
      x = x_clock - 0.035,
      y = y_clock,
      label = clock_label,
      color = epoch_organ
    ),
    hjust = 1,
    size = 2.35,
    fontface = "bold",
    show.legend = FALSE
  ) +
  geom_text(
    data = disease_bands,
    aes(
      x = -0.55,
      y = y_header,
      label = disease_label
    ),
    hjust = 0,
    vjust = 0,
    size = 3.35,
    fontface = "bold",
    color = "black"
  ) +
  geom_point(
    data = domain_nodes,
    aes(
      x = x_domain,
      y = y_domain
    ),
    shape = 21,
    fill = "grey95",
    color = "grey25",
    stroke = 0.95,
    size = 2.8,
    show.legend = FALSE
  ) +
  geom_text(
    data = domain_nodes,
    aes(
      x = x_domain + 0.035,
      y = y_domain,
      label = endpoint_domain
    ),
    hjust = 0,
    size = 2.75,
    fontface = "bold",
    color = "grey15"
  ) +
  geom_text(
    data = word_terms_positioned,
    aes(
      x = x_word,
      y = y_word,
      label = word,
      size = word_size,
      color = dominant_organ,
      angle = angle_word
    ),
    fontface = "bold",
    alpha = 0.92,
    show.legend = TRUE
  ) +
  annotate(
    "text",
    x = -0.34,
    y = plot_y_max - 0.35,
    label = "Disease EPOCH clocks",
    fontface = "bold",
    size = 3.7,
    hjust = 0.5
  ) +
  annotate(
    "text",
    x = X_DOMAIN + 0.10,
    y = plot_y_max - 0.35,
    label = "Disease endpoint domains",
    fontface = "bold",
    size = 3.7,
    hjust = 0.5
  ) +
  annotate(
    "text",
    x = (DOMAIN_BOX_XMIN + DOMAIN_BOX_XMAX) / 2,
    y = plot_y_max - 0.35,
    label = "Bonferroni-significant disease endpoints",
    fontface = "bold",
    size = 3.7,
    hjust = 0.5
  ) +
  scale_fill_manual(
    values = disease_band_colors,
    guide = "none"
  ) +
  scale_color_manual(
    values = organ_colors,
    name = "EPOCH organ/system",
    drop = FALSE
  ) +
  scale_linewidth_continuous(
    range = c(0.20, 2.60),
    name = "No. significant\nassociations"
  ) +
  scale_alpha_continuous(
    range = c(0.16, 0.80),
    guide = "none"
  ) +
  scale_size_identity(
    guide = "none"
  ) +
  coord_cartesian(
    xlim = c(-0.66, 5.65),
    ylim = c(plot_y_min, plot_y_max),
    clip = "off"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_void(base_size = 10.5) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(12, 125, 12, 120)
  )

# ============================================================
# 8. Save panels
# ============================================================

assoc_pdf <- file.path(
  OUTDIR,
  "Panel_G1_DE_epoch_organ_specific_endpoint_domain_map.pdf"
)

assoc_png <- file.path(
  OUTDIR,
  "Panel_G1_DE_epoch_organ_specific_endpoint_domain_map.png"
)

network_pdf <- file.path(
  OUTDIR,
  "Panel_G2_DE_epoch_clock_domain_network_all_solid_deduplicated_endpoint_word_boxes.pdf"
)

network_png <- file.path(
  OUTDIR,
  "Panel_G2_DE_epoch_clock_domain_network_all_solid_deduplicated_endpoint_word_boxes.png"
)

ggsave(assoc_pdf, p_assoc, width = 12.8, height = 5.2, device = cairo_pdf)
ggsave(assoc_png, p_assoc, width = 12.8, height = 5.2, dpi = 300)

ggsave(network_pdf, p_network, width = 22.0, height = 13.8, device = cairo_pdf)
ggsave(network_png, p_network, width = 22.0, height = 13.8, dpi = 300)

# ============================================================
# 9. Combined panel
# ============================================================

combined_pdf <- file.path(
  OUTDIR,
  "Panel_G_DE_epoch_disease_endpoint_network_all_solid_deduplicated_endpoint_word_boxes_combined.pdf"
)

combined_png <- file.path(
  OUTDIR,
  "Panel_G_DE_epoch_disease_endpoint_network_all_solid_deduplicated_endpoint_word_boxes_combined.png"
)

if (requireNamespace("patchwork", quietly = TRUE)) {
  
  combined_plot <- (
    p_assoc +
      ggtitle("G1  Organ-specific disease endpoint-domain associations")
  ) /
    (
      p_network +
        ggtitle("G2  Disease EPOCH clocks, endpoint domains, and deduplicated endpoint traits")
    )
  
  combined_plot <- combined_plot +
    patchwork::plot_layout(heights = c(0.68, 1.45)) +
    patchwork::plot_annotation(
      title = "Bonferroni-significant LDSC genetic correlations between disease EPOCH clocks and disease endpoints",
      subtitle = paste0(
        "FinnGen/PGC disease endpoints only; P < 0.05/527/47 = ",
        signif(BONF_P, 3),
        ". Solid curves connect disease EPOCH clocks to endpoint domains; right-side boxes show deduplicated endpoint traits as short disease names."
      ),
      theme = theme(
        plot.title = element_text(face = "bold", size = 15, hjust = 0),
        plot.subtitle = element_text(size = 10.5, hjust = 0)
      )
    )
  
  ggsave(combined_pdf, combined_plot, width = 22.0, height = 18.4, device = cairo_pdf)
  ggsave(combined_png, combined_plot, width = 22.0, height = 18.4, dpi = 300)
  
} else if (requireNamespace("cowplot", quietly = TRUE)) {
  
  combined_plot <- cowplot::plot_grid(
    p_assoc + ggtitle("G1  Organ-specific disease endpoint-domain associations"),
    p_network + ggtitle("G2  Disease EPOCH clocks, endpoint domains, and deduplicated endpoint traits"),
    ncol = 1,
    rel_heights = c(0.68, 1.45)
  )
  
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(
      paste0(
        "Bonferroni-significant LDSC genetic correlations between disease EPOCH clocks and disease endpoints\n",
        "FinnGen/PGC only; P < 0.05/527/47 = ",
        signif(BONF_P, 3)
      ),
      x = 0,
      hjust = 0,
      fontface = "bold",
      size = 12
    )
  
  combined_plot <- cowplot::plot_grid(
    title_grob,
    combined_plot,
    ncol = 1,
    rel_heights = c(0.08, 1)
  )
  
  ggsave(combined_pdf, combined_plot, width = 22.0, height = 18.4, device = cairo_pdf)
  ggsave(combined_png, combined_plot, width = 22.0, height = 18.4, dpi = 300)
  
} else {
  message("Neither patchwork nor cowplot is installed. Individual panels were saved, but combined panel was not created.")
}

# ============================================================
# 10. Summary tables
# ============================================================

summary_by_disease <- de_df %>%
  group_by(disease_label) %>%
  summarise(
    n_sig_rows = n(),
    n_clocks = n_distinct(clock_folder),
    n_targets = n_distinct(target_id),
    mean_rg = mean(gc_mean, na.rm = TRUE),
    median_rg = median(gc_mean, na.rm = TRUE),
    max_abs_rg = max(abs(gc_mean), na.rm = TRUE),
    min_p = min(P, na.rm = TRUE),
    .groups = "drop"
  )

summary_by_clock <- de_df %>%
  group_by(disease_label, clock_folder, epoch_organ, modality) %>%
  summarise(
    n_sig_rows = n(),
    n_endpoint_domains = n_distinct(endpoint_domain),
    n_targets = n_distinct(target_id),
    mean_rg = mean(gc_mean, na.rm = TRUE),
    max_abs_rg = max(abs(gc_mean), na.rm = TRUE),
    min_p = min(P, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(disease_label, desc(n_sig_rows))

summary_by_domain <- de_df %>%
  group_by(endpoint_domain) %>%
  summarise(
    n_sig_rows = n(),
    n_disease_clocks = n_distinct(clock_folder),
    n_epoch_organs = n_distinct(epoch_organ),
    n_targets = n_distinct(target_id),
    mean_rg = mean(gc_mean, na.rm = TRUE),
    median_rg = median(gc_mean, na.rm = TRUE),
    max_abs_rg = max(abs(gc_mean), na.rm = TRUE),
    min_p = min(P, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig_rows))

summary_by_organ_domain <- de_df %>%
  group_by(epoch_organ, endpoint_domain) %>%
  summarise(
    n_sig_rows = n(),
    n_disease_clocks = n_distinct(clock_folder),
    n_targets = n_distinct(target_id),
    mean_rg = mean(gc_mean, na.rm = TRUE),
    max_abs_rg = max(abs(gc_mean), na.rm = TRUE),
    min_p = min(P, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(epoch_organ, desc(n_sig_rows))

fwrite(
  summary_by_disease,
  file.path(OUTDIR, "LDSC_gc_47_disease_epoch_DE_Bonf_summary_by_disease.tsv"),
  sep = "\t"
)

fwrite(
  summary_by_clock,
  file.path(OUTDIR, "LDSC_gc_47_disease_epoch_DE_Bonf_summary_by_clock.tsv"),
  sep = "\t"
)

fwrite(
  summary_by_domain,
  file.path(OUTDIR, "LDSC_gc_47_disease_epoch_DE_Bonf_summary_by_endpoint_domain.tsv"),
  sep = "\t"
)

fwrite(
  summary_by_organ_domain,
  file.path(OUTDIR, "LDSC_gc_47_disease_epoch_DE_Bonf_summary_by_organ_and_endpoint_domain.tsv"),
  sep = "\t"
)

# ============================================================
# 11. Console output
# ============================================================

message("============================================================")
message("Finished disease endpoint LDSC organ-resolved plotting.")
message("Bonferroni threshold: ", signif(BONF_P, 4))
message("Significant rows: ", nrow(de_df))
message("Unique disease EPOCH clocks: ", dplyr::n_distinct(de_df$clock_folder))
message("Unique disease endpoint targets: ", dplyr::n_distinct(de_df$target_id))
message("")
message("Filtered significant data:")
message(sig_out)
message("")
message("Endpoint short-name mapping:")
message(name_map_out)
message("")
message("Deduplicated endpoint words:")
message(dedup_word_out)
message("")
message("Panel G1:")
message(assoc_pdf)
message(assoc_png)
message("")
message("Panel G2, deduplicated word boxes:")
message(network_pdf)
message(network_png)
message("")
message("Combined panel:")
message(combined_pdf)
message(combined_png)
message("")
message("Output directory:")
message(OUTDIR)
message("============================================================")