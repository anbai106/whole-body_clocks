import pandas as pd
#### here is to check what data columns we should use to develop the AD L'EPOCH
istaging_pickle_file = '/Users/hao/cubic-projects/ISTAGING/Pipelines/ISTAGING_Data_Consolidation_2020/v2.0/istaging.pkl.gz'
data = pd.read_pickle(istaging_pickle_file)
df_adni = data.loc[data['Study'].isin(['ADNI'])]
# df_adni_test = df_adni.iloc[:1]
# df_adni_test.to_csv("~/test.tsv", sep="\t", index=False)
df_adni.to_csv("/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/adni_istaging.tsv", sep="\t", index=False)
print(df_adni.columns.to_list())


df_adni["Date_parsed"] = pd.to_datetime(
    df_adni["Date"],
    errors="coerce"
)

df_adni["Age_numeric"] = pd.to_numeric(
    df_adni["Age"],
    errors="coerce"
)

df_adni["DX_Binary_clean"] = (
    df_adni["DX_Binary"]
    .astype("string")
    .str.strip()
)

# Treat common missing-value strings as missing
df_adni.loc[
    df_adni["DX_Binary_clean"].isin(
        ["", "NA", "NaN", "nan", "None", "null", "<NA>"]
    ),
    "DX_Binary_clean"
] = pd.NA

# Keep only rows with an available diagnosis
df_with_dx = df_adni.loc[
    df_adni["DX_Binary_clean"].notna()
].copy()

# Define ordering:
#   1. Earliest valid Date
#   2. If Date is unavailable, youngest Age
#
# Rows with dates are prioritized over rows without dates.
df_with_dx["date_missing"] = df_with_dx["Date_parsed"].isna()

df_with_dx = df_with_dx.sort_values(
    by=[
        "Study",
        "PTID",
        "date_missing",
        "Date_parsed",
        "Age_numeric"
    ],
    ascending=[
        True,
        True,
        True,
        True,
        True
    ],
    na_position="last"
)

# One baseline diagnosis row per participant within each study
df_baseline = (
    df_with_dx
    .groupby(
        ["Study", "PTID"],
        as_index=False,
        sort=False
    )
    .first()
)

# Restore the diagnosis column name for downstream analyses
df_baseline["DX_Binary"] = df_baseline["DX_Binary_clean"]

# Count unique baseline participants by study and diagnosis
baseline_dx_counts = (
    df_baseline
    .groupby(
        ["Study", "DX_Binary"],
        dropna=False
    )["PTID"]
    .nunique()
    .reset_index(name="n_unique_participants")
    .sort_values(
        ["Study", "DX_Binary"],
        na_position="last"
    )
)

print(baseline_dx_counts.to_string(index=False))


print('Stop...')