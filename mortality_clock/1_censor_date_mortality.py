import pandas as pd

death_xlsx = "/cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx"

df = pd.read_excel(death_xlsx)
death_date = pd.to_datetime(df["40000-0.0"], errors="coerce")

print("Number of non-missing death dates:", death_date.notna().sum())
print("Minimum death date:", death_date.min())
print("Maximum death date:", death_date.max())

print("\nDeaths by year:")
print(death_date.dt.year.value_counts().sort_index())

print("\nDeaths by month near the end:")
print(
    death_date.dropna()
    .dt.to_period("M")
    .value_counts()
    .sort_index()
    .tail(24)
)