python - <<'PY'
import pandas as pd

f = "/Users/hao/Dropbox/NHANES/output_dir/model2_nondisease_mortality_epoch/nhanes_model2_analysis_table.tsv.gz"
df = pd.read_csv(f, sep="\t")

print("Rows:", len(df))
print("Unique SEQN:", df["SEQN"].nunique())
print("Duplicated SEQN rows:", df["SEQN"].duplicated().sum())

dup = df[df["SEQN"].duplicated(keep=False)][["SEQN", "cycle"]].sort_values("SEQN")
print(dup.head(20))
PY