import numpy as np
import pandas as pd

am_path = 'AlphaMissense_hg38.tsv' # **replace with actual path
am_scores = (
    pd.read_csv(am_path, sep='\t', skiprows=3)
    .rename(columns={'am_pathogenicity': 'am_score'})
)
am_scores['#CHROM'] = am_scores['#CHROM'].astype(str).str.replace('chr', '')

merge_cols = ['#CHROM', 'POS', 'REF', 'ALT']
am_var = am_scores.loc[:, merge_cols + ['am_score']]

# --- Step 1: Deduplicate am_scores before merging ---
# Drop exact duplicates
am_var = am_var.drop_duplicates()

# Find inconsistent scores
dupes = am_var[am_var.duplicated(subset=merge_cols, keep=False)]

def inconsistent(group):
    return group['am_score'].nunique(dropna=False) > 1

bad_keys = (
    dupes.groupby(merge_cols)
    .filter(inconsistent)
    [merge_cols]
    .drop_duplicates()
)

if not bad_keys.empty:
    print(f"❌ Found {len(bad_keys)} inconsistent keys. Setting am_score=NaN")
    mask = am_var[merge_cols].apply(tuple, axis=1).isin(
        bad_keys.apply(tuple, axis=1)
    )
    am_var.loc[mask, 'am_score'] = np.nan
else:
    print("✅ No inconsistent am_scores")

# After this, each key maps to ≤1 am_score
am_var = am_var.drop_duplicates(subset=merge_cols, keep='first')


all_scores = pd.read_csv("./benchmark/all_scores_final.csv")
all_scores = all_scores.drop(columns=['am_score'], errors='ignore')

new_var = pd.merge(
    all_scores,
    am_var,
    on=merge_cols,
    how='left'  # keeps only rows in all_scores
)

print(f"benchmark len before = {len(all_scores)}")
print(f"benchmark len after  = {len(new_var)}")

new_var.to_csv('./benchmark/all_scores_am.csv', index=False)
print('alphamissense_merge.py done ✅')
