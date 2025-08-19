import numpy as np
import pandas as pd

pai_path = 'PrimateAI-3D.hg38.txt' # **replace with actual path
pai_scores = pd.read_csv(pai_path, sep='\t').rename(columns = {'chr' : '#CHROM', 'pos' : 'POS', 'non_flipped_ref' : 'REF', 'non_flipped_alt' : 'ALT', 'score_PAI3D' : 'pai_score'})
pai_scores['#CHROM'] = pai_scores['#CHROM'].astype(str).str.replace('chr', '')
print(pai_scores.columns.values)

pai_var = pai_scores.iloc[:, [0, 1, 2, 3, 8]]

all_scores = pd.read_csv("./benchmark/all_scores_stop_gain_am.csv")
all_scores.drop(columns = ['pai_score'], inplace=True)

merge_cols = pai_scores.columns[[0, 1, 2, 3]].tolist()


duplicate_keys = pai_scores.duplicated(subset=merge_cols, keep=False)
print(f"{duplicate_keys.sum()} duplicate key rows in pai_var")

pd.set_option('display.max_columns', None)
dupes = pai_scores[pai_scores.duplicated(subset=merge_cols, keep=False)]
print(dupes.sort_values(merge_cols).head(10))
pd.reset_option('display.max_columns')


new_var = pd.merge(all_scores,
    pai_var,
    on = merge_cols,
    how='left'
)

print(f"benchmark len = {len(new_var)}")
new_var.to_csv('./benchmark/all_scores_pai.csv', index = False)

print('primateai3d_merge.py done')
