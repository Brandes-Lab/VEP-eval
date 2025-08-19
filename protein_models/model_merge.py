import numpy as np
import pandas as pd

DNA_scores = pd.read_csv("./benchmark/all_scores_v6.csv")
DNA_scores.drop(columns = ['ESM1b_score', 'ESM1v-1_score', 'ESM2_score'], inplace=True)
print(DNA_scores.columns.values)

protein_scores = pd.read_csv("./protein_dataset/benchmark_scores.csv")
protein_scores = protein_scores[['#CHROM', 'POS', 'ID', 'ESM1b_score', 'ESM1v-1_score', 'ESM2_score']]

# merge_cols = protein_scores.columns[[0, 1, 2] + list(range(6, protein_scores.shape[1]))].tolist()
merge_cols = protein_scores.columns[[0, 1, 2]].tolist()

all_scores = pd.merge(DNA_scores,
    protein_scores,
    on = merge_cols,
    how='inner'
)

print(all_scores.shape[0])
print('Any NA?')
print(all_scores.isna().any())

all_scores.to_csv('./benchmark/all_scores_final.csv', index = False)

print('model_merge.py done')

