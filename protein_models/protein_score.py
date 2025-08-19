import numpy as np
import pandas as pd
from tqdm import tqdm

# Load data in chunks if needed
all_scores = pd.read_csv('./protein_dataset/ESM_score.csv')
benchmark = pd.read_csv('./benchmark/ClinVarBenchmark_PB_202504.csv')

# Preprocessing (unchanged)
benchmark = benchmark.rename(columns={'#CHROM': 'chrom', 'POS': 'pos', 'INFO': 'label'})
benchmark['chrom'] = benchmark['chrom'].astype(str).str.replace('chr', '')
benchmark['pos'] = benchmark['pos'].astype(int)
benchmark['ID'] = benchmark['ID'].astype(int)
# benchmark = benchmark[pd.to_numeric(benchmark['ClinVarName_AAPOS'], errors='coerce').notna()]
# benchmark['ClinVarName_AAPOS'] = benchmark['ClinVarName_AAPOS'].astype(int)

all_scores[['REF', 'POS', 'ALT']] = all_scores['mut_name'].str.extract(r'([A-Z*])(\d+)([A-Z*])')
all_scores['POS'] = all_scores['POS'].astype(int)

# Create merge keys and merge
benchmark['merge_key'] = (
    benchmark['ClinVarName_refseq_ids'].astype(str) + '_' +
    benchmark['ClinVarName_AAREF'].astype(str) +
    benchmark['ClinVarName_AAPOS'].apply(lambda x: str(int(x)) if pd.notna(x) else 'nan') +
    benchmark['ClinVarName_AAALT'].astype(str)
)

all_scores['merge_key'] = (
    all_scores['seq_id'].astype(str) + '_' +
    all_scores['REF'].astype(str) +
    all_scores['POS'].apply(lambda x: str(int(x)) if pd.notna(x) else 'nan') +
    all_scores['ALT'].astype(str)
)

models = ['ESM1b', 'ESM1v-1', 'ESM2']
result = benchmark.merge(
    all_scores[['merge_key'] + [f'{model}_score' for model in models]],
    on='merge_key',
    how='left'
).drop(columns = ['merge_key'])

# ========== PROPER DOWNSTREAM SCORING ==========
print("Finding true downstream worst scores...")

# Group variants by protein for efficient lookup
protein_variants = {
    seq_id: group[['POS'] + [f'{model}_score' for model in models]]
    for seq_id, group in all_scores.groupby('seq_id')
}

# Process stop codons in chunks
stop_mask = result['ClinVarName_AAALT'] == '*'
chunk_size = 5000  # Reduce if OOM occurs
n_chunks = (stop_mask.sum() // chunk_size) + 1

for i in tqdm(range(n_chunks), desc="Processing stop codons"):
    chunk = result[stop_mask].iloc[i*chunk_size : (i+1)*chunk_size]
    
    for idx, row in chunk.iterrows():
        refseq = row['ClinVarName_refseq_ids']
        stop_pos = row['ClinVarName_AAPOS']
        
        if refseq not in protein_variants:
            continue
            
        variants = protein_variants[refseq]
        downstream = variants[variants['POS'] >= stop_pos]
        
        if not downstream.empty:
            for model in models:
                worst_score = downstream[f'{model}_score'].min()
                result.at[idx, f'{model}_score'] = worst_score

result = result.rename(columns={'chrom' : '#CHROM', 'pos' : 'POS', 'label' : 'INFO'})
result = result.iloc[:, [0, 1, 2, -3, -2, -1] + list(range(3, result.shape[1] - 3))]

# Save results
result.to_csv(
    "./protein_dataset/benchmark_scores.csv", 
    index=False
)
print("Done")