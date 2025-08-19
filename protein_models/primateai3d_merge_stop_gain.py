import pandas as pd
from tqdm import tqdm
import mygene

# --------------------------
# 1. Load benchmark data
# --------------------------
benchmark_df = pd.read_csv("./benchmark/all_scores_pai.csv")
# --------------------------
# 2. Map RefSeq -> ENST using mygene
# --------------------------
mg = mygene.MyGeneInfo()
refseq_ids = benchmark_df["ClinVarName_refseq_ids"].dropna().unique().tolist()

results = mg.querymany(
    refseq_ids,
    scopes="refseq",
    fields="ensembl.transcript",
    species="human"
)

refseq_to_enst = {}
for res in results:
    refseq = res['query']
    enst = None
    if 'ensembl' in res:
        if isinstance(res['ensembl'], dict):
            enst = res['ensembl'].get('transcript')
        elif isinstance(res['ensembl'], list):
            # Take the first transcript from the first dict
            first_entry = res['ensembl'][0]
            enst = first_entry.get('transcript')
    refseq_to_enst[refseq] = enst

print(f"len of enst dict = {len(refseq_to_enst)}")
benchmark_df['enst_id'] = benchmark_df['ClinVarName_refseq_ids'].map(refseq_to_enst)

# --------------------------
# 3. Load PrimateAI3D scores
# --------------------------
pai3d_path = "PrimateAI-3D.hg38.txt"  # **replace with actual path
pai_cols = [
    '#CHROM', 'POS', 'REF', 'ALT', 'gene_name',
    'change_position_1based', 'ref_aa', 'alt_aa',
    'score_PAI3D', 'percentile_PAI3D'
]

pai_scores = pd.read_csv(pai3d_path, sep='\t', names=pai_cols, header=0)

# Rename for convenience
pai_scores = pai_scores.rename(columns={
    'gene_name': 'enst_id',
    'score_PAI3D': 'pai_score',
    'change_position_1based': 'position',
    'ref_aa': 'ref',
    'alt_aa': 'alt'
})

# --------------------------
# 4. Separate stop-gain vs non-stop-gain
# --------------------------
stop_gain_mask = benchmark_df['ClinVarName_AAALT'] == "*"
stop_gain_df = benchmark_df[stop_gain_mask].copy()
non_stop_gain_df = benchmark_df[~stop_gain_mask].copy()

# --------------------------
# 5. Process stop-gain variants: worst score after stop codon
# --------------------------
AA_vocab = ['K','R','H','E','D','N','Q','T','S','C','G',
            'A','V','L','I','M','P','Y','F','W']

# Group PrimateAI3D by ENST for fast lookup
pai_dict = {enst: subdf for enst, subdf in pai_scores.groupby("enst_id")}

stop_gain_results = []
counter = 0

for idx, row in tqdm(stop_gain_df.iterrows(), desc="Processing stop-gains"):
    ensts = row['enst_id']
    stop_pos = row['ClinVarName_AAPOS']

    if not ensts:
        stop_gain_results.append(None)
        continue

    if not isinstance(ensts, list):
        ensts = [ensts]

    worst_scores = []
    for enst in ensts:
        if enst not in pai_dict:
            continue
        subdf = pai_dict[enst]
        after = subdf[subdf['position'] >= stop_pos]
        if not after.empty:
            worst_scores.append(after['pai_score'].min())

    if worst_scores:
        stop_gain_results.append(min(worst_scores))  # worst overall
        counter += 1
    else:
        stop_gain_results.append(None)

stop_gain_df['pai_score'] = stop_gain_results
print(f"Total stop_gain processed = {counter}")

# --------------------------
# 6. Merge non-stop-gain variants by enst_id
# --------------------------
pai_scores['enst_id_single'] = pai_scores['enst_id'].apply(
    lambda x: sorted(x)[0] if isinstance(x, list) else x
)
non_stop_gain_df['enst_id_single'] = non_stop_gain_df['enst_id'].apply(
    lambda x: sorted(x)[0] if isinstance(x, list) else x
)

# Merge on the single ENST
non_stop_gain_merged = non_stop_gain_df.merge(
    pai_scores[['enst_id_single', 'pai_score']].drop_duplicates(subset='enst_id_single'),
    on='enst_id_single',
    how='left'
).drop(columns='enst_id_single')

# --------------------------
# 7. Combine stop-gain and non-stop-gain
# --------------------------
final_df = pd.concat([stop_gain_df, non_stop_gain_merged], ignore_index=True)

# --------------------------
# 8. Save final result
# --------------------------

print(f"benchmark len = {len(final_df)}")
final_df.to_csv("./benchmark/all_scores_stop_gain_pai.csv", index=False)
print("primateai3d_merge_stop_gain.py done")
