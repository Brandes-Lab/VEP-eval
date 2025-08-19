import pandas as pd
import re
import mygene
from tqdm import tqdm

# ============================================================
# 1. Load benchmark data
# ============================================================
benchmark_df = pd.read_csv("./benchmark/all_scores_am.csv")

# ============================================================
# 2. Map NM_id -> UniProt ID (using MyGene)
# ============================================================
def map_nm_to_uniprot(nm_ids):
    """Map RefSeq NM IDs to UniProt IDs via MyGene."""
    mg = mygene.MyGeneInfo()
    results = mg.querymany(
        nm_ids, scopes="refseq", fields="uniprot", species="human"
    )

    nm_to_uniprot = {}
    for res in results:
        refseq = res.get("query")
        uniprot = res.get("uniprot", {})
        uniprot_id = None

        if isinstance(uniprot, dict):
            # Prefer Swiss-Prot, else TrEMBL
            ids = uniprot.get("Swiss-Prot") or uniprot.get("TrEMBL")
            uniprot_id = ids[0] if isinstance(ids, list) else ids
        elif isinstance(uniprot, list):
            uniprot_id = uniprot[0]
        else:
            uniprot_id = uniprot

        nm_to_uniprot[refseq] = uniprot_id

    return nm_to_uniprot


nm_ids = benchmark_df["ClinVarName_refseq_ids"].dropna().unique().tolist()
nm_to_uniprot = map_nm_to_uniprot(nm_ids)

benchmark_df["uniprot_id"] = benchmark_df["ClinVarName_refseq_ids"].map(nm_to_uniprot)
needed_uniprots = set(benchmark_df["uniprot_id"].dropna())

print(f"✅ Finished mapping: {len(needed_uniprots)} UniProt IDs found.")

# ============================================================
# 3. Helper for parsing protein_variant
# ============================================================
AA_VOCAB = [
    "K","R","H","E","D","N","Q","T","S","C",
    "G","A","V","L","I","M","P","Y","F","W"
]

def parse_mut(protein_variant):
    """
    Parse protein variant like 'A123T' or 'W50*'.
    Returns (ref, pos, alt).
    """
    match = re.match(r"([A-Z])(\d+)([A-Z*])", str(protein_variant))
    if match:
        ref, pos, alt = match.groups()
        return ref, int(pos), alt
    return None, None, None

# ============================================================
# 4. Process AlphaMissense (in chunks)
# ============================================================
am_path = "AlphaMissense_hg38.tsv"  # **replace with actual path
protein_variants = {}  # dict: uniprot_id -> DataFrame of variants

for chunk in tqdm(
    pd.read_csv(am_path, sep="\t", skiprows=3, chunksize=500_000),
    desc="Reading AlphaMissense"
):
    chunk = chunk.rename(columns={"am_pathogenicity": "am_score"})
    if "uniprot_id" not in chunk.columns:
        continue

    # Keep only proteins we care about
    chunk = chunk[chunk["uniprot_id"].isin(needed_uniprots)]
    if chunk.empty:
        continue

    # Parse ref, pos, alt from protein_variant
    chunk[["ref", "pos", "alt"]] = chunk["protein_variant"].apply(
        lambda x: pd.Series(parse_mut(x))
    )

    # Group by UniProt and append
    for uni, subdf in chunk.groupby("uniprot_id"):
        if uni not in protein_variants:
            protein_variants[uni] = subdf
        else:
            protein_variants[uni] = pd.concat(
                [protein_variants[uni], subdf], ignore_index=True
            )

print(f"✅ Collected variants for {len(protein_variants)} proteins.")

# ============================================================
# 5. Fill stop_gain AlphaMissense scores
# ============================================================
stop_mask = benchmark_df["ClinVarName_AAALT"] == "*"
stop_indices = benchmark_df[stop_mask].index

counter = 0
log = {}

for idx in tqdm(stop_indices, desc="Processing stop_gains"):
    row = benchmark_df.loc[idx]
    uni = row["uniprot_id"]
    stop_pos = row["ClinVarName_AAPOS"]

    if uni not in protein_variants:
        continue

    after = protein_variants[uni][protein_variants[uni]["pos"] >= stop_pos]
    if after.empty:
        continue

    # Check full coverage of AA substitutions (optional)
    full_coverage = True
    for p in sorted(after["pos"].unique()):
        subdf = after.loc[after["pos"] == p]
        missing_aas = set(AA_VOCAB) - set(subdf["alt"])
        if missing_aas:
            wt_ref = subdf["ref"].iloc[0]
            log[uni] = f"{wt_ref}{p}{''.join(missing_aas)}"
            full_coverage = False
            break

    # Worst-case score (maximum)
    worst_score = after["am_score"].max()
    benchmark_df.at[idx, "am_score"] = worst_score
    counter += 1

# ============================================================
# 6. Report & Save
# ============================================================
print(f"✅ Number of stop_gain am_scores updated: {counter}")
print(f"🔍 First 5 log entries: {dict(list(log.items())[:5])}")

print(f"benchmark len = {len(benchmark_df)}")
benchmark_df.to_csv("./benchmark/all_scores_stop_gain_am.csv", index=False)
print("alphamissense_merge_stop_gain.py done.")
