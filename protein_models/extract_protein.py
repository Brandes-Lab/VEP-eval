import pandas as pd
import subprocess
from tqdm import tqdm
import glob
from Bio import Entrez, SeqIO
import xml.etree.ElementTree as ET
import time

# Set your email for NCBI API (required)
Entrez.email = "anonymous.user@example.com"

gene2acc_path = './protein_dataset/gene2accession' # **replace with actual path

# 1. Load your CSV and get unique NM IDs
df = pd.read_csv('./benchmark/ClinVarBenchmark_PB_202504.csv')
target_nm_ids = set(df['ClinVarName_refseq_ids'].dropna().unique())  # Only NM IDs we care about
print(f'Total targets = {len(target_nm_ids)}')

# 2. Build custom NM→NP mapping (first pass: gene2accession)
print("Building targeted NM→NP mapping (gene2accession)...")
nm_to_np = {}
with open(gene2acc_path) as f:
    total_lines = sum(1 for _ in f)

with open(gene2acc_path) as f:
    for line in tqdm(f, total=total_lines, desc="Scanning gene2accession"):
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) > 5:
            nm_id = parts[3]  # NM_XXXX.X
            np_id = parts[5]  # NP_XXXX.X
            if nm_id in target_nm_ids and np_id.startswith('NP_'):
                nm_to_np[nm_id] = np_id

# 3. Find missing NM IDs and fetch NP IDs via NCBI API
missing_nm_ids = target_nm_ids - set(nm_to_np.keys())
print(f"Found {len(missing_nm_ids)} unmapped NM IDs")

'''
def fetch_np_id(nm_id):
    try:
        handle = Entrez.esummary(db="nuccore", id=nm_id)
        record = Entrez.read(handle)
        handle.close()

        if record[0].get("Status") == "replaced":
            nm_id = record[0]["ReplacedBy"]

        # Always try linking to protein regardless of status
        handle = Entrez.elink(dbfrom="nuccore", db="protein", id=nm_id, retmode='xml')
        xml_data = handle.read()
        handle.close()

        root = ET.fromstring(xml_data)
        for link_set in root.findall('.//LinkSet'):
            for link_db in link_set.findall('.//LinkSetDb'):
                for link in link_db.findall('.//Link'):
                    protein_id = link.find('.//Id')
                    if protein_id is not None:
                        # Fetch protein accession
                        handle = Entrez.efetch(
                            db="protein",
                            id=protein_id.text,
                            rettype="acc",
                            retmode="text"
                        )
                        np_id = handle.read().strip()
                        handle.close()
                        if np_id.startswith("NP_"):
                            return np_id
    except Exception as e:
        print(f"Error fetching for {nm_id}: {e}")

    return None


# Batch-process missing NM IDs (with delay to avoid API rate limits)
for nm_id in tqdm(missing_nm_ids, desc="Querying NCBI for missing mappings"):
    np_id = fetch_np_id(nm_id)
    if np_id:
        nm_to_np[nm_id] = np_id
    time.sleep(0.4)  # NCBI limits to 3 requests/second

missing_nm_ids = target_nm_ids - set(nm_to_np.keys())
print(f"Found {len(missing_nm_ids)} unmapped NM IDs")

# read missing nm ids locally
local_nm = []
local_np = []
with open('missing_NM_id.txt', 'r+') as f:
    local_nm = [line.strip() for line in f]

assert set(missing_nm_ids) == set(local_nm)

with open('missing_NP_id.txt', 'r+') as f:
    local_np = [line.strip() for line in f]

local_nm_to_np = dict(zip(local_nm, local_np))
nm_to_np.update(local_nm_to_np)

print(f"current NM dict = {len(nm_to_np)} / {len(target_nm_ids)}")
'''


# 4. Fetch sequences for all NP IDs (from local FASTA or NCBI)
print("Fetching protein sequences...")

def get_sequence(np_id, fasta_files):
    """Search local FASTA files for NP sequence."""
    for fasta_file in fasta_files:
        with open(fasta_file) as f:
            record_id = None
            sequence = []
            for line in f:
                if line.startswith('>'):
                    if record_id == np_id:
                        return ''.join(sequence)
                    record_id = line[1:].split()[0]  # Get first word after >
                    sequence = []
                elif record_id == np_id:
                    sequence.append(line.strip())
            if record_id == np_id:
                return ''.join(sequence)
    return None

fasta_files = glob.glob('./protein_dataset/human.*.protein.faa')
nm_to_sequence = {}

# Try local FASTA files first
for nm_id, np_id in tqdm(nm_to_np.items(), desc="Processing sequences"):
    sequence = get_sequence(np_id, fasta_files)
    if sequence:
        nm_to_sequence[nm_id] = sequence
    else:
        missing_nm_ids.add(nm_id)

print(f"Found {len(missing_nm_ids)} NMs with no sequence.")


def fetch_protein_sequence_from_refseq(transcript_id):
    if transcript_id == 'NM_000557.1': # special case
        transcript_id = 'NM_000557.5'

    handle = Entrez.efetch(db="nucleotide", id=transcript_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    for feature in record.features:
        if feature.type == "CDS" and "translation" in feature.qualifiers:
            protein_id = feature.qualifiers.get("protein_id", ["N/A"])[0]
            amino_acid_seq = feature.qualifiers["translation"][0]
            return protein_id, amino_acid_seq

    return None, None  

ext_nm_to_sequence = {}
for nm_id in tqdm(missing_nm_ids):
    protein_id, aa_seq = fetch_protein_sequence_from_refseq(nm_id)
    if aa_seq is not None:
        ext_nm_to_sequence[nm_id] = aa_seq

nm_to_sequence.update(ext_nm_to_sequence)

print(f"current NP seq = {len(nm_to_sequence)} / {len(target_nm_ids)}")
assert len(nm_to_sequence) == len(target_nm_ids)


# 5. Write FASTA (all successful mappings)
print("Writing FASTA output...")
with open('./protein_dataset/Preprocessed_all_prot.fasta', 'w') as f:
    for nm_id, seq in tqdm(nm_to_sequence.items()):
        f.write(f">{nm_id}\n{seq}\n")

print(f"Successfully mapped {len(nm_to_sequence)}/{len(target_nm_ids)} NM IDs")