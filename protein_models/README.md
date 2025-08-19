ESM Models:

0. Clone the repository from https://github.com/ntranoslab/esm-variants
1. `python extract_protein.py` to get protein_dataset/Preprocessed_all_prot.fasta
2. ```python ../esm-variants/esm_score_missense_mutations.py --input-fasta-file protein_dataset/Preprocessed_all_prot.fasta --output-csv-file protein_dataset/ESM1b_score.csv --model-name esm1b_t33_650M_UR50S
   python ../esm-variants/esm_score_missense_mutations.py --input-fasta-file protein_dataset/Preprocessed_all_prot.fasta --output-csv-file protein_dataset/ESM1v_score.csv --model-name esm1v_t33_650M_UR90S_1
   python ../esm-variants/esm_score_missense_mutations.py --input-fasta-file protein_dataset/Preprocessed_all_prot.fasta --output-csv-file protein_dataset/ESM2_score.csv --model-name esm2_t33_650M_UR50D``` to get separate scores for ESM1b, ESM1v-1 and ESM2
3. `python protein_merge.py` to merge 3 ESM model scores to ./protein_dataset/ESM_score.csv
4. `python protein_score.py` to match scores according to the benchmark, and also compute scores for stop_gain
5. `python model_merge.py` to merge the result to other models



AlphaMissense:
1. get precomputed AlphaMissense scores AlphaMissense_hg38.tsv from https://alphamissense.hegelab.org/
2. ```python alphamissense_merge.py
   python alphamissense_merge_stop_gain.py``` to match scores according to the benchmark, and also compute scores for stop_gain



PrimateAI3D:
1. get precomputed PrimateAI3D scores PrimateAI-3D.hg38.txt from https://primateai3d.basespace.illumina.com/
2. ```python primateai3d_merge.py
   python primateai3d_merge_stop_gain.py``` to match scores according to the benchmark, and also compute scores for stop_gain
