ESM Models:

0\. Clone the repository from https://github.com/ntranoslab/esm-variants

1. `python extract\_protein.py` to get protein\_dataset/Preprocessed\_all\_prot.fasta

2. ```python ../esm-variants/esm\_score\_missense\_mutations.py --input-fasta-file protein\_dataset/Preprocessed\_all\_prot.fasta --output-csv-file protein\_dataset/ESM1b\_score.csv --model-name esm1b\_t33\_650M\_UR50S
   python ../esm-variants/esm\_score\_missense\_mutations.py --input-fasta-file protein\_dataset/Preprocessed\_all\_prot.fasta --output-csv-file protein\_dataset/ESM1v\_score.csv --model-name esm1v\_t33\_650M\_UR90S\_1
   python ../esm-variants/esm\_score\_missense\_mutations.py --input-fasta-file protein\_dataset/Preprocessed\_all\_prot.fasta --output-csv-file protein\_dataset/ESM2\_score.csv --model-name esm2\_t33\_650M\_UR50D``` to get separate scores for ESM1b, ESM1v-1 and ESM2
3. `python protein\_merge.py` to merge 3 ESM model scores to ./protein\_dataset/ESM\_score.csv

3. `python protein\_merge.py` to merge 3 ESM model scores to ./protein\_dataset/ESM\_score.csv
4. `python protein\_score.py` to match scores according to the benchmark, and also compute scores for stop\_gain
5. `python model\_merge.py` to merge the result to other models



AlphaMissense:

1. get precomputed AlphaMissense scores AlphaMissense\_hg38.tsv from https://alphamissense.hegelab.org/

```python alphamissense\_merge.py
   python alphamissense\_merge\_stop\_gain.py``` to match scores according to the benchmark, and also compute scores for stop\_gain



PrimateAI3D:

1. get precomputed PrimateAI3D scores PrimateAI-3D.hg38.txt from https://primateai3d.basespace.illumina.com/

2. ```python primateai3d\_merge.py
   python primateai3d\_merge\_stop\_gain.py``` to match scores according to the benchmark, and also compute scores for stop\_gain
