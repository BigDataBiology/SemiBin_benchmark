# SemiBin benchmark

CAMI_I.py: benchmark code for CAMI I datasets

CAMI_II.py: benchmark code for CAMI II datasets

Real.py: benchmark code for real datasets

### Evaluation

For CAMI I and CAMI II datasets, we used AMBER to evaluate the results.

CAMI I

```bash
python amber.py -g gsa_mapping.binning \
-l "Method" \
-p 1 \
-r unique_common.tsv \
-k "circular element" \
Method/bins.tsv \
-o output_dir/
```

CAMI II(for every sample) 

```bash
python amber.py -g gsa_mapping.binning \
-l "Method" \
Method/bins.tsv \
-o output_dir/
```

For real datasets, we used CheckM and GUNC to evaluate the binning performance.

CheckM
```bash
checkm lineage_wf -x  fa -t 48 SemiBin/output_recluster_bins SemiBin/checkm_output > SemiBin/checkm_result.txt
```

GUNC
```bash
mkdir GUNC
gunc run -d SemiBin/output_recluster_bins -o GUNC -r ~/gunc_path/gunc_db_2.0.4.dmnd -t 32
```

### Analysis

We used GTDB-TK to annotate the high-quality bins.

```bash
gtdbtk classify_wf --genome_dir high_quality_bins --out_dir gtdbtk_output -x 
fa --cpus 32
```

We used Mash to identify the overlap between SemiBin and Metabat2.
For two bins from SemiBin and Metabat2, we used
```bash
mash sketch SemiBin.fa
mash sketch Metabat2.fa
mash dist SemiBin.fa.msh Metabat2.fa.msh
```
