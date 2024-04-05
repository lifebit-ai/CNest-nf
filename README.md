# Nextflow pipeline for copy number estimation and analysis using CNest

Primary tool used in this pipeline - https://github.com/tf2/CNest

## Run locally

This pipeline can run in two modes - 
1. Running End to end for all the steps
2. Running individual steps

### 1. Running End to end for all the steps

```bash
nextflow run main.nf \
    --project "1kgp_subset_test" \
    --design "testdata/1kg_CEU_7_samples_subset.csv" \
    --bedgz "testdata/hg38.1kb.chr-19-22.baits.bed.gz" \
    --gender "testdata/gender_classification.txt" \
    --cov "testdata/hg38.1kb.chr-19-22-mean_coverage.txt"
```

or run from the test profile 

```bash
nextflow run main.nf -profile test,standard
```

### 2. Running individual steps


```bash
# Step 1 : Make index
nextflow run main.nf \
    --step 1 \
    --project "1kgp_subset_test" \
    --bedgz "testdata/hg38.1kb.chr-19-22.baits.bed.gz"

# Step 2 : Bait read count
nextflow run main.nf \
    --step 2 \
    --project "1kgp_subset_test" \
    --design "testdata/1kg_CEU_7_samples_subset.csv" \
    --indexb "results/1kgp_subset_test/index.bed"

# Step 3 : Gender QC
nextflow run main.nf \
    --step 3 \
    --project "1kgp_subset_test" \
    --bindir "results/1kgp_subset_test/bin" \
    --index_tab "results/1kgp_subset_test/index_tab.txt"

# Step 4 : logR-ratio calculation
nextflow run main.nf \
    --step 4 \
    --batch_size 4 \
    --target_size 4 \
    --project "1kgp_subset_test" \
    --bindir "results/1kgp_subset_test/bin" \
    --index_tab "results/1kgp_subset_test/index_tab.txt"

# Step 5 : log2 rbin file conversion
nextflow run main.nf \
    --step 5 \
    --batch_size 4 \
    --target_size 4 \
    --project "1kgp_subset_test" \
    --bindir "results/1kgp_subset_test/bin" \
    --cordir "results/1kgp_subset_test/cor" \
    --gender "testdata/gender_classification.txt" \
    --index_tab "results/1kgp_subset_test/index_tab.txt"

# Step 6 : HMM call
nextflow run main.nf \
    --step 6 \
    --batch_size 4 \
    --target_size 4 \
    --project "1kgp_subset_test" \
    --rbindir "results/1kgp_subset_test/rbin" \
    --cordir "results/1kgp_subset_test/cor" \
    --gender "testdata/gender_classification.txt" \
    --index_tab "results/1kgp_subset_test/index_tab.txt" \
    --cov "testdata/hg38.1kb.chr-19-22-mean_coverage.txt"
```
