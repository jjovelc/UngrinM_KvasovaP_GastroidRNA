# Data Processing

## 1. Quality inspection

```shell
  fastqc *
  conda activate multiqc
  multiqc .
```

## 2. Formatting Mus musculus transcriptome

```shell
  wget https://ftp.ensembl.org/pub/current/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
  unzip Mus_musculus.GRCm39.cdna.all.fa.gz
  sbatch format_index_salmon.slurm
```

## 3. Quantify transcripts abundance with Salmon

```shell
  sbatch salmon_quant.slurm
  bash extract_salmon-quant_results.sh
```

### 4. Generate dotplots from normalized data (TPMs)

```shell
  conda activate r_3.6.3
  Rscript create_dot_plots.R
```

### 5. Conduct differential expression and gene ontology analyses

Script extract_salmon-quant_results.sh produces tables `allSamples_counts.tsv` and `allSamples_tpms.tsv`. The former one was used for postprocessing with script DESeq2.R.
