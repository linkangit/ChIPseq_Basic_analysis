# Simple ChIP-seq Analysis Pipeline

A clean, easy-to-use Snakemake pipeline for ChIP-seq data analysis.

## Quick Setup

### 1. Required Software
Install these tools (or use conda/mamba):
```bash
conda install -c bioconda snakemake fastqc bowtie2 samtools macs2
```

### 2. Directory Structure
Organize your files like this:
```
project/
├── Snakefile
├── config.yaml
└── data/
    ├── input_sample.fastq.gz
    ├── chip_sample1.fastq.gz
    └── chip_sample2.fastq.gz
```

### 3. Configure Your Analysis
Edit `config.yaml`:
- Set your genome index path
- Define your sample names
- Adjust genome size if needed

### 4. Run the Pipeline
```bash
# Dry run to check everything
snakemake -n

# Run with 8 cores
snakemake -j 8

# Run on cluster (SLURM example)
snakemake -j 100 --cluster "sbatch -n 1 -t 60"
```

## Results

The pipeline creates:
- `results/qc/` - Quality control reports
- `results/bam/` - Aligned and sorted BAM files  
- `results/peaks/` - Peak calling results

## What It Does

1. **Quality Control**: FastQC on raw reads
2. **Alignment**: Bowtie2 alignment to genome
3. **Processing**: Sort and index BAM files
4. **Statistics**: Generate alignment stats
5. **Peak Calling**: MACS2 peak calling (ChIP vs Input)

## Customization

### Add More Samples
Just add them to the `samples` section in `config.yaml`:
```yaml
samples:
  input_sample: "input"
  chip_sample1: "H3K4me3"
  chip_sample2: "H3K27ac" 
  chip_sample3: "H3K4me1"  # Add this line
```

### Change Parameters
Modify the shell commands in the Snakefile. For example, to change MACS2 q-value:
```bash
-q 0.01  # Change to -q 0.05 for less stringent peaks
```

### Different File Names
If your files have different names, either:
1. Rename them to match the sample names, OR
2. Modify the input paths in the rules

## Common Issues

**"No such file"**: Check that your fastq.gz files are in the `data/` directory with the correct names.

**Bowtie2 index error**: Make sure the `genome_index` path in config.yaml points to your actual genome index files.

**MACS2 fails**: Ensure you have both ChIP and input samples, and they're properly aligned.

## Need Help?

This pipeline covers the basics. For more advanced features:
- Add quality trimming (Trimmomatic)
- Include duplicate removal (Picard)
- Add more QC steps (MultiQC)
- Implement peak annotation

The design is intentionally simple - modify as needed for your specific requirements!
