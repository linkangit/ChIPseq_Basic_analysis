# Simple ChIP-seq Analysis Pipeline (Until peak calling. Stay tuned for downstream analysis)

A clean, easy-to-use Snakemake pipeline for ChIP-seq data analysis. This pipeline takes raw sequencing reads and produces high-quality peak calls, handling all the standard processing steps automatically.

## What This Pipeline Does

**ChIP-seq (Chromatin Immunoprecipitation sequencing)** identifies where proteins bind to DNA genome-wide. This pipeline processes your raw sequencing data through these essential steps:

1. **Quality Assessment** - Checks if your sequencing data is good quality
2. **Read Alignment** - Maps your reads to the reference genome
3. **Data Processing** - Sorts and indexes alignment files for downstream analysis
4. **Peak Calling** - Identifies regions where your protein of interest binds to DNA
5. **Quality Control** - Provides statistics about your experiment's success

## Before You Start

### What You Need

**Data Requirements:**
- Raw sequencing files (`.fastq.gz` format) for your ChIP samples
- At least one input/control sample (essential for peak calling)
- A reference genome index (we'll help you get this)

**Computing Requirements:**
- Linux/Mac computer or cluster access
- At least 8GB RAM (16GB+ recommended)
- ~50GB free disk space per sample

## Detailed Setup Instructions

### Step 1: Install Required Software

**Option A: Using Conda/Mamba (Recommended)**
```bash
# Install conda if you don't have it
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Create environment for ChIP-seq analysis
conda create -n chipseq -c bioconda -c conda-forge snakemake fastqc bowtie2 samtools macs2
conda activate chipseq
```

**Option B: Manual Installation**
If conda isn't available, install each tool separately:
- [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Bowtie2](https://github.com/BenLangmead/bowtie2)
- [SAMtools](http://www.htslib.org/)
- [MACS2](https://github.com/macs3-project/MACS)

### Step 2: Get a Reference Genome

You need a bowtie2 genome index. Here are common options:

**For Human (hg38):**
```bash
mkdir -p genome
cd genome
wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
unzip GRCh38_noalt_as.zip
cd ..
```

**For Mouse (mm10):**
```bash
mkdir -p genome
cd genome  
wget https://genome-idx.s3.amazonaws.com/bt/mm10.zip
unzip mm10.zip
cd ..
```

**Or build your own:**
```bash
# Download genome FASTA
wget [your_genome.fa]
# Build index
bowtie2-build your_genome.fa genome_index
```

### Step 3: Organize Your Project

Create this exact directory structure:

```
my_chipseq_project/
├── Snakefile                    # The pipeline code
├── config.yaml                  # Your settings
├── data/                        # Put your sequencing files here
│   ├── input_sample.fastq.gz    # Your input/control sample
│   ├── H3K4me3_rep1.fastq.gz    # Your ChIP samples
│   ├── H3K4me3_rep2.fastq.gz
│   └── H3K27ac_rep1.fastq.gz
└── genome/                      # Your genome index files
    ├── genome_index.1.bt2
    ├── genome_index.2.bt2
    └── ... (other index files)
```

**Important File Naming Rules:**
- Use `.fastq.gz` extension (compressed FASTQ format)
- Avoid spaces and special characters in filenames
- Use consistent, descriptive names
- One input/control sample is required

### Step 4: Configure Your Analysis

Edit the `config.yaml` file with your specific settings:

```yaml
# Path to your genome index (without file extensions)
genome_index: "genome/GRCh38_noalt_as"  # or "genome/mm10" for mouse

# Genome size for MACS2 peak calling
genome_size: "hs"  # Options: "hs"=human, "mm"=mouse, or specific number like "2.7e9"

# Define your samples (must match your filename prefixes)
samples:
  input_sample: "input"           # This MUST be your control sample
  H3K4me3_rep1: "H3K4me3_r1"    # Change these to match your files
  H3K4me3_rep2: "H3K4me3_r2"
  H3K27ac_rep1: "H3K27ac_r1"
```

**Critical Notes:**
- The `input_sample` key must stay exactly as written - this is your control
- Sample names on the left should match your `.fastq.gz` filenames
- Values on the right are just labels for the analysis

### Step 5: Test Your Setup

Before running the full analysis, test everything:

```bash
# Check if all files are found and rules make sense
snakemake --dry-run --printshellcmds

# This should show you exactly what will run without actually doing it
```

If you see errors:
- **"No such file"**: Check your file paths and names
- **"Missing input files"**: Ensure your `.fastq.gz` files are in the `data/` directory
- **Config errors**: Double-check your `config.yaml` syntax

## Running the Pipeline

### Quick Start
```bash
# Run with 8 CPU cores (adjust based on your computer)
snakemake --cores 8

# With more detailed output
snakemake --cores 8 --printshellcmds
```

### For Computer Clusters

**SLURM cluster:**
```bash
snakemake --cores 100 --cluster "sbatch --time=2:00:00 --mem=8G --cpus-per-task=1"
```

**PBS cluster:**
```bash
snakemake --cores 50 --cluster "qsub -l walltime=2:00:00,mem=8gb"
```

### Monitoring Progress

**Check what's running:**
```bash
# See rule progress
snakemake --cores 8 --printshellcmds

# Create a visual workflow diagram
snakemake --dag | dot -Tpng > workflow.png
```

**Typical runtime expectations:**
- FastQC: 5-10 minutes per sample
- Alignment: 30-60 minutes per sample
- Peak calling: 10-30 minutes per comparison
- **Total time**: 2-4 hours for 4 samples

## Understanding Your Results

The pipeline creates organized output directories:

```
results/
├── qc/                          # Quality control reports
│   ├── input_sample_fastqc.html     # Open these in web browser
│   ├── H3K4me3_rep1_fastqc.html
│   ├── input_sample.flagstat        # Alignment statistics
│   └── H3K4me3_rep1.flagstat
├── bam/                         # Aligned sequencing data
│   ├── input_sample.bam             # Sorted alignment files
│   ├── input_sample.bam.bai         # Index files (needed for viewing)
│   └── H3K4me3_rep1.bam
└── peaks/                       # Peak calling results
    ├── H3K4me3_rep1_peaks.narrowPeak   # Main peak file (BED format)
    ├── H3K4me3_rep1_summits.bed        # Peak summit coordinates
    └── H3K27ac_rep1_peaks.narrowPeak
```

### Key Output Files Explained

**Quality Control (`results/qc/`):**
- `*_fastqc.html`: Open in web browser to see read quality, contamination, etc.
- `*.flagstat`: Text files with alignment statistics (% mapped reads, duplicates)

**Peak Files (`results/peaks/`):**
- `*_peaks.narrowPeak`: Main results! These are your protein binding sites
  - Columns: chromosome, start, end, peak_name, score, strand, fold_enrichment, -log10(pvalue), -log10(qvalue), summit_offset
- `*_summits.bed`: Exact peak summit coordinates (use for motif analysis)

### Evaluating Your Results

**Good Quality Indicators:**
- FastQC shows >90% high-quality bases
- >70% of reads align to genome (check `.flagstat` files)
- Peak files contain hundreds to thousands of peaks
- Clear enrichment in peak regions vs background

**Warning Signs:**
- <50% alignment rate (possible contamination or wrong genome)
- Very few peaks (<100) might indicate failed ChIP
- FastQC shows adapter contamination or low quality

## Advanced Usage

### Adding More Samples

Simply add them to your `config.yaml`:

```yaml
samples:
  input_sample: "input"
  H3K4me3_rep1: "H3K4me3_r1"
  H3K4me3_rep2: "H3K4me3_r2"      # Biological replicate
  H3K4me3_treat: "H3K4me3_treat"  # Treatment condition
  H3K27ac_ctrl: "H3K27ac_ctrl"    # Different histone mark
  CTCF_sample: "CTCF"             # Transcription factor
```

Each ChIP sample will be compared against your input sample automatically.

### Customizing Parameters

**Change peak calling stringency:**
Edit the MACS2 command in the Snakefile:
```bash
# More stringent (fewer, higher confidence peaks)
-q 0.001

# Less stringent (more peaks, some false positives)  
-q 0.05
```

**Adjust for different read lengths:**
MACS2 auto-detects, but you can specify:
```bash
--extsize 200  # For 50bp reads, extend to ~200bp fragments
```

### Working with Paired-End Data

If you have paired-end sequencing data:

1. **Merge R1 and R2 files** or modify the alignment rule:
```bash
# In the align rule, change:
bowtie2 -x {params.genome} -1 data/{sample}_R1.fastq.gz -2 data/{sample}_R2.fastq.gz
```

2. **Use BAMPE format** in MACS2:
```bash
# Change -f BAM to:
-f BAMPE
```

## Troubleshooting Common Issues

### "Command not found" errors
**Problem**: Software not installed or not in PATH  
**Solution**: 
```bash
# Check if tools are available
which snakemake bowtie2 samtools macs2

# If missing, reinstall conda environment
conda activate chipseq
conda install -c bioconda [missing_tool]
```

### "No such file or directory"
**Problem**: Input files not found  
**Solution**: 
```bash
# Check your files exist and are named correctly
ls -la data/
# Files must be named exactly as specified in config.yaml
```

### "Bowtie2 index not found"
**Problem**: Genome index path is wrong  
**Solution**:
```bash
# Check your index files exist
ls genome/
# Update config.yaml with correct path (no file extensions)
genome_index: "genome/your_actual_index_name"
```

### MACS2 fails with "Too few peaks"
**Problem**: Poor ChIP efficiency or very stringent parameters  
**Solutions**:
- Check if your input files are swapped
- Try less stringent q-value: `-q 0.05` or `-q 0.1`
- Verify your ChIP experiment worked (check raw data quality)

### Out of memory errors
**Problem**: Not enough RAM for large genomes/files  
**Solutions**:
- Reduce number of parallel jobs: `--cores 4` instead of `--cores 8`
- Use cluster mode to access more memory
- Split large FASTQ files if possible

### Very low alignment rates (<30%)
**Problem**: Wrong genome, contamination, or poor quality  
**Solutions**:
- Verify you're using the correct genome build (hg38 vs hg19, etc.)
- Check FastQC for contamination (high duplication, adapter content)
- Try aligning a subset of reads to test different genomes

## Next Steps After Peak Calling

This pipeline gives you the foundation - high-quality peak calls. Common downstream analyses include:

1. **Peak Annotation**: What genes are near your peaks?
2. **Motif Analysis**: What DNA sequences are enriched in peaks?
3. **Differential Binding**: How do peaks change between conditions?
4. **Functional Enrichment**: What pathways are associated with bound genes?
5. **Data Visualization**: Create genome browser tracks

Popular tools for these analyses:
- **ChIPseeker** (R package): Peak annotation
- **HOMER**: Motif discovery and annotation  
- **DiffBind** (R package): Differential binding analysis
- **GREAT**: Functional annotation of regulatory regions

## Getting Help

**First steps:**
1. Check the log files in your working directory for error messages
2. Run with `--dry-run` to test without executing
3. Verify your input files and config.yaml settings

**Still stuck?**
- Check the [Snakemake documentation](https://snakemake.readthedocs.io/)
- MACS2 help: `macs2 callpeak --help`
- Post on bioinformatics forums with specific error messages

**Want to extend this pipeline?**
The simple design makes it easy to add steps like:
- Read trimming (Trimmomatic)
- Duplicate removal (Picard)
- Multiple QC reports (MultiQC) 
- Peak visualization (deepTools)

Remember: This pipeline is designed to be a solid, understandable foundation. Feel free to modify it for your specific needs!
