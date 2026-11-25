# OrCA-seq Processing Workflow (Snakemake)

A Snakemake-based workflow for processing COI/28S/18S amplicon data locally!

## Version Info

- **Snakemake**: ≥7.0
- **Python**: 3.8+
- **Conda/Mamba**: Required
- **Tested on**:
   - MacBook Pro M3 (8 cores, 16GB RAM)
   - HPC interactive job (16 cores, 24GB RAM)

## Workflow Steps

1 & 2. **Barbell**: Reorient cDNA reads, quality trim, and remove indexing primers, plus extra filtering to remove reads with aberrant partial indexing primer segments
3. **Amplicon Sorter**: Cluster reads by sequence similarity and filter by size
4. **Primer Removal**: Strip 5'/3' primers from amplicon sequences
5. **Pybarrnap/Reorganize**: Extract rRNA genes or reorganize COIs
6. **BLAST**: BLAST search against NR v5 (25th May, 2025), either with blastn for rRNA barcodes or tblastx for COI barcodes

## Setup

### Prerequisites

- **macOS** with M3 or M-series chip
- **Conda/Mamba** installed (recommended: [Mambaforge](https://github.com/conda-forge/miniforge))
- **Snakemake** ≥7.0
- **Git** (to clone the OrCA-seq-processing repo)

### Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/srisarya/OrCA-seq-processing.git
   cd OrCA-seq-processing
   ```

2. **Create a Snakemake environment:**
   ```bash
   conda create -n snakemake-orca snakemake=7.32.4 -c bioconda -c conda-forge
   conda activate snakemake-orca
   ```

3. **Download amplicon_sorter.py (if not present):**
   ```bash
   curl -o scripts/auxiliary_code/amplicon_sorter.py \
     https://raw.githubusercontent.com/avierstr/amplicon_sorter/master/amplicon_sorter.py
   chmod +x scripts/auxiliary_code/amplicon_sorter.py
   ```

**Your directory structure should look like:**
   ```
   OrCA-seq-processing/
   ├── OrCAseq_processing.smk (the snakefile)
   ├── configs/
   │   └── config_for_dataset.yaml
   ├── envs/
   │   ├── barbell.yaml
   |   ├── barbell.post-deploy.sh # since barbell installation is more robust with rust through cargo
   │   ├── cutadapt.yaml
   │   ├── amplicon_sorter.yaml
   │   ├── blast.yaml
   │   └── pybarrnap.yaml
   ├── adapters_primers/
   ├── scripts/
   │   ├── auxiliary_code/
   │   │   ├── amplicon_sorter.py
   │   │   └── make_filters.py # generate barbell filters.txt from `barbell inspect` output.
   │   └── ...
   └── input_reads/        # Your FASTQ files here
   ```

## Usage

### Configuration Files
You will need a configuration file (`config.yaml`) in order to run the workflow.
This is so it knows:
   - What the dataset is and where it is
   - What dataset name to give
   - Whether you're using macOS or PC/windows
   - The adapters, configuration files and primers for demultiplexing and cleaning
   - Which size thresholds to cluster reads by
   - Which resources to provide (CPUs/threads)

In this study the dataset configs are in configs/ 

### Basic Workflow Execution

1. **Prepare input files:**
   - Place FASTQ files (`.fastq.gz`) in the `input_reads/` directory

2. **Dry-run the workflow with default config:**
   ```bash
   snakemake -s OrCAseq_processing.smk --use-conda --configfile configs/config_all_experimental_data.yaml --dry-run
   ```

3. **Run the protocol with 6 threads if the summary looks right**
   ```bash
   snakemake -s OrCAseq_processing.smk -c 6 --use-conda --configfile configs/config_all_experimental_data.yaml
   ```
   
- **Options:**
   - `-c 4`: Use 4 cores (adjust for your system; M3 can handle 4-8)
   - `--use-conda`: Automatically create/activate conda environments.
      - If you want to install tools in your PATH, then you can. See [Advanced Usage](#advanced-usage)
   - `--rerun-incomplete`: Re-run jobs if they fail partway through
   - once the workflow has completed, you can run it like so `snakemake -s OrCAseq_processing.smk --configfile configs/config_Lakesday1.yaml --report` to get a HTML report of runtimes and timestamps :)


## Output Structure (summary)

```
results_NTC_barbell_M13filt_3mm/                                                                                                                      
├── amplicon_sorted
│   └── <sample>
│       └── SP27_NNN_SP5_NNN
│           ├── COIs
│               ├── results.csv
│               ├── SP27_NNN_SP5_NNN_consensus_COIs.fasta
│               └── SP27_NNN_SP5_NNN_NTC_exp/
|               ...
│           └── rRNAs
│               ├── results.csv
│               ├── SP27_NNN_SP5_NNN_consensus_rRNAs.fasta
│               └── SP27_NNN_SP5_NNN_NTC_exp/
|               ...
├── barbell
│   └── <sample>
│       ├── annotate.log
│       ├── annotate.tsv
│       ├── filtered.tsv
│       ├── filters.txt
│       ├── m13_filter_summary.tsv
│       ├── patterns.txt
│       └── demuxed
|           ├── SP27_NNN_SP5_NNN_NTC_exp.fastq.gz
|           ...
├── BLAST
│   ├── <sample>
|       ├── SP27_NNN_SP5_NNN_<barcode>_blast_top5.tsv
|       ...
├── COI_gene
│   └── <sample>
│       ├── SP27_NNN_SP5_001_NNN_COI.fasta
├── primerless
│   └── <sample>
│       └── SP27_NNN_SP5_NNN
│           ├── COIs
│               ├── round1_amplicon_SP27_NNN_SP5_NNN.fasta
│               └── cleaned_amplicon_SP27_NNN_SP5_NNN.fasta
|               ...
│           └── rRNAs
│               ├── round1_amplicon_SP27_NNN_SP5_NNN.fasta
│               └── cleaned_amplicon_SP27_NNN_SP5_NNN.fasta
|               ...
├── rRNA_genes
│   └── <sample>
│       ├── SP27_NNN_SP5_001_NNN_18S.fa
│       ├── SP27_NNN_SP5_001_NNN_28S.fa
│       ...
├── failures
│   └── failures.tsv
├── logs
│   ├── <rule>_<sample>_<adapter_combo>.log
...
```
- `failures.tsv` is a tab-separated file of the stages and amplicons yielding no results OR failing
   - it's good to look at this as the workflow runs to check if something is wrong
- logs are per rule and also explain issues

## Performance on MacBook Pro M3

TBC!

## Memory Considerations

The MacBook Pro M3 with 16GB RAM is sufficient for:
- Single samples use ~2-4GB peak memory
- Recommend running with `-c 6` per job to avoid OOM

If memory is an issue:
- Reduce thread counts (`-c 2`)

If it's too slow:
- Increase thread counts if your laptop is more powerful and can handle it
- NOTE: don't bother upping threads for pychopper, since it's I/O bound. An informal test with 8 threads made it very slow.

## Advanced Usage

### Run only specific rules

```bash
# Only certain rules
snakemake {rule_name} -c 6 --use-conda

```

### Force re-run of failed steps

```bash
snakemake -c 6 --use-conda --rerun-incomplete
```

### Generate reports (fun)

```bash
snakemake --report report.html --use-conda
```

### Use all available cores (caution: may use >8GB RAM)

```bash
snakemake -c 8 --use-conda  # M3 has 8 cores total
```

### Run without conda (if tools already installed)

```bash
snakemake -c 4 --rerun-incomplete
```
(Assumes tools, etc. are in your PATH. Not advised, since the version of the tool will be static)

## Notes

- **Adapter sequences** must match your exact wet-lab protocol. Check `adapters_primers/` files.
- **Primer sequences** should be 5'→3' orientation. Use `seqkit seq -r` to reverse-complement if needed before using.
- **M13 sequences** are provided in both orientations so the search doesn't go for 'both' strands, but searches all options. 
- **Size filters** (min/max amplicon size) are optional but recommended for specificity.
- **Multiple dataset runs** can be done in separate directories using different configs.
- The workflow is **idempotent**: re-running with `--rerun-incomplete` will skip completed steps.

## Citation

This Snakemake workflow wraps the bash scripts from the main branch of:
https://github.com/srisarya/OrCA-seq-processing
Please cite the original repository if publishing results.

Snakemake was made by Johannes Köster et al., ;
Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33.

## License
Same as the OrCA-seq-processing repository main branch :)
