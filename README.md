# OrCA-seq Processing Workflow (Snakemake)

A Snakemake-based workflow for processing COI/28S/18S amplicon data locally!

## Version Info

- **Snakemake**: ≥7.0
- **Python**: 3.8+
- **Conda/Mamba**: Required
- **Tested on**: MacBook Pro M3 (8 cores, 16GB RAM)

## Workflow Steps

1. **Pychopper** (01): Reorient cDNA reads and quality trim
2. **Cutadapt (SP5/SP27)** (02): Two-round demultiplexing by adapter indices
3. **Amplicon Sorter** (03): Cluster reads by sequence similarity and filter by size
4. **Primer Removal** (04): Strip 5'/3' primers from amplicon sequences
5. **Pybarrnap/Reorganize** (05): Extract rRNA genes or reorganize COIs

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
   ├── Snakefile
   ├── config.yaml
   ├── configs/
   │   └── config_for_dataset.yaml
   ├── envs/
   │   ├── pychopper.yaml
   │   ├── cutadapt.yaml
   │   ├── amplicon_sorter.yaml
   │   └── pybarrnap.yaml
   ├── adapters_primers/
   ├── scripts/
   │   ├── auxiliary_code/
   │   │   └── amplicon_sorter.py
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
   snakemake -s OrCAseq_processing.smk --use-conda --configfile configs/config_Lakesday1.yaml --dry-run
   ```

3. **Run the protocol with 6 threads if the summary looks right**
   ```bash
   snakemake -s OrCAseq_processing.smk -c 6 --use-conda --configfile configs/config_Lakesday1.yaml
   ```
   
- **Options:**
   - `-c 4`: Use 4 cores (adjust for your system; M3 can handle 4-8)
   - `--use-conda`: Automatically create/activate conda environments.
      - If you want to install tools in your PATH, then you can. See [Advanced Usage](#advanced-usage)
   - `--rerun-incomplete`: Re-run jobs if they fail partway through
   - once the workflow has completed, you can run it like so `snakemake -s OrCAseq_processing.smk --configfile configs/config_Lakesday1.yaml --report` to get a HTML report of runtimes and timestamps :)


## Output Structure (summary)

```
results/
├── pychopped/
│   └── {sample}_pass.fastq.gz
├── demuxed/
│   ├── SP5/{sample}/{SP5_id}_{DATASET_NAME}.fastq.gz
│   └── SP27/{sample}/{combo}_{DATASET_NAME}.fastq.gz      # combo = SP27_xxx_SP5_yyy
├── amplicon_sorted/
│   └── {sample}/{combo}/
│       ├── rRNAs/{combo}_consensus_rRNAs.fasta
│       └── COIs/{combo}_consensus_COIs.fasta
├── primerless/
│   └── {sample}/{combo}/
│       ├── rRNAs/cleaned_amplicon_{combo}.fasta
│       └── COIs/cleaned_amplicon_{combo}.fasta
├── rRNA_genes/
│   └── {sample}/
│       ├── {combo}_18S.fa
│       └── {combo}_28S.fa
├── COI_gene/
│   └── {sample}/
│       └── {combo}_COI.fasta
├── failures/
    └── failures.tsv
└── logs/
    ├── pychopper_{sample}.log
    ├── cutadapt_sp5_{sample}.log
    ├── cutadapt_sp27_{sample}.log
    ├── amplicon_sorter_{sample}_{combo}.log
    ├── primer_removal_{sample}_{combo}_{amplicon_type}.log
    └── pybarrnap_{sample}_{combo}.log
```
- `failures.tsv` is a tab-separated file of the stages and amplicons yielding no results OR failing
   - it's good to look at this as the workflow runs to check if something is wrong
- logs are per rule and also explain issues

## Performance on MacBook Pro M3

Estimated runtime for ~2.6M reads in the raw dataset (Lakes Day 1 data, 4.4G, from the NHM training workshop in Ambleside)

| Step | Mean compute time across all jobs | Sum compute time for all jobs | Notes |
|------|------|-------|
| Pychopper | ~30 min | ~30 min | I/O bound, ~4 threads |
| Cutadapt SP5 | ~10 min | ~10 min | Fast demultiplexing |
| Cutadapt SP27 | ~5 min | ~5 min | Per-adapter loop |
| Amplicon Sorter | ~20 min | ~72h | *see below |
| Primer Removal | ~1 second| ~4 min | We use the round 1 trimming onl
| Pybarrnap | ~1 second | ~3 min | covariance model search for rRNAs |
| COI reorganisation | milliseconds | ~2 seconds | just moving files |

* After demultiplexing, you will have MANY samples to run! 
* While the amplicon_sorter rule runs reasonably fast, it has a lot to get through.
* With multithreading, it does not take the whole 72h!!!!
* I ran this with 6 cores and 2 cores per job, so it took 72/3=24h to finish amplicon sorting
* Total run time was 1 day, 56 minutes, 8.21 seconds

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
snakemake -c 6 --use-conda --rerun-incomplete --rerun-all
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
(Assumes `pychopper`, `cutadapt`, etc. are in your PATH)

## Notes

- **Adapter sequences** must match your exact wet-lab protocol. Check `adapters_primers/` files.
- **Primer sequences** in config should be 5'→3' orientation. Use `seqkit seq -r` to reverse-complement if needed before using.
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