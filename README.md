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

**AN IMPORTANT NOTE**
If you're using MacOS, as I am, then you'll need to add this line below line 18 of the `amplicon_sorter.py` code
```
multiprocessing.set_start_method("fork", force=True)
```
Since MacOS `spawns` threads by default, and this is the workaround I have found. 

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
   - The adapters, configuration files and primers for demultiplexing and cleaning
   - Which size thresholds to cluster reads by

In this study the dataset configs are in configs/ 

```yaml
# configs/config_mydataset.yaml
dataset_name: "my_dataset"
raw_reads: "input_reads/my_dataset"

# Adapter and primer file paths (relative to repo root)
m13_seqs: "adapters_primers/my_adapters_sequences.fa" # needed for pychopper
m13_config: "adapters_primers/my_adapters_config.txt" # needed for pychopper
fwd_adapters: "adapters_primers/forward_adapters.fa" # needed for demux
rvs_adapters: "adapters_primers/reverse_adapters_reverse_complemented.fa" # needed for demux

# Your specific primers/adapters (relative to repo root)
r1_primers: "adapters_primers/my_primers.fa"
r2_primers: null

# Amplicon types to process, each with its own size filter.
amplicon_types:
  amplicon1:
    min_size: <int>                    # Minimum rRNA length
    max_size: <int>                      # Maximum rRNA length
  amplicon2:
    min_size: <int>                       # Minimum COI length
    max_size: <int>                       # Maximum COI length
```

### Basic Workflow Execution

1. **Prepare input files:**
   - Place FASTQ files (`.fastq.gz`) in the `input_reads/` directory

2. **Generate a workflow visualization:**
   ```bash
   snakemake -s OrCAseq_processing.smk --dag | dot -Tpng > workflow.png
   ```

3. **Dry-run the workflow with default config:**
   ```bash
   snakemake -s OrCAseq_processing.smk -c 4 --use-conda --configfile configs/config_Lakesday1.yaml --dry-run
   ```

4. **Run the protocol if the summary looks right**
   ```bash
   snakemake -s OrCAseq_processing.smk -c 4 --use-conda --configfile configs/config_Lakesday1.yaml
   ```
   
- **Options:**
   - `-c 4`: Use 4 cores (adjust for your system; M3 can handle 4-8)
   - `--use-conda`: Automatically create/activate conda environments.
      - If you want to install tools in your PATH, then you can. See [Advanced Usage](#advanced-usage)
   - `--rerun-incomplete`: Re-run jobs if they fail partway through


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
└── logs/
    ├── pychopper_{sample}.log
    ├── cutadapt_sp5_{sample}.log
    ├── cutadapt_sp27_{sample}.log
    ├── amplicon_sorter_{sample}_{combo}.log
    ├── primer_removal_{sample}_{combo}_{amplicon_type}.log
    └── pybarrnap_{sample}_{combo}.log
```

## Performance on MacBook Pro M3

Estimated runtime for ~1M reads in the raw dataset:

| Step | Time | Notes |
|------|------|-------|
| Pychopper | ~60 min | I/O bound, ~4 threads |
| Cutadapt SP5 | ~20 min | Fast demultiplexing |
| Cutadapt SP27 | ~20-30 min | Per-adapter loop |
| Amplicon Sorter | ~10-20 min per sample | *see below |
| Primer Removal | ~5-10 min | Fast with cutadapt |
| Pybarrnap | ~10-15 min | covariance model search for rRNAs |
| COI reorganisation | ~1-5 min | just moving files |

* After demultiplexing, you will have MANY samples to run! 
* While the amplicon_sorter rule runs reasonably fast, it has a lot to get through.
* So, bear with!
* If you ran a 96-well plate of samples as we did, this will take anywhere between 16h to 23h (960 min - 1920 min)

## Memory Considerations

The MacBook Pro M3 with 16GB RAM is sufficient for:
- Single samples use ~2-4GB peak memory
- Recommend running with `-c 4` per job to avoid OOM

If memory is an issue:
- Reduce thread counts (`-c 2`)

If it's too slow:
- Increase thread counts (`-c 6`) or if your laptop is more powerful, up threads more
- NOTE: don't bother upping threads for pychopper, since it's I/O bound. An informal test with 8 threads made it very slow.

## Advanced Usage

### Run only specific rules

```bash
# Only certain rules
snakemake {rule_name} -c 4 --use-conda

```

### Force re-run of failed steps

```bash
snakemake -c 4 --use-conda --rerun-incomplete --rerun-all
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