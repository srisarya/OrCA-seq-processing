"""
OrCA-seq amplicon processing workflow for local MacBook Pro M3 execution
Covers: pychopper -> cutadapt (SP5 x SP27 demux) -> amplicon_sorter -> primer_removal -> pybarrnap/COI reorganization
"""
# ----------------------------------------
# Setup
# ----------------------------------------
import os
import re
import glob
from pathlib import Path
configfile: "config.yaml"
# ----------------------------------------
# Configuration
# ----------------------------------------
DATASET_NAME = config.get("dataset_name", "dataset")
RAW_READS = config.get("raw_reads", "raw_reads")
WORK_DIR = config.get("work_dir", "results")
os.makedirs(WORK_DIR, exist_ok=True)
# ----------------------------------------
# Helper: raw FASTQ path for a sample (single source of truth)
# ----------------------------------------
def sample_fastq_path(sample):
    """Return the raw FASTQ path for a sample (does not check existence)."""
    return Path(RAW_READS) / f"{sample}.fastq.gz"
# ----------------------------------------
# Discover samples
# ----------------------------------------
raw_dir = Path(RAW_READS)
if not raw_dir.exists():
    raise ValueError(f"Raw reads directory does not exist: {RAW_READS}")
SAMPLES = sorted(
    {
        path.name.removesuffix(".fastq.gz")
        for path in raw_dir.iterdir()
        if path.is_file() and path.name.endswith(".fastq.gz")
    }
)
if not SAMPLES:
    raise ValueError(
        f"No FASTQ files found in {RAW_READS}. Does it have the `fastq.gz` suffix?"
    )
print("Detected samples:")
for sample in SAMPLES:
    print(f"  {sample}")
# ----------------------------------------
# Wildcard constraints
# ----------------------------------------
wildcard_constraints:
    sample = "|".join(
        re.escape(s)
        for s in SAMPLES
    ),
    combo = r"[^/]+",
    amplicon_type = r"COIs|rRNAs"
# ----------------------------------------
# Parameters
# ----------------------------------------
PARAM_DEFAULTS = {
    "m13_seqs": "adapters_primers/M13_seqs_for_pychopper.fa",
    "m13_config": "adapters_primers/M13_config_for_pychopper.txt",
    "sp5_adapters": "adapters_primers/M13_amplicon_indices_forward.fa",
    "sp27_adapters": "adapters_primers/M13_amplicon_indices_reverse_rc.fa",
    "pychopper_q_score": 10,
    "cutadapt_threads": 4,
    "pychopper_threads": 4,
    "amplicon_sorter_threads": 4,
}

M13_SEQS = config.get("m13_seqs", PARAM_DEFAULTS["m13_seqs"])
M13_CONFIG = config.get("m13_config", PARAM_DEFAULTS["m13_config"])
SP5_ADAPTERS = config.get("sp5_adapters", PARAM_DEFAULTS["sp5_adapters"])
SP27_ADAPTERS = config.get("sp27_adapters", PARAM_DEFAULTS["sp27_adapters"])
PYCHOPPER_Q = config.get("pychopper_q_score", PARAM_DEFAULTS["pychopper_q_score"])
CUTADAPT_THREADS = config.get("cutadapt_threads", PARAM_DEFAULTS["cutadapt_threads"])
PYCHOPPER_THREADS = config.get("pychopper_threads", PARAM_DEFAULTS["pychopper_threads"])
AS_THREADS = config.get("amplicon_sorter_threads", PARAM_DEFAULTS["amplicon_sorter_threads"])

INVALID_SP27 = [
    "SP27_009",
    "SP27_010",
    "SP27_011",
    "SP27_012",
]

# ----------------------------------------
# Amplicon decision schema
# ----------------------------------------
_amplicon_types_cfg = config.get("amplicon_types")
if not _amplicon_types_cfg:
    raise ValueError(
        "config 'amplicon_types' must define at least one of 'COIs'/'rRNAs', "
        "each with optional min_size/max_size, e.g.:\n"
        "amplicon_types:\n  COIs:\n    min_size: 300\n    max_size: 900"
    )

_VALID_AMPLICON_TYPES = ("COIs", "rRNAs")
unknown = set(_amplicon_types_cfg) - set(_VALID_AMPLICON_TYPES)
if unknown:
    raise ValueError(
        f"Invalid amplicon_types entries {sorted(unknown)}; "
        f"must be one or more of {_VALID_AMPLICON_TYPES}"
    )

AMPLICON_TYPES = [t for t in _VALID_AMPLICON_TYPES if t in _amplicon_types_cfg]

AMPLICON_SIZE_LIMITS = {}
for _t in AMPLICON_TYPES:
    _entry = _amplicon_types_cfg[_t] or {}
    AMPLICON_SIZE_LIMITS[_t] = {
        "min_size": _entry.get("min_size"),
        "max_size": _entry.get("max_size"),
    }

# ----------------------------------------
# Primer files per amplicon type
# ----------------------------------------
_primers_cfg = config.get("primers")
if not _primers_cfg:
    raise ValueError(
        "config 'primers' must map each amplicon_type to a primer FASTA, e.g.:\n"
        "primers:\n  COIs: adapters_primers/COI_primers.fa\n"
        "  rRNAs: adapters_primers/RNA_primers.fa"
    )

missing_primers = set(AMPLICON_TYPES) - set(_primers_cfg)
if missing_primers:
    raise ValueError(
        f"config 'primers' is missing entries for amplicon_types "
        f"{sorted(missing_primers)}; every configured amplicon_type needs "
        f"a primer file"
    )

PRIMER_FILES = {t: _primers_cfg[t] for t in AMPLICON_TYPES}
# ----------------------------------------
# Helper: discover SP5 identifiers
# ----------------------------------------
def get_sp5_identifiers(wildcards):
    """
    Discover the SP5 identifiers produced by the SP5 checkpoint.
    This function is evaluated after cutadapt_sp5 has completed.
    """
    checkpoint_output = (
        checkpoints.cutadapt_sp5
        .get(sample=wildcards.sample)
        .output.demux_dir
    )
    ck_dir = str(checkpoint_output)
    suffix = f"_{DATASET_NAME}.fastq.gz"
    pattern = os.path.join(
        ck_dir,
        f"*{suffix}"
    )
    files = glob.glob(pattern)
    identifiers = []
    for filepath in files:
        basename = os.path.basename(filepath)
        if not basename.endswith(suffix):
            continue
        identifier = basename[:-len(suffix)]
        if identifier == "unknown":
            continue
        # Ignore the cutadapt JSON file if its name happens to match anything unexpectedly.
        if not filepath.endswith(".fastq.gz"):
            continue
        identifiers.append(identifier)
    return sorted(set(identifiers))
# ----------------------------------------
# Helper: discover SP27 combinations
# ----------------------------------------
def get_sp27_combos(wildcards):
    """
    Discover the concrete SP27/SP5 combinations produced by
    cutadapt_sp27 for a given sample.
    Example:
        SP27_001_SP5_003
        SP27_001_SP5_002
        SP27_006_SP5_008
    """
    checkpoint_output = (
        checkpoints.cutadapt_sp27
        .get(sample=wildcards.sample)
        .output.demux_dir
    )
    ck_dir = str(checkpoint_output)
    suffix = f"_{DATASET_NAME}.fastq.gz"
    files = glob.glob(
        os.path.join(
            ck_dir,
            f"*{suffix}"
        )
    )
    combos = []
    for filepath in files:
        basename = os.path.basename(filepath)
        if not basename.endswith(suffix):
            continue
        combo = basename[:-len(suffix)]
        if combo == "unknown":
            continue
        combos.append(combo)
    return sorted(set(combos))
# ----------------------------------------
# Helper: discover final targets
# ----------------------------------------
def get_final_targets(wildcards):
    """
    Return the actual final output files generated from the
    combinations discovered by cutadapt_sp27, for every amplicon
    type this run is configured for (AMPLICON_TYPES). When both
    COIs and rRNAs are listed in config `amplicon_types`, both
    sets of final outputs are requested for every combo.
    The dependency chain is:
        cutadapt_sp5
            ↓
        cutadapt_sp27
            ↓
        amplicon_sorter (once per amplicon_type)
            ↓
        primer_removal
            ↓
        pybarrnap_extract (rRNAs) / reorganize_cois (COIs)
            ↓
        final targets
    """
    targets = []
    for sample in SAMPLES:
        combos = get_sp27_combos(
            type(
                "Wildcards",
                (),
                {"sample": sample}
            )()
        )
        for combo in combos:
            if "rRNAs" in AMPLICON_TYPES:
                targets.append(
                    os.path.join(
                        WORK_DIR,
                        "rRNA_genes",
                        sample,
                        f"{combo}_18S.fa"
                    )
                )
            if "COIs" in AMPLICON_TYPES:
                targets.append(
                    os.path.join(
                        WORK_DIR,
                        "COI_gene",
                        sample,
                        f"{combo}_COI.fasta"
                    )
                )
    return targets
# ----------------------------------------
# Rule: all
# ----------------------------------------
rule all:
    input:
        get_final_targets
# ----------------------------------------
# 1: Pychopper
# ----------------------------------------
rule pychopper:
    input:
        fastq=lambda wc: str(sample_fastq_path(wc.sample))
    output:
        pass_fastq=(
            f"{WORK_DIR}/pychopped/"
            f"{{sample}}_pass.fastq"
        ),
        rescued_fastq=(
            f"{WORK_DIR}/pychopped/"
            f"{{sample}}_rescued.fastq"
        ),
        unclass_fastq=(
            f"{WORK_DIR}/pychopped/"
            f"{{sample}}_unclass.fastq"
        ),
        short_fastq=(
            f"{WORK_DIR}/pychopped/"
            f"{{sample}}_short.fastq"
        ),
        stats=(
            f"{WORK_DIR}/pychopped/"
            f"{{sample}}_stats.out"
        )
    threads:
        PYCHOPPER_THREADS
    conda:
        "envs/pychopper.yaml"
    log:
        f"{WORK_DIR}/logs/pychopper_{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.pass_fastq}")"
        mkdir -p "$(dirname "{log}")"
        pychopper \
            -b "{M13_SEQS}" \
            -c "{M13_CONFIG}" \
            -k LSK114 \
            -Q {PYCHOPPER_Q} \
            -w "{output.rescued_fastq}" \
            -u "{output.unclass_fastq}" \
            -l "{output.short_fastq}" \
            -S "{output.stats}" \
            -p \
            -t {threads} \
            -m edlib \
            "{input.fastq}" \
            > "{output.pass_fastq}" \
            2> "{log}"
        """
# ----------------------------------------
# Gzip pychopper output
# ----------------------------------------
rule gzip_pychopped:
    input:
        f"{WORK_DIR}/pychopped/{{sample}}_pass.fastq"
    output:
        f"{WORK_DIR}/pychopped/{{sample}}_pass.fastq.gz"
    threads:
        1
    shell:
        r"""
        set -euo pipefail
        gzip -c "{input}" > "{output}"
        rm -f "{input}"
        """
# ----------------------------------------
# 2a: Cutadapt SP5
# ----------------------------------------
checkpoint cutadapt_sp5:
    input:
        f"{WORK_DIR}/pychopped/{{sample}}_pass.fastq.gz"
    output:
        demux_dir=directory(
            f"{WORK_DIR}/demuxed/SP5/{{sample}}"
        )
    threads:
        CUTADAPT_THREADS
    conda:
        "envs/cutadapt.yaml"
    log:
        f"{WORK_DIR}/logs/cutadapt_sp5_{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{output.demux_dir}"
        mkdir -p "$(dirname "{log}")"
        cutadapt \
            --action=trim \
            -e 0.1 \
            -j {threads} \
            --rc \
            -g "file:{SP5_ADAPTERS}" \
            -o "{output.demux_dir}/{{name}}_{DATASET_NAME}.fastq.gz" \
            "{input}" \
            --json="{output.demux_dir}/cutadapt_SP5_{DATASET_NAME}.json" \
            2> "{log}"
        find "{output.demux_dir}" \
            -type f \
            -name "*unknown*" \
            -delete || true
        """
# ----------------------------------------
# 2b: Cutadapt SP27
# ----------------------------------------
checkpoint cutadapt_sp27:
    input:
        sp5_dir=(
            f"{WORK_DIR}/demuxed/SP5/{{sample}}"
        )
    output:
        demux_dir=directory(
            f"{WORK_DIR}/demuxed/SP27/{{sample}}"
        )
    params:
        identifiers=get_sp5_identifiers,
        invalid=INVALID_SP27
    threads:
        CUTADAPT_THREADS
    conda:
        "envs/cutadapt.yaml"
    log:
        f"{WORK_DIR}/logs/cutadapt_sp27_{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{output.demux_dir}"
        mkdir -p "$(dirname "{log}")"
        : > "{log}"
        # Process every SP5 identifier discovered by the cutadapt_sp5 checkpoint.
        for identifier in {params.identifiers}; do
            echo "Processing: ${{identifier}}" >> "{log}"
            input_file="{input.sp5_dir}/${{identifier}}_{DATASET_NAME}.fastq.gz"
            if [ ! -f "$input_file" ]; then
                echo "WARNING: missing $input_file" >> "{log}"
                continue
            fi
            cutadapt \
                --action=trim \
                -e 0.1 \
                -j {threads} \
                --rc \
                -a "file:{SP27_ADAPTERS}" \
                -o "{output.demux_dir}/{{name}}_${{identifier}}_{DATASET_NAME}.fastq.gz" \
                "$input_file" \
                --json="{output.demux_dir}/${{identifier}}_{DATASET_NAME}.json" \
                >> "{log}" 2>&1
        done
        # Remove unknown reads.
        find "{output.demux_dir}" \
            -type f \
            -name "*unknown*" \
            -delete || true
        # Remove invalid SP27 combinations.
        for bad in {params.invalid}; do
            find "{output.demux_dir}" \
                -type f \
                -name "${{bad}}_*_{DATASET_NAME}.fastq.gz" \
                -delete || true
        done
        """
# ----------------------------------------
# 3: AmpliconSorter
# ----------------------------------------
rule amplicon_sorter:
    """
    Cluster and sort amplicons by sequence similarity and size.
    One job runs per (sample, combo, amplicon_type) discovered/configured,
    each with its own size filter from AMPLICON_SIZE_LIMITS, so a single
    workflow run can produce both COI and rRNA outputs from the same
    demuxed reads when both are listed in config `amplicon_types`.
    """
    input:
        fastq=lambda wc: os.path.join(
            WORK_DIR,
            "demuxed",
            "SP27",
            wc.sample,
            f"{wc.combo}_{DATASET_NAME}.fastq.gz"
        ),
        sp27_checkpoint=lambda wc: (
            checkpoints.cutadapt_sp27
            .get(sample=wc.sample)
            .output.demux_dir
        )
    params:
        min_flag=lambda wc: (
            f"-min {AMPLICON_SIZE_LIMITS[wc.amplicon_type]['min_size']}"
            if AMPLICON_SIZE_LIMITS[wc.amplicon_type]["min_size"] is not None
            else ""
        ),
        max_flag=lambda wc: (
            f"-max {AMPLICON_SIZE_LIMITS[wc.amplicon_type]['max_size']}"
            if AMPLICON_SIZE_LIMITS[wc.amplicon_type]["max_size"] is not None
            else ""
        )
    output:
        consensus=(
            f"{WORK_DIR}/amplicon_sorted/"
            f"{{sample}}/{{combo}}/"
            f"{{amplicon_type}}/{{combo}}_consensus_{{amplicon_type}}.fasta"
        )
    threads:
        AS_THREADS
    conda:
        "envs/amplicon_sorter.yaml"
    log:
        f"{WORK_DIR}/logs/"
        f"amplicon_sorter_{{sample}}_{{combo}}_{{amplicon_type}}.log"
    shell:
        r"""
        set -euo pipefail
        
        # AmpliconSorter output directory -- one per (sample, combo, amplicon_type)
        # so that COI and rRNA runs for the same combo never share intermediates.
        # Example:
        # results/amplicon_sorted/
        #     Lakes_day1/
        #         SP27_003_SP5_001/
        #             COIs/
        #             rRNAs/
        
        outdir="$(dirname "{output.consensus}")"
        mkdir -p "$outdir"
        mkdir -p "$(dirname "{log}")"
        
        # Build AmpliconSorter command
        
        cmd=(
            python3
            scripts/auxiliary_code/amplicon_sorter.py
            -i "{input.fastq}"
            -o "$outdir"
            -ar
            -np {threads}
            {params.min_flag}
            {params.max_flag}
        )
        
        # Run AmpliconSorter
        
        "${{cmd[@]}}" \
            > "{log}" \
            2>&1
        
        # Check that AmpliconSorter produced its expected
        # consensus file.
        
        if [ ! -f "$outdir/consensusfile.fasta" ]; then
            echo \
                "ERROR: consensusfile.fasta not created" \
                >> "{log}"
            exit 1
        fi
        
        # Convert AmpliconSorter read-count notation.
        # Example: sequence(123) becomes sequence_readcount_123
        
        seqkit replace \
            -p '\\((\d+)\\)$' \
            -r '_readcount_$1' \
            "$outdir/consensusfile.fasta" \
            > "$outdir/temp.fa"
        
        # Replace AmpliconSorter group numbering with explicit
        # group identifiers.
        
        awk '
        BEGIN {{ counter = 1 }}
        /^>/ {{
            if (match($0, /_[0-9]+_[0-9]+_readcount/)) {{
                sub(/_[0-9]+_[0-9]+_readcount/, "_group" counter "_readcount")
                counter++
            }}
        }}
        {{ print }}
        ' "$outdir/temp.fa" \
            > "$outdir/classified.fasta"
        
        # Remove intermediate files.
        
        rm -f \
            "$outdir/temp.fa" \
            "$outdir/consensusfile.fasta"
        
        # AmpliconSorter clusters by similarity/size; it does not itself
        # classify sequences as COI vs rRNA. Each (sample, combo,
        # amplicon_type) job applies that type's own size filter
        # (-min/-max above) and writes straight to that type's own
        # directory, so classified.fasta just needs renaming in place.
        
        cp \
            "$outdir/classified.fasta" \
            "{output.consensus}"
        """
# ----------------------------------------
# 4: Primer removal
# ----------------------------------------
rule primer_removal:
    """
    Remove forward and reverse primers from consensus sequences using cutadapt.
    One job runs per (sample, combo, amplicon_type), using the primer file
    configured for that amplicon_type.
    """
    input:
        fasta=(
            f"{WORK_DIR}/amplicon_sorted/"
            f"{{sample}}/{{combo}}/"
            f"{{amplicon_type}}/{{combo}}_consensus_{{amplicon_type}}.fasta"
        ),
        primers=lambda wc: PRIMER_FILES[wc.amplicon_type]
    output:
        fasta=(
            f"{WORK_DIR}/primerless/"
            f"{{sample}}/{{combo}}/"
            f"{{amplicon_type}}/cleaned_amplicon_{{combo}}.fasta"
        )
    threads:
        2
    conda:
        "envs/cutadapt.yaml"
    log:
        f"{WORK_DIR}/logs/"
        f"primer_removal_{{sample}}_{{combo}}_{{amplicon_type}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.fasta}")"
        mkdir -p "$(dirname "{log}")"
        
        # Remove forward primers (5' end) from all sequences
        cutadapt \
            -g "file:{input.primers}" \
            -o "{output.fasta}" \
            "{input.fasta}" \
            2> "{log}"
        
        # Remove reverse primers (3' end) from all sequences
        # cutadapt will automatically reverse-complement for -a mode
        cutadapt \
            -a "file:{input.primers}" \
            -o "{output.fasta}.tmp" \
            "{output.fasta}" \
            >> "{log}" 2>&1
        
        mv "{output.fasta}.tmp" "{output.fasta}"
        
        # Report stats
        echo "Primer removal complete for {wildcards.sample}/{wildcards.combo}/{wildcards.amplicon_type}" \
            >> "{log}"
        """

# ----------------------------------------
# 5a: Extract rRNAs
# ----------------------------------------
rule pybarrnap_extract:
    """
    Extract 18S and 28S rRNA sequences using barrnap.
    """
    input:
        fasta=(
            f"{WORK_DIR}/primerless/"
            f"{{sample}}/{{combo}}/"
            f"rRNAs/cleaned_amplicon_{{combo}}.fasta"
        )
    output:
        fasta_18s=(
            f"{WORK_DIR}/rRNA_genes/"
            f"{{sample}}/{{combo}}_18S.fa"
        ),
        fasta_28s=(
            f"{WORK_DIR}/rRNA_genes/"
            f"{{sample}}/{{combo}}_28S.fa"
        )
    threads:
        2
    conda:
        "envs/pybarrnap.yaml"
    log:
        f"{WORK_DIR}/logs/"
        f"pybarrnap_{{sample}}_{{combo}}.log"
    shell:
        r"""
        set -euo pipefail
        outdir="$(dirname "{output.fasta_18s}")"
        mkdir -p "$outdir"
        mkdir -p "$(dirname "{log}")"
        temp_dir="$outdir/{wildcards.combo}_barrnap_temp"
        mkdir -p "$temp_dir"
        barrnap \
            -k euk \
            --incseq \
            -o "$temp_dir/{wildcards.combo}_euk.fa" \
            "{input.fasta}" \
            > "$temp_dir/{wildcards.combo}_euk.gff3" \
            2> "{log}"
        seqkit grep \
            -r \
            -p "18S_rRNA" \
            "$temp_dir/{wildcards.combo}_euk.fa" \
            > "{output.fasta_18s}" \
            2>> "{log}" || true
        seqkit grep \
            -r \
            -p "28S_rRNA" \
            "$temp_dir/{wildcards.combo}_euk.fa" \
            > "{output.fasta_28s}" \
            2>> "{log}" || true
        rm -rf "$temp_dir"
        """
# ----------------------------------------
# 5b: Reorganize COIs
# ----------------------------------------
rule reorganize_cois:
    """
    Copy the primer-cleaned COI FASTA into the final
    COI_gene output directory.
    """
    input:
        fasta=(
            f"{WORK_DIR}/primerless/"
            f"{{sample}}/{{combo}}/"
            f"COIs/cleaned_amplicon_{{combo}}.fasta"
        )
    output:
        fasta=(
            f"{WORK_DIR}/COI_gene/"
            f"{{sample}}/{{combo}}_COI.fasta"
        )
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.fasta}")"
        cp \
            "{input.fasta}" \
            "{output.fasta}"
        """
