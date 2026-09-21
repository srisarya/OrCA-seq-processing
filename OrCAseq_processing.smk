"""
OrCA-seq amplicon processing workflow for local MacBook Pro M3 execution
Covers: pychopper -> cutadapt (SP5 x SP27 demux) -> amplicon_sorter -> primer_removal -> barrnap/COI reorganization
"""
# ----------------------------------------
# Setup
# ----------------------------------------
import os
import re
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
    "primer_removal_threads": 1,
    "barrnap_threads": 1,
}

M13_SEQS = config.get("m13_seqs", PARAM_DEFAULTS["m13_seqs"])
M13_CONFIG = config.get("m13_config", PARAM_DEFAULTS["m13_config"])
SP5_ADAPTERS = config.get("sp5_adapters", PARAM_DEFAULTS["sp5_adapters"])
SP27_ADAPTERS = config.get("sp27_adapters", PARAM_DEFAULTS["sp27_adapters"])
PYCHOPPER_Q = config.get("pychopper_q_score", PARAM_DEFAULTS["pychopper_q_score"])
CUTADAPT_THREADS = config.get("cutadapt_threads", PARAM_DEFAULTS["cutadapt_threads"])
PYCHOPPER_THREADS = config.get("pychopper_threads", PARAM_DEFAULTS["pychopper_threads"])
AS_THREADS = config.get("amplicon_sorter_threads", PARAM_DEFAULTS["amplicon_sorter_threads"])
PRIMER_REMOVAL_THREADS = config.get("primer_removal_threads", PARAM_DEFAULTS["primer_removal_threads"])
BARRNAP_THREADS = config.get("barrnap_threads", PARAM_DEFAULTS["barrnap_threads"])  

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
# Static barcode combinations
# ----------------------------------------
def read_adapter_ids(path, prefix):
    """Read unique adapter IDs from FASTA headers and validate their prefix."""
    identifiers = []
    with open(path) as adapter_file:
        for line in adapter_file:
            if line.startswith(">"):
                identifier = line[1:].strip().split()[0]
                if not identifier.startswith(prefix):
                    raise ValueError(
                        f"Unexpected adapter ID '{identifier}' in {path}; "
                        f"expected prefix '{prefix}'"
                    )
                identifiers.append(identifier)
    if not identifiers or len(identifiers) != len(set(identifiers)):
        raise ValueError(f"Adapter FASTA must contain unique IDs: {path}")
    return sorted(identifiers)


SP5_IDENTIFIERS = read_adapter_ids(SP5_ADAPTERS, "SP5_")
SP27_IDENTIFIERS = read_adapter_ids(SP27_ADAPTERS, "SP27_")
VALID_SP27_IDENTIFIERS = [
    identifier
    for identifier in SP27_IDENTIFIERS
    if identifier not in INVALID_SP27
]
COMBOS = [
    f"{sp27}_{sp5}"
    for sp27 in VALID_SP27_IDENTIFIERS
    for sp5 in SP5_IDENTIFIERS
]
SP5_FASTQ_OUTPUTS = [
    f"{WORK_DIR}/demuxed/SP5/{{{{sample}}}}/{identifier}_{DATASET_NAME}.fastq.gz"
    for identifier in SP5_IDENTIFIERS
]
SP27_FASTQ_OUTPUTS = [
    f"{WORK_DIR}/demuxed/SP27/{{{{sample}}}}/{combo}_{DATASET_NAME}.fastq.gz"
    for combo in COMBOS
]
SP27_REPORT_OUTPUTS = [
    f"{WORK_DIR}/demuxed/SP27/{{{{sample}}}}/{identifier}_{DATASET_NAME}.json"
    for identifier in SP5_IDENTIFIERS
]

if not VALID_SP27_IDENTIFIERS:
    raise ValueError("No valid SP27 adapters remain after INVALID_SP27 filtering")

FINAL_TARGETS = []
for sample in SAMPLES:
    for combo in COMBOS:
        if "rRNAs" in AMPLICON_TYPES:
            FINAL_TARGETS.extend(
                [
                    os.path.join(WORK_DIR, "rRNA_genes", sample, f"{combo}_18S.fa"),
                    os.path.join(WORK_DIR, "rRNA_genes", sample, f"{combo}_28S.fa"),
                ]
            )
        if "COIs" in AMPLICON_TYPES:
            FINAL_TARGETS.append(
                os.path.join(WORK_DIR, "COI_gene", sample, f"{combo}_COI.fasta")
            )
# ----------------------------------------
# Rule: all
# ----------------------------------------
rule all:
    input:
        FINAL_TARGETS
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
rule cutadapt_sp5:
    input:
        f"{WORK_DIR}/pychopped/{{sample}}_pass.fastq.gz"
    output:
        fastqs=SP5_FASTQ_OUTPUTS,
        report=f"{WORK_DIR}/demuxed/SP5/{{sample}}/cutadapt_SP5_{DATASET_NAME}.json"
    threads:
        CUTADAPT_THREADS
    conda:
        "envs/cutadapt.yaml"
    log:
        f"{WORK_DIR}/logs/cutadapt_sp5_{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        outdir="$(dirname "{output.fastqs[0]}")"
        mkdir -p "$outdir"
        mkdir -p "$(dirname "{log}")"
        cutadapt \
            --action=trim \
            -e 0.1 \
            -j {threads} \
            --rc \
            -g "file:{SP5_ADAPTERS}" \
            -o "$outdir/{{name}}_{DATASET_NAME}.fastq.gz" \
            "{input}" \
            --json="{output.report}" \
            2> "{log}"
        find "$outdir" \
            -type f \
            -name "*unknown*" \
            -delete || true
        """
# ----------------------------------------
# 2b: Cutadapt SP27
# ----------------------------------------
rule cutadapt_sp27:
    input:
        sp5_fastqs=SP5_FASTQ_OUTPUTS
    output:
        fastqs=SP27_FASTQ_OUTPUTS,
        reports=SP27_REPORT_OUTPUTS
    params:
        identifiers=" ".join(SP5_IDENTIFIERS),
        invalid=" ".join(INVALID_SP27)
    threads:
        CUTADAPT_THREADS
    conda:
        "envs/cutadapt.yaml"
    log:
        f"{WORK_DIR}/logs/cutadapt_sp27_{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        outdir="$(dirname "{output.fastqs[0]}")"
        mkdir -p "$outdir"
        mkdir -p "$(dirname "{log}")"
        : > "{log}"
        for identifier in {params.identifiers}; do
            echo "Processing: ${{identifier}}" >> "{log}"
            input_file="$(dirname "{input.sp5_fastqs[0]}")/${{identifier}}_{DATASET_NAME}.fastq.gz"
            cutadapt \
                --action=trim \
                -e 0.1 \
                -j {threads} \
                --rc \
                -a "file:{SP27_ADAPTERS}" \
                -o "$outdir/{{name}}_${{identifier}}_{DATASET_NAME}.fastq.gz" \
                "$input_file" \
                --json="$outdir/${{identifier}}_{DATASET_NAME}.json" \
                >> "{log}" 2>&1
        done
        # Remove unknown reads.
        find "$outdir" \
            -type f \
            -name "*unknown*" \
            -delete || true
        for bad in {params.invalid}; do
            find "$outdir" \
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

        outdir="$(dirname "{output.consensus}")"
        mkdir -p "$outdir"
        mkdir -p "$(dirname "{log}")"
        mkdir -p "{WORK_DIR}/failures"

        cmd=(
            python
            scripts/auxiliary_code/amplicon_sorter.py
            -i "{input.fastq}"
            -o "$outdir"
            -ar
            -np {threads}
            {params.min_flag}
            {params.max_flag}
        )

        if ! "${{cmd[@]}}" > "{log}" 2>&1; then
            echo -e "{wildcards.sample}\t{wildcards.combo}\t{wildcards.amplicon_type}\tamplicon_sorter\t{log}" >> "{WORK_DIR}/failures/failures.tsv"
            touch "{output.consensus}"
            exit 0
        fi

        if [ ! -f "$outdir/consensusfile.fasta" ]; then
            echo -e "{wildcards.sample}\t{wildcards.combo}\t{wildcards.amplicon_type}\tamplicon_sorter\t{log}" >> "{WORK_DIR}/failures/failures.tsv"
            touch "{output.consensus}"
            exit 0
        fi

                awk -v combo="{wildcards.combo}" \
                    -v dataset="{DATASET_NAME}" \
                '
        BEGIN {{
            dataset_label = dataset
            sub(/^[^_]+_/, "", dataset_label)
        }}
        /^>/ {{
            prefix = ">consensus_" combo "_" dataset "_"
            if (index($0, prefix) == 1) {{
                suffix = substr($0, length(prefix) + 1)
                gsub(/[()]/, "_", suffix)
                split(suffix, fields, "_")
                if (fields[1] ~ /^[0-9]+$/ && fields[3] ~ /^[0-9]+$/) {{
                    $0 = ">consensus_" combo "_" dataset_label \
                        "_pass_group" fields[1] "_readcount_" fields[3]
                }}
            }}
        }}
        {{ print }}
          "$outdir/consensusfile.fasta" \
            > "$outdir/classified.fasta"

        rm -f "$outdir/consensusfile.fasta"

        cp "$outdir/classified.fasta" "{output.consensus}"
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
        PRIMER_REMOVAL_THREADS
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
            -j {threads} \
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

        # Drop any zero-length records left after primer trimming
        seqkit seq -m 1 "{output.fasta}" > "{output.fasta}.filtered"
        mv "{output.fasta}.filtered" "{output.fasta}"
        
        # Report stats
        echo "Primer removal complete for {wildcards.sample}/{wildcards.combo}/{wildcards.amplicon_type}" \
            >> "{log}"
        """

# ----------------------------------------
# 5a: Extract rRNAs
# ----------------------------------------
rule barrnap_extract:
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
        ),
        filtered_fasta=temp(
            f"{WORK_DIR}/rRNA_genes/"
            f"{{sample}}/{{combo}}_barrnap_nonempty.fa"
        ),
        barrnap_fasta=temp(
            f"{WORK_DIR}/rRNA_genes/"
            f"{{sample}}/{{combo}}_barrnap_euk.fa"
        ),
        barrnap_gff=temp(
            f"{WORK_DIR}/rRNA_genes/"
            f"{{sample}}/{{combo}}_barrnap_euk.gff3"
        )
    threads:
        BARRNAP_THREADS
    conda:
        "envs/barrnap.yaml"
    log:
        f"{WORK_DIR}/logs/"
        f"barrnap_{{sample}}_{{combo}}.log"
    shell:
        r"""
        set -euo pipefail
        outdir="$(dirname "{output.fasta_18s}")"
        mkdir -p "$outdir"
        mkdir -p "$(dirname "{log}")"

        seqkit seq -m 10 "{input.fasta}" > "{output.filtered_fasta}"
        if [ ! -s "{output.filtered_fasta}" ]; then
            echo "WARNING: {input.fasta} contains no FASTA entries at least 10 bp long, skipping barrnap" > "{log}"
            touch "{output.fasta_18s}" "{output.fasta_28s}" \
                "{output.barrnap_fasta}" "{output.barrnap_gff}"
            exit 0
        fi

        barrnap \
            -k euk \
            --incseq \
            --threads {threads} \
            -o "{output.barrnap_fasta}" \
            "{output.filtered_fasta}" \
            > "{output.barrnap_gff}" \
            2> "{log}"
        seqkit grep \
            -r \
            -p "18S_rRNA" \
            "{output.barrnap_fasta}" \
            > "{output.fasta_18s}" \
            2>> "{log}" || true
        seqkit grep \
            -r \
            -p "28S_rRNA" \
            "{output.barrnap_fasta}" \
            > "{output.fasta_28s}" \
            2>> "{log}" || true
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
