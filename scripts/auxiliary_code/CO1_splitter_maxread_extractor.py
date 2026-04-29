#!/usr/bin/env python3
"""
SLURM script to extract reads with highest readcount from FASTA files,
separated by sequence length. Extracts longest reads (>600bp) to 'moorea'
and shortest reads (<350bp) to 'sauron' output files.
Written by claude.ai, Sonnet 4.5, accessed 12/01/25.
"""

#SBATCH --job-name=extract_by_size
#SBATCH --output=extract_size_%j.out
#SBATCH --error=extract_size_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

import os
import sys
import re
from pathlib import Path
from datetime import datetime
from typing import Tuple, Optional, List


def get_readcount(header: str) -> int:
    """Extract readcount value from FASTA header."""
    match = re.search(r'readcount_(\d+)', header)
    return int(match.group(1)) if match else 0


def read_fasta(filepath: Path) -> List[Tuple[str, str]]:
    """
    Read FASTA file and return list of (header, sequence) tuples.
    
    Args:
        filepath: Path to FASTA file
        
    Returns:
        List of tuples containing (header, sequence)
    """
    entries = []
    current_header = None
    current_sequence = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous entry
                if current_header is not None:
                    entries.append((current_header, ''.join(current_sequence)))
                # Start new entry
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)
        
        # Save last entry
        if current_header is not None:
            entries.append((current_header, ''.join(current_sequence)))
    
    return entries


def categorize_by_length(entries: List[Tuple[str, str]]) -> Tuple[List, List]:
    """
    Categorize entries by sequence length.
    
    Args:
        entries: List of (header, sequence) tuples
        
    Returns:
        Tuple of (moorea_entries, sauron_entries)
        moorea: sequences >= 600bp
        sauron: sequences < 350bp
    """
    moorea = []  # >= 600bp
    sauron = []  # < 350bp
    
    for header, sequence in entries:
        seq_length = len(sequence)
        if seq_length >= 600:
            moorea.append((header, sequence, seq_length))
        elif seq_length < 350:
            sauron.append((header, sequence, seq_length))
        # Sequences between 350-599bp are not included in either category
    
    return moorea, sauron


def find_max_readcount_entry(entries: List[Tuple[str, str, int]]) -> Tuple[Optional[str], Optional[str], int, int]:
    """
    Find entry with maximum readcount from categorized entries.
    
    Args:
        entries: List of (header, sequence, length) tuples
        
    Returns:
        Tuple of (max_header, max_sequence, max_readcount, sequence_length)
    """
    max_readcount = 0
    max_header = None
    max_sequence = None
    max_length = 0
    
    for header, sequence, length in entries:
        readcount = get_readcount(header)
        if readcount > max_readcount:
            max_readcount = readcount
            max_header = header
            max_sequence = sequence
            max_length = length
    
    return max_header, max_sequence, max_readcount, max_length


def process_fasta_file(fa_file: Path, output_files: dict, log_file) -> dict:
    """
    Process a single FASTA file: extract max readcount entries by size category.
    
    Args:
        fa_file: Path to FASTA file
        output_files: Dictionary of open file handles for moorea and sauron
        log_file: File handle for logging
        
    Returns:
        Dictionary with processing statistics
    """
    subdir_name = fa_file.parent.name
    fa_basename = fa_file.name
    
    log_msg = f"Processing: {subdir_name}/{fa_basename}"
    print(log_msg)
    log_file.write(log_msg + "\n")
    log_file.flush()
    
    stats = {"moorea": False, "sauron": False}
    
    # Check if file exists and is not empty
    if not fa_file.exists() or fa_file.stat().st_size == 0:
        log_msg = "  Skipping empty or non-existent file"
        print(log_msg)
        log_file.write(log_msg + "\n")
        return stats
    
    try:
        # Read all entries
        entries = read_fasta(fa_file)
        
        if not entries:
            log_msg = "  No reads found in file"
            print(log_msg)
            log_file.write(log_msg + "\n")
            return stats
        
        # Categorize by length
        moorea_entries, sauron_entries = categorize_by_length(entries)
        
        log_msg = f"  Found {len(moorea_entries)} reads >=600bp and {len(sauron_entries)} reads <350bp"
        print(log_msg)
        log_file.write(log_msg + "\n")
        
        # Process moorea (>=600bp) entries
        if moorea_entries:
            max_header, max_sequence, max_readcount, seq_len = find_max_readcount_entry(moorea_entries)
            if max_header:
                output_files["moorea"].write(f"{max_header}\n")
                output_files["moorea"].write(f"{max_sequence}\n")
                output_files["moorea"].flush()
                
                log_msg = f"  MOOREA: Extracted read with readcount={max_readcount}, length={seq_len}bp"
                print(log_msg)
                log_file.write(log_msg + "\n")
                stats["moorea"] = True
        else:
            log_msg = "  MOOREA: No reads >=600bp found"
            print(log_msg)
            log_file.write(log_msg + "\n")
        
        # Process sauron (<350bp) entries
        if sauron_entries:
            max_header, max_sequence, max_readcount, seq_len = find_max_readcount_entry(sauron_entries)
            if max_header:
                output_files["sauron"].write(f"{max_header}\n")
                output_files["sauron"].write(f"{max_sequence}\n")
                output_files["sauron"].flush()
                
                log_msg = f"  SAURON: Extracted read with readcount={max_readcount}, length={seq_len}bp"
                print(log_msg)
                log_file.write(log_msg + "\n")
                stats["sauron"] = True
        else:
            log_msg = "  SAURON: No reads <350bp found"
            print(log_msg)
            log_file.write(log_msg + "\n")
        
        return stats
        
    except Exception as e:
        log_msg = f"  Error processing file: {str(e)}"
        print(log_msg)
        log_file.write(log_msg + "\n")
        return stats


def main():
    """Main processing function."""
    # Get base directory from command line or use current directory
    base_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
    
    # Create output directory
    output_dir = base_dir / "extracted_by_size"
    output_dir.mkdir(exist_ok=True)
    
    # Output files for moorea (>=600bp) and sauron (<350bp)
    output_moorea = output_dir / "moorea.fa"
    output_sauron = output_dir / "sauron.fa"
    
    # Open log file
    log_file_path = output_dir / "extraction_log.txt"
    
    with open(log_file_path, 'w') as log_file, \
         open(output_moorea, 'w') as f_moorea, \
         open(output_sauron, 'w') as f_sauron:
        
        output_files = {"moorea": f_moorea, "sauron": f_sauron}
        
        start_time = datetime.now()
        log_msg = f"Starting extraction at {start_time}"
        print(log_msg)
        log_file.write(log_msg + "\n")
        
        log_msg = f"Searching for .fa and .fasta files in {base_dir}"
        print(log_msg)
        log_file.write(log_msg + "\n")
        
        log_msg = f"Output directory: {output_dir}"
        print(log_msg)
        log_file.write(log_msg + "\n")
        
        log_msg = f"Moorea (>=600bp) reads will be written to: {output_moorea.name}"
        print(log_msg)
        log_file.write(log_msg + "\n")
        
        log_msg = f"Sauron (<350bp) reads will be written to: {output_sauron.name}"
        print(log_msg)
        log_file.write(log_msg + "\n")
        
        log_file.write("-" * 60 + "\n")
        print("-" * 60)
        
        # Process all .fa and .fasta files in subdirectories
        fa_count = 0
        count_moorea = 0
        count_sauron = 0
        
        for subdir in sorted(base_dir.iterdir()):
            if subdir.is_dir() and subdir != output_dir:
                # Look for both .fa and .fasta files
                for fa_file in sorted(list(subdir.glob("*.fa")) + list(subdir.glob("*.fasta"))):
                    fa_count += 1
                    stats = process_fasta_file(fa_file, output_files, log_file)
                    if stats["moorea"]:
                        count_moorea += 1
                    if stats["sauron"]:
                        count_sauron += 1
        
        # Final summary
        end_time = datetime.now()
        duration = end_time - start_time
        
        print("-" * 60)
        log_file.write("-" * 60 + "\n")
        
        log_msg = f"Processing complete at {end_time}"
        print(log_msg)
        log_file.write(log_msg + "\n")
        
        log_msg = f"Duration: {duration}"
        print(log_msg)
        log_file.write(log_msg + "\n")
        
        log_msg = f"Total files processed: {fa_count}"
        print(log_msg)
        log_file.write(log_msg + "\n")
        
        log_msg = f"Moorea reads extracted (>=600bp): {count_moorea}"
        print(log_msg)
        log_file.write(log_msg + "\n")
        
        log_msg = f"Sauron reads extracted (<350bp): {count_sauron}"
        print(log_msg)
        log_file.write(log_msg + "\n")
        
        log_msg = f"Results saved to: {output_dir}"
        print(log_msg)
        log_file.write(log_msg + "\n")


if __name__ == "__main__":
    main()
