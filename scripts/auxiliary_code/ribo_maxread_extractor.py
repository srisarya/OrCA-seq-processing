#!/usr/bin/env python3
"""
SLURM script to extract reads with highest readcount from FASTA files.
Processes subdirectories containing .fa files and extracts the read with
the maximum readcount value from each file (original files remain unchanged).
Written using Claude, from Anthropic. Generated 15/12/2025 using Sonnet 4.5.
"""

#SBATCH --job-name=extract_max_readcount
#SBATCH --output=extract_max_%j.out
#SBATCH --error=extract_max_%j.err
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
from typing import Tuple, Optional


def get_readcount(header: str) -> int:
    """Extract readcount value from FASTA header."""
    match = re.search(r'readcount_(\d+)', header)
    return int(match.group(1)) if match else 0


def read_fasta(filepath: Path) -> list:
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


def write_fasta(filepath: Path, header: str, sequence: str):
    """Write a single FASTA entry to file."""
    with open(filepath, 'w') as f:
        f.write(f"{header}\n")
        f.write(f"{sequence}\n")


def find_max_readcount_entry(entries: list) -> Tuple[Optional[str], Optional[str], int]:
    """
    Find entry with maximum readcount.
    
    Args:
        entries: List of (header, sequence) tuples
        
    Returns:
        Tuple of (max_header, max_sequence, max_readcount)
    """
    max_readcount = 0
    max_header = None
    max_sequence = None
    
    for header, sequence in entries:
        readcount = get_readcount(header)
        if readcount > max_readcount:
            max_readcount = readcount
            max_header = header
            max_sequence = sequence
    
    return max_header, max_sequence, max_readcount


def process_fasta_file(fa_file: Path, output_files: dict, log_file) -> bool:
    """
    Process a single FASTA file: extract max readcount entry and append to consolidated file.
    
    Args:
        fa_file: Path to FASTA file
        output_files: Dictionary of open file handles for 18S and 28S
        log_file: File handle for logging
        
    Returns:
        True if successful, False otherwise
    """
    subdir_name = fa_file.parent.name
    fa_basename = fa_file.name
    
    log_msg = f"Processing: {subdir_name}/{fa_basename}"
    print(log_msg)
    log_file.write(log_msg + "\n")
    log_file.flush()
    
    # Check if file exists and is not empty
    if not fa_file.exists() or fa_file.stat().st_size == 0:
        log_msg = "  Skipping empty or non-existent file"
        print(log_msg)
        log_file.write(log_msg + "\n")
        return False
    
    # Determine if this is 18S or 28S based on filename
    if "18S" in fa_basename or "18s" in fa_basename:
        rna_type = "18S"
    elif "28S" in fa_basename or "28s" in fa_basename:
        rna_type = "28S"
    else:
        log_msg = f"  Warning: Could not determine RNA type (18S/28S) from filename, skipping"
        print(log_msg)
        log_file.write(log_msg + "\n")
        return False
    
    try:
        # Read all entries
        entries = read_fasta(fa_file)
        
        if not entries:
            log_msg = "  No reads found in file"
            print(log_msg)
            log_file.write(log_msg + "\n")
            return False
        
        # Find entry with max readcount
        max_header, max_sequence, max_readcount = find_max_readcount_entry(entries)
        
        if max_header is None:
            log_msg = "  No valid readcount found"
            print(log_msg)
            log_file.write(log_msg + "\n")
            return False
        
        # Append extracted entry to appropriate consolidated file
        output_files[rna_type].write(f"{max_header}\n")
        output_files[rna_type].write(f"{max_sequence}\n")
        output_files[rna_type].flush()
        
        log_msg = f"  Extracted read with readcount={max_readcount} to {rna_type}_max_readcount.fa"
        print(log_msg)
        log_file.write(log_msg + "\n")
        
        return True
        
    except Exception as e:
        log_msg = f"  Error processing file: {str(e)}"
        print(log_msg)
        log_file.write(log_msg + "\n")
        return False


def main():
    """Main processing function."""
    # Get base directory from command line or use current directory
    base_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
    
    # Create output directory
    output_dir = base_dir / "extracted_max_readcount"
    output_dir.mkdir(exist_ok=True)
    
    # Output files for 18S and 28S
    output_18S = output_dir / "18S_max_readcount.fa"
    output_28S = output_dir / "28S_max_readcount.fa"
    
    # Open log file
    log_file_path = output_dir / "extraction_log.txt"
    
    with open(log_file_path, 'w') as log_file, \
         open(output_18S, 'w') as f_18S, \
         open(output_28S, 'w') as f_28S:
        
        output_files = {"18S": f_18S, "28S": f_28S}
        
        start_time = datetime.now()
        log_msg = f"Starting extraction at {start_time}"
        print(log_msg)
        log_file.write(log_msg + "\n")
        
        log_msg = f"Searching for .fa files in {base_dir}"
        print(log_msg)
        log_file.write(log_msg + "\n")
        
        log_msg = f"Output directory: {output_dir}"
        print(log_msg)
        log_file.write(log_msg + "\n")
        
        log_msg = f"18S reads will be written to: {output_18S.name}"
        print(log_msg)
        log_file.write(log_msg + "\n")
        
        log_msg = f"28S reads will be written to: {output_28S.name}"
        print(log_msg)
        log_file.write(log_msg + "\n")
        
        log_file.write("-" * 60 + "\n")
        print("-" * 60)
        
        # Process all .fa files in subdirectories
        fa_count = 0
        success_count = 0
        count_18S = 0
        count_28S = 0
        
        for subdir in sorted(base_dir.iterdir()):
            if subdir.is_dir() and subdir != output_dir:
                for fa_file in sorted(subdir.glob("*.fa")):
                    fa_count += 1
                    if process_fasta_file(fa_file, output_files, log_file):
                        success_count += 1
                        # Count by type
                        if "18S" in fa_file.name or "18s" in fa_file.name:
                            count_18S += 1
                        elif "28S" in fa_file.name or "28s" in fa_file.name:
                            count_28S += 1
        
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
        
        log_msg = f"Successfully processed: {success_count}"
        print(log_msg)
        log_file.write(log_msg + "\n")
        
        log_msg = f"18S reads extracted: {count_18S}"
        print(log_msg)
        log_file.write(log_msg + "\n")
        
        log_msg = f"28S reads extracted: {count_28S}"
        print(log_msg)
        log_file.write(log_msg + "\n")
        
        log_msg = f"Results saved to: {output_dir}"
        print(log_msg)
        log_file.write(log_msg + "\n")


if __name__ == "__main__":
    main()
