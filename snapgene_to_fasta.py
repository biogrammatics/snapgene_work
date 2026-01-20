#!/usr/bin/env python3
"""
Extract DNA sequences from SnapGene .dna files and write to FASTA format.

Designed to run via cron to maintain a FASTA database for BLAST searching.

Usage:
    ./snapgene_to_fasta.py [OPTIONS]

Options:
    --source DIR      Source directory to scan for .dna files
                      (default: /Users/studio/Library/CloudStorage/Dropbox/SnapGene)
    --output FILE     Output FASTA file path
                      (default: ./snapgene_sequences.fasta)
    --min-length N    Minimum sequence length to include (default: 200)
                      Use 0 to include all sequences
    --line-width N    Characters per line in FASTA output (default: 60)
    --quiet           Suppress progress output
    --verbose         Show detailed progress
    --stats           Print statistics at the end

Output format:
    >1_filename_without_extension
    ATGCATGC...

Each entry is numbered sequentially and the filename (with spaces replaced
by underscores and .dna extension removed) is appended to ensure unique IDs.
"""

import struct
import sys
import argparse
from pathlib import Path
from datetime import datetime


def extract_sequence_from_dna(filepath: Path) -> tuple:
    """
    Extract DNA sequence from a SnapGene .dna file.

    Returns:
        Tuple of (sequence: str, is_circular: bool, error: str or None)
    """
    try:
        with open(filepath, 'rb') as f:
            data = f.read()
    except Exception as e:
        return None, False, f"Could not read file: {e}"

    # Parse packets looking for Type 0 (DNA sequence)
    pos = 0
    while pos < len(data):
        if pos + 5 > len(data):
            break

        packet_type = data[pos]
        length = struct.unpack('>I', data[pos+1:pos+5])[0]

        if pos + 5 + length > len(data):
            break

        packet_data = data[pos+5:pos+5+length]

        if packet_type == 0:  # DNA sequence packet
            if len(packet_data) < 1:
                return None, False, "Empty sequence packet"

            flags = packet_data[0]
            is_circular = bool(flags & 0x01)

            try:
                sequence = packet_data[1:].decode('ascii')
                return sequence.upper(), is_circular, None
            except UnicodeDecodeError as e:
                return None, False, f"Could not decode sequence: {e}"

        pos += 5 + length

    return None, False, "No sequence packet found"


def sanitize_name(filename: str) -> str:
    """
    Convert filename to FASTA-safe identifier.

    - Removes .dna extension
    - Replaces spaces with underscores
    - Removes or replaces other problematic characters
    - Converts Unicode characters to ASCII equivalents

    Note: Truncation is handled separately to account for numeric prefix.
    """
    # Remove .dna extension
    name = filename
    if name.lower().endswith('.dna'):
        name = name[:-4]

    # Replace spaces with underscores
    name = name.replace(' ', '_')

    # Replace other problematic characters
    # FASTA headers should avoid: > | ; : and whitespace
    for char in ['>','|', ';', ':', '\t', '\n', '\r']:
        name = name.replace(char, '_')

    # Replace common Unicode characters with ASCII equivalents
    unicode_replacements = {
        '–': '-',   # en-dash
        '—': '-',   # em-dash
        'Δ': 'D',   # Greek delta
        '±': '+',   # plus-minus
        '°': 'deg', # degree
        'µ': 'u',   # micro
        '×': 'x',   # multiplication
        ''': "'",   # smart quote
        ''': "'",   # smart quote
        '"': '"',   # smart quote
        '"': '"',   # smart quote
    }
    for unicode_char, ascii_char in unicode_replacements.items():
        name = name.replace(unicode_char, ascii_char)

    return name


def truncate_id(full_id: str, max_length: int = 50) -> str:
    """
    Truncate FASTA ID to max_length (default 50 for BLAST compatibility).

    If truncation needed, uses first 47 chars + "..." format.
    """
    if len(full_id) <= max_length:
        return full_id
    return full_id[:max_length - 3] + "..."


def format_sequence(sequence: str, line_width: int = 60) -> str:
    """Format sequence with specified line width."""
    lines = []
    for i in range(0, len(sequence), line_width):
        lines.append(sequence[i:i+line_width])
    return '\n'.join(lines)


def find_dna_files(source_dir: Path) -> list:
    """Recursively find all .dna files in source directory."""
    return sorted(source_dir.rglob("*.dna"))


def main():
    parser = argparse.ArgumentParser(
        description="Extract DNA sequences from SnapGene files to FASTA format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage with defaults
    %(prog)s

    # Specify output file
    %(prog)s --output /path/to/sequences.fasta

    # Different source directory
    %(prog)s --source /path/to/dna/files --output sequences.fasta

    # Quiet mode for cron
    %(prog)s --quiet --output /var/data/blast/sequences.fasta

    # Verbose mode with stats
    %(prog)s --verbose --stats
"""
    )

    parser.add_argument(
        '--source', '-s',
        type=Path,
        default=Path("/Users/studio/Library/CloudStorage/Dropbox/SnapGene"),
        help="Source directory to scan for .dna files"
    )

    parser.add_argument(
        '--output', '-o',
        type=Path,
        default=Path("./snapgene_sequences.fasta"),
        help="Output FASTA file path"
    )

    parser.add_argument(
        '--min-length', '-m',
        type=int,
        default=200,
        help="Minimum sequence length to include (default: 200, use 0 for all)"
    )

    parser.add_argument(
        '--line-width', '-w',
        type=int,
        default=60,
        help="Characters per line in FASTA output (default: 60)"
    )

    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help="Suppress progress output"
    )

    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help="Show detailed progress"
    )

    parser.add_argument(
        '--stats',
        action='store_true',
        help="Print statistics at the end"
    )

    args = parser.parse_args()

    # Validate source directory
    if not args.source.exists():
        print(f"Error: Source directory does not exist: {args.source}", file=sys.stderr)
        sys.exit(1)

    if not args.source.is_dir():
        print(f"Error: Source path is not a directory: {args.source}", file=sys.stderr)
        sys.exit(1)

    # Find all .dna files
    if not args.quiet:
        print(f"Scanning {args.source} for .dna files...")

    dna_files = find_dna_files(args.source)
    total_files = len(dna_files)

    if total_files == 0:
        print("No .dna files found", file=sys.stderr)
        sys.exit(1)

    if not args.quiet:
        print(f"Found {total_files} .dna files")

    # Process files and write FASTA
    success_count = 0
    error_count = 0
    skipped_short = 0
    total_bases = 0
    errors = []

    # Ensure output directory exists
    args.output.parent.mkdir(parents=True, exist_ok=True)

    start_time = datetime.now()

    with open(args.output, 'w') as fasta:
        for idx, filepath in enumerate(dna_files, start=1):
            # Progress indicator
            if args.verbose:
                print(f"  [{idx}/{total_files}] {filepath.name}")
            elif not args.quiet and idx % 500 == 0:
                print(f"  Processed {idx}/{total_files} files...")

            # Extract sequence
            sequence, is_circular, error = extract_sequence_from_dna(filepath)

            if error:
                error_count += 1
                errors.append((filepath, error))
                if args.verbose:
                    print(f"    Error: {error}")
                continue

            if not sequence:
                error_count += 1
                errors.append((filepath, "Empty sequence"))
                continue

            # Skip sequences shorter than minimum length
            if len(sequence) < args.min_length:
                skipped_short += 1
                if args.verbose:
                    print(f"    Skipped: {len(sequence)} bp < {args.min_length} bp minimum")
                continue

            # Build FASTA header
            # Format: >NUMBER_FILENAME (truncated to 50 chars for BLAST compatibility)
            safe_name = sanitize_name(filepath.name)
            full_id = f"{idx}_{safe_name}"
            truncated_id = truncate_id(full_id, max_length=50)
            header = f">{truncated_id}"

            # Write entry
            fasta.write(f"{header}\n")
            fasta.write(format_sequence(sequence, args.line_width))
            fasta.write("\n")

            success_count += 1
            total_bases += len(sequence)

    elapsed = datetime.now() - start_time

    # Summary
    if not args.quiet:
        print(f"\nComplete: {success_count} sequences written to {args.output}")
        if skipped_short > 0:
            print(f"Skipped: {skipped_short} sequences shorter than {args.min_length} bp")
        if error_count > 0:
            print(f"Errors: {error_count} files could not be processed")

    # Statistics
    if args.stats or args.verbose:
        print(f"\n{'='*50}")
        print("Statistics:")
        print(f"  Source directory: {args.source}")
        print(f"  Output file: {args.output}")
        print(f"  Minimum length filter: {args.min_length} bp")
        print(f"  Total .dna files found: {total_files}")
        print(f"  Successfully extracted: {success_count}")
        print(f"  Skipped (too short): {skipped_short}")
        print(f"  Errors: {error_count}")
        print(f"  Total bases: {total_bases:,}")
        if success_count > 0:
            print(f"  Average sequence length: {total_bases // success_count:,} bp")
        print(f"  Processing time: {elapsed.total_seconds():.1f} seconds")
        print(f"  Output file size: {args.output.stat().st_size:,} bytes")

    # Show errors if verbose
    if args.verbose and errors:
        print(f"\n{'='*50}")
        print("Errors encountered:")
        for filepath, error in errors[:20]:  # Show first 20
            print(f"  {filepath.name}: {error}")
        if len(errors) > 20:
            print(f"  ... and {len(errors) - 20} more")

    # Exit with error code if there were failures
    if error_count > 0 and success_count == 0:
        sys.exit(1)

    sys.exit(0)


if __name__ == '__main__':
    main()
