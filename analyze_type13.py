#!/usr/bin/env python3
"""
Analyze Type 13 packets across many SnapGene .dna files to understand the structure.
"""

import struct
import os
import sys
from collections import defaultdict, Counter
from pathlib import Path


def extract_packets(filepath):
    """Extract all packets from a .dna file."""
    try:
        with open(filepath, 'rb') as f:
            data = f.read()
    except Exception as e:
        return None, str(e)

    packets = {}
    pos = 0
    while pos < len(data):
        if pos + 5 > len(data):
            break
        ptype = data[pos]
        length = struct.unpack('>I', data[pos+1:pos+5])[0]
        if pos + 5 + length > len(data):
            break
        pdata = data[pos+5:pos+5+length]

        if ptype not in packets:
            packets[ptype] = []
        packets[ptype].append(pdata)

        pos += 5 + length

    return packets, None


def analyze_type13(data):
    """Analyze a Type 13 packet and extract key fields."""
    if len(data) < 200:
        return None

    result = {
        'length': len(data),
        'header': data[0:8].hex(),
    }

    # Find THO marker
    tho_pos = data.find(b'THO')
    result['tho_position'] = tho_pos

    # Extract enzyme set name (at offset 17-18)
    if len(data) > 82:
        name_len = data[17]
        if name_len > 0 and name_len < 64:
            try:
                name = data[18:18+name_len].decode('ascii', errors='replace')
                result['enzyme_set_name'] = name
            except:
                result['enzyme_set_name'] = None
        else:
            result['enzyme_set_name'] = None

    # Extract bytes around known settings area (using THO as anchor if found)
    if tho_pos > 0:
        # Key bytes relative to THO
        for offset in range(-10, 60):
            pos = tho_pos + offset
            if 0 <= pos < len(data):
                result[f'tho{offset:+d}'] = data[pos]
    else:
        # Absolute positions
        for pos in [126, 127, 128, 129, 134, 135, 136, 137, 138, 139,
                    153, 154, 155, 156, 157, 158, 159, 160, 169, 173,
                    176, 177, 178, 179, 180, 181, 182, 183, 184, 185,
                    186, 187, 188, 189, 190]:
            if pos < len(data):
                result[f'byte_{pos}'] = data[pos]

    return result


def main():
    base_path = "/Users/tom/Library/CloudStorage/Dropbox/SnapGene/BioGrammatics"

    # Find all .dna files
    dna_files = list(Path(base_path).rglob("*.dna"))
    print(f"Found {len(dna_files)} .dna files")

    # Analyze Type 13 packets
    type13_data = []
    files_with_type13 = 0
    files_without_type13 = 0
    errors = 0

    # Track byte value distributions
    byte_values = defaultdict(Counter)  # byte_position -> Counter of values

    # Track THO positions
    tho_positions = Counter()

    # Track packet lengths
    lengths = Counter()

    # Track enzyme set names
    enzyme_sets = Counter()

    for i, filepath in enumerate(dna_files):
        if i % 500 == 0:
            print(f"Processing {i}/{len(dna_files)}...", file=sys.stderr)

        packets, err = extract_packets(filepath)
        if err:
            errors += 1
            continue

        if 13 not in packets:
            files_without_type13 += 1
            continue

        files_with_type13 += 1

        for pdata in packets[13]:
            analysis = analyze_type13(pdata)
            if analysis:
                type13_data.append(analysis)
                lengths[analysis['length']] += 1

                if analysis.get('tho_position') is not None:
                    tho_positions[analysis['tho_position']] += 1

                if analysis.get('enzyme_set_name'):
                    enzyme_sets[analysis['enzyme_set_name']] += 1

                # Track byte values at each position (relative to THO)
                for key, value in analysis.items():
                    if key.startswith('tho') and isinstance(value, int):
                        byte_values[key][value] += 1

    print()
    print("=" * 70)
    print("ANALYSIS RESULTS")
    print("=" * 70)
    print()
    print(f"Files analyzed: {len(dna_files)}")
    print(f"Files with Type 13: {files_with_type13}")
    print(f"Files without Type 13: {files_without_type13}")
    print(f"Errors: {errors}")
    print()

    print("=" * 70)
    print("TYPE 13 PACKET LENGTHS")
    print("=" * 70)
    for length, count in lengths.most_common(20):
        print(f"  {length} bytes: {count} files")
    print()

    print("=" * 70)
    print("THO MARKER POSITIONS")
    print("=" * 70)
    for pos, count in tho_positions.most_common(10):
        print(f"  Offset {pos}: {count} files")
    print()

    print("=" * 70)
    print("ENZYME SET NAMES (Top 30)")
    print("=" * 70)
    for name, count in enzyme_sets.most_common(30):
        print(f"  {count:5d}: {name}")
    print()

    print("=" * 70)
    print("BYTE VALUE ANALYSIS (relative to THO)")
    print("=" * 70)
    print()
    print("Bytes with HIGH variance (likely settings):")
    print("-" * 50)

    # Find bytes with high variance (multiple different values)
    variable_bytes = []
    constant_bytes = []

    def sort_key(x):
        # Handle keys like 'tho+5', 'tho-3', 'byte_126'
        if x.startswith('tho'):
            offset_str = x[3:]  # Remove 'tho' prefix
            if offset_str.startswith('+'):
                return int(offset_str[1:])
            elif offset_str.startswith('-'):
                return int(offset_str)
            else:
                return 0
        elif x.startswith('byte_'):
            return int(x[5:])
        return 0

    for key in sorted(byte_values.keys(), key=sort_key):
        values = byte_values[key]
        if len(values) > 1:
            variable_bytes.append((key, values))
        else:
            constant_bytes.append((key, values))

    for key, values in variable_bytes:
        total = sum(values.values())
        print(f"  {key}:")
        for val, count in values.most_common(5):
            pct = count / total * 100
            char = chr(val) if 32 <= val < 127 else '.'
            print(f"    0x{val:02x} ({val:3d}) '{char}': {count:5d} ({pct:5.1f}%)")
        print()

    print()
    print("Bytes that are CONSTANT (same value in all files):")
    print("-" * 50)
    for key, values in constant_bytes:
        val = list(values.keys())[0]
        count = values[val]
        char = chr(val) if 32 <= val < 127 else '.'
        print(f"  {key}: 0x{val:02x} ({val:3d}) '{char}' - {count} files")


if __name__ == '__main__':
    main()
