#!/usr/bin/env python3
"""Analyze features across many .dna files, grouping by DNA sequence."""

import struct
import xml.etree.ElementTree as ET
from pathlib import Path
from collections import defaultdict, Counter
import hashlib

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

def get_sequence(packets):
    """Extract DNA sequence from Type 0 packet."""
    if 0 not in packets:
        return None
    seq_data = packets[0][0]
    return seq_data[1:].decode('ascii', errors='replace').upper()

def parse_features(packets):
    """Parse features from Type 10 packet."""
    if 10 not in packets:
        return []
    
    try:
        xml_data = packets[10][0].decode('utf-8')
        root = ET.fromstring(xml_data)
    except Exception:
        return []
    
    features = []
    for feat in root.findall('.//Feature'):
        name = feat.get('name', '')
        feat_type = feat.get('type', '')
        directionality = feat.get('directionality', '1')
        
        # Get segments (ranges)
        segments = []
        for seg in feat.findall('.//Segment'):
            range_str = seg.get('range', '')
            if range_str:
                segments.append(range_str)
        
        features.append({
            'name': name,
            'type': feat_type,
            'directionality': directionality,
            'segments': segments
        })
    
    return features

def extract_feature_sequence(sequence, segments, directionality):
    """Extract the DNA sequence for a feature based on its segments."""
    if not sequence or not segments:
        return None
    
    # Parse segments and extract sequence
    seq_parts = []
    for seg in segments:
        try:
            if '..' in seg:
                start, end = seg.split('..')
                start = int(start)
                end = int(end)
            elif '-' in seg:
                start, end = seg.split('-')
                start = int(start)
                end = int(end)
            else:
                continue
            
            # SnapGene uses 1-based coordinates
            # Handle circular sequences (start > end means wrapping)
            if start <= end:
                seq_parts.append(sequence[start-1:end])
            else:
                # Wrapping around circular sequence
                seq_parts.append(sequence[start-1:] + sequence[:end])
        except (ValueError, IndexError):
            continue
    
    if not seq_parts:
        return None
    
    feat_seq = ''.join(seq_parts)
    
    # Reverse complement if on reverse strand
    if directionality == '2':
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        feat_seq = ''.join(complement.get(b, b) for b in reversed(feat_seq))
    
    return feat_seq

# Main analysis
base_path = "/Users/tom/Library/CloudStorage/Dropbox/SnapGene/BioGrammatics"
dna_files = list(Path(base_path).rglob("*.dna"))

print(f"Analyzing {len(dna_files)} .dna files...")

# Group features by sequence
# Key: sequence hash -> { 'sequence': str, 'names': Counter, 'types': Counter, 'count': int, 'files': list }
features_by_seq = defaultdict(lambda: {
    'sequence': None,
    'names': Counter(),
    'types': Counter(),
    'count': 0,
    'files': [],
    'length': 0
})

file_count = 0
feature_count = 0

for i, filepath in enumerate(dna_files):
    if i % 500 == 0:
        print(f"Processing {i}/{len(dna_files)}...")
    
    packets, err = extract_packets(filepath)
    if err:
        continue
    
    sequence = get_sequence(packets)
    if not sequence:
        continue
    
    features = parse_features(packets)
    if not features:
        continue
    
    file_count += 1
    
    for feat in features:
        feat_seq = extract_feature_sequence(
            sequence, 
            feat['segments'], 
            feat['directionality']
        )
        
        if not feat_seq or len(feat_seq) < 10:  # Skip very short features
            continue
        
        feature_count += 1
        
        # Use hash of sequence as key (for efficiency)
        seq_hash = hashlib.md5(feat_seq.encode()).hexdigest()
        
        entry = features_by_seq[seq_hash]
        entry['sequence'] = feat_seq
        entry['names'][feat['name']] += 1
        entry['types'][feat['type']] += 1
        entry['count'] += 1
        entry['length'] = len(feat_seq)
        if len(entry['files']) < 5:  # Keep up to 5 example files
            entry['files'].append(filepath.name)

print(f"\nAnalyzed {file_count} files with features")
print(f"Found {feature_count} total features")
print(f"Found {len(features_by_seq)} unique feature sequences")

# Sort by count and get top 50
top_features = sorted(features_by_seq.values(), key=lambda x: x['count'], reverse=True)[:50]

print("\n" + "=" * 80)
print("TOP 50 FEATURES BY OCCURRENCE")
print("=" * 80)

for i, feat in enumerate(top_features, 1):
    # Get most common name
    top_name = feat['names'].most_common(1)[0][0]
    top_type = feat['types'].most_common(1)[0][0]
    name_variants = len(feat['names'])
    
    print(f"\n{i}. {top_name} ({top_type})")
    print(f"   Count: {feat['count']} occurrences in files")
    print(f"   Length: {feat['length']} bp")
    
    if name_variants > 1:
        other_names = [n for n, c in feat['names'].most_common(5) if n != top_name]
        if other_names:
            print(f"   Also known as: {', '.join(other_names[:3])}")
    
    # Show first 60 bp of sequence
    seq_preview = feat['sequence'][:60]
    if len(feat['sequence']) > 60:
        seq_preview += "..."
    print(f"   Sequence: {seq_preview}")

# Save full results to file
print("\n" + "=" * 80)
print("Saving full results to feature_analysis.txt...")

with open('/Users/tom/Claude/snapgene_scripts/feature_analysis.txt', 'w') as f:
    f.write("TOP 50 FEATURES BY OCCURRENCE\n")
    f.write("=" * 80 + "\n\n")
    
    for i, feat in enumerate(top_features, 1):
        top_name = feat['names'].most_common(1)[0][0]
        top_type = feat['types'].most_common(1)[0][0]
        
        f.write(f"{i}. {top_name} ({top_type})\n")
        f.write(f"   Count: {feat['count']}\n")
        f.write(f"   Length: {feat['length']} bp\n")
        f.write(f"   Names: {dict(feat['names'].most_common(10))}\n")
        f.write(f"   Types: {dict(feat['types'])}\n")
        f.write(f"   Example files: {feat['files'][:3]}\n")
        f.write(f"   Sequence:\n")
        # Write sequence in 80-char lines
        seq = feat['sequence']
        for j in range(0, len(seq), 80):
            f.write(f"   {seq[j:j+80]}\n")
        f.write("\n")

print("Done!")
