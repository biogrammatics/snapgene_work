#!/usr/bin/env python3
"""
Convert Twist synthesis output to SnapGene .dna files with history.

This script takes:
1. CSV data from Twist (insert sequence, final sequence, backbone name, coordinates)
2. A backbone .dna file (your custom vector)
3. Secretion signal definitions

And produces a SnapGene .dna file with:
- The final sequence
- Features from the backbone (adjusted for length changes)
- Split CDS annotation for the ORF (signal peptide + mature protein)
- History showing the backbone vector and synthetic fragment as inputs
"""

import struct
import lzma
import re
import os
from datetime import datetime, timezone
from dataclasses import dataclass
from typing import Optional
import uuid


# Codon table for translation (standard genetic code)
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def translate(dna_seq: str) -> str:
    """Translate DNA sequence to protein."""
    protein = []
    seq = dna_seq.upper().replace('U', 'T')
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa = CODON_TABLE.get(codon, 'X')
        protein.append(aa)
    return ''.join(protein)


def calc_mw(protein_seq: str) -> float:
    """Calculate molecular weight of protein sequence."""
    aa_weights = {
        'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.15,
        'E': 147.13, 'Q': 146.15, 'G': 75.07, 'H': 155.16, 'I': 131.17,
        'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
        'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15,
        '*': 0, 'X': 110.0
    }
    mass = sum(aa_weights.get(aa, 110.0) for aa in protein_seq)
    mass -= 18.015 * (len(protein_seq.replace('*', '')) - 1)
    return round(mass, 2)


@dataclass
class SecretionSignal:
    """Definition of a secretion signal peptide."""
    name: str
    aa_sequence: str
    color: str = "#0000ff"


# Common yeast secretion signals
SECRETION_SIGNALS = [
    SecretionSignal("alpha factor", "MRFPSIFTAVLFAASSALA"),
    SecretionSignal("alpha factor prepro", "MRFPSIFTAVLFAASSALAAPVNTTTEDETAQIPAEAVIGYSDLEGDFDVAVLPFSNSTNNGLLFINTTIASIAAKEEGVSLEKREAEA"),
]


def make_packet(packet_type: int, data: bytes) -> bytes:
    """Create a packet with type byte and big-endian length."""
    return bytes([packet_type]) + struct.pack('>I', len(data)) + data


def extract_packet(data: bytes, packet_type: int) -> Optional[bytes]:
    """Extract a specific packet type from a .dna file."""
    pos = 0
    while pos < len(data):
        if pos + 5 > len(data):
            break
        ptype = data[pos]
        length = struct.unpack('>I', data[pos+1:pos+5])[0]
        if pos + 5 + length > len(data):
            break
        if ptype == packet_type:
            return data[pos+5:pos+5+length]
        pos += 5 + length
    return None


def extract_all_packets(data: bytes) -> list:
    """Extract all packets from a .dna file."""
    packets = []
    pos = 0
    while pos < len(data):
        if pos + 5 > len(data):
            break
        ptype = data[pos]
        length = struct.unpack('>I', data[pos+1:pos+5])[0]
        if pos + 5 + length > len(data):
            break
        packets.append((ptype, data[pos+5:pos+5+length]))
        pos += 5 + length
    return packets


def parse_features_xml(xml_str: str) -> list:
    """Parse features from XML into a list of feature dicts."""
    features = []
    feature_pattern = re.compile(r'<Feature([^>]*)>(.*?)</Feature>', re.DOTALL)
    segment_pattern = re.compile(r'<Segment([^>]*)/?>')

    for match in feature_pattern.finditer(xml_str):
        attrs_str = match.group(1)
        content = match.group(2)

        feature = {'attrs': {}, 'segments': [], 'qualifiers': []}
        for attr_match in re.finditer(r'(\w+)="([^"]*)"', attrs_str):
            feature['attrs'][attr_match.group(1)] = attr_match.group(2)

        for seg_match in segment_pattern.finditer(content):
            seg_attrs = {}
            for attr_match in re.finditer(r'(\w+)="([^"]*)"', seg_match.group(1)):
                seg_attrs[attr_match.group(1)] = attr_match.group(2)
            feature['segments'].append(seg_attrs)

        for q_match in re.finditer(r'<Q name="([^"]+)"><V ([^/]+)/></Q>', content):
            q_name = q_match.group(1)
            q_val_str = q_match.group(2)
            if 'int="' in q_val_str:
                q_val = int(re.search(r'int="(\d+)"', q_val_str).group(1))
            elif 'text="' in q_val_str:
                q_val = re.search(r'text="([^"]*)"', q_val_str).group(1)
            else:
                q_val = q_val_str
            feature['qualifiers'].append((q_name, q_val))

        features.append(feature)

    return features


def extract_features_xml_content(xml_str: str) -> str:
    """Extract just the Feature elements from a Features XML string."""
    features_content = ""
    feature_pattern = re.compile(r'<Feature[^>]*>.*?</Feature>', re.DOTALL)
    for match in feature_pattern.finditer(xml_str):
        features_content += match.group(0)
    return features_content


def adjust_feature_positions(feature: dict, insert_start: int, insert_end: int,
                            old_length: int, new_length: int) -> Optional[dict]:
    """
    Adjust feature positions after a replacement.
    """
    delta = new_length - old_length
    new_feature = {'attrs': feature['attrs'].copy(), 'segments': [], 'qualifiers': feature['qualifiers'][:]}

    for seg in feature['segments']:
        if 'range' not in seg:
            new_feature['segments'].append(seg.copy())
            continue

        range_match = re.match(r'(\d+)-(\d+)', seg['range'])
        if not range_match:
            new_feature['segments'].append(seg.copy())
            continue

        start = int(range_match.group(1))
        end = int(range_match.group(2))

        if start >= insert_start and end <= insert_end:
            continue

        if start <= insert_end and end >= insert_start:
            continue

        new_seg = seg.copy()
        if start > insert_end:
            new_start = start + delta
            new_end = end + delta
            new_seg['range'] = f"{new_start}-{new_end}"

        new_feature['segments'].append(new_seg)

    if not new_feature['segments']:
        return None

    return new_feature


def build_features_xml(features: list, next_valid_id: int = None) -> str:
    """Build Features XML from list of feature dicts."""
    if next_valid_id is None:
        next_valid_id = len(features)

    xml_parts = [f'<?xml version="1.0"?><Features nextValidID="{next_valid_id}">']

    for i, f in enumerate(features):
        attrs = ' '.join(f'{k}="{v}"' for k, v in f['attrs'].items())
        if 'recentID' not in f['attrs']:
            attrs = f'recentID="{i}" ' + attrs
        xml_parts.append(f'<Feature {attrs}>')

        for seg in f['segments']:
            seg_attrs = ' '.join(f'{k}="{v}"' for k, v in seg.items())
            xml_parts.append(f'<Segment {seg_attrs}/>')

        for q_name, q_val in f.get('qualifiers', []):
            if isinstance(q_val, int):
                xml_parts.append(f'<Q name="{q_name}"><V int="{q_val}"/></Q>')
            else:
                xml_parts.append(f'<Q name="{q_name}"><V text="{q_val}"/></Q>')

        xml_parts.append('</Feature>')

    xml_parts.append('</Features>')
    return ''.join(xml_parts)


def create_orf_feature(name: str, start: int, end: int, sequence: str,
                       signal_peptide: Optional[SecretionSignal] = None,
                       directionality: int = 1) -> dict:
    """Create a CDS feature for an ORF, optionally with signal peptide."""
    protein = translate(sequence)
    mw = calc_mw(protein)

    feature = {
        'attrs': {
            'name': name,
            'directionality': str(directionality),
            'type': 'CDS',
            'translationMW': str(mw),
            'allowSegmentOverlaps': '0',
            'consecutiveTranslationNumbering': '1',
        },
        'segments': [],
        'qualifiers': [
            ('codon_start', 1),
            ('transl_table', 1),
            ('translation', protein),
        ]
    }

    if signal_peptide:
        signal_len_aa = len(signal_peptide.aa_sequence)
        signal_len_nt = signal_len_aa * 3

        if protein.startswith(signal_peptide.aa_sequence):
            signal_end = start + signal_len_nt - 1
            mature_start = signal_end + 1

            feature['attrs']['cleavageArrows'] = str(signal_end)
            feature['segments'] = [
                {'name': 'signal peptide', 'range': f'{start}-{signal_end}',
                 'color': signal_peptide.color, 'type': 'standard', 'translated': '1'},
                {'range': f'{mature_start}-{end}', 'color': '#00bfff',
                 'type': 'standard', 'translated': '1'},
            ]
            signal_aa = protein[:signal_len_aa]
            mature_aa = protein[signal_len_aa:]
            feature['qualifiers'] = [
                ('codon_start', 1),
                ('transl_table', 1),
                ('translation', f'{signal_aa},{mature_aa}'),
            ]
        else:
            feature['segments'] = [
                {'range': f'{start}-{end}', 'color': '#00bfff',
                 'type': 'standard', 'translated': '1'},
            ]
    else:
        feature['segments'] = [
            {'range': f'{start}-{end}', 'color': '#00bfff',
             'type': 'standard', 'translated': '1'},
        ]

    return feature


def build_history_tree_with_fragment(
    final_name: str,
    final_seq_len: int,
    is_circular: bool,
    backbone_name: str,
    backbone_seq_len: int,
    backbone_features_xml: str,
    fragment_name: str,
    fragment_seq_len: int,
    fragment_features_xml: str,
    replace_start: int,
    replace_end: int,
    insert_length: int,
) -> str:
    """
    Build a History Tree XML showing a replacement with backbone and synthetic fragment.

    Structure:
    - Root node: Final construct (operation="replace")
      - Child 1: Backbone vector with features (resurrectable="1", operation="invalid")
      - Child 2: Synthetic fragment with features (resurrectable="1", operation="invalid")

    Args:
        final_name: Name of the final construct
        final_seq_len: Length of final sequence
        is_circular: Whether the final sequence is circular
        backbone_name: Name of the backbone vector
        backbone_seq_len: Length of backbone sequence
        backbone_features_xml: Features XML from the backbone
        fragment_name: Name of the synthetic fragment (e.g., "Twist Synthetic Fragment")
        fragment_seq_len: Length of the synthetic fragment
        fragment_features_xml: Features XML for the fragment (ORF annotation)
        replace_start: Start of replaced region in backbone (1-based)
        replace_end: End of replaced region in backbone (1-based)
        insert_length: Length of the new insert
    """
    circular = "1" if is_circular else "0"

    # Extract just the Feature elements
    backbone_features = extract_features_xml_content(backbone_features_xml)
    fragment_features = extract_features_xml_content(fragment_features_xml)

    # Calculate the insert end position in the final construct
    insert_end_final = replace_start + insert_length - 1

    history = f'''<?xml version="1.0" encoding="UTF-8"?><HistoryTree><Node name="{final_name}" type="DNA" seqLen="{final_seq_len}" strandedness="double" ID="2" circular="{circular}" operation="replace"><HistoryColors><StrandColors type="Product"><TopStrand><ColorRange range="{replace_start}..{insert_end_final}" colors="red"/></TopStrand><BottomStrand><ColorRange range="{replace_start}..{insert_end_final}" colors="red"/></BottomStrand></StrandColors></HistoryColors><InputSummary manipulation="replace" val1="{replace_start}" val2="{replace_end}"/><InputSummary manipulation="overlapAndInsert" val1="0" val2="{insert_length - 1}"/><Node name="{backbone_name}" type="DNA" seqLen="{backbone_seq_len}" strandedness="double" ID="0" circular="{circular}" resurrectable="1" operation="invalid"><Features>{backbone_features}</Features><HistoryColors><StrandColors type="Input"><TopStrand><ColorRange range="{replace_start}..{replace_end}" colors="white"/></TopStrand><BottomStrand><ColorRange range="{replace_start}..{replace_end}" colors="white"/></BottomStrand></StrandColors></HistoryColors></Node><Node name="{fragment_name}" type="DNA" seqLen="{fragment_seq_len}" strandedness="double" ID="1" circular="0" upstreamModification="Unmodified" downstreamModification="Unmodified" resurrectable="1" operation="invalid"><Features>{fragment_features}</Features><HistoryColors><StrandColors type="Input"><TopStrand><ColorRange range="0..{fragment_seq_len - 1}" colors="red"/></TopStrand><BottomStrand><ColorRange range="0..{fragment_seq_len - 1}" colors="red"/></BottomStrand></StrandColors></HistoryColors></Node></Node></HistoryTree>'''

    return history


def make_display_settings(
    enzyme_set_name: str = "Unique 6+ Cutters",
    features_on_circle: bool = True,
    show_enzymes: bool = True,
    circular_view: bool = True
) -> bytes:
    """Create Type 13 display settings packet data."""
    data = bytearray(345)

    data[0] = 0x00
    data[1] = 0x00
    data[2] = 0x01 if features_on_circle else 0x00
    data[3] = 0x00
    data[4] = 0x01
    data[5:7] = b'\x00\x00'
    data[7] = 0x4b

    name_bytes = enzyme_set_name.encode('ascii')[:63]
    data[17] = len(name_bytes)
    data[18:18+len(name_bytes)] = name_bytes

    data[126] = 0x19
    data[127] = 0xff
    data[128] = 0x00
    data[129] = 0x64

    data[134] = 0x00
    data[135] = ord('T') if circular_view else ord('C')
    data[136:138] = b'HO'
    data[138] = 0x00
    data[139] = 0xff
    data[140] = 0xfe

    data[153] = 0x01
    data[154] = 0x01
    data[155] = 0x01 if show_enzymes else 0x00
    data[156] = 0x01
    data[157] = 0x00
    data[158] = 0x01
    data[159] = 0x00
    data[160] = 0x45
    data[169] = 0x01

    data[172] = 0x00
    data[173] = 0x01
    data[174] = 0x00
    data[175] = 0x00
    data[176] = 0x01
    data[177:181] = b'\xff\xff\xff\xff'
    data[181] = 0x01
    data[182] = 0x59
    data[183] = 0x01
    data[184] = 0xf4
    data[185] = 0x01
    data[186] = 0x01
    data[187] = 0x3f
    data[188] = 0x00
    data[189] = 0x50

    return bytes(data)


def generate_twist_construct(
    final_sequence: str,
    backbone_path: str,
    insert_start: int,
    insert_end: int,
    construct_name: str,
    orf_name: str = "ORF",
    fragment_name: str = None,
    signal_peptide: Optional[SecretionSignal] = None,
    output_path: str = None,
) -> str:
    """
    Generate a SnapGene .dna file for a Twist construct.

    The history will show:
    - The backbone vector with its original features
    - The synthetic fragment from Twist with the ORF annotation
    - The final construct as the result of replacing the stuffer region

    Args:
        final_sequence: The complete final DNA sequence
        backbone_path: Path to the backbone .dna file
        insert_start: Start position of replaced region in backbone (1-based)
        insert_end: End position of replaced region in backbone (1-based)
        construct_name: Name for the final construct
        orf_name: Name for the ORF feature
        fragment_name: Name for the synthetic fragment (default: "Twist {orf_name}")
        signal_peptide: Optional secretion signal to detect
        output_path: Output file path (default: construct_name.dna)
    """
    if output_path is None:
        output_path = f"{construct_name}.dna"
    if fragment_name is None:
        fragment_name = f"Twist {orf_name}"

    # Load backbone
    with open(backbone_path, 'rb') as f:
        backbone_data = f.read()

    seq_packet = extract_packet(backbone_data, 0)
    is_circular = bool(seq_packet[0] & 0x01)
    backbone_seq = seq_packet[1:].decode('ascii')
    backbone_name = os.path.basename(backbone_path)

    # Parse backbone features
    features_xml_raw = extract_packet(backbone_data, 10)
    if features_xml_raw:
        features_xml = features_xml_raw.decode('utf-8')
        backbone_features = parse_features_xml(features_xml)
    else:
        backbone_features = []
        features_xml = '<?xml version="1.0"?><Features nextValidID="0"/>'

    # Calculate insert properties
    old_insert_length = insert_end - insert_start + 1
    prefix = backbone_seq[:insert_start-1]
    suffix = backbone_seq[insert_end:]

    # Verify flanking regions match
    if not (final_sequence.upper().startswith(prefix.upper()) and
            final_sequence.upper().endswith(suffix.upper())):
        print(f"Warning: Final sequence doesn't match expected backbone flanking regions")

    insert_sequence = final_sequence[insert_start-1:len(final_sequence)-len(suffix)]
    new_insert_length = len(insert_sequence)

    print(f"Backbone: {backbone_name} ({len(backbone_seq)} bp)")
    print(f"Replaced region: {insert_start}-{insert_end} ({old_insert_length} bp)")
    print(f"New insert: {new_insert_length} bp")
    print(f"Final sequence: {len(final_sequence)} bp")

    # Adjust backbone features for the length change
    adjusted_features = []
    for f in backbone_features:
        adjusted = adjust_feature_positions(f, insert_start, insert_end,
                                           old_insert_length, new_insert_length)
        if adjusted:
            adjusted_features.append(adjusted)

    # Create ORF feature for the insert
    orf_end = insert_start + new_insert_length - 1

    orf_feature = create_orf_feature(
        name=orf_name,
        start=insert_start,
        end=orf_end,
        sequence=insert_sequence,
        signal_peptide=signal_peptide,
    )

    # Add ORF feature at the beginning
    adjusted_features.insert(0, orf_feature)

    # Build final features XML
    final_features_xml = build_features_xml(adjusted_features)

    # Create fragment feature XML (just the ORF, with positions relative to fragment)
    fragment_orf = create_orf_feature(
        name=orf_name,
        start=1,
        end=new_insert_length,
        sequence=insert_sequence,
        signal_peptide=signal_peptide,
    )
    fragment_features_xml = build_features_xml([fragment_orf])

    # Build history tree with backbone and fragment as separate inputs
    history_xml = build_history_tree_with_fragment(
        final_name=construct_name,
        final_seq_len=len(final_sequence),
        is_circular=is_circular,
        backbone_name=backbone_name,
        backbone_seq_len=len(backbone_seq),
        backbone_features_xml=features_xml,
        fragment_name=fragment_name,
        fragment_seq_len=new_insert_length,
        fragment_features_xml=fragment_features_xml,
        replace_start=insert_start,
        replace_end=insert_end,
        insert_length=new_insert_length,
    )

    # Compress history
    history_compressed = lzma.compress(history_xml.encode('utf-8'))

    # Build packets
    packets = []

    # Type 9: File header
    header_data = b'SnapGene\x00\x01\x00\x0f\x00\x14'
    packets.append(make_packet(9, header_data))

    # Type 0: DNA sequence
    flags = 0x01 if is_circular else 0x00
    seq_data = bytes([flags]) + final_sequence.upper().encode('ascii')
    packets.append(make_packet(0, seq_data))

    # Type 7: History tree (XZ compressed)
    packets.append(make_packet(7, history_compressed))

    # Type 8: Additional Sequence Properties
    additional_props = extract_packet(backbone_data, 8)
    if additional_props:
        packets.append(make_packet(8, additional_props))

    # Type 10: Features
    packets.append(make_packet(10, final_features_xml.encode('utf-8')))

    # Type 5: Primers
    primers = extract_packet(backbone_data, 5)
    if primers:
        packets.append(make_packet(5, primers))
    else:
        primers_xml = '<?xml version="1.0"?><Primers nextValidID="0"><HybridizationParams minContinuousMatchLen="15" allowMismatch="1" minMeltingTemperature="40" showAdditionalFivePrimeMatches="1" minimumFivePrimeAnnealing="15"/></Primers>'
        packets.append(make_packet(5, primers_xml.encode('utf-8')))

    # Type 6: Notes
    now = datetime.now(timezone.utc)
    notes_xml = f'''<Notes>
<UUID>{str(uuid.uuid4())}</UUID>
<Type>Synthetic</Type>
<ConfirmedExperimentally>0</ConfirmedExperimentally>
<Created UTC="{now.strftime('%H:%M:%S')}">{now.strftime('%Y.%m.%d')}</Created>
<LastModified UTC="{now.strftime('%H:%M:%S')}">{now.strftime('%Y.%m.%d')}</LastModified>
<SequenceClass>UNA</SequenceClass>
<TransformedInto>unspecified</TransformedInto>
<CustomMapLabel>{construct_name}</CustomMapLabel>
</Notes>'''
    packets.append(make_packet(6, notes_xml.encode('utf-8')))

    # Type 13: Display settings
    display_settings = make_display_settings(
        features_on_circle=True,
        show_enzymes=True,
        circular_view=is_circular,
    )
    packets.append(make_packet(13, display_settings))

    # Type 28: Enzyme visibilities
    enzyme_vis_xml = '<?xml version="1.0"?><EnzymeVisibilities vals=""/>'
    packets.append(make_packet(28, enzyme_vis_xml.encode('utf-8')))

    # Write file
    with open(output_path, 'wb') as f:
        for packet in packets:
            f.write(packet)

    print(f"Written {output_path} ({sum(len(p) for p in packets)} bytes)")
    print(f"History shows: {backbone_name} + {fragment_name} -> {construct_name}")
    return output_path


def main():
    """Example usage."""
    print("Twist to SnapGene Converter")
    print("=" * 40)
    print()
    print("Usage:")
    print("  from twist_to_snapgene import generate_twist_construct, SecretionSignal")
    print()
    print("  # Define your secretion signal (optional)")
    print("  alpha_factor = SecretionSignal('alpha factor', 'MRFPSIFTAVLFAASSALA')")
    print()
    print("  # Generate the construct")
    print("  generate_twist_construct(")
    print("      final_sequence='ATGC...',  # Full final sequence from Twist")
    print("      backbone_path='my_vector.dna',")
    print("      insert_start=100,  # Where insert begins in backbone")
    print("      insert_end=200,    # Where stuffer ends in backbone")
    print("      construct_name='My Construct',")
    print("      orf_name='Client Protein',")
    print("      fragment_name='Twist Client Protein',  # Name for synthetic fragment")
    print("      signal_peptide=alpha_factor,  # Optional")
    print("  )")
    print()
    print("The history will show:")
    print("  - Backbone vector (with all original features)")
    print("  - Synthetic fragment from Twist (with ORF annotation)")
    print("  - Final construct as the result")


if __name__ == '__main__':
    main()
