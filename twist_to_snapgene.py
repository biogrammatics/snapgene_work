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
import html
import hashlib
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
    name: str           # Full name for main CDS segment (e.g., "Sc AGA2 signal peptide")
    aa_sequence: str    # Amino acid sequence
    short_name: str = None  # Short name for separate feature (e.g., "s5")
    color: str = "#3366ff"

    def __post_init__(self):
        if self.short_name is None:
            self.short_name = self.name.split()[0]


# Common yeast secretion signals
SECRETION_SIGNALS = {
    "alpha_factor": SecretionSignal("alpha factor", "MRFPSIFTAVLFAASSALA", "alpha factor"),
    "alpha_factor_prepro": SecretionSignal("alpha factor prepro", "MRFPSIFTAVLFAASSALAAPVNTTTEDETAQIPAEAVIGYSDLEGDFDVAVLPFSNSTNNGLLFINTTIASIAAKEEGVSLEKREAEA", "alpha factor prepro"),
    "s5": SecretionSignal("Sc AGA2 signal peptide", "MQLLRCFSIFSVIASVLA", "s5"),
}

# His tag sequence
HIS_TAG_6X = "HHHHHH"


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
    # Match <Feature with space to avoid matching <Features
    feature_pattern = re.compile(r'<Feature(\s[^>]*)>(.*?)</Feature>', re.DOTALL)
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


def extract_features_xml_content(xml_str: str, strip_recent_id: bool = False) -> str:
    """Extract just the Feature elements from a Features XML string.

    Args:
        xml_str: The full Features XML string
        strip_recent_id: If True, remove recentID attributes (for history nodes)
    """
    features_content = ""
    # Match <Feature with space or > to avoid matching <Features
    feature_pattern = re.compile(r'<Feature(?:\s[^>]*)?>.*?</Feature>', re.DOTALL)
    for match in feature_pattern.finditer(xml_str):
        feature_xml = match.group(0)
        if strip_recent_id:
            # Remove recentID attribute for history nodes
            feature_xml = re.sub(r'\s*recentID="[^"]*"', '', feature_xml)
        features_content += feature_xml
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
        attrs = ' '.join(f'{k}="{html.escape(str(v), quote=True)}"' for k, v in f['attrs'].items())
        if 'recentID' not in f['attrs']:
            attrs = f'recentID="{i}" ' + attrs
        xml_parts.append(f'<Feature {attrs}>')

        for seg in f['segments']:
            seg_attrs = ' '.join(f'{k}="{html.escape(str(v), quote=True)}"' for k, v in seg.items())
            xml_parts.append(f'<Segment {seg_attrs}/>')

        for q_name, q_val in f.get('qualifiers', []):
            if isinstance(q_val, int):
                xml_parts.append(f'<Q name="{q_name}"><V int="{q_val}"/></Q>')
            else:
                # Escape XML special characters in text values
                escaped_val = html.escape(str(q_val), quote=True)
                xml_parts.append(f'<Q name="{q_name}"><V text="{escaped_val}"/></Q>')

        xml_parts.append('</Feature>')

    xml_parts.append('</Features>')
    return ''.join(xml_parts)


def create_orf_feature(name: str, start: int, end: int, sequence: str,
                       signal_peptide: Optional[SecretionSignal] = None,
                       mature_name: str = None,
                       directionality: int = 1) -> dict:
    """Create a CDS feature for an ORF, optionally with signal peptide.

    Args:
        name: Name of the full CDS feature
        start: Start position (1-based)
        end: End position (1-based)
        sequence: DNA sequence
        signal_peptide: Optional SecretionSignal to detect at N-terminus
        mature_name: Name for the mature protein segment (default: "secreted {name}")
        directionality: 1 for forward, 2 for reverse
    """
    protein = translate(sequence)
    mw = calc_mw(protein)

    if mature_name is None:
        mature_name = f"secreted {name}"

    feature = {
        'attrs': {
            'name': name,
            'directionality': str(directionality),
            'type': 'CDS',
            'translationMW': str(mw),
            'allowSegmentOverlaps': '0',
            'consecutiveTranslationNumbering': '1',
            'hitsStopCodon': '1',
        },
        'segments': [],
        'qualifiers': [
            ('codon_start', 1),
            ('gene', f'<html><body>{name}</body></html>'),
            ('product', f'<html><body>{name}</body></html>'),
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

            feature['segments'] = [
                {'name': signal_peptide.name, 'range': f'{start}-{signal_end}',
                 'color': signal_peptide.color, 'type': 'standard', 'translated': '1'},
                {'name': mature_name, 'range': f'{mature_start}-{end}',
                 'color': '#00bfff', 'type': 'standard', 'translated': '1'},
            ]
            signal_aa = protein[:signal_len_aa]
            mature_aa = protein[signal_len_aa:]
            feature['qualifiers'] = [
                ('codon_start', 1),
                ('gene', f'<html><body>{name}</body></html>'),
                ('product', f'<html><body>{name}</body></html>'),
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


def create_signal_peptide_feature(name: str, start: int, end: int, sequence: str,
                                   orf_end: int, color: str = "#3366ff") -> dict:
    """Create a separate CDS feature for a signal peptide annotation."""
    protein = translate(sequence)
    mw = calc_mw(protein)

    return {
        'attrs': {
            'name': name,
            'directionality': '1',
            'type': 'CDS',
            'translationMW': str(mw),
            'allowSegmentOverlaps': '0',
            'consecutiveTranslationNumbering': '1',
            'maxRunOn': str(orf_end),
            'maxFusedRunOn': str(orf_end),
            'detectionMode': 'exactProteinMatch',
        },
        'segments': [
            {'range': f'{start}-{end}', 'color': color, 'type': 'standard', 'translated': '1'},
        ],
        'qualifiers': [
            ('codon_start', 1),
            ('transl_table', 1),
            ('translation', protein),
        ]
    }


def create_his_tag_feature(start: int, end: int, orf_end: int) -> dict:
    """Create a CDS feature for a 6xHis tag annotation."""
    return {
        'attrs': {
            'name': '6xHis',
            'directionality': '1',
            'type': 'CDS',
            'translationMW': '840.86',
            'allowSegmentOverlaps': '0',
            'consecutiveTranslationNumbering': '1',
            'maxRunOn': str(orf_end),
            'maxFusedRunOn': str(orf_end),
            'detectionMode': 'exactProteinMatch',
        },
        'segments': [
            {'range': f'{start}-{end}', 'color': '#cc99b2', 'type': 'standard', 'translated': '1'},
        ],
        'qualifiers': [
            ('codon_start', 1),
            ('product', '<html><body>6xHis affinity tag</body></html>'),
            ('transl_table', 1),
            ('translation', 'HHHHHH'),
        ]
    }


def build_external_data_from_file(external_data_path: str) -> bytes:
    """
    Load pre-built EXTERNAL data from a binary file.

    This is used when we have a known-working EXTERNAL data structure
    that we want to reuse (e.g., extracted from a manually-created file).
    """
    with open(external_data_path, 'rb') as f:
        return f.read()


def build_external_data(
    old_sequence: str,
    backbone_packets: list,
    external_template_path: str = None,
) -> bytes:
    """
    Build the EXTERNAL reference data for a Type 11 manipulation packet.

    The EXTERNAL data contains:
    1. A 112-byte header with sequence metadata
    2. XZ-compressed packets (features, notes, etc.)

    Args:
        old_sequence: The sequence being stored (the backbone stuffer region)
        backbone_packets: List of (type, data) tuples from the backbone file
        external_template_path: Optional path to a pre-built EXTERNAL data file
    """
    # If we have a template file, use it directly
    if external_template_path and os.path.exists(external_template_path):
        return build_external_data_from_file(external_template_path)

    seq_len = len(old_sequence)

    # Build the 112-byte header
    # Copy the exact structure from a working manual file, only changing seq_len
    # Header format reverse-engineered from pJAG v2 s5-bLee-Hi.dna
    header = bytearray(112)

    # First 32 bytes: metadata (little-endian except seq_len which is big-endian)
    struct.pack_into('<I', header, 0, 31)         # 0x1f
    struct.pack_into('<I', header, 4, 102)        # 0x66
    struct.pack_into('<I', header, 8, 256)        # 0x100
    struct.pack_into('<I', header, 12, 23808)     # 0x5d00
    struct.pack_into('<I', header, 16, 142337)    # 0x22c01
    struct.pack_into('<I', header, 20, 65536)     # 0x10000
    struct.pack_into('<I', header, 24, 16777216)  # 0x1000000
    struct.pack_into('>I', header, 28, seq_len)   # sequence length (big-endian)

    # Bytes 32-111: Copy exact bytes from working file (80 bytes)
    # These appear to be binary data, possibly checksum or sequence encoding
    # Using the exact pattern from the working manual file
    manual_header_tail = bytes.fromhex(
        '6311f306a82d8469b41badf91beaa5b3b1d6966a7a6a82d7fd590a1ab67a0ada'
        'd559492ea846ead6699656679974f1e8ab6ab3162bff4e9bb131a556146a56d6'
        'e549998495e82e0f30bb251e00001694'
    )
    header[32:112] = manual_header_tail

    # Build the inner packet data - only include specific packet types
    # Based on reverse engineering: 16, 17, 27, 8, 10, 5, 6, 14, 28
    allowed_types = {16, 17, 27, 8, 10, 5, 6, 14, 28}
    inner_packets = []
    for ptype, pdata in backbone_packets:
        if ptype not in allowed_types:
            continue
        inner_packets.append(make_packet(ptype, pdata))

    inner_data = b''.join(inner_packets)

    # Compress the inner packet data
    inner_compressed = lzma.compress(inner_data)

    return bytes(header) + inner_compressed


def extract_type11_from_dna_file(dna_path: str) -> tuple:
    """
    Extract the XZ stream and EXTERNAL data from an existing .dna file's Type 11 packet.

    Returns:
        Tuple of (xz_stream_bytes, external_data_bytes, original_node_id)
    """
    with open(dna_path, 'rb') as f:
        data = f.read()

    for ptype, pdata in extract_all_packets(data):
        if ptype == 11:
            node_id = struct.unpack('>I', pdata[:4])[0]
            xz_data = pdata[9:]

            # Find XZ stream end (YZ marker)
            xz_end = None
            for i in range(len(xz_data) - 2):
                if xz_data[i:i+2] == b'YZ':
                    xz_end = i + 2
                    break

            if xz_end is None:
                raise ValueError(f"Could not find XZ end marker in {dna_path}")

            xz_stream = xz_data[:xz_end]
            external_data = xz_data[xz_end:]

            return xz_stream, external_data, node_id

    raise ValueError(f"No Type 11 packet found in {dna_path}")


def build_manipulation_packet(
    node_id: int,
    replace_start: int,
    replace_end: int,
    new_sequence: str,
    old_sequence: str,
    backbone_packets: list = None,
    template_dna_path: str = None,
) -> bytes:
    """
    Build a Type 11 manipulation packet for a replace operation.

    This packet is required for history nodes to be clickable/restorable.
    It contains the undo/redo actions to navigate between the parent and child states.

    Args:
        node_id: The history node ID this manipulation corresponds to (child node)
        replace_start: Start of replaced region (1-based, will be converted to 0-based)
        replace_end: End of replaced region (1-based)
        new_sequence: The new sequence inserted at replace_start (for Redo INSERT)
        old_sequence: The original sequence that was replaced (for Undo INSERT)
        backbone_packets: List of (type, data) tuples from backbone for EXTERNAL reference
        template_dna_path: Path to existing .dna file to extract XZ stream and EXTERNAL data from
    """
    # Build header (9 bytes):
    # Bytes 0-3: Node ID (big-endian u32)
    # Bytes 4-8: Flags (observed: 0x1d000002d4)
    header = struct.pack('>I', node_id) + bytes([0x1d, 0x00, 0x00, 0x02, 0xd4])

    # If we have a template file, use its XZ stream and EXTERNAL data
    # This ensures compatibility since Python's lzma produces different output than SnapGene's
    if template_dna_path and os.path.exists(template_dna_path):
        xz_stream, external_data, _ = extract_type11_from_dna_file(template_dna_path)
        return header + xz_stream + external_data

    # Fallback: build our own (may not work for resurrection)
    # Convert to 0-based positions (SnapGene Type 11 uses 0-based like history ColorRange)
    replace_start_0 = replace_start - 1

    # The new sequence end position (0-based, inclusive)
    new_end_0 = replace_start_0 + len(new_sequence) - 1
    old_end_0 = replace_start_0 + len(old_sequence) - 1

    # Build the Manipulation XML
    manipulation_xml = f'''<?xml version="1.0" encoding="UTF-8"?><Manipulation><AdditionalSequenceProperties><UpstreamStickiness>0</UpstreamStickiness><DownstreamStickiness>0</DownstreamStickiness><UpstreamModification>FivePrimePhosphorylated</UpstreamModification><DownstreamModification>FivePrimePhosphorylated</DownstreamModification></AdditionalSequenceProperties><Redo><Action type="REMOVE" range="{replace_start_0}-{old_end_0}"/><Action type="INSERT" position="{replace_start_0}"><Residues type="GENERIC" residues="{new_sequence}"/></Action></Redo><Undo><Action type="REMOVE" range="{replace_start_0}-{new_end_0}"/><Action type="INSERT" position="{replace_start_0}"><Residues type="EXTERNAL" length="{len(old_sequence)}" ID="0"/></Action></Undo></Manipulation>
'''

    # Compress the manipulation XML with LZMA
    compressed_xml = lzma.compress(manipulation_xml.encode('utf-8'))

    # Build the EXTERNAL data for the Undo action
    if backbone_packets is not None:
        external_data = build_external_data(old_sequence, backbone_packets)
    else:
        external_data = b''

    return header + compressed_xml + external_data


def build_history_tree(
    final_name: str,
    final_seq_len: int,
    is_circular: bool,
    backbone_name: str,
    backbone_seq_len: int,
    replace_start: int,
    replace_end: int,
    insert_length: int,
    existing_history: str = None,
    backbone_features_xml: str = None,
) -> tuple:
    """
    Build a History Tree XML showing a replacement operation.

    If existing_history is provided, wraps it as a child of the new replace operation.
    This preserves the backbone's existing history (e.g., "change origin" steps).

    Structure:
    - Root node: Final construct (operation="replace")
      - Child: The existing backbone history tree (with features added to root)

    Args:
        final_name: Name of the final construct (typically same as backbone)
        final_seq_len: Length of final sequence
        is_circular: Whether the final sequence is circular
        backbone_name: Name of the backbone vector
        replace_start: Start of replaced region in backbone (1-based)
        replace_end: End of replaced region in backbone (1-based)
        insert_length: Length of the new insert
        existing_history: The existing history XML from the backbone (optional)
        backbone_features_xml: The backbone's Type 10 features XML (needed for child node)
    """
    circular = "1" if is_circular else "0"

    # Calculate the insert end position in the final construct
    insert_end_final = replace_start + insert_length - 1

    # SnapGene history uses 0-based positions for val1/val2 in InputSummary
    replace_start_0 = replace_start - 1  # Convert to 0-based for InputSummary

    # Get backbone info from existing history if available
    if existing_history:
        match = re.search(r'<HistoryTree>(.*)</HistoryTree>', existing_history, re.DOTALL)
        if match:
            inner_history = match.group(1)
            existing_ids = [int(m) for m in re.findall(r'ID="(\d+)"', inner_history)]
            # Use ID one higher than any existing ID
            new_root_id = max(existing_ids) + 1 if existing_ids else 3
            child_id = new_root_id - 1
        else:
            new_root_id = 3
            child_id = 2
    else:
        new_root_id = 3
        child_id = 2

    # Build features content for child node (backbone features)
    features_content = ""
    if backbone_features_xml:
        backbone_features = extract_features_xml_content(backbone_features_xml, strip_recent_id=True)
        if backbone_features:
            features_content = f"<Features>{backbone_features}</Features>"

    # Build HistoryColors for child node (Input colors showing replaced region)
    child_history_colors = f'<HistoryColors><StrandColors type="Input"><TopStrand><ColorRange range="{replace_start_0}..{replace_end}" colors="white"/></TopStrand><BottomStrand><ColorRange range="{replace_start_0}..{replace_end}" colors="white"/></BottomStrand></StrandColors></HistoryColors>'

    # Create child node with resurrectable="2" to make it clickable/openable
    child_content = f'<Node name="{backbone_name}" type="DNA" seqLen="{backbone_seq_len}" strandedness="double" ID="{child_id}" circular="{circular}" resurrectable="2">{features_content}{child_history_colors}</Node>'

    history = f'''<?xml version="1.0" encoding="UTF-8"?><HistoryTree><Node name="{final_name}" type="DNA" seqLen="{final_seq_len}" strandedness="double" ID="{new_root_id}" circular="{circular}" operation="replace"><HistoryColors><StrandColors type="Product"><TopStrand><ColorRange range="{replace_start_0}..{insert_end_final}" colors="red"/></TopStrand><BottomStrand><ColorRange range="{replace_start_0}..{insert_end_final}" colors="red"/></BottomStrand></StrandColors></HistoryColors><InputSummary manipulation="replace" val1="{replace_start_0}" val2="{replace_end}"/>{child_content}</Node></HistoryTree>'''

    # Return both history XML and child_id (needed for Type 11 manipulation packet)
    return history, child_id


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
    signal_peptide: Optional[SecretionSignal] = None,
    output_path: str = None,
    type11_template_path: str = None,
) -> str:
    """
    Generate a SnapGene .dna file for a Twist construct.

    The history will show:
    - The backbone vector with its original features (parent node)
    - The replaced region highlighted
    - No separate fragment node (matches SnapGene's native behavior)

    The final construct will have:
    - Adjusted backbone features (shifted for insert length change)
    - New ORF feature for the inserted sequence

    Args:
        final_sequence: The complete final DNA sequence
        backbone_path: Path to the backbone .dna file
        insert_start: Start position of replaced region in backbone (1-based)
        insert_end: End position of replaced region in backbone (1-based)
        construct_name: Name for the final construct
        orf_name: Name for the ORF feature
        signal_peptide: Optional secretion signal to detect
        output_path: Output file path (default: construct_name.dna)
        type11_template_path: Path to existing .dna file to use as Type 11 template
                             (required for working history resurrection)
    """
    if output_path is None:
        output_path = f"{construct_name}.dna"

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

    # Check if we need to extend ORF to include a split stop codon
    # If the insert ends with partial stop codon (T, TA) and suffix completes it
    orf_sequence = insert_sequence
    suffix_start_in_final = len(prefix) + new_insert_length
    if suffix_start_in_final < len(final_sequence):
        # Check for split stop codons: TAA, TAG, TGA
        insert_end_2nt = insert_sequence[-2:].upper() if len(insert_sequence) >= 2 else ""
        insert_end_1nt = insert_sequence[-1:].upper() if len(insert_sequence) >= 1 else ""
        suffix_start_2nt = final_sequence[suffix_start_in_final:suffix_start_in_final+2].upper()
        suffix_start_1nt = final_sequence[suffix_start_in_final:suffix_start_in_final+1].upper()

        # Check for TA + A = TAA, TA + G = TAG, TG + A = TGA
        if insert_end_2nt == "TA" and suffix_start_1nt in ("A", "G"):
            orf_end += 1
            orf_sequence = insert_sequence + final_sequence[suffix_start_in_final:suffix_start_in_final+1]
        elif insert_end_2nt == "TG" and suffix_start_1nt == "A":
            orf_end += 1
            orf_sequence = insert_sequence + final_sequence[suffix_start_in_final:suffix_start_in_final+1]
        # Check for T + AA, T + AG, T + GA
        elif insert_end_1nt == "T" and suffix_start_2nt in ("AA", "AG", "GA"):
            orf_end += 2
            orf_sequence = insert_sequence + final_sequence[suffix_start_in_final:suffix_start_in_final+2]

    protein = translate(orf_sequence)

    orf_feature = create_orf_feature(
        name=orf_name,
        start=insert_start,
        end=orf_end,
        sequence=orf_sequence,
        signal_peptide=signal_peptide,
    )

    # Add ORF feature at the beginning of the feature list
    adjusted_features.insert(0, orf_feature)

    # Add separate signal peptide and His tag features if applicable
    if signal_peptide and protein.startswith(signal_peptide.aa_sequence):
        signal_len_nt = len(signal_peptide.aa_sequence) * 3
        signal_end = insert_start + signal_len_nt - 1

        # Add separate signal peptide feature
        sp_feature = create_signal_peptide_feature(
            name=f"{signal_peptide.short_name} signal peptide",
            start=insert_start,
            end=signal_end,
            sequence=insert_sequence[:signal_len_nt],
            orf_end=orf_end,
            color=signal_peptide.color,
        )
        adjusted_features.append(sp_feature)

    # Check for His tag at C-terminus
    if protein.rstrip('*').endswith(HIS_TAG_6X):
        his_len_nt = len(HIS_TAG_6X) * 3  # 18 nt
        # Account for stop codon if present
        stop_offset = 3 if protein.endswith('*') else 0
        his_end = orf_end - stop_offset
        his_start = his_end - his_len_nt + 1

        his_feature = create_his_tag_feature(
            start=his_start,
            end=his_end,
            orf_end=orf_end,
        )
        adjusted_features.append(his_feature)

    # Build final features XML
    final_features_xml = build_features_xml(adjusted_features)

    # Extract existing backbone history to preserve it
    backbone_history_compressed = extract_packet(backbone_data, 7)
    existing_history = None
    if backbone_history_compressed:
        existing_history = lzma.decompress(backbone_history_compressed).decode('utf-8')

    # Build history tree - wraps existing backbone history with our replace operation
    # The final construct uses the backbone name (SnapGene convention)
    history_xml, child_node_id = build_history_tree(
        final_name=backbone_name,  # SnapGene keeps the backbone name
        final_seq_len=len(final_sequence),
        is_circular=is_circular,
        backbone_name=backbone_name,
        backbone_seq_len=len(backbone_seq),
        replace_start=insert_start,
        replace_end=insert_end,
        insert_length=new_insert_length,
        existing_history=existing_history,
        backbone_features_xml=features_xml,  # Pass backbone features for child node
    )

    # Compress history
    history_compressed = lzma.compress(history_xml.encode('utf-8'))

    # Build Type 11 manipulation packet for the replace operation
    # This links to the child node and is required for it to be clickable/restorable
    #
    # IMPORTANT: Type 11 ranges are END-INCLUSIVE in 0-based coordinates
    # The REMOVE range will be "(insert_start-1)-(insert_end)" in 0-based
    # Example: insert_start=940, insert_end=1238 -> range "939-1238"
    # This covers 1238 - 939 + 1 = 300 positions (0-based indices 939..1238 inclusive)
    #
    # old_length = (insert_end) - (insert_start - 1) + 1 = insert_end - insert_start + 2
    # new_length = final_len - backbone_len + old_length
    type11_old_length = insert_end - (insert_start - 1) + 1  # = insert_end - insert_start + 2
    type11_new_length = len(final_sequence) - len(backbone_seq) + type11_old_length
    type11_new_sequence = final_sequence[insert_start-1:insert_start-1+type11_new_length]
    type11_old_sequence = backbone_seq[insert_start-1:insert_start-1+type11_old_length]

    # Get all backbone packets for the EXTERNAL reference in Type 11
    backbone_packets = extract_all_packets(backbone_data)

    manipulation_data = build_manipulation_packet(
        node_id=child_node_id,
        replace_start=insert_start,
        replace_end=insert_end,
        new_sequence=type11_new_sequence,
        old_sequence=type11_old_sequence,
        backbone_packets=backbone_packets,
        template_dna_path=type11_template_path,
    )

    # Build new file by copying ALL packets from backbone, replacing only what's needed
    # This preserves enzyme cache, display settings, and other UI state
    packets = []

    # Type 0 replacement: New sequence
    # Flags: bit 0 = circular, bit 1 = double-stranded
    # Must set double-stranded bit (0x02) or enzymes will be disabled
    flags = 0x03 if is_circular else 0x02  # Always double-stranded
    new_seq_data = bytes([flags]) + final_sequence.upper().encode('ascii')

    # Type 6 replacement: Updated notes with new name
    now = datetime.now(timezone.utc)
    new_notes_xml = f'''<Notes>
<UUID>{str(uuid.uuid4())}</UUID>
<Type>Synthetic</Type>
<ConfirmedExperimentally>0</ConfirmedExperimentally>
<Created UTC="{now.strftime('%H:%M:%S')}">{now.strftime('%Y.%m.%d')}</Created>
<LastModified UTC="{now.strftime('%H:%M:%S')}">{now.strftime('%Y.%m.%d')}</LastModified>
<SequenceClass>UNA</SequenceClass>
<TransformedInto>unspecified</TransformedInto>
<CustomMapLabel>{construct_name}</CustomMapLabel>
</Notes>'''

    # Iterate through all backbone packets and copy/replace as needed
    for packet_type, packet_data in extract_all_packets(backbone_data):
        if packet_type == 0:
            # Replace sequence
            packets.append(make_packet(0, new_seq_data))
        elif packet_type == 7:
            # Replace history with our new history
            packets.append(make_packet(7, history_compressed))
        elif packet_type == 10:
            # Replace features with adjusted features
            packets.append(make_packet(10, final_features_xml.encode('utf-8')))
        elif packet_type == 6:
            # Replace notes with new metadata
            packets.append(make_packet(6, new_notes_xml.encode('utf-8')))
        elif packet_type == 17:
            # Replace with empty AlignableSequences (keeping the packet seems important)
            empty_alignable = b'<AlignableSequences trimStringency="Medium"/>'
            packets.append(make_packet(17, empty_alignable))
        elif packet_type == 27:
            # Skip BAM alignment data - not relevant for new construct
            pass
        elif packet_type == 11:
            # Skip backbone's Type 11 manipulation packets - we'll add our own
            # The backbone's packets reference old node IDs and sequence positions
            pass
        else:
            # Copy all other packets unchanged (3, 5, 8, 9, 13, 14, 16, 28, 35, etc.)
            packets.append(make_packet(packet_type, packet_data))

    # Add our new Type 11 manipulation packet (for the replace operation)
    packets.append(make_packet(11, manipulation_data))

    # Write file
    with open(output_path, 'wb') as f:
        for packet in packets:
            f.write(packet)

    print(f"Written {output_path} ({sum(len(p) for p in packets)} bytes)")
    print(f"History shows: {backbone_name} with replaced region {insert_start}-{insert_end}")
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
    print("      insert_start=940,  # Where insert begins in backbone (1-based)")
    print("      insert_end=1238,   # Where stuffer ends in backbone (1-based)")
    print("      construct_name='My Construct',")
    print("      orf_name='Client Protein',")
    print("      signal_peptide=alpha_factor,  # Optional")
    print("  )")
    print()
    print("The history will show:")
    print("  - Backbone vector as parent (with all original features)")
    print("  - Replaced region highlighted in white (input) / red (product)")
    print("  - ORF features are added to the final construct, not shown in history")


if __name__ == '__main__':
    main()
