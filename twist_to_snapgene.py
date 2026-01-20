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
import argparse
import json
from datetime import datetime, timezone
from dataclasses import dataclass
from typing import Optional, Union
from pathlib import Path
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


# =============================================================================
# Configuration Loading
# =============================================================================

# Path to external configuration file
CONFIG_FILE = Path(__file__).parent / "signals_and_tags.json"


def load_config() -> dict:
    """Load signals and tags from external JSON configuration file."""
    if not CONFIG_FILE.exists():
        raise FileNotFoundError(
            f"Configuration file not found: {CONFIG_FILE}\n"
            "This file should contain secretion_signals and tags definitions."
        )
    with open(CONFIG_FILE, 'r') as f:
        return json.load(f)


def load_secretion_signals() -> dict:
    """
    Load secretion signals from configuration file.

    Returns dict of signal_key -> SecretionSignal objects.
    Signals are sorted by sequence length (longest first) for proper detection.
    """
    config = load_config()
    signals = {}

    for key, data in config.get("secretion_signals", {}).items():
        if key.startswith("_"):  # Skip comment fields
            continue
        signals[key] = SecretionSignal(
            name=data["name"],
            aa_sequence=data["sequence"],
            short_name=data.get("short_name", key),
            color=data.get("color", "#3366ff"),
        )

    # Sort by sequence length (longest first) for proper detection
    sorted_signals = dict(
        sorted(signals.items(), key=lambda x: -len(x[1].aa_sequence))
    )
    return sorted_signals


def load_tags() -> dict:
    """Load affinity/detection tags from configuration file."""
    config = load_config()
    tags = {}

    for key, data in config.get("tags", {}).items():
        if key.startswith("_"):
            continue
        tags[key] = data["sequence"]

    return tags


def load_naming_rules() -> dict:
    """Load naming rules from configuration file."""
    config = load_config()
    rules = config.get("naming_rules", {})

    # Flatten the nested structure for easy lookup
    flat_rules = {}
    for category in ["signals", "codon_variants", "tags", "proteins"]:
        if category in rules:
            flat_rules.update(rules[category])

    return flat_rules


# Load configuration at module import time
SECRETION_SIGNALS = load_secretion_signals()
TAGS = load_tags()

# Backwards compatibility
HIS_TAG_6X = TAGS.get("6xHis", "HHHHHH")


@dataclass
class VectorConfig:
    """Configuration for a backbone vector."""
    name: str                    # Display name (e.g., "pJAG v2")
    backbone_path: str           # Path to backbone .dna file (relative to DATA_DIR)
    insert_start: int            # Start of insertion site (1-based)
    insert_end: int              # End of insertion site (1-based)
    template_path: str = None    # Path to Type 11 template .dna file (relative to DATA_DIR)
    default_signal: str = None   # Default secretion signal key (from SECRETION_SIGNALS)


# Default data directory for vector files
# Can be overridden by setting SNAPGENE_DATA_DIR environment variable
DATA_DIR = Path(os.environ.get("SNAPGENE_DATA_DIR", "~/snapgene_data")).expanduser()


def resolve_data_path(relative_path: str) -> Path:
    """Resolve a path relative to DATA_DIR, or return as-is if absolute."""
    path = Path(relative_path)
    if path.is_absolute():
        return path
    return DATA_DIR / path


# Registered vectors
# Add new vectors here with their insertion site coordinates and template files
VECTORS: dict[str, VectorConfig] = {
    "pJAG_v2": VectorConfig(
        name="pJAG v2",
        backbone_path="pJAG v2.dna",
        insert_start=1715,
        insert_end=2014,
        template_path="pJAG v2 s5-bLee-Hi.dna",
        default_signal="s5",
    ),
    # pSAN vector - needs configuration once template file is available
    # Uncomment and fill in insert_start, insert_end, and template_path
    # "pSAN": VectorConfig(
    #     name="pSAN",
    #     backbone_path="pSAN.dna",
    #     insert_start=1,     # TODO: Find BsaI insertion site start
    #     insert_end=1,       # TODO: Find BsaI insertion site end
    #     template_path="pSAN_template.dna",  # Manually created in SnapGene
    #     default_signal=None,  # No default signal (varies per construct)
    # ),
}


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


def build_manipulation_packet(node_id: int, template_dna_path: str) -> bytes:
    """
    Build a Type 11 manipulation packet for a replace operation.

    This packet is required for history nodes to be clickable/restorable.
    It contains the undo/redo actions to navigate between the parent and child states.

    NOTE: A template file (created manually in SnapGene) is REQUIRED because
    Python's lzma library produces XZ streams that SnapGene cannot correctly
    parse for resurrection.

    Args:
        node_id: The history node ID this manipulation corresponds to (child node)
        template_dna_path: Path to existing .dna file to extract XZ stream and EXTERNAL data from
    """
    if not template_dna_path or not os.path.exists(template_dna_path):
        raise ValueError(
            "A Type 11 template file is required for working history resurrection. "
            "Create a template by manually performing the same cloning operation in SnapGene."
        )

    # Build header (9 bytes):
    # Bytes 0-3: Node ID (big-endian u32)
    # Bytes 4-8: Flags (observed: 0x1d000002d4)
    header = struct.pack('>I', node_id) + bytes([0x1d, 0x00, 0x00, 0x02, 0xd4])

    # Extract XZ stream and EXTERNAL data from template
    # This ensures compatibility since Python's lzma produces different output than SnapGene's
    xz_stream, external_data, _ = extract_type11_from_dna_file(template_dna_path)
    return header + xz_stream + external_data


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
    manipulation_data = build_manipulation_packet(
        node_id=child_node_id,
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


def generate_from_vector(
    vector_name: str,
    final_sequence: str,
    construct_name: str,
    orf_name: str = "ORF",
    signal_peptide: Union[str, SecretionSignal, None] = None,
    output_path: str = None,
) -> str:
    """
    Generate a SnapGene .dna file using a registered vector configuration.

    This is the recommended high-level function for most use cases.

    Args:
        vector_name: Key from VECTORS registry (e.g., "pJAG_v2")
        final_sequence: The complete final DNA sequence
        construct_name: Name for the final construct
        orf_name: Name for the ORF feature
        signal_peptide: Secretion signal - can be:
            - None: Use vector's default signal (if any)
            - str: Key from SECRETION_SIGNALS (e.g., "s5", "alpha_factor")
            - SecretionSignal: Custom signal peptide object
        output_path: Output file path (default: construct_name.dna)

    Returns:
        Path to the generated .dna file

    Example:
        generate_from_vector(
            vector_name="pJAG_v2",
            final_sequence="ATGCAG...",
            construct_name="My Protein",
            orf_name="MyProtein",
        )
    """
    if vector_name not in VECTORS:
        available = ", ".join(VECTORS.keys())
        raise ValueError(f"Unknown vector '{vector_name}'. Available: {available}")

    vector = VECTORS[vector_name]

    # Resolve signal peptide
    if signal_peptide is None and vector.default_signal:
        signal_peptide = vector.default_signal

    if isinstance(signal_peptide, str):
        if signal_peptide not in SECRETION_SIGNALS:
            available = ", ".join(SECRETION_SIGNALS.keys())
            raise ValueError(f"Unknown signal peptide '{signal_peptide}'. Available: {available}")
        signal_peptide = SECRETION_SIGNALS[signal_peptide]

    # Resolve paths
    backbone_path = resolve_data_path(vector.backbone_path)
    template_path = resolve_data_path(vector.template_path) if vector.template_path else None

    return generate_twist_construct(
        final_sequence=final_sequence,
        backbone_path=str(backbone_path),
        insert_start=vector.insert_start,
        insert_end=vector.insert_end,
        construct_name=construct_name,
        orf_name=orf_name,
        signal_peptide=signal_peptide,
        output_path=output_path,
        type11_template_path=str(template_path) if template_path else None,
    )


# =============================================================================
# CSV Import and ORF Detection
# =============================================================================

@dataclass
class TwistConstruct:
    """Parsed construct from Twist CSV."""
    name: str                   # Twist name (e.g., "pSAN-s3bPhi_Ea1h-6xH")
    vector_name: str            # Vector name from CSV (e.g., "pSAN")
    insert_sequence: str        # Insert DNA sequence
    construct_sequence: str     # Full construct DNA sequence
    insert_length: int
    construct_length: int
    status: str                 # "delivered" or "failed"
    shipping_est: str           # Estimated shipping date

    # Parsed components (filled by analyze_construct)
    signal_type: str = None     # Detected signal type (e.g., "s1", "s3", "s5", None)
    codon_variant: str = None   # Codon variant (e.g., "b", "c")
    protein_name: str = None    # Protein name extracted from Twist name
    c_terminal_tag: str = None  # C-terminal tag (e.g., "6xH")
    display_name: str = None    # Human-readable name
    orf_name: str = None        # Name for ORF feature


# Load naming rules from configuration file
# These convert Twist-compatible names (no spaces) to human-readable display names
NAMING_RULES = load_naming_rules()


def detect_signal_peptide(protein_sequence: str) -> Optional[tuple]:
    """
    Detect which signal peptide is present at the N-terminus.

    Returns:
        Tuple of (signal_key, SecretionSignal) or None if no match
    """
    for key, signal in SECRETION_SIGNALS.items():
        if protein_sequence.startswith(signal.aa_sequence):
            return (key, signal)
    return None


def detect_his_tag(protein_sequence: str) -> bool:
    """Check if protein has C-terminal 6xHis tag."""
    return protein_sequence.rstrip('*').endswith(HIS_TAG_6X)


def validate_orf(dna_sequence: str) -> tuple:
    """
    Validate that sequence is a valid ORF.

    Returns:
        Tuple of (is_valid, protein_sequence, error_message)
    """
    seq = dna_sequence.upper()

    # Check for ATG start
    if not seq.startswith('ATG'):
        return (False, None, f"ORF does not start with ATG (starts with {seq[:3]})")

    # Translate
    protein = translate(seq)

    # Check for stop codon
    if '*' not in protein:
        return (False, protein, "ORF does not contain a stop codon")

    # Check stop is at end (allowing trailing sequence from vector)
    stop_pos = protein.index('*')
    if stop_pos < len(protein) - 1:
        # There's sequence after the stop - this is OK for inserts
        protein = protein[:stop_pos + 1]

    return (True, protein, None)


def parse_twist_name(name: str) -> dict:
    """
    Parse a Twist construct name into components.

    Examples:
        "pSAN-s3bPhi_Ea1h-6xH" -> {
            'vector': 'pSAN',
            'signal': 's3',
            'codon': 'b',
            'protein': 'Phi_Ea1h',
            'tag': '6xH'
        }
        "pSANna-bProA" -> {
            'vector': 'pSAN',
            'signal': None,
            'codon': 'b',
            'protein': 'ProA',
            'tag': None
        }
    """
    result = {
        'vector': None,
        'signal': None,
        'codon': None,
        'protein': None,
        'tag': None,
    }

    # Split by hyphens
    parts = name.split('-')
    if not parts:
        return result

    # First part is vector (may include 'na' for no signal)
    vector_part = parts[0]
    if vector_part.endswith('na'):
        result['vector'] = vector_part[:-2]
        result['signal'] = None
    else:
        result['vector'] = vector_part

    if len(parts) < 2:
        return result

    # Second part contains: signal code + codon variant + protein name
    middle = parts[1]

    # Extract signal type (s1, s2, s3, s5, etc.)
    signal_match = re.match(r'^(s\d+)([bc]?)(.*)$', middle)
    if signal_match:
        result['signal'] = signal_match.group(1)
        result['codon'] = signal_match.group(2) or 'b'
        result['protein'] = signal_match.group(3)
    else:
        # No signal prefix, might be "na" case or just codon+protein
        codon_match = re.match(r'^([bc])(.*)$', middle)
        if codon_match:
            result['codon'] = codon_match.group(1)
            result['protein'] = codon_match.group(2)
        else:
            result['protein'] = middle

    # Third part (if present) is C-terminal tag
    if len(parts) >= 3:
        result['tag'] = parts[2]

    return result


def twist_name_to_display(name: str, parsed: dict = None) -> str:
    """
    Convert Twist-compatible name to human-readable display name.

    Args:
        name: Twist construct name
        parsed: Optional pre-parsed name dict
    """
    if parsed is None:
        parsed = parse_twist_name(name)

    parts = []

    # Protein name (apply naming rules)
    protein = parsed.get('protein', '')
    if protein:
        # Apply protein name replacements
        display_protein = NAMING_RULES.get(protein, protein.replace('_', ' '))
        parts.append(display_protein)

    # Signal peptide suffix (only if notable)
    signal = parsed.get('signal')
    if signal and signal != 's5':  # s5 is common, don't suffix
        signal_name = NAMING_RULES.get(signal, signal)
        if signal_name:
            parts.append(f"({signal_name})")

    # Codon variant suffix
    codon = parsed.get('codon')
    if codon and codon != 'b':
        codon_suffix = NAMING_RULES.get(codon, f" {codon}")
        parts.append(codon_suffix.strip())

    # C-terminal tag
    tag = parsed.get('tag')
    if tag:
        tag_name = NAMING_RULES.get(tag, tag)
        parts.append(f"-{tag_name}")

    return ' '.join(parts).strip()


def analyze_construct(construct: TwistConstruct) -> TwistConstruct:
    """
    Analyze a Twist construct to detect signal peptides, tags, and validate ORF.

    Modifies construct in place and returns it.
    """
    # Parse the Twist name
    parsed = parse_twist_name(construct.name)
    construct.codon_variant = parsed.get('codon')
    construct.protein_name = parsed.get('protein')
    construct.c_terminal_tag = parsed.get('tag')

    # Validate the insert sequence as an ORF
    is_valid, protein, error = validate_orf(construct.insert_sequence)
    if not is_valid:
        print(f"Warning: {construct.name}: {error}")
        return construct

    # Detect signal peptide from sequence (more reliable than name parsing)
    signal_result = detect_signal_peptide(protein)
    if signal_result:
        signal_key, _ = signal_result
        construct.signal_type = signal_key
    else:
        # Check if name indicates signal that we didn't detect
        name_signal = parsed.get('signal')
        if name_signal:
            print(f"Warning: {construct.name}: Name indicates {name_signal} signal but not detected in sequence")
        construct.signal_type = None

    # Generate display name
    construct.display_name = twist_name_to_display(construct.name, parsed)

    # Generate ORF name (protein name with tag)
    orf_parts = [parsed.get('protein', 'ORF').replace('_', ' ')]
    if construct.c_terminal_tag:
        tag_name = NAMING_RULES.get(construct.c_terminal_tag, construct.c_terminal_tag)
        orf_parts.append(f"-{tag_name}")
    construct.orf_name = ''.join(orf_parts)

    return construct


def read_twist_csv(csv_path: str, skip_failed: bool = True) -> list:
    """
    Read constructs from a Twist CSV file.

    Args:
        csv_path: Path to the CSV file
        skip_failed: If True, skip constructs with status "failed"

    Returns:
        List of TwistConstruct objects
    """
    import csv

    constructs = []
    with open(csv_path, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            status = row.get('Step', '').lower()
            if skip_failed and status == 'failed':
                continue

            construct = TwistConstruct(
                name=row['Name'],
                vector_name=row['Vector name'],
                insert_sequence=row['Insert sequence'].upper(),
                construct_sequence=row['Construct sequence'].upper(),
                insert_length=int(row['Insert length']),
                construct_length=int(row['Construct length']),
                status=status,
                shipping_est=row.get('Shipping est', ''),
            )

            # Analyze the construct
            analyze_construct(construct)
            constructs.append(construct)

    return constructs


def process_twist_csv(
    csv_path: str,
    output_dir: str = None,
    skip_failed: bool = True,
    dry_run: bool = False,
) -> list:
    """
    Process a Twist CSV file and generate SnapGene .dna files.

    Args:
        csv_path: Path to the Twist CSV file
        output_dir: Directory for output files (default: same as CSV)
        skip_failed: Skip constructs marked as "failed"
        dry_run: If True, analyze but don't generate files

    Returns:
        List of (construct, output_path or None) tuples
    """
    if output_dir is None:
        output_dir = os.path.dirname(csv_path) or '.'

    constructs = read_twist_csv(csv_path, skip_failed=skip_failed)
    results = []

    print(f"Processing {len(constructs)} constructs from {csv_path}")
    print("=" * 60)

    for construct in constructs:
        print(f"\n{construct.name}:")
        print(f"  Vector: {construct.vector_name}")
        print(f"  Signal: {construct.signal_type or 'none'}")
        print(f"  Protein: {construct.protein_name}")
        print(f"  Tag: {construct.c_terminal_tag or 'none'}")
        print(f"  Display name: {construct.display_name}")
        print(f"  ORF name: {construct.orf_name}")
        print(f"  Insert length: {construct.insert_length} bp")

        if dry_run:
            results.append((construct, None))
            continue

        # Look up vector configuration
        # Try both exact match and normalized name
        vector_key = construct.vector_name.replace(' ', '_')
        if vector_key not in VECTORS:
            print(f"  ERROR: Vector '{construct.vector_name}' not in VECTORS registry")
            results.append((construct, None))
            continue

        vector = VECTORS[vector_key]

        # Resolve signal peptide
        signal = None
        if construct.signal_type and construct.signal_type in SECRETION_SIGNALS:
            signal = SECRETION_SIGNALS[construct.signal_type]

        # Generate output path
        output_filename = f"{construct.display_name or construct.name}.dna"
        output_path = os.path.join(output_dir, output_filename)

        try:
            generate_from_vector(
                vector_name=vector_key,
                final_sequence=construct.construct_sequence,
                construct_name=construct.display_name or construct.name,
                orf_name=construct.orf_name or "ORF",
                signal_peptide=signal,
                output_path=output_path,
            )
            print(f"  Generated: {output_path}")
            results.append((construct, output_path))
        except Exception as e:
            print(f"  ERROR: {e}")
            results.append((construct, None))

    return results


def list_vectors():
    """Print available vectors and their configurations."""
    print("Available vectors:")
    print("-" * 60)
    for key, vec in VECTORS.items():
        print(f"  {key}:")
        print(f"    Name: {vec.name}")
        print(f"    Backbone: {vec.backbone_path}")
        print(f"    Insert site: {vec.insert_start}-{vec.insert_end}")
        if vec.template_path:
            print(f"    Template: {vec.template_path}")
        if vec.default_signal:
            print(f"    Default signal: {vec.default_signal}")
        print()


def list_signals():
    """Print available secretion signals."""
    print(f"Configuration file: {CONFIG_FILE}")
    print()
    print("Available secretion signals:")
    print("-" * 60)
    for key, sig in SECRETION_SIGNALS.items():
        print(f"  {key}:")
        print(f"    Name: {sig.name}")
        print(f"    Sequence: {sig.aa_sequence}")
        print()
    print("Available tags:")
    print("-" * 60)
    for key, seq in TAGS.items():
        print(f"  {key}: {seq}")
    print()


def main():
    """Command-line interface for twist_to_snapgene."""
    parser = argparse.ArgumentParser(
        description="Convert Twist synthesis output to SnapGene .dna files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process a Twist CSV file (dry run first to check)
  %(prog)s --csv "Lytos FTE3.csv" --dry-run

  # Process CSV and generate .dna files
  %(prog)s --csv "Lytos FTE3.csv" --output-dir ./output

  # Generate a single construct using a registered vector
  %(prog)s --vector pJAG_v2 --sequence ATGCAG... --name "My Protein"

  # Read sequence from file
  %(prog)s --vector pJAG_v2 --sequence-file insert.txt --name "My Protein"

  # List available vectors
  %(prog)s --list-vectors

  # List available secretion signals
  %(prog)s --list-signals

Environment:
  SNAPGENE_DATA_DIR  Directory containing vector and template files
                     (default: ~/snapgene_data)
""",
    )

    # Info commands
    parser.add_argument("--list-vectors", action="store_true",
                        help="List available vector configurations")
    parser.add_argument("--list-signals", action="store_true",
                        help="List available secretion signals")

    # CSV import mode
    parser.add_argument("--csv", metavar="FILE",
                        help="Process a Twist CSV file")
    parser.add_argument("--output-dir", metavar="DIR",
                        help="Output directory for generated .dna files")
    parser.add_argument("--dry-run", action="store_true",
                        help="Analyze CSV without generating files")
    parser.add_argument("--include-failed", action="store_true",
                        help="Include failed constructs from CSV")

    # Single construct mode
    parser.add_argument("--vector", "-v", metavar="NAME",
                        help=f"Vector name (available: {', '.join(VECTORS.keys())})")
    parser.add_argument("--sequence", "-s", metavar="SEQ",
                        help="Final DNA sequence")
    parser.add_argument("--sequence-file", "-f", metavar="FILE",
                        help="File containing the final DNA sequence")
    parser.add_argument("--name", "-n", metavar="NAME",
                        help="Construct name")
    parser.add_argument("--orf-name", metavar="NAME", default="ORF",
                        help="ORF feature name (default: ORF)")
    parser.add_argument("--signal", metavar="NAME",
                        help=f"Secretion signal (available: {', '.join(SECRETION_SIGNALS.keys())})")
    parser.add_argument("--output", "-o", metavar="FILE",
                        help="Output file path (default: <name>.dna)")

    args = parser.parse_args()

    # Handle info commands
    if args.list_vectors:
        list_vectors()
        return

    if args.list_signals:
        list_signals()
        return

    # CSV import mode
    if args.csv:
        results = process_twist_csv(
            csv_path=args.csv,
            output_dir=args.output_dir,
            skip_failed=not args.include_failed,
            dry_run=args.dry_run,
        )
        # Print summary
        success = sum(1 for _, path in results if path is not None)
        failed = sum(1 for _, path in results if path is None)
        print(f"\n{'=' * 60}")
        if args.dry_run:
            print(f"Dry run complete: {len(results)} constructs analyzed")
        else:
            print(f"Generated: {success} files")
            if failed > 0:
                print(f"Failed: {failed} constructs")
        return

    # Single construct mode - validate required arguments
    if not args.vector:
        parser.error("--vector is required for generation")
    if not args.name:
        parser.error("--name is required for generation")
    if not args.sequence and not args.sequence_file:
        parser.error("Either --sequence or --sequence-file is required")

    # Get sequence
    if args.sequence_file:
        with open(args.sequence_file, 'r') as f:
            sequence = f.read()
    else:
        sequence = args.sequence

    # Clean sequence (remove whitespace, newlines)
    sequence = ''.join(sequence.split()).upper()

    # Generate
    try:
        output_path = generate_from_vector(
            vector_name=args.vector,
            final_sequence=sequence,
            construct_name=args.name,
            orf_name=args.orf_name,
            signal_peptide=args.signal,
            output_path=args.output,
        )
        print(f"Successfully generated: {output_path}")
    except Exception as e:
        print(f"Error: {e}")
        raise SystemExit(1)


if __name__ == '__main__':
    main()
