#!/usr/bin/env python3
"""
Generate a SnapGene .dna file.

Packet types:
- Type 9: File header (required, must be first)
- Type 0: DNA sequence (required)
- Type 5: Primers XML
- Type 6: Notes XML
- Type 8: Additional Sequence Properties XML
- Type 10: Features XML
- Type 13: Display settings (binary)
"""

import struct
from datetime import datetime, timezone


def make_packet(packet_type: int, data: bytes) -> bytes:
    """Create a packet with type byte and big-endian length."""
    return bytes([packet_type]) + struct.pack('>I', len(data)) + data


def extract_packet(data: bytes, packet_type: int) -> bytes:
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


def make_display_settings(
    enzyme_set_name: str = "Unsaved Enzyme Set",
    features_on_circle: bool = True,
    show_enzymes: bool = True,
    circular_view: bool = True
) -> bytes:
    """
    Create Type 13 display settings packet data.

    Args:
        enzyme_set_name: Name of the enzyme set to display
        features_on_circle: If True, features display on the circle line.
                           If False, features display inside the circle.
        show_enzymes: If True, restriction enzymes are displayed.
        circular_view: If True, opens in circular map view. If False, linear view.
    """
    data = bytearray(345)

    # Header bytes - KEY DISPLAY SETTINGS
    data[0] = 0x00
    data[1] = 0x00
    data[2] = 0x01 if features_on_circle else 0x00  # Features on circle
    data[3] = 0x00
    data[4] = 0x01
    data[5] = 0x00
    data[6] = 0x00
    data[7] = 0x4b  # 'K'

    # Enzyme set name length and name
    name_bytes = enzyme_set_name.encode('ascii')[:63]  # Max 63 chars
    data[17] = len(name_bytes)
    data[18:18+len(name_bytes)] = name_bytes

    # Pattern at offset 126-127
    data[126] = 0x19
    data[127] = 0xff

    # Display settings section
    data[128] = 0x00
    data[129] = 0x64  # 100 decimal - possibly zoom percentage

    # View mode marker at position 134-137
    # 'T' = circular view (THO), 'C' = linear view (CHO)
    data[134] = 0x00
    data[135] = ord('T') if circular_view else ord('C')
    data[136:138] = b'HO'

    # Bytes after marker (from original pJAG)
    data[138] = 0x00
    data[139] = 0xff
    data[140] = 0xfe

    # Flag section (matched to pJAG v2.dna pattern)
    data[153] = 0x01
    data[154] = 0x01
    data[155] = 0x01 if show_enzymes else 0x00  # Show enzymes
    data[156] = 0x01
    data[157] = 0x00
    data[158] = 0x01
    data[159] = 0x00
    data[160] = 0x45  # 'E'

    # More settings
    data[169] = 0x01

    # Constant pattern section (from pJAG v2.dna - verified working)
    # Note: byte 172 must be 0x00, pattern starts at byte 173
    data[172] = 0x00
    data[173] = 0x01
    data[174] = 0x00
    data[175] = 0x00
    data[176] = 0x01
    data[177:181] = b'\xff\xff\xff\xff'
    data[181] = 0x01
    data[182] = 0x59  # 'Y'
    data[183] = 0x01
    data[184] = 0xf4
    data[185] = 0x01
    data[186] = 0x01
    data[187] = 0x3f  # '?'
    data[188] = 0x00
    data[189] = 0x50  # 'P'

    return bytes(data)


def generate_dna_file(
    sequence: str,
    is_circular: bool = True,
    features_xml: str = None,
    notes_xml: str = None,
    primers_xml: str = None,
    additional_props_xml: str = None,
    custom_enzymes_xml: str = None,
    enzyme_visibilities_xml: str = None,
    display_settings: bytes = None,
    features_on_circle: bool = True,
    show_enzymes: bool = True,
    enzyme_set_name: str = "Unsaved Enzyme Set",
    enzyme_names: list = None,
    output_path: str = "test.dna"
):
    """Generate a .dna file."""

    packets = []

    # Type 9: File header (MUST BE FIRST)
    header_data = b'SnapGene\x00\x01\x00\x0f\x00\x14'
    packets.append(make_packet(9, header_data))

    # Type 0: DNA sequence
    flags = 0x01 if is_circular else 0x00
    seq_data = bytes([flags]) + sequence.upper().encode('ascii')
    packets.append(make_packet(0, seq_data))

    # Type 8: Additional Sequence Properties (before features)
    if additional_props_xml:
        packets.append(make_packet(8, additional_props_xml.encode('utf-8')))

    # Type 10: Features XML
    if features_xml is None:
        features_xml = '<?xml version="1.0"?><Features nextValidID="0"/>'
    packets.append(make_packet(10, features_xml.encode('utf-8')))

    # Type 5: Primers XML
    if primers_xml is None:
        primers_xml = '<?xml version="1.0"?><Primers nextValidID="0"><HybridizationParams minContinuousMatchLen="15" allowMismatch="1" minMeltingTemperature="40" showAdditionalFivePrimeMatches="1" minimumFivePrimeAnnealing="15"/></Primers>'
    packets.append(make_packet(5, primers_xml.encode('utf-8')))

    # Type 6: Notes XML
    if notes_xml is None:
        now = datetime.now(timezone.utc)
        notes_xml = f'''<Notes>
<UUID>00000000-0000-0000-0000-000000000000</UUID>
<Type>Synthetic</Type>
<ConfirmedExperimentally>0</ConfirmedExperimentally>
<Created UTC="{now.strftime('%H:%M:%S')}">{now.strftime('%Y.%m.%d')}</Created>
<LastModified UTC="{now.strftime('%H:%M:%S')}">{now.strftime('%Y.%m.%d')}</LastModified>
<SequenceClass>UNA</SequenceClass>
<TransformedInto>unspecified</TransformedInto>
</Notes>'''
    packets.append(make_packet(6, notes_xml.encode('utf-8')))

    # Type 13: Display settings
    if display_settings is None:
        display_settings = make_display_settings(
            enzyme_set_name=enzyme_set_name,
            features_on_circle=features_on_circle,
            show_enzymes=show_enzymes,
            circular_view=is_circular
        )
    packets.append(make_packet(13, display_settings))

    # Type 28: Enzyme Visibilities (before Type 14)
    if enzyme_visibilities_xml is None:
        enzyme_visibilities_xml = '<?xml version="1.0"?><EnzymeVisibilities vals=""/>'
    packets.append(make_packet(28, enzyme_visibilities_xml.encode('utf-8')))

    # Type 14: Custom Enzyme Sets
    if custom_enzymes_xml is None and enzyme_names:
        enzymes_str = ",".join(enzyme_names)
        custom_enzymes_xml = f'<?xml version="1.0"?><CustomEnzymeSets><CustomEnzymeSet type="2" name="{enzyme_set_name}" enzymeNames="{enzymes_str}"/></CustomEnzymeSets>'
    if custom_enzymes_xml:
        packets.append(make_packet(14, custom_enzymes_xml.encode('utf-8')))

    # Write file
    with open(output_path, 'wb') as f:
        for packet in packets:
            f.write(packet)

    print(f"Written {output_path} ({sum(len(p) for p in packets)} bytes)")
    return output_path


def main():
    """Copy from original file, using exact XML content."""

    with open('pJAG v2.dna', 'rb') as f:
        orig_data = f.read()

    # Extract sequence from Type 0
    seq_packet = extract_packet(orig_data, 0)
    is_circular = bool(seq_packet[0] & 0x01)
    sequence = seq_packet[1:].decode('ascii')

    # Extract exact XML from original file
    features_xml = extract_packet(orig_data, 10).decode('utf-8')
    notes_xml = extract_packet(orig_data, 6).decode('utf-8')
    primers_xml = extract_packet(orig_data, 5).decode('utf-8')
    additional_props_xml = extract_packet(orig_data, 8).decode('utf-8')

    print(f"Extracted from original:")
    print(f"  Sequence: {len(sequence)} bp, circular={is_circular}")
    print(f"  Features XML: {len(features_xml)} bytes")
    print(f"  Notes XML: {len(notes_xml)} bytes")
    print(f"  Primers XML: {len(primers_xml)} bytes")
    print(f"  Additional Props XML: {len(additional_props_xml)} bytes")
    print()

    # Generate file with exact same XML, features on circle
    generate_dna_file(
        sequence=sequence,
        is_circular=is_circular,
        features_xml=features_xml,
        notes_xml=notes_xml,
        primers_xml=primers_xml,
        additional_props_xml=additional_props_xml,
        features_on_circle=True,  # Display features ON the circle line
        show_enzymes=True,
        enzyme_names=["PmeI"],  # Show only PmeI
        output_path='test.dna'
    )

    # Verify
    print("\nVerifying generated file...")
    from parse_snapgene import parse_snapgene
    test_result = parse_snapgene('test.dna')
    print(f"  Sequence length: {len(test_result['sequence'])} bp")
    print(f"  Circular: {test_result['is_circular']}")
    print(f"  Features: {len(test_result['features'])}")
    print(f"  Features on circle: True")


if __name__ == '__main__':
    main()
