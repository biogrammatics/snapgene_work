#!/usr/bin/env python3
"""
Enhanced SnapGene parser - extract sequence and features properly
"""

import struct
import xml.etree.ElementTree as ET

def parse_snapgene(filepath):
    """Parse a SnapGene .dna file"""

    with open(filepath, 'rb') as f:
        data = f.read()

    pos = 0
    result = {
        'sequence': None,
        'is_circular': False,
        'features': [],
        'notes': {},
        'primers': []
    }

    while pos < len(data):
        if pos + 5 > len(data):
            break

        packet_type = data[pos]
        length = struct.unpack('>I', data[pos+1:pos+5])[0]

        if pos + 5 + length > len(data):
            break

        pdata = data[pos+5:pos+5+length]

        if packet_type == 0:  # DNA sequence header + sequence
            flags = pdata[0]
            result['is_circular'] = bool(flags & 0x01)
            try:
                result['sequence'] = pdata[1:].decode('ascii')
            except:
                pass

        elif packet_type == 10:  # Features (XML)
            try:
                xml_str = pdata.decode('utf-8')
                root = ET.fromstring(xml_str)
                for feature in root.iter('Feature'):
                    feat_info = {
                        'type': feature.get('type'),
                        'name': feature.get('name'),
                        'directionality': feature.get('directionality'),
                        'segments': [],
                        'qualifiers': {}
                    }
                    for seg in feature.findall('.//Segment'):
                        feat_info['segments'].append({
                            'range': seg.get('range'),
                            'color': seg.get('color'),
                            'type': seg.get('type')
                        })
                    for qual in feature.findall('.//Q'):
                        qname = qual.get('name')
                        qval = qual.findtext('V')
                        if qval:
                            feat_info['qualifiers'][qname] = qval
                    result['features'].append(feat_info)
            except Exception as e:
                print(f"Error parsing features: {e}")

        elif packet_type == 6:  # Notes (XML)
            try:
                xml_str = pdata.decode('utf-8')
                root = ET.fromstring(xml_str)
                for child in root:
                    result['notes'][child.tag] = child.text
            except:
                pass

        pos += 5 + length

    return result

def main():
    import sys
    filepath = sys.argv[1] if len(sys.argv) > 1 else 'input.dna'
    result = parse_snapgene(filepath)

    print("=" * 70)
    print(f"PLASMID ANALYSIS: {filepath}")
    print("=" * 70)

    print(f"\n📋 BASIC INFORMATION")
    print(f"   Sequence Length: {len(result['sequence']) if result['sequence'] else 0} bp")
    print(f"   Topology: {'Circular' if result['is_circular'] else 'Linear'}")
    print(f"   Type: {result['notes'].get('Type', 'Unknown')}")

    print(f"\n🧬 FEATURES ({len(result['features'])} total)")
    print("-" * 70)

    for feat in result['features']:
        name = feat.get('name', 'Unnamed')
        ftype = feat.get('type', 'unknown')
        direction = feat.get('directionality', '')
        dir_symbol = '→' if direction == '1' else '←' if direction == '2' else '↔'

        locations = [seg['range'] for seg in feat.get('segments', []) if seg.get('range')]
        loc_str = ', '.join(locations) if locations else 'unknown'

        print(f"  {dir_symbol} {name} ({ftype}): {loc_str}")

    if result['sequence']:
        seq = result['sequence'].upper()
        gc = (seq.count('G') + seq.count('C')) / len(seq) * 100
        print(f"\n📊 GC Content: {gc:.1f}%")

if __name__ == '__main__':
    main()
    