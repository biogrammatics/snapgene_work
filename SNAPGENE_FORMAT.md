# SnapGene .dna File Format Specification

## Overview

SnapGene .dna files are binary files consisting of a series of **packets**. Each packet has a 5-byte header followed by variable-length data.

## Packet Structure

```
+--------+------------------+------------------+
| Type   | Length (4 bytes) | Data (variable)  |
| 1 byte | Big-endian u32   | Length bytes     |
+--------+------------------+------------------+
```

## Packet Types

### Type 9: File Header (Required, First Packet)
- **Purpose**: Identifies the file as SnapGene format
- **Content**:
  - Bytes 0-7: ASCII "SnapGene" magic string
  - Bytes 8-13: Version/format info (e.g., `0001000f0014`)

### Type 0: DNA Sequence (Required)
- **Purpose**: Contains the actual DNA sequence
- **Content**:
  - Byte 0: Flags
    - Bit 0 (0x01): Circular if set, linear if clear
  - Bytes 1+: ASCII DNA sequence (A, T, G, C, N, etc.)

### Type 3: Restriction Enzyme Database
- **Purpose**: List of all restriction enzyme recognition sequences
- **Content**:
  - Byte 0: Flags (usually 0x01)
  - Bytes 1-4: Unknown (possibly count or offset)
  - Bytes 5+: Comma-separated enzyme recognition sequences
    - Uses IUPAC ambiguity codes (N, W, S, R, Y, K, M, B, D, H, V)
    - Example: `GACNNNNNNGTC,CCAGGCCTGG,GATCCGGATC,...`

### Type 5: Primers XML
- **Purpose**: Primer definitions and hybridization parameters
- **Content**: XML document
```xml
<?xml version="1.0"?>
<Primers nextValidID="0">
  <HybridizationParams minContinuousMatchLen="15" allowMismatch="1"
    minMeltingTemperature="40" showAdditionalFivePrimeMatches="1"
    minimumFivePrimeAnnealing="15"/>
</Primers>
```

### Type 6: Notes XML
- **Purpose**: Plasmid metadata
- **Content**: XML document
```xml
<Notes>
  <UUID>7a9e4685-36b1-41af-a7f5-52a4681f4318</UUID>
  <Type>Synthetic</Type>
  <ConfirmedExperimentally>0</ConfirmedExperimentally>
  <Created UTC="23:49:54">2025.8.12</Created>
  <LastModified UTC="9:10:36">2025.9.15</LastModified>
  <SequenceClass>UNA</SequenceClass>
  <TransformedInto>unspecified</TransformedInto>
</Notes>
```

### Type 7: History Tree (XZ Compressed)
- **Purpose**: Undo/redo history with full feature snapshots
- **Content**: XZ-compressed XML document containing complete edit history
- **Compression**: LZMA/XZ (`\xfd7zXZ\x00` signature)

### Type 8: Additional Sequence Properties XML
- **Purpose**: Sequence end modifications
- **Content**: XML document
```xml
<AdditionalSequenceProperties>
  <UpstreamStickiness>0</UpstreamStickiness>
  <DownstreamStickiness>0</DownstreamStickiness>
  <UpstreamModification>FivePrimePhosphorylated</UpstreamModification>
  <DownstreamModification>FivePrimePhosphorylated</DownstreamModification>
</AdditionalSequenceProperties>
```

### Type 10: Features XML
- **Purpose**: Sequence annotations (genes, promoters, etc.)
- **Content**: XML document
```xml
<?xml version="1.0"?>
<Features nextValidID="9">
  <Feature recentID="0" name="G418 Resistance" directionality="1"
    type="CDS" translationMW="30717.01" hitsStopCodon="1"
    originalSequence="ATG..." isFavorite="1" detectionMode="exactProteinMatch">
    <Segment range="1980-2789" color="#00bfff" type="standard" translated="1"/>
    <Q name="codon_start"><V int="1"/></Q>
    <Q name="translation"><V text="MGKEKTHVS..."/></Q>
  </Feature>
  <!-- More features... -->
</Features>
```

**Feature attributes:**
- `directionality`: 1 = forward, 2 = reverse
- `type`: CDS, promoter, terminator, rep_origin, misc_feature, etc.
- `swappedSegmentNumbering`: Present for reverse-strand features
- `cleavageArrows`: Position for signal peptide cleavage

### Type 11: Manipulation Data (XZ Compressed with Header)
- **Purpose**: Undo/redo manipulation records
- **Structure**:
  - Bytes 0-3: Sequence index (big-endian u32)
  - Bytes 4-5: Unknown
  - Bytes 6-7: Unknown
  - Byte 8: Unknown
  - Bytes 9+: XZ-compressed XML
- **Content**: Manipulation actions (e.g., sequence shifts)

### Type 13: Display Settings (Binary)
- **Purpose**: View settings including enzyme set and feature display options
- **Structure**: Fixed-size binary structure (345 bytes)
- **Anchor**: "THO" marker at offset 134-136 (in most files)

#### Structure Overview

| Offset | Size | Description |
|--------|------|-------------|
| 0 | 1 | Unknown (values: 0x00, 0x01, 0x06) |
| 1 | 1 | Unknown (usually 0x00) |
| **2** | 1 | **Features on circle**: `0x01` = on circle, `0x00` = inside circle |
| 3 | 1 | Always 0x00 |
| 4 | 1 | Always 0x01 |
| 5-6 | 2 | Always 0x00 |
| 7 | 1 | Usually 0x4b ('K') |
| 8-16 | 9 | Zeros |
| 17 | 1 | Length of enzyme set name |
| 18-81 | 64 | Enzyme set name (null-padded ASCII string) |
| 82-125 | 44 | Unknown |
| 126-127 | 2 | Pattern `\x19\xff` (constant) |
| 128 | 1 | Unknown (constant `\x00`) |
| 129 | 1 | Value 100 (`\x64`) - possibly zoom percentage |
| 130-133 | 4 | Zeros |
| 134-137 | 4 | **View mode marker** (see HO Marker Variants below) |

#### Key Display Settings

| Offset | Values | Description |
|--------|--------|-------------|
| **2** | `0x01` / `0x00` | **Features on circle**: `0x01` = on circle line, `0x00` = inside circle |
| **135** | `T` / `C` / `\x00` | **View mode prefix**: `T` = circular view, `C` = linear view |
| **155** | `0x01` / `0x00` | **Show enzymes**: `0x01` = show, `0x00` = hide |

#### Other Settings (less certain)

| Offset | Values | Description |
|--------|--------|-------------|
| 153 | 0x01 (95%), 0x00 (5%) | Unknown display flag |
| 154 | 0x01 (86%), 0x00 (14%) | Unknown display flag |
| 156 | 0x01 / 0x00 | Unknown flag |
| 157 | 0x01 (53%), 0x00 (47%) | Unknown display option |
| 159 | 'E' (96%), 'I' (4%) | Unknown character flag |
| 168 | 0x01 (72%), 0x00 (28%) | Unknown display option |
| 169 | 0x01 (53%), 0x00 (47%) | Unknown display option |

#### Constant Byte Patterns (THO-relative)

| THO Offset | Abs Offset | Value | Description |
|------------|------------|-------|-------------|
| THO+38 | 172 | 0x01 | Unknown |
| THO+41 | 175 | 0x01 | Unknown |
| THO+42-45 | 176-179 | 0xffffffff | Possibly -1 or max value |
| THO+46 | 180 | 0x01 | Unknown |
| THO+47 | 181 | 'Y' (0x59) | Unknown marker |
| THO+48 | 182 | 0x01 | Unknown |
| THO+49 | 183 | 0xf4 (244) | Unknown |
| THO+50 | 184 | 0x01 | Unknown |
| THO+51 | 185 | 0x01 | Unknown |
| THO+52 | 186 | '?' (0x3f) | Unknown marker |
| THO+54 | 188 | 'P' (0x50) | Unknown marker |

**Example enzyme set names:**
- "Unsaved Enzyme Set"
- "Unique 6+ Cutters"
- "Unique Cutters"
- "6+ Cutters"
- "None"

#### HO Marker Variants

The marker at position 134-137 controls the **default view mode**:

| Bytes at 134-137 | Marker | Default View | Count |
|------------------|--------|--------------|-------|
| `\x00THO` | THO | **Circular view** | 337 |
| `\x00CHO` | CHO | **Linear view** | 396 |
| `\x00\x00HO` | HO | Default | 172 |

**Important**: Despite 'C' correlating 85% with circular topology in the analyzed files, 'C' (CHO) triggers **linear view** while 'T' (THO) triggers **circular view**. Use 'T' prefix for circular sequences that should display in circular view.

**Analysis notes** (from 928 files with Type 13):
- All Type 13 packets are exactly 345 bytes
- The "HO" portion is always at bytes 136-137
- Byte 135 controls view mode: 'T' = circular view, 'C' = linear view, 0x00 = default

### Type 14: Custom Enzyme Sets XML
- **Purpose**: User-defined enzyme sets with specific enzymes
- **Content**: XML document
```xml
<?xml version="1.0"?>
<CustomEnzymeSets>
  <CustomEnzymeSet type="2" name="Unsaved Enzyme Set" enzymeNames="PmeI"/>
</CustomEnzymeSets>
```

### Type 16: Unknown (Contains "BASE" marker)
- **Purpose**: Possibly base quality or modification data
- **Structure**: Binary with embedded "BASE" marker and zlib-compressed data

### Type 17: Alignable Sequences XML
- **Purpose**: Sequence alignment/trace data references
- **Content**: XML document
```xml
<AlignableSequences trimStringency="Medium">
  <Sequence name="pJAG-v2_consensus" ID="0" isTrace="1"
    sortOrder="0" trimmedRange="0..4840"/>
</AlignableSequences>
```

### Type 27: BAM Alignment Data (GZIP Compressed)
- **Purpose**: Sequencing trace alignment data
- **Content**: GZIP-compressed BAM format data
- **Signature**: `\x1f\x8b` (GZIP magic)

### Type 28: Enzyme Visibilities XML
- **Purpose**: Which enzymes are currently visible
- **Content**: XML document
```xml
<?xml version="1.0"?>
<EnzymeVisibilities vals=""/>
```

### Type 35: Unknown (XZ Compressed)
- **Purpose**: Unknown, decompresses to `[]` in this file
- **Content**: XZ-compressed data

## Packet Order (Observed)

1. Type 9: Header
2. Type 0: DNA Sequence
3. Type 3: Enzyme Database
4. Type 11: Manipulation Data (multiple)
5. Type 16: BASE data
6. Type 7: History Tree
7. Type 17: Alignable Sequences
8. Type 27: BAM Data
9. Type 8: Additional Sequence Properties
10. Type 10: Features
11. Type 5: Primers
12. Type 6: Notes
13. Type 13: Enzyme Display Settings
14. Type 35: Unknown
15. Type 28: Enzyme Visibilities
16. Type 14: Custom Enzyme Sets

## Compression Formats Used

- **XZ/LZMA**: Signature `\xfd7zXZ\x00` - Used for history, manipulations
- **GZIP**: Signature `\x1f\x8b` - Used for BAM alignment data
- **Zlib**: Used within some binary packets

## Minimum Required Packets for Valid File

For generating a basic .dna file:
1. **Type 9**: File header with "SnapGene" magic
2. **Type 0**: DNA sequence with topology flag
3. **Type 10**: Features XML (can be empty `<Features nextValidID="0"/>`)
4. **Type 6**: Notes XML (basic metadata)

Optional but recommended:
- **Type 5**: Primers XML
- **Type 8**: Sequence properties
- **Type 13**: Display settings (required for features to display on circle instead of inside)
- **Type 14**: Custom enzyme sets (for restriction site visualization)
- **Type 28**: Enzyme visibilities
