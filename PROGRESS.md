# Twist to SnapGene Converter - Progress Notes

## Current Status (2025-01-19)

The `twist_to_snapgene.py` script generates SnapGene .dna files from Twist synthesis data.

### Working Features

1. **Enzyme toggle** - Now visible and functional
   - Fixed by setting double-stranded bit (0x02) in Type 0 sequence flags
   - Flags byte must be 0x03 for circular double-stranded DNA (not 0x01)
   - Also required keeping Type 17 (AlignableSequences) packet

2. **Features** - All backbone features are preserved and adjusted for insert length
   - ORF feature with split CDS (signal peptide + mature protein segments)
   - Separate signal peptide annotation (e.g., "s5 signal peptide")
   - 6xHis tag annotation when detected at C-terminus

3. **History structure** - Shows replace operation with parent/child nodes
   - Root node: Final construct with `operation="replace"`
   - Child node: Backbone vector with features and HistoryColors

### Known Issue: History Child Node Not Clickable

The child node in the history tree is not clickable/openable in SnapGene.

**What we've tried:**
- Adding `resurrectable="2"` attribute - didn't work
- Removing `operation` attribute from child - didn't work
- Matching exact structure from manually-created working file - didn't work
- Even when user manually deleted grandchild nodes in SnapGene, child still wasn't clickable

**Working example (pJAG v2 s5-bLee-Hi.dna) structure:**
```xml
<Node name="pJAG v2.dna" ... ID="5" operation="replace">
  <Node name="pJAG v2.dna" ... ID="4" resurrectable="2">
    <Features>...</Features>
    <HistoryColors>...</HistoryColors>
  </Node>
</Node>
```

**Our generated structure (test_denovo.dna):**
```xml
<Node name="pJAG v2.dna" ... ID="3" operation="replace">
  <Node name="pJAG v2.dna" ... ID="2">
    <Features>...</Features>
    <HistoryColors>...</HistoryColors>
  </Node>
</Node>
```

The structures look similar but behavior differs. May be related to:
- How the file was originally created (manual cloning vs script generation)
- Some other packet or metadata we're not aware of
- Internal SnapGene state not stored in the file

### Key Technical Details

**Type 0 (Sequence) flags byte:**
- Bit 0: Circular topology
- Bit 1: Double-stranded (REQUIRED for enzyme toggle)
- Use 0x03 for circular double-stranded DNA

**Packets we copy from backbone:**
- Type 3, 11, 16: Enzyme cache data
- Type 13: Display settings
- Type 14, 28: Custom enzyme sets and visibilities
- Type 35: Unknown (XZ compressed)

**Packets we replace:**
- Type 0: New sequence with correct flags
- Type 6: Notes with new name/UUID
- Type 7: History tree (XZ compressed)
- Type 10: Features XML
- Type 17: Empty AlignableSequences

**Packets we skip:**
- Type 27: BAM alignment data

### Files

- `twist_to_snapgene.py` - Main converter script
- `pJAG v2.dna` - Backbone vector (4841 bp)
- `pJAG v2 s5-bLee-Hi.dna` - Manually created working example
- `test_denovo.dna` - Generated test output
- `SNAPGENE_FORMAT.md` - Binary format documentation
