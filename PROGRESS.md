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
   - Child node: Backbone vector with features, HistoryColors, and `resurrectable="2"`

4. **Type 11 Manipulation packet** - Required for clickable history nodes
   - Each restorable history node needs a corresponding Type 11 packet
   - Type 11 contains undo/redo actions linking node states
   - The packet's node_id (bytes 0-3 of header) must match the child node's ID

### Type 11 Manipulation Packet Structure

The Type 11 packet is critical for making history nodes clickable/restorable:

```
Header (9 bytes):
  Bytes 0-3: Node ID (big-endian u32) - must match history node ID
  Bytes 4-8: Unknown flags (observed: 0x1d000002d4)

Content (XZ compressed XML):
<Manipulation>
  <AdditionalSequenceProperties>...</AdditionalSequenceProperties>
  <Redo>
    <Action type="REMOVE" range="start-end"/>
    <Action type="INSERT" position="start">
      <Residues type="GENERIC" residues="SEQUENCE"/>
    </Action>
  </Redo>
  <Undo>
    <Action type="REMOVE" range="start-new_end"/>
    <Action type="INSERT" position="start">
      <Residues type="EXTERNAL" length="old_length" ID="0"/>
    </Action>
  </Undo>
</Manipulation>
```

- **Redo**: Actions to transform from child state to root state
- **Undo**: Actions to transform from root state back to child state
- `GENERIC` residues contain the actual sequence data
- `EXTERNAL` residues reference externally stored sequences

### Key Technical Details

**Type 0 (Sequence) flags byte:**
- Bit 0: Circular topology
- Bit 1: Double-stranded (REQUIRED for enzyme toggle)
- Use 0x03 for circular double-stranded DNA

**Packets we copy from backbone:**
- Type 3, 16: Enzyme cache data
- Type 13: Display settings
- Type 14, 28: Custom enzyme sets and visibilities
- Type 35: Unknown (XZ compressed)

**Packets we replace:**
- Type 0: New sequence with correct flags
- Type 6: Notes with new name/UUID
- Type 7: History tree (XZ compressed)
- Type 10: Features XML
- Type 11: New manipulation packet for replace operation
- Type 17: Empty AlignableSequences

**Packets we skip:**
- Type 27: BAM alignment data
- Backbone's Type 11 packets (they reference old node IDs)

### Files

- `twist_to_snapgene.py` - Main converter script
- `pJAG v2.dna` - Backbone vector (4841 bp)
- `pJAG v2 s5-bLee-Hi.dna` - Manually created working example
- `test_denovo.dna` - Generated test output
- `SNAPGENE_FORMAT.md` - Binary format documentation

### Testing Status

**Working** - Child node is clickable and resurrection correctly restores the backbone with all features.

### Critical Limitation: Type 11 Template Requirement

**Python's lzma library produces XZ streams that SnapGene cannot correctly parse for resurrection.**

Even with identical XML content, the compressed output differs from SnapGene's native compression, causing resurrection to fail (wrong sequence, garbage appended, missing features).

**Solution**: Use a template file created manually in SnapGene.

The `type11_template_path` parameter in `generate_twist_construct()` extracts the XZ stream and EXTERNAL data from an existing .dna file and uses it verbatim (only updating the node ID in the header).

**Template requirements**:
- Must be a .dna file created by manually performing the same operation in SnapGene
- The operation must be at the same insertion site (e.g., replacing the 300 bp Stuffer/AarI region)
- The template defines what backbone will be restored on resurrection
- The template is **NOT** generic - it's tied to the specific operation and backbone

**Example usage**:
```python
generate_twist_construct(
    final_sequence=my_sequence,
    backbone_path="pJAG v2.dna",
    insert_start=1715,
    insert_end=2014,
    construct_name="My Construct",
    orf_name="MyProtein",
    signal_peptide=SecretionSignal.S5,
    output_path="output.dna",
    type11_template_path="pJAG v2 s5-bLee-Hi.dna"  # Manually created template
)
```

**Without a template**: History will show the replace operation, but clicking the child node will not properly restore the backbone sequence and features.

### Files

- `twist_to_snapgene.py` - Main converter script
- `pJAG v2.dna` - Backbone vector (4841 bp)
- `pJAG v2 s5-bLee-Hi.dna` - Working manual example (use as template)
- `test_denovo.dna` - Generated test output
- `SNAPGENE_FORMAT.md` - Binary format documentation
