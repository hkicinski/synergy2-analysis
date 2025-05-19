# C. glabrata GFF3 Modification

This directory contains a modified GFF3 file for C. glabrata where exon coordinates have been updated to match their corresponding CDS coordinates. This eliminates UTRs (untranslated regions) from the annotation. This discordance between the Exon-feature coordinate and the CDS-feature coordinate induced a frameshift mutation in the internalized translation performed by the *SYNERGY2* algorithm. As a result, the coordinates were adjusted. The original Candida Genome Database **CBS 138** genome was preserved.

## Files

- `C_glabrata_CBS138_current_features_modified.gff3`: The modified GFF3 file where exon coordinates match CDS coordinates
- Original GFF3 file is located at: `/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/data/feature_inputs/cgla_utr/C_glabrata_CBS138_current_features.gff3`

## Modification Script

The Python script used for modification is located at:
`/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/scripts/modify_gff3.py`

### Usage

```
python modify_gff3.py input.gff3 output.gff3
```

### What the script does

1. Reads the input GFF3 file
2. For each gene in the file, finds its CDS feature and corresponding exon feature
3. Updates the exon's start/end coordinates to match the CDS coordinates
4. Keeps all other data (IDs, attributes, etc.) unchanged
5. Writes the modified GFF3 to the output file

### Statistics

- Total exons in original file: 5,934
- Total modified exons: 5,060 out of 5,858 exons
- All CAGL gene IDs are preserved
- 419 exons already had identical coordinates as their CDS features -
	- The script only modified exons whose coordinates differed from their
  	corresponding CDS coordinates. For these exons, the script finds the CDS
  	but doesn't modify them because they already match exactly.
- 334 exons had no corresponding CDS features - The primary reason here
  is that RNA genes (those with IDs ending in 'r' instead of 'g') don't
  have CDS features since they're non-coding. The script couldn't find CDS
  coordinates to apply to these exons.
- Other unmapped exons (45) - These likely include special cases where
  the gene ID parsing didn't correctly match between exon and CDS.


## Purpose

This modification ensures that exon coordinates exactly match CDS coordinates, which is important for:
1. Eliminating UTRs from gene models
2. Consistent coordinate mapping between transcript and protein sequences
3. Improved compatibility with tools that expect exon and CDS coordinates to match