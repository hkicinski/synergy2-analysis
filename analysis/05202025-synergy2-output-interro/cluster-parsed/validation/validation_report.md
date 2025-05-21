# SYNERGY2 Parsing Validation Report

Generated: 2025-05-20 14:28:14.442629

## File Existence Check

- final_clusters.txt: ✓
- parsed_clusters.csv: ✓
- genes_by_cluster.csv: ✓

## Independent Original File Analysis

Total lines in original file: 9563

Original file metrics:
- Total clusters: 9563
- Total genes: 23602
- Average genes per cluster: 2.47

Sample lines from original file:
```
Cluster1 (taxa: 1, genes: 1)	gene-YNCM0023C
```
```
Cluster2 (taxa: 4, genes: 4)	orf19.6237 CAGL0J11242g gene-YNL180C gene:KLLA0_D08327g
```
```
Cluster3 (taxa: 1, genes: 1)	gene:KLLA0_F22671g
```
```
Cluster4 (taxa: 1, genes: 1)	gene-YNCL0036C
```
```
Cluster5 (taxa: 1, genes: 1)	gene-YNCJ0015W
```

## Parsed Data Analysis

Parsed data metrics:
- Total clusters: 9563
- Total genes: 23602
- Average genes per cluster: 2.47

## Count Validation

- Cluster count match: ✓ (Original: 9563, Parsed: 9563)
- Gene count match: ✓ (Original: 23602, Parsed: 23602)

## Cluster Size Distribution

### Size Distribution Comparison

| Size | Original Count | Parsed Count |
|------|---------------|-------------|
| 1 | 4259 | 4259 |
| 2 | 726 | 726 |
| 3 | 960 | 960 |
| 4 | 3249 | 3249 |
| 5 | 225 | 225 |
| 6 | 125 | 125 |
| 7 | 14 | 14 |
| 8 | 3 | 3 |
| 9 | 2 | 2 |

**Single-gene clusters**:
- Original: 4259 (44.54%)
- Parsed: 4259 (44.54%)

![Size Distribution Comparison](validation/size_distribution_comparison.png)

## Gene-to-Cluster Ratio

- Original ratio: 2.47 genes per cluster
- Parsed ratio: 2.47 genes per cluster
- Ratio match: ✓

## Validation Verdict

The parsed results match the expected patterns from the original data.

This validation was performed by directly comparing the parsed results against an independent analysis of the original file.
