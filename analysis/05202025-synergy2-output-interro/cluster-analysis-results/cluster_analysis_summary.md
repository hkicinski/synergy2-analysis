# SYNERGY2 Cluster Analysis Summary

## Dataset Overview

- Total genes analyzed: 23602
- Total clusters: 9563
- Average genes per cluster: 2.47

## Species Distribution

| Species | Gene Count | Percentage |
|---------|------------|------------|
| S_cerevisiae | 6413 | 27.17% |
| C_albicans | 6194 | 26.24% |
| C_glabrata | 5569 | 23.60% |
| K_lactis | 5076 | 21.51% |
| unknown | 350 | 1.48% |

## Cluster Classifications

| Classification | Count | Percentage |
|----------------|-------|------------|
| single_copy_ortholog | 3176 | 33.21% |
| C_albicans_specific | 1889 | 19.75% |
| partial_ortholog | 1726 | 18.05% |
| S_cerevisiae_specific | 1102 | 11.52% |
| C_glabrata_specific | 619 | 6.47% |
| K_lactis_specific | 435 | 4.55% |
| multi_copy_ortholog | 346 | 3.62% |
| other | 270 | 2.82% |

## Single-Copy Orthologs

Found 3176 clusters with exactly one gene from each species.
These represent 33.21% of all clusters.

## Ortholog Conservation

| Species | Percentage with Orthologs |
|---------|-------------------------|
| C_albicans | 69.50% |
| C_glabrata | 89.12% |
| K_lactis | 91.45% |
| S_cerevisiae | 83.02% |

## Gene Family Expansions

### C_albicans

- 1 clusters with gene expansions
- Top 5 expansions:
  - Cluster7547: 3 genes

### C_glabrata

- 5 clusters with gene expansions
- Top 5 expansions:
  - Cluster3662: 3 genes
  - Cluster6139: 3 genes
  - Cluster6220: 3 genes
  - Cluster7510: 3 genes
  - Cluster7547: 3 genes

### S_cerevisiae

- 5 clusters with gene expansions
- Top 5 expansions:
  - Cluster3662: 3 genes
  - Cluster427: 3 genes
  - Cluster7421: 3 genes
  - Cluster7510: 3 genes
  - Cluster9216: 3 genes


## Evolutionary Insights

- **Core genome**: 3176 single-copy orthologs (33.21%) represent the core conserved genome.
- **Species-specific genes**: 4045 clusters (42.30%) contain genes unique to one species.
- **Gene family expansions**: 346 clusters (3.62%) show evidence of gene duplication while maintaining representation in all species.
- **Partial orthologs**: 1726 clusters (18.05%) show gene loss in one or more species.
