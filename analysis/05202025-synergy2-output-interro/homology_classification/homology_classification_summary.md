# SYNERGY2 Homology Classification Summary

## Overall Classification

| Category | Gene Count | Percentage |
|----------|------------|------------|
| pure_ortholog | 12704 | 90.30% |
| ohnolog | 1024 | 7.28% |
| paralog | 341 | 2.42% |

## Classification by Species

| Species | Pure Orthologs | Ohnologs | Paralogs | Total |
|---------|---------------|----------|----------|-------|
| C_albicans | 3176.0 | 0.0 | 131.0 | 3307.0 |
| C_glabrata | 3176.0 | 424.0 | 61.0 | 3661.0 |
| K_lactis | 3176.0 | 0.0 | 72.0 | 3248.0 |
| S_cerevisiae | 3176.0 | 600.0 | 73.0 | 3849.0 |
| unknown | 0.0 | 0.0 | 4.0 | 4.0 |

## Ohnolog Evidence

| Evidence Type | Count | Percentage |
|---------------|-------|------------|
| WGD_2:1_pattern | 1024 | 100.00% |

### Sample Ohnolog Pairs

| WGD Species | Ohnolog Pair | Non-WGD Orthologs |
|-------------|-------------|-------------------|
| S_cerevisiae | gene-YBR169C, gene-YPL106C | {'C_albicans': 'orf19.2435', 'K_lactis': 'gene:KLLA0_E24597g'} |
| S_cerevisiae | gene-YJL098W, gene-YKR028W | {'C_albicans': 'orf19.5160', 'K_lactis': 'gene:KLLA0_F14124g'} |
| C_glabrata | CAGL0C01199g, CAGL0F07865g | {'C_albicans': 'orf19.391', 'K_lactis': 'gene:KLLA0_A04169g'} |
| S_cerevisiae | gene-YLR228C, gene-YDR213W | {'C_albicans': 'orf19.391', 'K_lactis': 'gene:KLLA0_A04169g'} |
| S_cerevisiae | gene-YDR263C, gene-YOR033C | {'C_albicans': 'orf19.926', 'K_lactis': 'gene:KLLA0_E16743g'} |
| S_cerevisiae | gene-YGL135W, gene-YPL220W | {'C_albicans': 'orf19.3465', 'K_lactis': 'gene:KLLA0_B02002g'} |
| C_glabrata | CAGL0I04312g, CAGL0J01870g | {'C_albicans': 'orf19.7089', 'K_lactis': 'gene:KLLA0_A03157g'} |
| S_cerevisiae | gene-YLR287C-A, gene-YOR182C | {'C_albicans': 'orf19.4375.1', 'K_lactis': 'gene:KLLA0_C04809g'} |
| S_cerevisiae | gene-YGL225W, gene-YER039C | {'C_albicans': 'orf19.1232', 'K_lactis': 'gene:KLLA0_A01364g'} |
| C_glabrata | CAGL0D01782g, CAGL0K02101g | {'C_albicans': 'orf19.6760', 'K_lactis': 'gene:KLLA0_F12452g'} |

## Conclusion

This analysis identified three distinct classes of homology relationships in the SYNERGY2 output:

1. **Pure Orthologs**: 12704 genes (90.30%) with one-to-one orthologous relationships across species
2. **Ohnologs (WGD-derived paralogs)**: 1024 genes (7.28%) derived from whole genome duplication
3. **Other Paralogs**: 341 genes (2.42%) resulting from small-scale duplication events

These classifications provide insights into the evolutionary history of genes across the analyzed yeast species, highlighting both conservation (orthologs) and diversification (paralogs and ohnologs).
