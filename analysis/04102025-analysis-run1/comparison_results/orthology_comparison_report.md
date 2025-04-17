# SYNERGY2 vs YGOB/CGOB Orthology Comparison

**Date:** 2025-04-15

## Summary

This report compares orthology predictions between SYNERGY2 and YGOB/CGOB databases.

## Comparison Results

| Species Pair | YGOB Pairs | SYNERGY2 Pairs | True Positives | False Positives | False Negatives | Precision | Recall | F1 Score |
|-------------|------------|---------------|---------------|----------------|-----------------|-----------|--------|----------|
| C.albicans - S.cerevisiae | 4011 | 4343 | 3623 | 720 | 388 | 0.8342 | 0.9033 | 0.8674 |
| C.glabrata - K.lactis | 4494 | 1740 | 1509 | 231 | 2985 | 0.8672 | 0.3358 | 0.4841 |
| C.glabrata - S.cerevisiae | 4601 | 2191 | 1549 | 642 | 3052 | 0.7070 | 0.3367 | 0.4561 |
| K.lactis - S.cerevisiae | 4565 | 4700 | 4040 | 660 | 525 | 0.8596 | 0.8850 | 0.8721 |
| C.albicans - C.glabrata | 3602 | 1542 | 1205 | 337 | 2397 | 0.7815 | 0.3345 | 0.4685 |
| C.albicans - K.lactis | 3605 | 3988 | 3249 | 739 | 356 | 0.8147 | 0.9012 | 0.8558 |

## Interpretation

### Precision
Measures the fraction of SYNERGY2 ortholog predictions that are also found in YGOB/CGOB.
Higher precision indicates fewer false positive ortholog predictions.

### Recall
Measures the fraction of YGOB/CGOB orthologs that were correctly identified by SYNERGY2.
Higher recall indicates fewer false negative ortholog predictions.

### F1 Score
The harmonic mean of precision and recall, providing a single metric for overall performance.

## Conclusions

Determining whether SYNERGY2 results are reliable depends on the specific research question and threshold of accuracy required.
Generally, an F1 score above 0.7 is considered good, and above 0.8 is excellent for orthology prediction.

### Visualizations

Please refer to the generated PNG files for visual representation of these metrics.
