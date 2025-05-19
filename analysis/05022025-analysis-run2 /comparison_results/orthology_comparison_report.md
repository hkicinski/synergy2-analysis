# SYNERGY2 vs YGOB/CGOB Orthology Comparison

**Date:** 2025-05-02

## Summary

This report compares orthology predictions between SYNERGY2 and YGOB/CGOB databases.

## Comparison Results

| Species Pair | YGOB Pairs | SYNERGY2 Pairs | True Positives | False Positives | False Negatives | Precision | Recall | F1 Score |
|-------------|------------|---------------|---------------|----------------|-----------------|-----------|--------|----------|
| C.albicans - S.cerevisiae | 4011 | 4362 | 3636 | 726 | 375 | 0.8336 | 0.9065 | 0.8685 |
| C.glabrata - K.lactis | 4494 | 4413 | 3972 | 441 | 522 | 0.9001 | 0.8838 | 0.8919 |
| C.glabrata - S.cerevisiae | 4601 | 5324 | 4102 | 1222 | 499 | 0.7705 | 0.8915 | 0.8266 |
| K.lactis - S.cerevisiae | 4565 | 4709 | 4037 | 672 | 528 | 0.8573 | 0.8843 | 0.8706 |
| C.albicans - C.glabrata | 3602 | 4056 | 3261 | 795 | 341 | 0.8040 | 0.9053 | 0.8517 |
| C.albicans - K.lactis | 3605 | 3996 | 3253 | 743 | 352 | 0.8141 | 0.9024 | 0.8559 |

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
