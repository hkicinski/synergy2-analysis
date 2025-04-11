# Synergy2 Yeast Ortholog Analysis 

**Author:** Hubert Kicinski 

**Affiliation:** Dr. Bin Z He Lab @ The University of Iowa 

**Date:** April 10 2025 

**Contact:** hkicinski@uiowa.edu

## Overview 
This directory contains the results of running SYNERGY2 on four yeast species, performed on April 10, 2025: 

- *Candida albicans*
- *Candida glabrata*
- *Kluyveromyces lactis*
- *Saccharomyces cerevisiae*

The initial setup and configuration steps can be found in the `setup_README.md` file. 

## 1. Preparation of Analysis directories 
For this analysis, the following directory was created: 
`mkdir -p /space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/synergy_runs/yeast_analysis`

## 2. Data Preparation 
## 2.1 Data Catalog Creation 
As with the *E. coli* benchmark, a data catalog file was created to specify the locations of genome sequences (in FASTA format) and annotations (in GFF3 format) for the 4 species. 

The contents of the `data_catalog.txt` are as follows: 

```
//
Genome  Candida_albicans
Sequence        /space/hkicinski/Documents/E009-hk-Homology_SYNERGY/data/feature_inputs/calb/C_albicans_SC5314_version_A21-s02-m09-r10_chromosomes.fasta
Annotation      /space/hkicinski/Documents/E009-hk-Homology_SYNERGY/data/feature_inputs/calb/C_albicans_SC5314_version_A21-s02-m09-r10_features.gff3
//
Genome  Candida_glabrata
Sequence        /space/hkicinski/Documents/E009-hk-Homology_SYNERGY/data/feature_inputs/cgla/C_glabrata_CBS138_current_chromosomes.fasta
Annotation      /space/hkicinski/Documents/E009-hk-Homology_SYNERGY/data/feature_inputs/cgla/C_glabrata_CBS138_current_features.gff3
//
Genome  Kluyveromyces_lactis
Sequence        /space/hkicinski/Documents/E009-hk-Homology_SYNERGY/data/feature_inputs/klac_fixed/klac_fixed.fasta
Annotation      /space/hkicinski/Documents/E009-hk-Homology_SYNERGY/data/feature_inputs/klac_fixed/Kluyveromyces_lactis_gca_000002515.ASM251v1.60.gff3
//
Genome  Saccharomyces_cerevisiae
Sequence        /space/hkicinski/Documents/E009-hk-Homology_SYNERGY/data/feature_inputs/scer/GCF_000146045.2_R64_genomic.fasta
Annotation      /space/hkicinski/Documents/E009-hk-Homology_SYNERGY/data/feature_inputs/scer/S-cer-R64-s288c.gff3
//
```

These are the absolute path locations of each species. The `Genome` identifier refers to the specific identifier used in the Newick phylogenetic tree. 
## 2.2 Phylogenetic Tree Preparation 
Before running the analysis, a Newick-formatted phylogenetic tree file for the four yeast species was generated via a ShinyApp module. This tree was obtained using [Treehouse](https://doi.org/10.1186/s13104-019-4581-6), a Shiny-App tool for extracting subtrees from large phylogenies (Steenwyk & Rokas, 2019). The contents of this tree are as follows: 
```
(Candida_albicans:1.106029968,(Kluyveromyces_lactis:0.6057108267,(Saccharomyces_cerevisiae:0.3736038106,Candida_glabrata:0.4806222863):0.206574918):0.595976674);
``` 

For this analysis, the Newick tree was stored in the /data/ directory. Specifically: 

`/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/data/feature_inputs/03282025-hk-4sps-tree-branchl.nwk`

## 2.3 *K. lactis* Fix 
The *K. lactis* headers needed modification to match the Gff3 file. To do so, I created a modified FASTA file with simplified headers to match the A, B, C... chromosome identifiers listed in the original FASTA. 

```bash
sed -e 's/^>CR[0-9]*\.[0-9]* .* chromosome \([A-F]\) .*/>\1/' \
/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/data/feature_inputs/klac/GCA_000002515.1_ASM251v1_genomic.fasta \
> /space/hkicinski/Documents/E009-hk-Homology_SYNERGY/data/feature_inputs/klac_fixed/klac_fixed.fasta
```

The matching headers became:

```
>A
>B
>C
>D
>E
>F
```

# 3. Running Synergy2 
## 3.1 Setting Up the Enviroment 
To run the algorithm, the conda environment created, and detailed in the `setup_README.md`, was activated. All dependencies were exported; additionally, the WorkFlow environment was sourced for scheduling use. 

```bash
conda activate synergy
source /space/hkicinski/Documents/E009-hk-Homology_SYNERGY/dependencies/workflow/exec_env.bash.token
export PATH=/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/dependencies/synergy2/Synergy2-1.1/bin:$PATH
```
And the rest of the workflow steps were performed, as detailed in the `setup_README.md`.

# 4. Post-Processing Error and Resolution 
## 4.1 Error Encountered 
After the SYNERGY2 workflow was completed, an error was encountered during the post-processing stage when running 

```bash
/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/dependencies/synergy2/Synergy2-1.1/bin/ClusterPostProcessing.py genomes/ nodes/root/locus_mappings.pkl 4
```
Specifically, the error was 

```bash
total genes: 23602
Traceback (most recent call last):
  File "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/dependencies/synergy2/Synergy2-1.1/bin/ClusterPostProcessing.py", line 171, in <module>
    main(sys.argv[1:])
  File "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/dependencies/synergy2/Synergy2-1.1/bin/ClusterPostProcessing.py", line 117, in main
    tNum = int(t)
ValueError: invalid literal for int() with base 10: 'CAGL0H03861g'
```
After looking at the `ClusterPostProcessing.py` file, I realized that the original script was designed for the *E. coli* benchmark dataset where gene IDs are numeric. As we know, yeast species used gene IDs containing letters. This led to an error in processing where the script was attempting to convert these string Gene IDs into an integer via `int(t)`. 

## 4.3 Solution Implementation
1. I created a modified version of the post-processing script: 
```bash
cp /space/hkicinski/Documents/E009-hk-Homology_SYNERGY/dependencies/synergy2/Synergy2-1.1/bin/ClusterPostProcessing.py \
/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/dependencies/synergy2/Synergy2-1.1/bin/ClusterPostProcessing_yeast.py
```
2. Then, I modified the script to handle string IDs instead of integers by changing this original piece of the `ClusterPostProcessing.py` file: 

```bash
for t in transcripts:
    tNum = int(t)
    for s in transcripts:
        if s==t:
            continue
        sNum = int(s)
        tlist = [sNum,tNum]
        tlist.sort()
        tup = (tlist[0],tlist[1])
        mypairs.append(tup)
        pairs.add(tup)
```

Modified to: 

```bash
for t in transcripts:
    for s in transcripts:
        if s==t:
            continue
        tlist = [s, t]
        tlist.sort()  # Still sort them, but as strings
        tup = (tlist[0], tlist[1])
        mypairs.append(tup)
        pairs.add(tup)
```
Thereafter, I ran the modified script: 

```bash
/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/dependencies/synergy2/Synergy2-1.1/bin/ClusterPostProcessing_yeast.py genomes/ nodes/root/locus_mappings.pkl 4
```

# 5. Key results 
- Total genes: 23,602 
- Ortholog pairs identified: 19,092
- Single-copy core genes (scc): 1,219 
- Multi-copy core genes (mcc): 177
- Auxiliary genes (aux): 3,744
- Orphan genes: 7,390 
- Non-orphan clusters: 5,140

# 6. Conclusion 
Now that this has been completed, I will be analyzing the data, which will be documented in a separate `README.md` file. This detailed the execution of the algorithm on 04/10/2025 by Hubert Kicinski for the GRE Lab at the University of Iowa. 

# 7. References
Steenwyk, J. L., & Rokas, A. (2019). Treehouse: a user-friendly application to obtain subtrees from large phylogenies. *BMC Research Notes, 12*, 541. https://doi.org/10.1186/s13104-019-4581-6

Synergy 2.0. Available at: https://synergytwo.sourceforge.net/. Accessed April 10, 2025.
