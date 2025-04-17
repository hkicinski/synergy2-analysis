# Synergy2 _E. coli_ Benchmark Dataset: Setup and Execution Guide 

**Author:** Hubert Kicinski 

**Affiliation:** Dr. Bin Z He Lab @ The University of Iowa 

**Date:** April 2025 

**Contact:** hkicinski@uiowa.edu

This README details instructions for downloading, installing, and running Synergy2 with the benchmark *E. coli* dataset. While this guide will focus on benchmark execution, the general protocol of execution to customized input data remains the same. 

## Overview 
Synergy2 is a computational tool that resolves gene ancestry by leveraging both homology and synteny in the context of a species phylogeny. Unlike simple ortholog clustering tools, Synergy2 reconstructs the complete evolutionary history of all genes across multiple genomes, identifying orthologs, paralogs, and even ohnologs through its phylogenetic resolution. The algorithm assembles "orthogroups"--sets of genes that descended from a single ancestral gene from their last common ancestral species--and creates gene trees that describe the history of species, duplication, and gene-loss events. 

This Tool is applicable to both prokaryotic and eukaryotic datasets of varying diversity. This guide focuses on using the provided benchmark dataset of 5 *E. coli* genomes to verify the installation and understand the basic workflow of the algorithm. Our current project focuses on analyzing the homology relationship between four yeast species (*C. albicans*, *C. glabrata*, *K. lactis*, and *S. cerevisiae*), but the workflow remains the same--even with customized input data. 

## Installation 

### Directory Structure 

Our project uses the following directory structure, which can be adapted for customized analyses and needs: 

```
/project_root/
├── data/                      
│   ├── genomes/              
│   ├── annotations/           
│   └── trees/                 
├── dependencies/              
│   ├── synergy2/              
│   ├── miniconda2/            
│   │   └── envs/synergy/     
│   └── workflow/              
├── analysis/                  
│   └── synergy_runs/          
└── output/                    
    ├── tables/                
    └── trees/                
```
### Directory Purpose 

- **data**/: contains all input data files for analysis 
  - **genomes**/: Genomic sequence files in FASTA format 
  - **annotations**/: Gene Annotation files in GFF3 format 
  - **trees**/: Phylogenetic trees in Newick format 
- **dependencies**/:Contains all software installations and environments 
  - **synergy2/**: Synergy2 software installation 
  - **miniconda2/**: Python environment manager 
    - **envs/synergy/**: Conda enviromment with all dependencies 
  - **workflow/**: Workflow package installation for job scheduling 

The remaining directories will become populated upon the execution of the `INSTALL.py` file for the SYNERGY2 algorithm. 

### Conda Environment Setup 

Since Synergy2 was developed in 2013, it requires specific older versions of dependencies that may conflict with modern software. At the time of writing this README document, the year is 2025; 12 years of software development strides have occured. Using a dedicated conda environment, as a result, is *strongly* recommended. 

#### Dependencies 

| **Software**       | **Version**          | **Installation Method** | **Purpose**                          |
|--------------------|----------------------|--------------------------|--------------------------------------|
| Python             | 2.7.18               | conda                    | Core language                        |
| NumPy              | 1.16.6               | conda                    | Numerical computing                  |
| SciPy              | 1.2.1                | conda                    | Scientific computing                 |
| NetworkX           | 1.8.1                | pip                      | Graph-based computations             |
| Java               | 1.8.0_412 (OpenJDK)  | conda                    | Required for WorkFlow                |
| BLAST-legacy       | 2.2.26               | conda (bioconda)         | Sequence alignment                   |
| MUSCLE             | 3.8.31               | manual install           | Multiple sequence alignment          |
| QuickTree          | 2.5                  | manual install           | Phylogenetic tree construction       |
| TIGR WorkFlow      | 3.2.0                | manual install           | Analysis workflow management         |

#### Setting up the Conda Environment 

To manage these dependencies, a dedicated conda environment is strongly recommended to avoid the, almost unavoidable, conflicts with modern software. 

```BASH
# Creates a conda environment directory if this doesn't exist yet 
mkdir -p dependencies/miniconda2 

# Download and install Miniconda (using the Python 2.7 version listed above)
wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh -b -p $PWD/dependencies/miniconda2
source dependencies/miniconda/bin/activate

# Create a dedicated environment for Synergy2
conda create -n synergy python-2.7 -y 
conda actiavte synergy

# Install the Python dependencies listed above
conda install -y numpy=1.16.6 scipy=1.2.1
pip install networkx==1.8.1

# And Finally, while you are doing downloading, Install Java and BLAST
conda install -y -c conda-forge openjdk=1.8.0.412
conda install -y -c bioconda blast-legacy=2.2.26
```
For all other dependencies, it is recommended to create a dedicated bin directory for manually installed tools: 

```bash
mkdir -p dependencies/bin
```

#### Installing MUSCLE 
One tool required for sequence alignment is MUSCLE; this is the MUlitple Sequence Comparison by Log-Expectation. This is the main tool, used in the SYNERGY2 algorithm, for multiple sequence alignment between the homologs of tested species. To download, the following steps were performed and are recommended: 

```bash
# Download and install MUSCLE
wget https://drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz
tar -xzf muscle3.8.31_i86linux64.tar.gz
mv muscle3.8.31_i86linux64 dependencies/bin/muscle
chmod +x dependencies/bin/muscle

# Add to PATH
export PATH=$PWD/dependencies/bin:$PATH

# Verify installation
muscle -version  # Should show MUSCLE v3.8.31
```
#### Installing QuickTree
This tool is used for the phylogenetic tree construction in the SYNERGY2 algorithm. 

```bash
# Download and install QuickTree
wget https://github.com/khowe/quicktree/archive/v2.5.tar.gz
tar -xzf v2.5.tar.gz
cd quicktree-2.5
make
cp quicktree ../dependencies/bin/
cd ..

# Verify installation
dependencies/bin/quicktree -v  # Should show QuickTree version 2.5
```

#### Installing TIGR WorkFlow 
The primary workflow management system used to control the execution of the analysis pipeline is the WorkFlow system. WorkFlow manages the entire computational workflow, coordinating multiple computational steps that must occur in discrete and specified orders. For each `.py` file, WorkFlow also ensures dependency tracking--avoiding the initiation of dependencies outside of the `.py` file in execution. To download, perform the following recommended steps: 

```bash
# Download WorkFlow
wget http://sourceforge.net/projects/tigr-workflow/files/latest/download -O workflow.tar.gz
tar -xzf workflow.tar.gz
mv workflow dependencies/workflow
cd dependencies/workflow

# Install WorkFlow
./configure --prefix=$PWD
make
make install
```
#### Installing SYNERGY2 
Now that all dependecies are installed, the remainder is the installation of the actual algorithm itself: Synergy2. For this installation, it is recommended to install the algorithm in its own directory: 

```bash 
mkdir dependencies/synergy2
```
 For installation, the following steps were performed and are recommended: 

```bash
# Download Synergy2
wget https://sourceforge.net/projects/synergy2/files/latest/download -O synergy2.tar.gz
tar -xzf synergy2.tar.gz
cd Synergy2

# Run the installation script
python INSTALL.py
```
The execution of `Install.py` will prompt for the full paths to the dependencies. When prompted, provide the full paths to: 

- QuickTree executable (`path/to/dependencies/bin/quicktree`)
- MUSCLE executable (`path/to/dependencies/bin/muscle`)
- BLAST executable (`path/to/dependencies/miniconda2/envs/synergy/bin`)

Thereafter, parameters for BLAST will be asked. Standard parameters were used, as recommanded by the `SYNERGY2-README.md`:

- Number of Cores for BLAST = 4 
- E-Value Threshold = 0.01

Now, to run SYNERGY2, you must first

```bash 
conda activate synergy
export PATH=$PWD/dependencies/bin:$PATH
source /path/to/workflow/exec_env.tcsh
```
This will set up the environment variables needed for the TIGR WorkFlow system to function properly. The conda activation ensures that all of the contained dependencies in your command-line environment. 

## Running the _E. coli_ Benchmark Dataset 

After successfully installing SYNERGY2 and its dependencies, the next step is to run the benchmark dataset! This section details how to configure and execute Synergy2 on the example _E. coli_ dataset that is made available upon the `INSTALL.py` execution.

### Preparing the Benchmark Dataset

The benchmark dataset consists of 5 _E. coli_ genomes with their genomic sequences and annotations. These files are included in the Synergy2 distribution under the `example/` directory. 

First: 

```bash
# Navigate to the Synergy2 example directory
cd dependencies/synergy2/example/

# Verify the contents
ls -la
# You should see directories for each E. coli genome, a data_catalog.txt file, and a tree.nwk file
```

Next, we must understand what the example directory contains: 

- Genome directories (Esch_coli_H296/, Esch_coli_H378_V1/, etc.)
- Each genome directory contains:
  - Genome sequence (.genome files)
  - Structural annotations (.annotation.gff3 files)
- data_catalog.txt: File listing paths to all genome sequences and annotations
- tree.nwk: Newick-formatted species tree for the 5 E. coli genomes

### Setting Up the Data Catalog

Before running SYNERGY2, we need to ensure the `data_catalog.txt` file contains the correct paths to all genome sequences and annotations. Since this is dependent on your installation directory, this `.txt` file must be updated. 

First, open the file with your favorite editor. For this guide, `vim` was used. Once in the example subdirectory, under your Synergy2 directory, open the file using the text editor of your choice: 

```bash 
vim data_catalog.txt
```

Next, manually update each path to reflect your absolute paths to the sequence and annotation files for the _E. coli_ data. 

```bash 
//
Genome  Esch_coli_H296
Sequence        /your/actual/path/to/Esch_coli_H296/Esch_coli_H296.genome
Annotation      /your/actual/path/to/Esch_coli_H296/Esch_coli_H296_PRODIGAL_2.annotation.gff3
//
```
Then, save the edits and exit the `data_catalog.txt` file. 

### Creating a Working Directory for the Benchmark Run 

Next, create a dedicated directory for the benchmark analysis run: 

```bash 
# Create and navigate to a working directory
mkdir -p analysis/synergy_runs/ecoli_benchmark
cd analysis/synergy_runs/ecoli_benchmark
```
### Running Synergy2 Initialization
Now, the next step is to initialize the Synergy2 workflow...

If you haven't already, first: 

```bash 
conda activate synergy
export PATH=$PWD/dependencies/bin:$PATH
source path/to/workflow/exec_env.tcsh
``` 
As this will initialize all of the dependencies in the conda environment, and initialize the WorkFlow scheduler. 

Next, you have to run the initialization script: 

```bash
../../bin/WF_runSynergy2.py -r /full/path/to/dependencies/synergy2/example/data_catalog.txt -w /full/path/to/analysis/synergy_runs/ecoli_benchmark -t /full/path/to/dependencies/synergy2/example/tree.nwk -f ecoli_benchmark
```

This will initialize the Synergy2 worjflow with: 

- **-r**: Path to the data catalog
- **-w**: Path to the working directory
- **-t**: Path to the species tree
- **-f**: Name prefix for the workflow files

After running this command, you'll see an output indicating that you need to execute commands in ```genomes/needed_extractions.cmd.txt```. 

```bash 
# Execute the extraction commands
bash genomes/needed_extractions.cmd.txt

# After the commands have finished, re-run the initialization script
../../bin/WF_runSynergy2.py \
  -r /full/path/to/dependencies/synergy2/example/data_catalog.txt \
  -w /full/path/to/analysis/synergy_runs/ecoli_benchmark \
  -t /full/path/to/dependencies/synergy2/example/tree.nwk \
  -f ecoli_benchmark
  ```

this second run will create the following WorkFlow files: 

- `ecoli_benchmark.xml`: WorkFlow template file 
- `ecoli_benchmark.ini`: Configuration file
- `ecoli_instance.xml`: Instance file 

### Launching the Workflow 
To be able to monitor the Synergy2 steps, use the `RunWorkflow` command: 

```bash 
RunWorkflow -t ecoli_benchmark.xml -c ecoli_benchmark.ini -i ecoli_instance.xml &
```
Finally, you have to launch the main set of SYNERGY2 jobs via the `instance_commands.txt` file. 

```bash 
# launch the jobs
bash instance_commands.txt
```
### Synergy2 Monitoring 
You can track the progress of the algorithm in 2 ways. Note: the more stable option is the Terminal-based monitoring option. 

#### GUI Monitoring 

```bash
MonitorWorkflow -i ecoli_instance.xml &
```
This assumes that X11 forwarding is enabled to run graphical applications on a Linux/Unix machine. This uses the X Window System (for GUI initialization) over the standard `SSH`.

#### Terminal-based Monitoring 
What was more reliable and persistently used was the Terminal-based monitoring system. 

In specific: 

```bash 
# Monitor with updates every 30 seconds
../../bin/WF_monitor.py ecoli_instance.xml 30
```
As the workflow continues, several steps will be achieved. 