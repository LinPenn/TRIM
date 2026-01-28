# TRIM: Targeted Reduction of Interaction Magnitude via Interface-guided Mutation Design

## Software and Environment Dependencies

This pipeline integrates multiple external software tools and requires
a specific conda environment for execution.

### External Software

Please download and install the following software packages manually
according to their official documentation:

- **FoldX (v5.1)**  
  Apply for and download the latest FoldX 5.1 release from the official website:  
  https://foldxsuite.crg.eu/

- **Rosetta (v3.14)**  
  Download Rosetta 3.14 from the RosettaCommons official website:  
  https://rosettacommons.org/

After installation, please ensure that the installation paths of FoldX
and Rosetta are correctly specified in the pipeline script.

---

### Conda Environment

This pipeline requires a dedicated conda environment.

An automated environment setup script is provided in the working directory.
To create the required conda environment, run:

```
conda env create -f environment.yml
```

## Usage

Before running the pipeline, please open `pipeline.sh` and update the following paths
according to your local environment:

- The installation directories of the required software/tools
- The file path to the protein complex structure used as input
- Define the ligand chain (modified chain) and receptor chain.

After configuring the script, run the following commands to execute the pipeline:

```bash
chmod +x pipeline.sh
./pipeline.sh
