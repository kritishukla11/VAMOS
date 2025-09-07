# VAMOS

Developed by Kriti Shukla, kritis@unc.edu

# What this does:
Using this workflow, you can obtain whole proteome structural information from Alphafold, create dense variant clusters, identify which clusters are associated with a given pathway, and run statistical analysis on these clusters to identify likelihood of association.

# Instructions to run:
Data : All input data is publicly available on DepMap, TCGA, Alphafold, and MSigDB.

Versions: Originally ran on Python 3.8.8, Conda 22.9.0

# Known issues:
This workflow uses the kneed package for automatic finding of the epsilon parameter during density based clustering. If your input data does not have a clear "elbow" or "knee," this automatic process will fail and a manual epsilon parameter will need to be assigned for that protein. 

# Citations:
Please cite the following paper: Shukla K, Idanwekhai K, Naradikian M, Ting S, Schoenberger SP, Brunk E. Machine Learning of Three-Dimensional Protein Structures to Predict the Functional Impacts of Genome Variation. J Chem Inf Model. 2024 Jul 8;64(13):5328–43.


## Examples

Example datasets and tutorials are available in the `examples/` directory.

- `examples/tutorial_quickstart.ipynb`: Quick installation and usage.
- `examples/tutorial_full_pipeline.ipynb`: Full pipeline demo using toy data.
