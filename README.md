# README
#### This directory contains all relevant code and results used in "Neuronal Identities Derived by Misexpression of the POU IV Sensory Determinant in a Proto-Vertebrate". Chacha et al., PNAS 2022.

![POU IV Misexpression](https://relevant-pou4-data.s3.us-east-2.amazonaws.com/fig5.png)

### Code
Is present in "original_rmarkdowns" directory.
All code was written by Prakriti Paul Chacha.

**1. "final_pou4_analysis.pipeline.Rmd"**
- Identifies all relevant WT and OE Epidermal Cell Types.
- Calculates Differentially Expressed Genes (DEGs) for each.
- Creates two "Transcriptomic Profiles" (250 and 300, refer Methods).
- Uses 250 Transcriptomic Profile to find Average Expression Levels for all WT and OE Cell Types.
    - These are inputs to Linear Model Analyses.
- Also finds human orthologs for OE DEGs.
    - These are inputs to OE DEG Analysis.

**Note:** When this project started, this code was written to be compatible with Seurat version 2.3.4.

**2. "final_pou4_revision_code.Rmd"**
- Contains additional supplemental code for new visualizations and analyses.

**Note:** When manuscript was under revision, code was written to be compatible with Seurat version 4.0.4. in order to accommodate certain visualization suggestions, namely the projection of OE subclusters using UMAP instead of TSNE embeddings (Fig 2c).

### Analysis Output Directories:
1. "DEGs"
2. "LM_inputs" 
3. "OE_DEG_Analysis"
4. "MEME_TOMTOM" 
Their contents are described in individual README files. 

### Respective directories containing images for knitted markdown files:
1. "final_pou4_analysis_pipeline_files"
2. "final_pou4_revision_code_files"

### Data
- Complete Data (raw and processed) for Control and Pou4 Misexpressed embryos is available using GEO Accession Number GSE192645.
  Please quote in any manuscript discussing the data.
- pou4_OE directory
	* Raw data from both embryos (not provided in GEO) was aggregated using the cellranger aggr function and csv file "pou4_agg.csv", which is also present in this directory.
	* Note: cellranger version 3.1.0 was used.
	```
	cellranger aggr --id=pou4_OE --csv=pou4_agg.csv
	```	   
	* All analyses were conducted using raw files in this directory.
	* Its Amazon S3 Bucket object URL is present in "original_rmarkdowns/final_pou4_analysis.pipeline.Rmd"

### Files 

**"ANISEED-Cirobu-GeneName-3bestBlastHitHuman.rnames"** contains khids and their human orthologs (compiled from the Aniseed Database). This file was kindly provided by Christelle Dantec of the Lemaire Lab, Institute of Biological Sciences, French National Centre for Scientific Research.
