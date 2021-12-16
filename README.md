# README
##### This directory contains all relevant code and results used in "Novel Neuronal Identities Derived by Misexpression of the POU IV Sensory Determinant in a Proto-Vertebrate". Chacha et. al, PNAS 2022

Note: Data can be downloaded from **(Insert Link here)** 

### Code
**(1) "final_pou4_analysis.pipeline.Rmd"**
Routine: 
- Identifies all relevant WT and OE Epidermal Cell Types.
- Calculates Differentially Expressed Genes (DEGs) for each.
- Creates two "Transcriptomic Profiles" (250 and 300, refer Methods).
- Uses 250 Transcriptomic Profile to find Average Expression Levels for all WT and OE Cell Types.
    - These are inputs to Linear Model Analyses.
- Also finds human orthologs for OE DEGs.
    - These are inputs to OE DEG Analysis.

**Note:** When this project started, this code was written to be compatible with Seurat version 2.3.4.

**(2) "final_pou4_revision_code.Rmd"**
- Contains additional supplemental code for new visualizations and analyses.

**Note:** When manuscript was under revision, code was written to be compatible with Seurat version 4.0.4. in order to accomodate certain visualization suggestions, namely the projection of OE subclusters using UMAP instead of TSNE embeddings (Fig 2c).

### "pou4_data" directory
Only contains the raw CellRanger Count Files "matrix.mtx", "barcodes.tsv" and "genes.tsv", the necessary inputs to perform the subsequent analyses. 
Rest of raw and processed data can be accessed from [Insert GEO Link]

### "pou4_revision_code_relevant_objects" directory
Contains "oe_epi_seurat" and "wt_epi_cns9_seurat" Seurat objects that are used in "final_pou4_revision_code.Rmd" 

### Contents of subdirectories "DEGs", "LM_inputs", "OE_DEG_Analysis", and "MEME_TOMTOM" are described in individual README files. 

### Files 

**"ANISEED-Cirobu-GeneName-3bestBlastHitHuman.rnames"** contains khids and their human orthologs (compiled from the Aniseed Database). This file was kindly provided by Christelle Dantec of the Lemaire Lab, Institute of Biological Sciences, French National Centre for Scientific Research.
