# README
#### MEME/TOMTOM directory

### Files
(1) **"manual_oe_intervals.txt"**
This contains the exact intervals I used to get putative 5' cis-regulatory regions of Qualifying OE DEGs for subsequent MEME and TOMTOM analyses.

### Note:
Up to 2kb sequence upstream the transcriptional start site was analyzed. If the sequence overlapped a gene body of a neighboring gene, it was rejected (Unqualified OE DEG). This rendered a total of 33 sequences. (Qualifying OE DEG)

(2) **"manual_oe50_fasta.txt"**
- Fasta file that contains sequences of above intervals. 
- This was used as input to MEME.

(3) **"forkhead.csv"**
- Contains khids and human orthologs of Forkhead proteins that are significantly upregulated in OE Cells compared to WT Epidermal Cells.
- They are also present in Supplemental List 1.

### "MEME_classic_anr_5/15_3" subdirectory.
- Contains inputs and outputs of MEME  and subsequent TOMTOM analyses.
- Refer to README within subdirectory for details
