Transcript quantifications: ENCODE K562 total RNA-seq
+ "Promoter" = transcription start site for a transcript

Ensembl annotations: Homo_sapiens.GRCh38.115.gtf.gz

1. Extract anchors from Aiden loops
2. Compute arithmetic mean of transcript quantifications
3. Join with Ensembl transcript annotations
4. Label transcript quantifications with TSS and TES
4. Drop nulls

Outputs:

01_extract_anchors.py: cleaned_loops.parquet
02_extract_anchors.py: anchors.parquet
03_extract_ensembl_annots.py: gene_annot_ensembl.parquet, transcript_annot_ensembl.parquet, exon_annot_ensembl.parquet, 
04_extract_transcript_quant.py: transcript_quant_annot.parquet
05_extract_enhancers.py: enhancers.parquet