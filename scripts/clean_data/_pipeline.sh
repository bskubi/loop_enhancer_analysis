mkdir -p output/data
rm raw
ln -s ../../raw .
python 01_extract_loops.py
python 02_extract_anchors.py
python 03_extract_transcript_annots.py
python 04_extract_transcript_quant.py
python 05_extract_enhancers.py
python 06_extract_peaks.py
python 07_anchor_neighbors.py
python 08_anchor_categories.py
python 09_loop_categories.py
python 10_fen1_fads_locus_b_loops.py
cd output/browser_tracks
ln -s ../../raw/ENCODE_Enhancers/FADS1:Active.hg38.CRISPRi_HCRFF.bed
ln -s ../../raw/ENCODE_Enhancers/FADS2:Active.hg38.CRISPRi_HCRFF.bed
ln -s ../../raw/ENCODE_Enhancers/FADS3:Active.hg38.CRISPRi_HCRFF.bed
ln -s ../../raw/ENCODE_Enhancers/FADS3:Repressive.hg38.CRISPRi_HCRFF.bed
ln -s ../../raw/ENCODE_Enhancers/FEN1:Active.hg38.CRISPRi_HCRFF.bed