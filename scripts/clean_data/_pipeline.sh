mkdir -p output/data
rm raw
ln -s ../../raw .
python 01_extract_loops.py
python 02_extract_anchors.py
python 03_extract_ensembl_annots.py
python 04_extract_transcript_quant.py
python 05_extract_enhancers.py
python 06_extract_peaks.py
python 07_anchor_neighbors.py
python 08_anchor_categories.py
python 09_loop_categories.py
