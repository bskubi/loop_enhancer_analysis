mkdir -p input/data
mkdir -p output/data
rm raw
ln -s ../../raw
rm -f input/data/*
cd input/data/
ln -s ../../../clean_data/output/data/loop_categories.parquet
ln -s ../../raw/TF_ChIP/sig_pval/* ./
cd ../..