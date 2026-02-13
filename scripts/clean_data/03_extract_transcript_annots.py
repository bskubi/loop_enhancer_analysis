#%%
import polars as pl
import duckdb

sql = f"""

SELECT
    chrom, 
    "start", 
    "end", 
    strand, 
    feature, 
    regexp_extract(attribute, 'gene_name \"(.*?)\";', 1) AS gene_name,
    regexp_extract(attribute, 'gene_id \"(.*?)\\.[0-9]+\";', 1) AS gene_id,
    regexp_extract(attribute, 'gene_id \"(.*?\\.[0-9]+)\";', 1) AS gene_id_full,
    regexp_extract(attribute, 'transcript_id \"(.*?)\\.[0-9]+\";', 1) AS transcript_id,
    regexp_extract(attribute, 'transcript_id \"(.*?\\.[0-9]+)\";', 1) AS transcript_id_full,
    regexp_extract(attribute, 'transcript_name \"(.*?)\";', 1) AS transcript_name,
    regexp_extract(attribute, 'gene_type \"(.*?)\";', 1) AS gene_type,
    regexp_extract(attribute, 'transcript_type \"(.*?)\";', 1) AS transcript_type,
    TRY_CAST(regexp_extract(attribute, 'transcript_support_level \"(.*?)\";', 1) AS INTEGER) AS transcript_support_level,
    TRY_CAST(regexp_extract(attribute, 'exon_number (.*?);', 1) AS INTEGER) AS exon_number
    
FROM
read_csv(
    'raw/Genes/gencode.v49.annotation.gtf.gz',
    columns = {{'chrom': 'VARCHAR', 'source': 'VARCHAR', 'feature': 'VARCHAR', 'start': INTEGER, 'end': INTEGER, 'score': VARCHAR, 'strand': VARCHAR, 'frame': INTEGER, 'attribute': VARCHAR}},
    skip=5
)
"""
with duckdb.connect() as conn:
    ensembl_annots = conn.sql(sql)
    
    gene_annots = (
        ensembl_annots
        .filter("feature = 'gene'")
        .project("* EXCLUDE (exon_number, transcript_id, transcript_id_full, transcript_name, transcript_type, transcript_support_level)")
    )
    print("gene_annots")
    print(gene_annots)
    gene_annots.write_parquet("output/data/gene_annot.parquet")
    transcript_annots = (
        ensembl_annots
        .filter("feature = 'transcript'")
        .project("* EXCLUDE (exon_number)")
    )
    print("transcript_annots")
    print(transcript_annots)
    transcript_annots.write_parquet("output/data/transcript_annot.parquet")
    exon_annots = (
        ensembl_annots
        .filter("feature = 'exon'")
    )
    print("exon_annots")
    print(exon_annots)
    exon_annots.write_parquet("output/data/exon_annot.parquet")
    gene_annotation_types = ensembl_annots.project("feature").distinct().pl()

gene_annotation_types
# %%
