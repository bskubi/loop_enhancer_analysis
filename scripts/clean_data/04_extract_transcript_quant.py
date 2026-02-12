#%%
import polars as pl
import duckdb

# ENCODE K562 total RNA-seq
# 1. Compute arithmetic mean expression
# 2. Join with Ensembl transcript annotations
# 3. Drop nulls

def get_replicate(path, rep: int):
    df = duckdb.sql(f"SELECT * FROM '{path}'").pl()
    return (
        df
        .rename({
            col: f"{col}{rep}"
            for col in df.columns if col not in ["transcript_id", "gene_id"]
        })
        .with_columns(
            gene_id = pl.col.gene_id.str.extract(r"(ENSG[0-9]+)+\.[0-9]", 1),
            transcript_id = pl.col.transcript_id.str.extract(r"(ENST[0-9]+)+\.[0-9]", 1)
        )
    )

rep1 = get_replicate("raw/RNA/quant/transcript_quant_rep1_ENCFF190NFH.tsv", 1)
rep2 = get_replicate("raw/RNA/quant/transcript_quant_rep2_ENCFF461FLA.tsv", 2)
transcript_annot = (
    pl.read_parquet("output/data/transcript_annot.parquet")
    .drop("gene_id")
)
with pl.Config(set_tbl_cols=-1):
    print(transcript_annot.select("transcript_id").sort("transcript_id"))
    print(rep2.filter(pl.col.transcript_id.str.starts_with("ENST")).select("transcript_id").sort("transcript_id"))

df = (
    rep1.join(rep2, on="transcript_id")
    .with_columns(
        TPM_arithm = (pl.col.TPM1 + pl.col.TPM2)/2,
        pme_TPM_arithm = (pl.col.pme_TPM1 + pl.col.pme_TPM2)/2,    )
    .join(transcript_annot, on = "transcript_id")
    # Label transcription start site (TSS) and end site (TES)
    .with_columns(
        tss = pl.when(pl.col.strand == pl.lit("+"))
            .then(pl.col.start)
            .otherwise(pl.col.end)
        ,
        tes = (pl.when(pl.col.strand == pl.lit("+"))
            .then(pl.col.end)
            .otherwise(pl.col.start))
    )

    # .drop_nulls()
    .sort("chrom", "start", "end")
    .rename({
        "gene_id":"gencode_gene_id", 
        "gene_id_full":"gencode_gene_id_full",
        "transcript_id": "gencode_transcript_id",
        "transcript_id_full":"gencode_transcript_id_full",
    })
    .drop("transcript_support_level")
    .with_row_index("transcript_id")
)

df.write_parquet("output/data/transcript_quant.parquet")
with pl.Config(set_tbl_cols=-1):
    print(df)
print(df.columns)
print(f"Of original {len(rep1)} (rep1) and {len(rep2)} (rep2) transcripts, got {len(df)} transcripts after integrating")
# %%
