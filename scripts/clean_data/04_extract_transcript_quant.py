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
    )

rep1 = get_replicate("input/data/transcript_quant_rep1_ENCFF190NFH.tsv", 1)
rep2 = get_replicate("input/data/transcript_quant_rep2_ENCFF461FLA.tsv", 2)
transcript_annot_ensembl = (
    pl.read_parquet("output/data/transcript_annot_ensembl.parquet")
    .drop("gene_id")
)

df = (
    rep1.join(rep2, on="transcript_id")
    .with_columns(
        TPM_arithm = (pl.col.TPM1 + pl.col.TPM2)/2,
        pme_TPM_arithm = (pl.col.pme_TPM1 + pl.col.pme_TPM2)/2,
        
        # Transcript quants have transcript_id format like 'ENST00000475007.5'.
        # Ensembl gene annotations are labeled like 'ENST00000475007'.
        # Extract transcript_id_base from transcript quants to enable join.
        transcript_id_base = pl.col.transcript_id.str.extract("(ENST[\d]+).(\d)", 1)
    )
    .join(transcript_annot_ensembl, left_on="transcript_id_base", right_on="transcript_id")

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

    .drop_nulls()
    .sort("chrom", "start", "end")
    .rename({"gene_id":"ensembl_gene_id", "transcript_id": "ensembl_transcript_id"})

    .with_row_index("transcript_id")
)
df.write_parquet("output/data/transcript_quant.parquet")
print(df)
# %%
