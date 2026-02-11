#%%
import duckdb
import polars as pl

# 1. Extract CTCF, RAD21, SMC3 peaks to parquet
# 2. Get union set of RAD21 and SMC3 peaks as "cohesin" peaks

def extract_tf(name, path):
    df = (
        duckdb.read_csv(
            path,
            names=["chrom", "start", "end"]
        )
        .pl()
        .with_columns(
            name=pl.lit(name),
            center = (pl.col.start + pl.col.end)//2
        )        
        .select("chrom", "start", "center", "end")

        .sort("chrom", "start", "end")
        .with_row_index(f"{name}_id")
    )
    print(name)
    print(df)
    df.write_parquet(f"output/data/{name}.parquet")

    return df

def combine_tfs(name, *dfs):
    shared_cols = set(dfs[0].columns)
    for df in dfs:
        shared_cols = shared_cols.intersection(df.columns)
    dfs = [df.select(shared_cols) for df in dfs]
    df = (
        pl.concat(dfs, how="vertical")
        .with_row_index(f"{name}_id")
    )
    print(name)
    print(df)
    df.write_parquet(f"output/data/{name}.parquet")

ctcf = extract_tf("CTCF", "raw/TF_ChIP/pooled_cons_IDR_peaks/CTCF_ENCFF901CBP.bed.gz")
rad21 = extract_tf("RAD21", "raw/TF_ChIP/pooled_cons_IDR_peaks/RAD21_ENCFF439DYW.bed.gz")
smc3 = extract_tf("SMC3", "raw/TF_ChIP/pooled_cons_IDR_peaks/SMC3_ENCFF289LLT.bed.gz")
combine_tfs("cohesin", rad21, smc3)

# %%
