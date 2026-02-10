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
        .select("chrom", "start", "end")
        .with_columns(
            name=pl.lit(name)
        )        
        .sort("chrom", "start", "end")
        .with_row_index("peak_id")
    )
    print(name)
    print(df)
    df.write_parquet(f"output/data/{name}.parquet")

    return df

def combine_tfs(name, *dfs):
    df = pl.concat(dfs)
    print(name)
    print(df)
    df.write_parquet(f"output/data/{name}.parquet")

ctcf = extract_tf("CTCF", "input/data/CTCF_ENCFF901CBP.bed.gz")
rad21 = extract_tf("RAD21", "input/data/RAD21_ENCFF439DYW.bed.gz")
smc3 = extract_tf("SMC3", "input/data/SMC3_ENCFF289LLT.bed.gz")
combine_tfs("cohesin", rad21, smc3)

# %%
