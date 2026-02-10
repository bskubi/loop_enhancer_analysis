#%%
import duckdb
import polars as pl

# Kuei consensus enhancers
# Add center position
# Mark Group = "*:Repressive" as "silencer" = True
# Add index
enhancers = (
    duckdb.read_csv("input/data/region.annotation.fcc_starrmpra.group.2025.08.22.tsv")
    .pl()
    .rename({
        "Chrom":"chrom",
        "ChromStart":"start",
        "ChromEnd":"end",
    })
    .with_columns(
        center = (pl.col.start + pl.col.end)//2,
        silencer = pl.col.Group.str.extract("(.*):(.*)",2) == pl.lit("Repressive")
    )
    .drop("Region", "Group")
    .sort("chrom", "start", "end")
    .with_row_index("enhancer_id")
    
)
print("enhancers")
enhancers.write_parquet("output/data/enhancers.parquet")
enhancers
# %%
