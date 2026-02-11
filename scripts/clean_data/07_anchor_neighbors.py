#%%
import polars as pl
from typing import Callable, Dict, Any
from functools import partial
from common.genomic_polars import partition_overlaps, nearest

# Label anchors with distance to nearest promoter, enhancer, cohesin peak
anchors = pl.read_parquet("output/data/anchors.parquet")
enhancers = (
    pl.read_parquet("output/data/enhancers.parquet")
    .select("enhancer_id", "chrom", "center")
)
transcripts = (
    pl.read_parquet("output/data/transcript_quant.parquet")
    .select("transcript_id", "chrom", "tss")
)
cohesin = (
    pl.read_parquet("output/data/cohesin.parquet")
    .select("cohesin_id", "chrom", "center")
)

nearest_partial = partial(
    nearest, 
    from_cols = ["center"], 
    query_kw = {"k": 1}
)
genomic_nearest = partial(
    partition_overlaps,
    func = nearest_partial,
    partition_by = ["chrom"]
)
anchor_neighbors = genomic_nearest(
    df1=anchors, 
    df2=enhancers, 
    kw={
        "to_cols": ["center"], 
        "dist_col": "nn_dist_enh", 
        "to_df_col": "nn_{col}",
        "join_to_df_cols": ["enhancer_id"]
    }
)
anchor_neighbors = genomic_nearest(
    df1=anchor_neighbors, 
    df2=transcripts, 
    kw={
        "to_cols": ["tss"], 
        "dist_col": "nn_dist_tss", 
        "to_df_col": "nn_{col}",
        "join_to_df_cols": ["transcript_id"]
    }
)
anchor_neighbors = genomic_nearest(
    df1=anchor_neighbors, 
    df2=cohesin, 
    kw={
        "to_cols": ["center"], 
        "dist_col": "nn_dist_coh", 
        "to_df_col": "nn_{col}",
        "join_to_df_cols": ["cohesin_id"]
    }
)

anchor_neighbors.write_parquet("output/data/anchor_neighbors.parquet")
print("anchor_neighbors")
print(anchor_neighbors)
# %%
