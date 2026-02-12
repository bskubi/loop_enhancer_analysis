#%%
import polars as pl
from params import BPThresholds, ProxDistAbbr

anchors = (
    pl.read_parquet("output/data/anchor_neighbors.parquet")
    .with_columns(
        nn_enh = (
            pl.when(pl.col.nn_dist_enh <= BPThresholds.enhancer_is_proximal)
            .then(pl.lit(ProxDistAbbr.prox))
            .otherwise(pl.lit(ProxDistAbbr.dist))
        ),
        nn_pro = (
            pl.when(pl.col.nn_dist_tss <= BPThresholds.promoter_is_proximal)
            .then(pl.lit(ProxDistAbbr.prox))
            .otherwise(pl.lit(ProxDistAbbr.dist))
        ),
        nn_coh = (
            pl.when(pl.col.nn_dist_coh <= BPThresholds.cohesin_is_proximal)
            .then(pl.lit(ProxDistAbbr.prox))
            .otherwise(pl.lit(ProxDistAbbr.dist))
        )
    )
    .select("anchor", "anchor_id", "loop_id", "chrom", "start", "end", "nn_enh", "nn_pro", "nn_coh")
)
print("anchor_categories")
print(anchors)
anchors.write_parquet("output/data/anchor_categories.parquet")
# %%
