#%%
import polars as pl
import duckdb
from common.nearest_pl import nearest_pl
# Label anchors with distance to nearest promoter, enhancer, cohesin peak

anchors = pl.read_parquet("output/data/anchors.parquet")
enhancers = pl.read_parquet("output/data/enhancers.parquet")
transcripts = pl.read_parquet("output/data/transcript_quant.parquet")

nearest_pl(anchors, enhancers, ["chrom"], ["center"], ["center"], "nearest_enhancer")

# df = (
#     enhancers
#     .join(
#         nearest_pl(anchors, enhancers, ["chrom"], ["center"], ["center"], "nearest_enhancer"),
#         left_on = "anchor_id",
#         right_on = "index_0"
#     )
#     .join(
#         nearest_pl(anchors, transcripts, ["chrom"], ["center"], ["tss"], "nearest_tss"),
#         left_on = "anchor_id",
#         right_on = "index_0"
#     )
# )
# df.write_parquet("output/data/anchors_with_distances.parquet")
# %%
