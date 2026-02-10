"""
Extract and label anchors from loops
"""

import polars as pl

# 1. Extract anchor position and loop ID
# 2. Rename anchor positions as: chrom, start, center, end
# 3. Label the anchor number
# 4. Concatenate A1 and A2 and add an anchor index

a1 = (
    pl.read_parquet("output/data/loops.parquet")
    .select("loop_id", "chrom1", "start1", "center1", "end1")
    .rename({"chrom1": "chrom", "start1": "start", "center1": "center", "end1": "end"})
    .with_columns(anchor=pl.lit(1))
)

a2 = (
    pl.read_parquet("output/data/loops.parquet")
    .select("loop_id", "chrom2", "start2", "center2", "end2")
    .rename({"chrom2": "chrom", "start2": "start", "center2": "center", "end2": "end"})
    .with_columns(anchor=pl.lit(2))
)


anchors = (
    pl.concat([a1, a2])
    .with_row_index("anchor_id")
)
print("anchors")
print(anchors)
anchors.write_parquet("output/data/anchors.parquet")