import polars as pl
from params import ProxDistAbbr, LoopCohesinCategoryAbbr, LoopRegulatoryCategoryAbbr

anchors = pl.read_parquet("output/data/anchor_categories.parquet")
a1 = (
    anchors.filter(pl.col.anchor == pl.lit(1))
    .drop("anchor")
    .rename({
        col: col + "1"
        for col in anchors.columns 
        if col not in ["loop_id", "anchor"]
    })
)
a2 = (
    anchors.filter(pl.col.anchor == pl.lit(2))
    .drop("anchor")
    .rename({
        col: col + "2"
        for col in anchors.columns 
        if col not in ["loop_id", "anchor"]
    })    
)

loop_categories = a1.join(a2, on="loop_id")
loop_categories = (
    loop_categories.with_columns(
        cohesin = LoopCohesinCategoryAbbr.label(loop_categories, "nn_coh1", "nn_coh2"),
        regulation = LoopRegulatoryCategoryAbbr.label(loop_categories, "nn_pro1", "nn_enh1", "nn_pro2", "nn_enh2"),
    )
    .select("loop_id", "chrom1", "start1", "end1", "chrom2", "start2", "end2", "cohesin", "regulation")
)

print("loop_categories")
with pl.Config(set_tbl_cols=-1):
    print(loop_categories)
loop_categories.write_parquet("output/data/loop_categories.parquet")