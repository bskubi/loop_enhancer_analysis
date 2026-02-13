#%%
import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.axis import Axis
import matplotlib.pyplot as plt
from matplotlib.legend import Legend
from matplotlib.patches import Patch
import polars as pl
import seaborn as sns
import numpy as np

loop_strength_expected_expr = pl.max_horizontal(
    pl.col.expectedBL,
    pl.col.expectedDonut,
    pl.col.expectedH,
    pl.col.expectedV
)
loop_strength_expr = pl.col.observed / loop_strength_expected_expr

loops = (
    pl.read_parquet("input/data/loop_categories.parquet")
    .join(
        (
            pl.read_parquet("input/data/loops.parquet")
            .with_columns(
                loop_strength = loop_strength_expr
            )
            .select("loop_id", "loop_strength")            
        ),
        on = "loop_id"
    )
)
# 80mm x 59mm
fig = plt.figure(figsize=(3.149,2.32), layout="constrained")
mpl.rcParams.update({"font.size":7})
grey = "#DEE2E6"
blue = "#007FFF"
ax = sns.violinplot(
    loops.filter(pl.col.regulation == pl.lit("B")),
    x="cohesin",
    y="loop_strength",
    order=["D", "H", "I"],
    log_scale=True,
    facecolor=grey,
    edgecolor="black"
)
sns.despine()
ax.set_xticks(ax.get_xticks(), ["Dependent\n(Peak +/+)", "Hemi-Independent\n(Peak +/-)", "Independent\n(Peak -/-)"], fontsize=7)
ax.set_yticks([1,2,5,10,20], ["1", "2", "5", "10", "20"], fontsize=7)
ax.tick_params(axis="y", which="minor", length=0)
ax.set_ylabel("Loop Strength (O/E)", fontsize=7)
ax.set_xlabel("Cohesin Dependence of Loops", fontsize=7)
fig.suptitle("Bridging P-E Loop Strength Is Cohesin-Invariant", fontsize=7)
fig.savefig("output/panels/metaloop_violinplot/metaloop_violinplot.png", dpi=300)
fig.savefig("output/panels/metaloop_violinplot/metaloop_violinplot.svg", dpi=300)
# %%
