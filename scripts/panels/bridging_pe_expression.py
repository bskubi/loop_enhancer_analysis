#%%
import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.axis import Axis
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from matplotlib.legend import Legend
from matplotlib.patches import Patch
import polars as pl
import seaborn as sns
import numpy as np
import statsmodels.formula.api as smf

loops = pl.read_parquet("input/data/transcript_pe_loop_counts.parquet")
enhancers = pl.read_parquet("input/data/enhancers.parquet")
transcripts = pl.read_parquet("input/data/transcript_quant.parquet")
x_levels = 4

min_expr = (
    transcripts
    .filter(pl.col.pme_TPM_arithm > 0)
    ["pme_TPM_arithm"]
    .min()
)

locus_resolution = 5_000_000

enhancer_count_by_locus = (
    enhancers
    .with_columns(
        locus = pl.col.chrom + ":" + ((pl.col.start+pl.col.end)//(2 * locus_resolution) * 1_000_000).cast(pl.String()) + "MB"
    )
    .group_by("locus")
    .agg(locus_enhancer_count = pl.len())
)

bridging_pe_expression = (
    loops
    .join(transcripts, on = "transcript_id")
    .with_columns(
        distal_enhancer_interaction = pl.col.bridging_pe_loop_count > 0,
        locus = pl.col.chrom + ":" + ((pl.col.tss // locus_resolution * 1_000_000).cast(pl.String())) + "MB",
        expr = pl.max_horizontal(pl.col.pme_TPM_arithm, min_expr).log10(),
        distal_enhancer_interaction_count = (
            pl.when(pl.col.bridging_pe_loop_count == 0)
            .then(pl.lit("0"))
            .when(pl.col.bridging_pe_loop_count == 1)
            .then(pl.lit("1"))
            .when(pl.col.bridging_pe_loop_count == 2)
            .then(pl.lit("2"))
            .when(pl.col.bridging_pe_loop_count == 3)
            .then(pl.lit("3"))
            .otherwise(pl.lit("4+"))
        )
    )
    .join(enhancer_count_by_locus, on = "locus")
    .select(
        "transcript_id", 
        "chrom", 
        "tss", 
        "bridging_pe_loop_count", 
        "loop_count",
        "locus_enhancer_count",
        "locus",
        "distal_enhancer_interaction", 
        "distal_enhancer_interaction_count",
        "expr"
    )
)


expr_normalized_by_loop_count = (
    bridging_pe_expression
    .group_by("loop_count")
    .agg(expr_loops = pl.col.expr.mean())
)

expr_normalized_by_locus_enhancer_count = (
    bridging_pe_expression
    .group_by("locus_enhancer_count")
    .agg(expr_enhancers = pl.col.expr.mean())
)

expr_normalized_by_locus = (
    bridging_pe_expression
    .group_by("locus")
    .agg(expr_locus = pl.col.expr.mean())
)

bridging_pe_expression = (
    bridging_pe_expression
    .join(expr_normalized_by_loop_count, on = "loop_count")
    .join(expr_normalized_by_locus_enhancer_count, on = "locus_enhancer_count")
    .join(expr_normalized_by_locus, on = "locus")
    
)

model = smf.ols(formula = "expr ~ expr_loops + expr_locus + expr_enhancers", data = bridging_pe_expression.to_pandas())
results = model.fit()

bridging_pe_expression = (
    bridging_pe_expression
    .with_columns(expr_norm = pl.Series(results.resid))
)

palette_map = {
    "0": "black",
    "1": "#007FFF",
    "2": "#0072BB",
    "3": "#00B7EB",
    "4+": "#0047AB"
}

fig = plt.figure(figsize=(3.15,1.38), layout="constrained")
mpl.rcParams.update({"font.size":7})
ax = sns.ecdfplot(
    bridging_pe_expression.filter(pl.col.distal_enhancer_interaction_count == pl.lit("0")),
    hue = "distal_enhancer_interaction_count",
    x = "expr_norm",
    hue_order = ["0", "1", "2", "3", "4+"],
    linewidth=1,
    palette = palette_map,
    ls="--",
    legend = False,
    stat="percent"
)
ax = sns.ecdfplot(
    bridging_pe_expression.filter(pl.col.distal_enhancer_interaction_count != pl.lit("0")),
    hue = "distal_enhancer_interaction_count",
    x = "expr_norm",
    hue_order = ["0", "1", "2", "3", "4+"],
    linewidth=1,
    palette = palette_map,
    legend = False,
    stat="percent"
)
ax.set_ylim(0,100)
ax.set_ylabel("%", fontsize=7)
ax.set_xscale("symlog", linthresh=2)
legend_elements = [
    Line2D([0], [0], color=palette_map["0"], lw=1, linestyle=':', label='0'),
    Line2D([0], [0], color=palette_map["1"], lw=1, linestyle='-', label='1'),
    Line2D([0], [0], color=palette_map["2"], lw=1, linestyle='-', label='2'),
    Line2D([0], [0], color=palette_map["3"], lw=1, linestyle='-', label='3'),
    Line2D([0], [0], color=palette_map["4+"], lw=1, linestyle='-', label='4+')
]
ax.legend(handles=legend_elements, title="P-E Loops", loc="lower right", frameon=False, fontsize=7, ncol=2)
ax.set_xlabel("log(TPM) Residuals")
ax.set_xticks([-2, -1, 0, 1, 2, 3, 4], ["-2", "-1", "0", "1", "2", "3", "4"], fontsize=7)
fig.suptitle("Bridging P-E Loops: A Binary On-Switch For Expression", fontsize=7)
sns.despine()
fig.savefig("output/panels/bridging_pe_expression/bridging_pe_expression.png", dpi=300)
fig.savefig("output/panels/bridging_pe_expression/bridging_pe_expression.svg", dpi=300)
# %%
