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

anchors = pl.read_parquet("input/data/anchor_neighbors.parquet")

# 80mm x 35mm
fig = plt.figure(figsize=(3.14,1.38), layout="constrained")
fig.suptitle("Anchors Are Bimodally Positioned\nWith Respect to CREs and Cohesin", fontsize=7)
fig.supxlabel("Distance From Loop Anchor Center To Nearest Element", fontsize=7)
axes = fig.subplots(1,3, sharey=True)
mpl.rcParams.update({"font.size":7})
tss_ax: Axes = axes[0]
enh_ax: Axes = axes[1]
coh_ax: Axes = axes[2]

#sns.histplot(x = anchors["nn_dist_enh"], stat="percent", bins=100, log_scale=True)
x_max = max([anchors["nn_dist_tss"].max(), anchors["nn_dist_enh"].max(), anchors["nn_dist_coh"].max()])
tss_bins_max = anchors["nn_dist_tss"].log10().max()
enh_bins_max = anchors["nn_dist_enh"].log10().max()
coh_bins_max = anchors["nn_dist_coh"].log10().max()
bins_max = max([enh_bins_max, tss_bins_max, coh_bins_max])
bins = np.logspace(1, bins_max, 100)
grey = "#DEE2E6"
blue = "#007FFF"
sns.histplot(
    x=anchors["nn_dist_tss"],
    color=blue,
    bins=bins,
    common_norm=False,
    fill=False,
    element="poly",
    line_kws={"linewidth":.1},
    stat="percent",
    ax=tss_ax
)
sns.histplot(
    x=anchors["nn_dist_enh"],
    color=blue,
    bins=bins,
    common_norm=False,
    fill=False,
    element="poly",
    line_kws={"linewidth":.1},
    stat="percent",
    ax=enh_ax
)
sns.histplot(
    x=anchors["nn_dist_coh"],
    color=blue,
    bins=bins,
    common_norm=False,
    fill=False,
    element="poly",
    line_kws={"linewidth":.1},
    stat="percent",
    ax=coh_ax
)
proximal_color = grey
distal_color = "#FFFFFF"
tss_ax.axvspan(xmin=0, xmax=1000, color=proximal_color)
tss_ax.axvspan(xmin=1000, xmax=x_max, color=distal_color)
tss_ax.set_xlabel("TSS", fontsize=7)
enh_ax.axvspan(xmin=0, xmax=1000, color=proximal_color)
enh_ax.axvspan(xmin=1000, xmax=x_max, color=distal_color)
enh_ax.set_xlabel("Enhancer", fontsize=7)
coh_ax.axvspan(xmin=0, xmax=2900, color=proximal_color)
coh_ax.axvspan(xmin=2900, xmax=x_max, color=distal_color)
coh_ax.set_xlabel("Cohesin Peak", fontsize=7)

fig.legend(
    handles=[
        Patch(label="Proximal", facecolor=proximal_color, edgecolor="grey", linewidth=1),
        Patch(label="Distal", facecolor="#FFFFFF", edgecolor="grey", linewidth=1)
    ],
    handlelength=1,
    loc="upper right",
    bbox_to_anchor=(1.02,1.02),
    frameon=False,
    fontsize=7
)

for ax in axes:
    ax: Axes
    ax.set_xscale("log")
    ax.spines[["top", "right"]].set_visible(False)
    ticks = [1000, 1000000]
    tick_labels = ["1kb", "1mb"]
    ax.set_xticks(ticks, tick_labels, fontsize=7)
    ax.set_yticks([1, 3, 5], ["1%", "3%", "5%"], fontsize=7)
    ax.tick_params(which="minor", bottom=False)
    ax.set_facecolor("white")
    ax.set_ylabel("")


fig.savefig("output/panels/loop_types/loop_types.png", dpi=300, transparent=True)
fig.savefig("output/panels/loop_types/loop_types.svg", dpi=300, transparent=True)
# %%
