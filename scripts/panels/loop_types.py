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
fig = plt.figure(figsize=(6,2), layout="constrained")
fig.suptitle("Anchors Are Bimodally Positioned With Respect to CREs and Cohesin")
fig.supxlabel("Distance (kb) From Loop Anchor Center To Nearest Element")
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
tss_ax.set_xlabel("TSS")
enh_ax.axvspan(xmin=0, xmax=1000, color=proximal_color)
enh_ax.axvspan(xmin=1000, xmax=x_max, color=distal_color)
enh_ax.set_xlabel("Enhancer Center")
coh_ax.axvspan(xmin=0, xmax=2900, color=proximal_color)
coh_ax.axvspan(xmin=2900, xmax=x_max, color=distal_color)
coh_ax.set_xlabel("Cohesin Peak Center")

coh_ax.legend(
    handles=[
        Patch(label="Proximal", facecolor=proximal_color, edgecolor="grey", linewidth=1),
        Patch(label="Distal", facecolor="#FFFFFF", edgecolor="grey", linewidth=1)
    ],
    loc="upper right",
    frameon=False
)

for ax in axes:
    ax: Axes
    ax.set_xscale("log")
    ax.spines[["top", "right"]].set_visible(False)
    ticks = [10, 100, 1000, 10000, 100000,1000000]
    tick_labels = [".01", ".1", "1", "10", "100", "1000"]
    ax.set_xticks(ticks, tick_labels)
    ax.tick_params(which="minor", bottom=False)
    ax.set_facecolor("white")

fig.savefig("output/panels/anchor_to_element_distance/anchor_to_element_distance.png", dpi=300, transparent=True)
fig.savefig("output/panels/anchor_to_element_distance/anchor_to_element_distance.svg", bbox_inches="tight", transparent=True)
# %%
