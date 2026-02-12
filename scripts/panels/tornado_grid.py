#%%
from dataclasses import dataclass
#%%
import pyBigWig
import polars as pl
from collections import defaultdict
import numpy as np

loops = (
    pl.read_parquet("input/data/loop_categories.parquet")
    .filter(pl.col.regulation.is_in(["B"]))
)

tfs = {
    "RAD21": "input/data/RAD21_ENCFF994GBG.bigWig",
    "SMC3": "input/data/SMC3_ENCFF596CNE.bigWig",
    "CTCF": "input/data/CTCF_ENCFF336UPT.bigWig"
}

tf_data = defaultdict(list)
for name, path in tfs.items():
    file = pyBigWig.open(path)
    for category, loop_cat in loops.partition_by("regulation", "cohesin", as_dict=True).items():
        loop_cat = loop_cat.sample(792, seed=42)
        for loop in loop_cat.iter_rows(named=True):
            center1 = (loop["start1"]+loop["end1"])//2
            center2 = (loop["start2"]+loop["end2"])//2
            a1 = file.stats(loop["chrom1"], center1-2000, center1+2000, nBins=20)
            a2 = file.stats(loop["chrom2"], center2-2000, center2+2000, nBins=20)
            k1 = (*category, name, "a1")
            k2 = (*category, name, "a2")
            tf_data[k1].append(a1)
            tf_data[k2].append(a2)
for k, a in tf_data.items():
    tf_data[k] = np.array(a)


#%%
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

nrows = 3
ncols = 8
width = ncols * .5
height = nrows * 2

# 1. Pyplot binds backend and IPython display hooks
fig = plt.figure(figsize=(width, height))
mpl.rcParams["font.size"]=7

fig.suptitle("Distribution of Cohesin and CTCF At Loop Anchors")


ax_dict = fig.subplot_mosaic(
"""
11.22.33
ae.im.qu
bf.jn.rv
cg.ko.sw
dh.lp.tx
""",
gridspec_kw={
    "width_ratios": [
        1,1, .5, 1,1, .5, 1,1
    ],
    "height_ratios": [.1,1,2,2,2]
}
)

mosaic = """
11.22.33
ae.im.qu
bf.jn.rv
cg.ko.sw
dh.lp.tx
""".strip()
rows = [row.split(".") for row in mosaic.split("\n")]
tf_names = ["RAD21", "SMC3", "CTCF"]
coh_names = ["D", "H", "I"]
data_extent = [-2000, 2000, 0, 20]
mean_keys = "aeimqu"

@dataclass
class ColorRange:
    vmin: int = None
    vmax: int = None

color_ranges = {}
for tf_name in tf_names:
    color_ranges[tf_name] = ColorRange()
    for key, data in tf_data.items():
        if tf_name in key:
            crange = color_ranges[tf_name]
            mean = np.nanmean(data, axis=0)
            crange.vmin = 0
            nanmax = np.nanmax(mean)
            crange.vmax = nanmax if crange.vmax is None else max(nanmax, crange.vmax)
            color_ranges[tf_name] = crange
            

mean_max = {tf_name: 0 for tf_name in tf_names}

for i, (row, coh_name) in enumerate(zip(rows[2:], coh_names)):
    mean_keys_iter = iter("aeimqu")

    for j, (tf_name, a1_a2) in enumerate(zip(tf_names, row)):
        a1, a2 = a1_a2
        a1_data = tf_data[("B", coh_name, tf_name, "a1")]
        a2_data = tf_data[("B", coh_name, tf_name, "a2")]
        crange = color_ranges[tf_name]
        vmin = crange.vmin
        vmax = crange.vmax
        im = ax_dict[a1].imshow(a1_data, extent=data_extent, cmap="cividis", interpolation="nearest", aspect="auto", vmin=vmin, vmax=vmax)
        if i == 0:
        # Get the bounding box of the target axes in figure coordinates
            pos = ax_dict[a1].get_position()
            
            # Define colorbar dimensions relative to that box
            cbar_width = 0.015
            pad = 0.005 # Space between axes and colorbar
            
            # Calculate [left, bottom, width, height]
            # x = axes_left_edge - cbar_width - padding
            cax_rect = [pos.x0 - cbar_width - pad, pos.y0, cbar_width, pos.height]
            
            cax = fig.add_axes(cax_rect)
            cb = fig.colorbar(im, cax=cax, ticks=[0, int(vmax)])
            # Format ticks to appear on the left side
            cax.yaxis.set_ticks_position("left")
            cax.yaxis.set_label_position("left")
        im = ax_dict[a2].imshow(a2_data, extent=data_extent, cmap="cividis", interpolation="nearest", aspect="auto", vmin=vmin, vmax=vmax)
        
        mean_ax1 = ax_dict[next(mean_keys_iter)]
        x = np.linspace(-2000, 2000, a1_data.shape[1])
        y = np.nanmean(a1_data, axis=0)
        mean_ax1.plot(x, y, label=coh_name)
        mean_max[tf_name] = max(mean_max[tf_name], y.max())
        mean_ax1.set_yticks([0, int(vmax)])

        mean_ax2 = ax_dict[next(mean_keys_iter)]
        x = np.linspace(-2000, 2000, a2_data.shape[1])
        y = np.nanmean(a2_data, axis=0)
        mean_ax2.plot(x, y)
        mean_ax2.set_yticks([])
        mean_max[tf_name] = max(mean_max[tf_name], y.max())

for i, (row, coh_name) in enumerate(zip(rows[2:], coh_names)):
    mean_keys_iter = iter("aeimqu")

    for j, (tf_name, a1_a2) in enumerate(zip(tf_names, row)):
        mean_ax1 = ax_dict[next(mean_keys_iter)]
        mean_ax1.set_ylim(0, mean_max[tf_name])
        mean_ax2 = ax_dict[next(mean_keys_iter)]
        mean_ax2.set_ylim(0, mean_max[tf_name])


for k, ax in ax_dict.items():
    if k not in mean_keys:
        ax.yaxis.set_ticks([])
    if k not in "aeimquAE":
        for key, spine in ax.spines.items():
            spine.set_visible(False)
    elif k in "aiq":
        ax.spines[["top", "right"]].set_visible(False)
    elif k in "emu":
        ax.spines[["top", "right", "left"]].set_visible(False)
    if k not in "dhlptxDH":
        ax.xaxis.set_ticks([])
    else:
        ax.xaxis.set_ticks([0], ["+/- 2kb"])
    ax.set_xlim(-2000,2000)
    ax.tick_params(axis="x", length=0)

ax_dict["a"].set_ylabel("Sig Pval\n(Mean)")
ax_dict["a"].yaxis.set_label_coords(-.5, .5)
ax_dict["b"].set_ylabel("D", rotation=0)
ax_dict["b"].yaxis.set_label_coords(-.5, .5)
ax_dict["c"].set_ylabel("H", rotation=0)
ax_dict["c"].yaxis.set_label_coords(-.5, .5)
ax_dict["d"].set_ylabel("I", rotation=0)
ax_dict["d"].yaxis.set_label_coords(-.5, .5)

handles, labels = ax_dict["a"].get_legend_handles_labels()

fig.legend(
    handles, 
    labels, 
    loc="upper left",
    bbox_to_anchor=(0.3, .95),
    ncol=3, 
    frameon=False,
    fontsize=7
)

for a in "aiq":
    ax_dict[a].set_xlabel("A1")
    ax_dict[a].xaxis.set_label_position("top")
for a in "emu":
    ax_dict[a].set_xlabel("A2")
    ax_dict[a].xaxis.set_label_position("top")

for a, tf in zip("123", ["RAD21", "SMC3", "CTCF"]):
    ax_dict[a].set_title(tf)
    ax_dict[a].set_axis_off()

fig.savefig("output/panels/tornado_grid/tornado_grid.png", dpi=300)
# %%
