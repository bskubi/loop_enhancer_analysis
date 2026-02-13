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

# 80mm x 60mm
fig = plt.figure(figsize=(3.15, 2.36))
mpl.rcParams["font.size"]=7

fig.suptitle("Distribution of Cohesin and CTCF At Loop Anchors", fontsize=7)


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
    "height_ratios": [.1,1,1,1,1]
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
class VminVmax:
    vmin: int = None
    vmax: int = None

# Get the min/max value for the color range and ylim for the mean lineplot
# For vmax, use 
vmin_vmax = {}
for tf_name in tf_names:
    vmin_vmax[tf_name] = VminVmax()
    for key, data in tf_data.items():
        if tf_name in key:
            tf_vmin_vmax = vmin_vmax[tf_name]
            mean = np.nanmean(data, axis=0)
            tf_vmin_vmax.vmin = 0
            nanmax = np.nanmax(mean)
            tf_vmin_vmax.vmax = nanmax if tf_vmin_vmax.vmax is None else max(nanmax, tf_vmin_vmax.vmax)
            vmin_vmax[tf_name] = tf_vmin_vmax


for i, (row, coh_name) in enumerate(zip(rows[2:], coh_names)):
    mean_keys_iter = iter("aeimqu")

    for j, (tf_name, a1_a2) in enumerate(zip(tf_names, row)):
        # Get the per-anchor data
        a1, a2 = a1_a2
        a1_data = tf_data[("B", coh_name, tf_name, "a1")]
        a2_data = tf_data[("B", coh_name, tf_name, "a2")]

        # Get the per-TF vmin/vmax (also used as the mean lineplot ylim)
        tf_vmin_vmax = vmin_vmax[tf_name]
        vmin = tf_vmin_vmax.vmin
        vmax = tf_vmin_vmax.vmax

        # Construct the A1 heatmap
        im = ax_dict[a1].imshow(
            a1_data, 
            extent=data_extent, 
            cmap="gray_r", 
            interpolation="nearest", 
            aspect="auto", 
            vmin=vmin, 
            vmax=vmax,
            rasterized=True,
            interpolation_stage='rgba'
        )

        # Construct the A2 heatmap
        im = ax_dict[a2].imshow(
            a2_data, 
            extent=data_extent, 
            cmap="gray_r", 
            interpolation="nearest", 
            aspect="auto", 
            vmin=vmin, 
            vmax=vmax,
            rasterized=True,
            interpolation_stage='rgba'
        )

        # Extract the colormap
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
            cb.outline.set_visible(False)
            # Format ticks to appear on the left side
            cax.yaxis.set_ticks_position("left")
            cax.yaxis.set_label_position("left")
            cax.tick_params(axis="y", which="both", pad=-.07)

        line_styles = {
            "D": "-",
            "H": "--",
            "I": ":"
        }
        
        # Construct the A1 mean lineplot
        mean_ax1 = ax_dict[next(mean_keys_iter)]
        x = np.linspace(-2000, 2000, a1_data.shape[1])
        y = np.nanmean(a1_data, axis=0)
        mean_ax1.plot(
            x, 
            y, 
            label=coh_name, 
            ls=line_styles[coh_name],
            color="black",
            linewidth=1
        )
        mean_ax1.set_ylim(vmin, vmax)
        mean_ax1.set_yticks([vmin, int(vmax)])
        mean_ax1.tick_params(axis="y", which="both", pad=-.07)

        # Construct the A1 mean lineplot
        mean_ax2 = ax_dict[next(mean_keys_iter)]
        x = np.linspace(-2000, 2000, a2_data.shape[1])
        y = np.nanmean(a2_data, axis=0)
        mean_ax2.plot(
            x, 
            y, 
            ls=line_styles[coh_name],
            color="black",
            linewidth=1
        )
        mean_ax2.set_yticks([])
        mean_ax2.set_ylim(vmin, vmax)





for k, ax in ax_dict.items():
    if k not in mean_keys:
        ax.yaxis.set_ticks([])
    if k not in "aeimquAE":
        # Leave spines for heatmaps
        pass
    elif k in "aiq":
        # Remove top/right spines for A1 mean plots
        ax.spines[["top", "right"]].set_visible(False)
    elif k in "emu":
        # Only bottom spines for A2 mean plots as they share a y axis
        ax.spines[["top", "right", "left"]].set_visible(False)
    if k not in "dhlptxDH":
        ax.xaxis.set_ticks([])
    else:
        ax.xaxis.set_ticks([0], ["+/- 2kb"])
    ax.set_xlim(-2000,2000)
    ax.tick_params(axis="x", length=0)

ax_dict["a"].set_ylabel("Mean\nPval", fontsize=7)
ax_dict["a"].yaxis.set_label_coords(-.3, .5)
ax_dict["b"].set_ylabel("Dep.", fontsize=7)
ax_dict["b"].yaxis.set_label_coords(-.3, .5)
ax_dict["c"].set_ylabel("Hemi-Ind.", fontsize=7)
ax_dict["c"].yaxis.set_label_coords(-.3, .5)
ax_dict["d"].set_ylabel("Ind.", fontsize=7)
ax_dict["d"].yaxis.set_label_coords(-.3, .5)

handles, labels = ax_dict["a"].get_legend_handles_labels()
labels[labels.index("D")] = "Dep."
labels[labels.index("H")] = "Hemi-Ind."
labels[labels.index("I")] = "Ind."

fig.legend(
    handles, 
    labels, 
    loc="upper left",
    # bbox_to_anchor=(-.2, 1.01),
    ncol=3, 
    frameon=False,
    fontsize=7
)

for a in "aiq":
    ax_dict[a].set_xlabel("A1", fontsize=7)
    ax_dict[a].xaxis.set_label_position("top")
for a in "emu":
    ax_dict[a].set_xlabel("A2", fontsize=7)
    ax_dict[a].xaxis.set_label_position("top")

for a, tf in zip("123", ["RAD21", "SMC3", "CTCF"]):
    ax_dict[a].set_title(tf, y=.95, fontsize=7)
    # ax_dict[a].plot([0.1, 0.9], [.5, .5], color='black', transform=ax_dict[a].transAxes, linewidth=.5)
    ax_dict[a].set_axis_off()

fig.savefig("output/panels/tornado_grid/tornado_grid.png", dpi=300, pad_inches=0, transparent=True)
fig.savefig("output/panels/tornado_grid/tornado_grid.svg", dpi=300, pad_inches=0, transparent=True)
# %%
