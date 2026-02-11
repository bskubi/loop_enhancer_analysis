#%%
%load_ext autoreload
%autoreload 2
#%%
from typing import List
import pickle
import time
import polars as pl
import numpy as np
import hictkpy
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
from params import LoopRegulatoryCategoryAbbr as Reg, MetaloopHeatmapPanel as Panel, seed, Aesthetics as Aes


#%%




######################################################
#   Helper methods to extract matrices for pileup
######################################################
def in_bounds(size, start, end):
    "Check that start and end are in chrom bounds"
    return start >= 0 and end < size

def as_ucsc(chrom, start, end):
    "Format as UCSC string"
    return "{}:{}-{}".format(chrom, start, end)

loops_loaded = 0

def extract_matrices(
        loops: pl.DataFrame, 
        hic: hictkpy.File, 
        flank_bp: int, 
        normalization: str
    ) -> np.typing.NDArray:
    """
    Create individual numpy matrices as basis for pileup
    """
    global loops_loaded
    resolution = hic.resolution()
    loops = (
        loops.with_columns(
            center1 = (pl.col.start1 + pl.col.end1)//2,
            center2 = (pl.col.start2 + pl.col.end2)//2,
        )
        .with_columns(
            center_snap1 = ((pl.col.center1 / resolution).round() * resolution).cast(pl.Int64),
            center_snap2 = ((pl.col.center2 / resolution).round() * resolution).cast(pl.Int64)
        )
        .with_columns(
            flank_start1 = pl.col.center_snap1 - flank_bp,
            flank_end1 = pl.col.center_snap1 + flank_bp,
            flank_start2 = pl.col.center_snap2 - flank_bp,
            flank_end2 = pl.col.center_snap2 + flank_bp
        )
    )
    matrices = []
    chromsizes = hic.chromosomes()
    
    for loop in loops.iter_rows(named=True):
        # Extract anchor position
        chrom1, start1, end1 = loop["chrom1"], loop["flank_start1"], loop["flank_end1"]
        chrom2, start2, end2 = loop["chrom2"], loop["flank_start2"], loop["flank_end2"]

        # Ensure loop anchor flanks are on chromosome
        if not in_bounds(chromsizes[chrom1], start1, end1):
            continue
        if not in_bounds(chromsizes[chrom2], start2, end2):
            continue

        # Format for retrieval from contact matrix
        range1 = as_ucsc(chrom1, start1, end1)
        range2 = as_ucsc(chrom2, start2, end2)

        # Extract contact matrix
        matrix = hic.fetch(
            range1, 
            range2, 
            normalization
        ).to_numpy()
        matrices.append(matrix)

        loops_loaded += 1
        if loops_loaded % 100 == 0:
            print("Loops loaded:", loops_loaded)

    return np.array(matrices)
#%%
###################################################
# Open input data
###################################################

# Get the loop categories to analyze
loops = (
    pl.read_parquet("input/data/loop_categories.parquet")
    .filter(
        pl.col.regulation.is_in(Panel.use_reg_categories),
        pl.col.end2 - pl.col.start1 > Panel.min_loop_distance
    )
)

#%%

# Open matrix to extract matrices from
hic = hictkpy.File(
    "raw/ContactMatrix/inter.hic", 
    resolution=Panel.resolution,
    matrix_type=Panel.matrix_type
)
#%%
######################################################
#   Generate pileup
######################################################

# Generate pileups for each category
matrices = {}
print(Panel.loop_sample_size_per_category)
loops_loaded=0
for category, df in loops.partition_by(["regulation", "cohesin"], as_dict=True).items():
    # Get uniform-size individual loop matrices going into the pileup
    flank_bp = Panel.resolution * Panel.flank_bins
    cat_matrices = extract_matrices(
        df.sample(Panel.loop_sample_size_per_category, seed=seed),
        hic,
        flank_bp,
        Panel.normalization
    )
    assert cat_matrices.shape[0] == Panel.loop_sample_size_per_category
    matrices[category] = cat_matrices

with open(f"output/data/metaloop_heatmap_matrices_{Panel.matrix_type}.pickle", "wb") as matrix_file:
    pickle.dump(matrices, file=matrix_file)
#%%

with open(f"output/data/metaloop_heatmap_matrices_{Panel.matrix_type}.pickle", "rb") as matrix_file:
    matrices = pickle.load(matrix_file)
pileups = {}
vmin, vmax = None, None
for category, cat_matrices in matrices.items():
    # Aggregate matrices into pileup
    pileup = Panel.aggregation_method(
        cat_matrices,
        axis=0
    )
    pileups[category] = pileup

    nanmin = np.nanmin(pileup[pileup>0])
    vmin = nanmin if vmin is None else min(vmin, nanmin)
    nanmax = np.nanmax(pileup[pileup>0])
    vmax = nanmax if vmax is None else max(vmax, nanmax)



#%%
######################################################
#   Plot panel
######################################################
# Produce panel
Aes.apply()

# Ensure using latest version (interactive mode)
flank_bp = Panel.resolution * Panel.flank_bins

# Get names of categories to analyze
regulation_categories = Panel.use_reg_categories
cohesin_categories = Panel.use_coh_categories

unit_size = Panel.panel_columns * Aes.col_w_mm * Aes.mm2inch / n_cols * 0.35
cbar_width_ratio = 0.15

width = unit_size * (n_cols + cbar_width_ratio)
height = unit_size * n_rows

# 2. Initialize unified figure
fig = plt.figure(figsize=(width, height))

# 3. Create unified GridSpec: Columns for heatmaps + 1 column for colorbar
gs = fig.add_gridspec(n_rows, n_cols + 1, width_ratios=[1] * n_cols + [cbar_width_ratio], wspace=0.05, hspace=0.05)

# 4. Initialize axes array to mimic plt.subplots() return structure
axes = []
for r in range(n_rows):
    row_axes = []
    for c in range(n_cols):
        row_axes.append(fig.add_subplot(gs[r, c]))
    axes.append(row_axes)

if n_cols == 1 or n_rows == 1:
    axes = [ax for row in axes for ax in row] # flatten

# 5. Nest sub-gridspec in final column for half-height colorbar
gs_side = gs[:, -1].subgridspec(3, 1, height_ratios=[1, 2, 1])
cax = fig.add_subplot(gs_side[1, 0])

fig.suptitle("Bridging P-E Loop Strengt\nIs Cohesin-Independent", y=1.02)

# Plot pileup heatmaps on uniform color scale
norm = LogNorm(vmin=vmin, vmax=vmax)

# Ensure that row and column position corresponds to Panel spec
panel_data = []
for reg in Panel.use_reg_categories:
    for coh in Panel.use_coh_categories:
        category = (reg, coh)
        pileup = pileups[category]
        panel_data.append((category, pileup))

for i, (category, pileup) in enumerate(panel_data):
    # Get row, column, and Axes object
    row = i % n_rows
    col = i // n_rows

    if n_cols == 1 or n_rows == 1:
        ax = axes[i]
    else:
        ax = axes[row][col]

    # Plot heatmap
    ax = sns.heatmap(
        pileup, 
        norm=norm, 
        ax = ax,
        cbar=False
    )
    ax.set_aspect('equal')
    ax.tick_params(axis='both', length=0)
    ax.tick_params(axis='x', labelrotation=0)
    ax.tick_params(axis='y', labelrotation=90)
    
    # Label axes
    if row == 0:
        # X label = cohesin status
        ax.set_xlabel(category[0])
        ax.xaxis.set_label_position('top')
    if col == 0:
        # Y label = regulatory status
        ax.set_ylabel(category[1], rotation=0, labelpad=10)

        # Plot +/- N kb in center of y axis (left column only)
        ax.set_yticks([pileup.shape[0]/2], [f"±{flank_bp/1000:.0f} kb"])
    else:
        ax.set_yticks([])
    
    if row == len(axes)-1:
        # Plot +/- N kb in center of x axis (bottom row only)
        ax.set_xticks([pileup.shape[0]/2], [f"±{flank_bp/1000:.0f} kb"])
    else:
        ax.set_xticks([])

mappable = ax.collections[0]
side_fig.colorbar(mappable, cax=cax)

plt.savefig("output/panels/metaloop_heatmap/metaloop_heatmap.png", dpi=300)

# %%
