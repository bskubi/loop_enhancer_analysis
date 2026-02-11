import polars as pl
import numpy as np
import matplotlib.pyplot as plt

seed = 42

class BPThresholds:
    enhancer_is_proximal = 1000
    promoter_is_proximal = 1000
    cohesin_is_proximal = 2900

class ProxDistAbbr:
    prox = "P"
    dist = "D"
    prox_pl = pl.lit("P")
    dist_pl = pl.lit("D")

class LoopCohesinCategoryAbbr:
    prox_prox = "D" # Dependent
    prox_dist = "H" # Hemi-independent
    dist_dist = "I" # Independent
    prox_prox_pl = pl.lit("D")
    prox_dist_pl = pl.lit("H")
    dist_dist_pl = pl.lit("I")

    @classmethod
    def label(cls, df: pl.DataFrame, col1, col2) -> pl.Expr:
        col1_pl = pl.col(col1)
        col2_pl = pl.col(col2)
        prox_pl = ProxDistAbbr.prox_pl
        dist_pl = ProxDistAbbr.dist_pl
        return (
            pl.when(col1_pl == prox_pl, col2_pl == prox_pl)
            .then(cls.prox_prox_pl)
            .when(col1_pl == prox_pl, col2_pl == dist_pl)
            .then(cls.prox_dist_pl)
            .when(col1_pl == dist_pl, col2_pl == prox_pl)
            .then(cls.prox_dist_pl)
            .when(col1_pl == dist_pl, col2_pl == dist_pl)
            .then(cls.dist_dist_pl)
        )

class LoopRegulatoryCategoryAbbr:
    bridging_pe = "B"
    nonregulatory = "NR"
    other = "O"
    bridging_pe_pl = pl.lit("B")
    nonregulatory_pl = pl.lit("NR")
    other_pl = pl.lit("O")

    @classmethod
    def label(cls, df: pl.DataFrame, nn_pro1, nn_enh1, nn_pro2, nn_enh2) -> pl.Expr:
        nn_pro1_pl = pl.col(nn_pro1)
        nn_enh1_pl = pl.col(nn_enh1)
        nn_pro2_pl = pl.col(nn_pro2)
        nn_enh2_pl = pl.col(nn_enh2)
        prox_pl = ProxDistAbbr.prox_pl
        dist_pl = ProxDistAbbr.dist_pl
        return (
            pl.when(nn_pro1_pl == prox_pl, nn_enh1_pl == dist_pl, nn_enh2_pl == prox_pl)
            .then(cls.bridging_pe_pl)
            .when(nn_pro2_pl == prox_pl, nn_enh2_pl == dist_pl, nn_enh1_pl == prox_pl)
            .then(cls.bridging_pe_pl)
            .when(nn_pro1_pl == dist_pl, nn_enh1_pl == dist_pl, nn_pro2_pl == dist_pl, nn_enh2_pl == dist_pl)
            .then(cls.nonregulatory_pl)
            .otherwise(cls.other_pl)
        )

class MetaloopHeatmapPanel:
    # Analysis
    resolution = 2000
    flank_bins = 10
    matrix_type = "oe"
    normalization = "RU"
    aggregation_method = np.nanmedian

    # All categories have 792+ loops (BI has exactly 792)
    loop_sample_size_per_category = 792
    min_loop_distance = 50_000
    use_reg_categories = [
        LoopRegulatoryCategoryAbbr.bridging_pe, 
        # LoopRegulatoryCategoryAbbr.nonregulatory
    ]
    use_coh_categories = [
        LoopCohesinCategoryAbbr.dist_dist,
        LoopCohesinCategoryAbbr.prox_dist,
        LoopCohesinCategoryAbbr.prox_prox,
    ]

    # Aesthetics
    panel_columns = 1
    page_depth_frac = .25

class Aesthetics:
    mm2inch = .03937
    col_w_mm = 90
    page_depth_mm = 170
    page_width_mm = 180

    col = 2
    font_size = 7

    @classmethod
    def apply(cls):
        plt.rcParams["font.size"] = cls.font_size
