import polars as pl
from scipy.spatial import KDTree
import numpy as np
from typing import Tuple

def nearest(from_pos, to_pos, tree_kws: dict = None, query_kws: dict = None) -> Tuple[np.ndarray, np.ndarray]:
    tree = KDTree( to_pos, **(tree_kws or {}) )
    dist, idx = tree.query( from_pos, **(query_kws or {}) )
    if len(dist.shape) == 1:
        return (
            dist.reshape(-1, 1),
            idx.reshape(-1, 1)
        )
    return dist, idx

def nearest_pl(
        from_df: pl.DataFrame, 
        to_df: pl.DataFrame,
        partition_by: list[str],
        from_cols: list[str], 
        to_cols: list[str], 
        nearest_prefix: str = "nearest", 
        index_prefix: str = "index",
        tree_kws: dict = None, 
        query_kws: dict = None,
        order_col = "__nearest_pl_index__"
    ) -> pl.DataFrame:
    from_partitioned = (
        from_df
        .with_row_index(order_col)
        .partition_by(*partition_by, as_dict=True)
    )
    to_partitioned = to_df.partition_by(*partition_by, as_dict=True)
    shared_keys = set(from_partitioned.keys()).intersection(to_partitioned.keys())
    all_results = []
    index_offset = 0
    for key in shared_keys:
        from_block = from_partitioned[key]
        to_block = to_partitioned[key]
        print(from_block)
        from_pos = from_block.select(*from_cols)
        to_pos = to_block.select(*to_cols)
        dist, idx = nearest( from_pos, to_pos, tree_kws, query_kws )
        idx += index_offset
        index_offset += to_block.shape[0]

        nearest_dict = {
            f"{nearest_prefix}_{i}": pl.Series(dist_i)
            for i, dist_i in enumerate(dist.T)
        } if nearest_prefix else {}
        
        index_dict = {
            f"{index_prefix}_{j}": pl.Series(idx_j)
            for j, idx_j in enumerate(idx.T)
        } if index_prefix else {}

        result_dict = {
            **{order_col: from_block[order_col]},
            **nearest_dict,
            **index_dict
        }

        all_results.append( pl.DataFrame(result_dict) )
    return (
        pl.concat(all_results)
        .sort(order_col)
        .drop(order_col)
    )
