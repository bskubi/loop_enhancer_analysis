import polars as pl
from scipy.spatial import KDTree
import numpy as np
from typing import Tuple

import polars as pl
import duckdb
from typing import Callable, Dict, Any
from scipy.spatial import KDTree
from functools import partial
from warnings import warn

def partition_overlaps(
        func: Callable, 
        df1: pl.DataFrame, 
        df2: pl.DataFrame, 
        partition_by: list[str], 
        intersection: bool = True,
        concat: bool = True,
        a = None, 
        kw = None
    ):
    df1_p = df1.partition_by(partition_by, as_dict=True)
    df2_p = df2.partition_by(partition_by, as_dict=True)
    key_union = sorted(list(set(df1_p.keys()).union(set(df2_p.keys()))))
    results = []
    for key in key_union:
        if (key in df1_p and key in df2_p) or not intersection:
            results.append(
                func(
                    df1_p.get(key), 
                    df2_p.get(key), 
                    *(a or []), 
                    **(kw or {})
                )
            )
    if concat:
        return pl.concat(results)
    else:
        return results

def nearest(
        from_df: pl.DataFrame, 
        to_df: pl.DataFrame, 
        from_cols: list[str], 
        to_cols: list[str], 
        query_kw: Dict[str, Any] = None,
        join_dist: bool = True,
        dist_col = r"dist{i}",
        join_to_df: bool = True,
        join_to_df_cols: list[str] = None,
        to_df_col = r"{col}{i}",
        temp_idx = "__nearest_index__"
    ):
    assert temp_idx not in from_df.columns, (
        f"temp_idx is '{temp_idx}', found in from_df.columns. "
        "Must be new temporary col name."
    )
    assert temp_idx not in to_df.columns, (
        f"temp_idx is '{temp_idx}', found in to_df.columns. "
        "Must be new temporary col name."
    )
    assert join_dist or join_to_df, (
        "At least one of 'join_dist' or 'join_to_df' must be True"
    )
    

    # Detect distances to rows in 'to_pos'
    to_pos = to_df.select(*to_cols).to_numpy()
    tree = KDTree(to_pos)

    # Get distances and indices for every row in 'from_df' 
    from_pos = from_df.select(*from_cols).to_numpy()
    dist, idx = tree.query(from_pos, **(query_kw or {}))

    # Homogenize dist and idx to arrays of shape (k, row_count)
    if len(dist.shape) == 2:
        # If k=2+, dist and idx initially have shape (row_count,k)
        dist = dist.T
        idx = idx.T
    else:
        # If k=1, dist and idx initially have shape (row_count,)
        dist = dist[None, :]
        idx = idx[None, :]

    # Safety checks
    assert dist.shape[1] == from_df.shape[0]
    assert idx.shape[1] == from_df.shape[0]

    # Number of nearest neighbors queried
    k = dist.shape[0]   
    
    # Build result df
    result_df = from_df.clone()

    if join_dist:
        # Add distances to result_df
        dist_dict = {dist_col.format(i=i+1): dist[i] for i in range(k)}
        result_df = pl.concat(
            [result_df, pl.DataFrame(dist_dict)], 
            how="horizontal"
        )

    if join_to_df:
        # 1. Temporarily add to_df indices to result_df 
        idx_dict = {f"{temp_idx}{i+1}": idx[i] for i in range(k)}
        result_df = pl.concat([result_df, pl.DataFrame(idx_dict)], how="horizontal")

        # 2. Iteratively join with 'to_df' on nearest-neighbor i
        to_df_idx = (
            to_df
            .select(*(join_to_df_cols or to_df.columns))
            .with_row_index(temp_idx)
        )
        for i, idx_col in enumerate(idx_dict, start=1):
            
            # 2.1 Rename columns in 'to_df' for current nearest neighbor
            to_df_idx_i = (
                to_df_idx.rename({
                    col: to_df_col.format(col=col, i=i)
                    for col in to_df_idx.columns
                    if col != temp_idx
                })
            )

            # 2.2 Join 'result_df' with renamed 'to_df' and drop idx col
            result_df = (
                result_df.join(
                    to_df_idx_i, 
                    left_on=idx_col, 
                    right_on=temp_idx
                )
                .drop(idx_col)
            )

    return result_df