"""
Convert raw Erez Aiden K562 micro-C loop calls to parquet.

- Select position and O/E columns
- Add center and id columns
"""
import polars as pl

# New labels for all the columns in the original bedpe
loops_header = [
    'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name', 'score', 'strand1', 'strand2', 
    'color', 'observed', 'expectedBL', 'expectedDonut', 'expectedH', 'expectedV', 
    'fdrBL', 'fdrDonut', 'fdrH', 'fdrV', 'numCollapsed', 'centroid1', 'centroid2', 
    'radius', 'coarse_start_1', 'coarse_end_1', 'coarse_start_2', 'coarse_end_2', 
    'highRes_start_1', 'highRes_end_1', 'highRes_start_2', 'highRes_end_2', 
    'upstream_start_1', 'upstream_end_1', 'downstream_start_2', 'downstream_end_2', 
    'localX', 'localY', 'localObserved', 'localPval', 'localPeakZ', 'local_xwidth1', 
    'local_xwidth2', 'local_ywidth1', 'local_ywidth2', 'localPeakID'
]


loops = (
    # Load in the bedpe
    pl.read_csv(
        "raw/Loops/localizedList_primary_10.bedpe",
        separator="\t",
        comment_prefix="#",
        has_header=False,
        new_columns=loops_header
    )

    # Add center columns
    .with_columns(
        center1 = (pl.col.start1+pl.col.end1)//2,
        center2 = (pl.col.start2+pl.col.end2)//2,
    )

    # Add row index
    .with_row_index('loop_id')

    # Select id, position, and O/E columns
    .select(
        'loop_id', # Unique loop ID
        'chrom1', 'start1', 'center1', 'end1', 'chrom2', 'start2', 'center2', 'end2', # Position
        "observed", "expectedBL", "expectedDonut", "expectedH", "expectedV" # Loop strength
    )
)
print("loops")
print(loops)
loops.write_parquet("output/data/loops.parquet")