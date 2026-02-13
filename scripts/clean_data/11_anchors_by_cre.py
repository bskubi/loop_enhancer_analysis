#%%
import polars as pl
import duckdb
from common import genomic_polars

anchors = pl.read_parquet("output/data/anchors.parquet")
enhancers = pl.read_parquet("output/data/enhancers.parquet")
transcripts = pl.read_parquet("output/data/transcript_quant.parquet")
anchor_neighbors = pl.read_parquet("output/data/anchor_neighbors.parquet")
with duckdb.connect() as conn:
    transcripts_enhancer_distal = conn.execute(
"""
SELECT *
FROM transcripts AS T
ANTI JOIN enhancers AS E
ON
    T.chrom = E.chrom
    AND ABS(T.tss - (E.start+E.end)//2) < 1000
"""
    ).pl()
    # Get all anchors that are distal (>1kb) from the nearest enhancer
    anchors_enhancer_distal = conn.execute(
"""
SELECT *
FROM anchors AS A
ANTI JOIN enhancers AS E
ON
    A.chrom = E.chrom
    AND ABS((A.start+A.end)//2 - (E.start+E.end)//2) < 1000
"""
    ).pl()
    anchors_enhancer_distal.write_parquet("output/data/anchors_enhancer_distal.parquet")

    # Get all combinations of transcripts and proximal (<1kb) enhancer-distal anchors
    # We expect that transcripts that don't have a local enhancer are more
    # dependent on looping interactions to bring an enhancer nearby.
    anchors_enhancer_distal_by_transcript = conn.execute(
"""
SELECT
    T.transcript_id AS transcript_id,
    A.anchor_id,
    A.loop_id
FROM transcripts_enhancer_distal AS T
JOIN anchors_enhancer_distal AS A
ON
    T.chrom = A.chrom
    AND ABS(T.tss - (A.start+A.end)//2) < 1000
"""
    ).pl()
    anchors_enhancer_distal_by_transcript.write_parquet("output/data/anchors_enhancer_distal_by_transcript.parquet")
    with pl.Config(set_tbl_cols=-1):
        print(anchors_enhancer_distal_by_transcript)

    # Get all combinations of enhancers and anchors
    anchors_by_enhancer = conn.execute(
"""
SELECT DISTINCT
    A.anchor_id,
    A.loop_id
FROM enhancers AS E
JOIN anchors AS A
ON
    E.chrom = A.chrom
    AND ABS((E.start+E.end)//2 - (A.start+A.end)//2) < 1000
"""
    ).pl()
    anchors_by_enhancer.write_parquet("output/data/anchors_by_enhancer.parquet")
    with pl.Config(set_tbl_cols=-1):
        print(anchors_by_enhancer)

    # Identify all transcripts and enhancers that are linked via a loop
    # where the transcript anchor does not have a proximal enhancer.
    bridging_pe_loops = conn.execute(
"""
SELECT
    T.loop_id,
    T.transcript_id

FROM anchors_enhancer_distal_by_transcript AS T
JOIN anchors_by_enhancer AS A
ON
    T.loop_id = A.loop_id
    AND T.anchor_id != A.anchor_id
"""
    ).pl()
    bridging_pe_loops.write_parquet("output/data/bridging_pe_loops.parquet")
    with pl.Config(set_tbl_cols=-1):
        print(bridging_pe_loops)
    
    # Identify all loops that transcripts lacking a nearby enhancer
    # participate in.
    transcript_loops = conn.execute(
"""
SELECT
    T.loop_id,
    T.transcript_id

FROM anchors_enhancer_distal_by_transcript AS T
JOIN anchors AS A
ON
    T.loop_id = A.loop_id
    AND T.anchor_id != A.anchor_id
"""
    ).pl()
    transcript_loops.write_parquet("output/data/transcript_loops.parquet")
    with pl.Config(set_tbl_cols=-1):
        print(transcript_loops)

    # For ALL transcripts participating in at least one loop (bridging PE or not)
    # count the number of loops and the number of bridging PE loops
    # they participate in.
    transcript_pe_loop_counts = conn.execute(
"""
WITH
PE_count AS (
    SELECT
        transcript_id,
        COUNT(*) AS bridging_pe_loop_count,
    FROM bridging_pe_loops
    GROUP BY transcript_id
),
L_count AS (
    SELECT
        transcript_id,
        COUNT(*) AS loop_count,
    FROM transcript_loops
    GROUP BY transcript_id
)
SELECT
    L_count.transcript_id,
    COALESCE(PE_count.bridging_pe_loop_count, 0) AS bridging_pe_loop_count,
    L_count.loop_count
FROM L_count
LEFT JOIN PE_count
    ON L_count.transcript_id = PE_count.transcript_id
"""
    ).pl()
    transcript_pe_loop_counts.write_parquet("output/data/transcript_pe_loop_counts.parquet")
    with pl.Config(set_tbl_cols=-1):
        print(transcript_pe_loop_counts)


# %%
