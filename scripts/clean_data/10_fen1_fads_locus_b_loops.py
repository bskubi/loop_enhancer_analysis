#%%
import polars as pl
import duckdb

#%%
flank = 30000
chrom = "chr11"
region_start = bridging_pe_loops["start1"].min() - flank
region_end = bridging_pe_loops["end1"].max() + flank
print(f"chr11:{region_start}-{region_end}")
gene_list = ["FEN1", "FADS1", "FADS2", "FADS3"]
threshold_bp = 1000
#%%

# For computational efficiency,
# filter for anchors and enhancers in FEN1/FADS1/2/3 locus (chr11:61-63 MB)
anchors = (
    pl.read_parquet("output/data/anchors.parquet")
    .filter(
        pl.col.chrom == pl.lit(chrom),
        pl.col.start >= region_start,
        pl.col.end <= region_end
    )
    .with_columns(
        center = (pl.col.start + pl.col.end)//2
    )
)
enhancers = (
    pl.read_parquet("output/data/enhancers.parquet")
    .filter(
        pl.col.chrom == pl.lit(chrom),
        pl.col.start >= region_start,
        pl.col.end < region_end
    )
    .with_columns(
        center = (pl.col.start + pl.col.end)//2
    )
)

transcripts = (
    pl.read_parquet("output/data/transcript_quant.parquet")
    .filter(
        pl.col.gene_name.is_in(gene_list),
        pl.col.transcript_biotype == pl.lit("protein_coding")
    )
)

# There are many nearby TSSs, so we will:
# 1. Get all TSS-enhancers connected by loops 
# 2. Select for unique loops for each gene name

promoter_anchors = duckdb.sql(
f"""
SELECT
    anchors.chrom,
    anchors.start,
    anchors.end,
    transcripts.gene_name,
    anchors.anchor_id,
    anchors.loop_id
FROM anchors
JOIN transcripts
ON abs(anchors.center - transcripts.tss) < {threshold_bp}
AND transcripts.gene_name IN {gene_list}
"""
).pl()

enhancer_anchors = duckdb.sql(
f"""
SELECT
    anchors.chrom,
    anchors.start,
    anchors.end,
    enhancers.silencer,
    anchors.anchor_id,
    anchors.loop_id
FROM anchors
JOIN enhancers
ON abs(anchors.center - enhancers.center) < {threshold_bp}
"""
).pl()
#%%
# Write the mass screen enhancers as a bed file (CRISPRi enhancers are already linked)
mass_screen_enhancers_path = "output/browser_tracks/mass_screen_enhancers.bed"
(
    enhancers
    .select("chrom", "start", "end")
    .write_csv(
        mass_screen_enhancers_path, 
        include_header=False, 
        separator="\t",

    )
)

mass_screen_enhancers_header = f"""
track type=bed name="Enhancer Screen" description=" " interactDirectional=true maxHeightPixels=20 visibility=dense
browser position {chrom}:{region_start}-{region_end}
"""


with open(mass_screen_enhancers_path) as file:
    mass_screen_enhancers_text = mass_screen_enhancers_header + file.read()
with open(mass_screen_enhancers_path, "w") as file:
    file.write(mass_screen_enhancers_text)


#%%

# None of the FEN1/FADS1-3 bridging PE loops are to silencers, so leave that out.
bridging_pe_loops = duckdb.sql(
"""
SELECT DISTINCT
    pa.chrom AS chrom1,
    LEAST(pa.start, ea.start) AS start1,
    LEAST(pa.end, ea.end) AS end1,
    pa.chrom AS chrom2,
    GREATEST(pa.start, ea.start) AS start2,
    GREATEST(pa.end, ea.end) AS end2,
    pa.gene_name
FROM promoter_anchors AS pa
JOIN enhancer_anchors AS ea
ON
    pa.loop_id = ea.loop_id
    AND pa.anchor_id != ea.anchor_id
"""
).pl()

#%%
for gene_name, gene_b_loops_bedpe in bridging_pe_loops.partition_by("gene_name", as_dict = True).items():
    print(gene_b_loops_bedpe)

# %%
# Extract the needed subset of the bigwigs
import subprocess
from pathlib import Path
def extract_bigwig(chrom, start, end, chromsizes, input, outpfx):
    bedgraph = f"output/browser_tracks/{outpfx}.bg"
    command1 = f"bigWigToBedGraph -chrom={chrom} -start={start} -end={end} {input} {bedgraph}"
    subprocess.run(command1, shell=True)
    bigwig = f"output/browser_tracks/{outpfx}.bw"
    command2 = f"bedGraphToBigWig {bedgraph} {chromsizes} {bigwig}"
    subprocess.run(command2, shell=True)
    Path(bedgraph).unlink()
    return bigwig

chromsizes = "raw/Reference/hg38.chrom.sizes"

extract_bigwig(chrom, region_start, region_end, chromsizes, "raw/TF_ChIP/sig_pval/RAD21_ENCFF994GBG.bigWig", "RAD21")
extract_bigwig(chrom, region_start, region_end, chromsizes, "raw/TF_ChIP/sig_pval/SMC3_ENCFF596CNE.bigWig", "SMC3")
extract_bigwig(chrom, region_start, region_end, chromsizes, "raw/TF_ChIP/sig_pval/CTCF_ENCFF336UPT.bigWig", "CTCF")

plus1 = extract_bigwig(chrom, region_start, region_end, chromsizes, "raw/RNA/bigwig/plus_strand_signal_of_unique_reads_rep1_ENCFF312ZLI.bigWig", "RNA_plus_rep1")
plus2 = extract_bigwig(chrom, region_start, region_end, chromsizes, "raw/RNA/bigwig/plus_strand_signal_of_unique_reads_rep2_ENCFF272TZC.bigWig", "RNA_plus_rep2")
minus1 = extract_bigwig(chrom, region_start, region_end, chromsizes, "raw/RNA/bigwig/minus_strand_signal_of_unique_reads_rep1_ENCFF530FJG.bigWig", "RNA_minus_rep1")
minus2 = extract_bigwig(chrom, region_start, region_end, chromsizes, "raw/RNA/bigwig/minus_strand_signal_of_unique_reads_rep2_ENCFF123ORL.bigWig", "RNA_minus_rep2")

def average_bigwig(bw1, bw2, outpfx):
    output = f"output/browser_tracks/{outpfx}.bw"
    subprocess.run(f"bigwigCompare -b1 {bw1} -b2 {bw2} --operation mean -o {output}", shell=True)
    Path(bw1).unlink()
    Path(bw2).unlink()
    return output

average_bigwig(plus1, plus2, "RNA_plus")
average_bigwig(minus1, minus2, "RNA_minus")

# %%
