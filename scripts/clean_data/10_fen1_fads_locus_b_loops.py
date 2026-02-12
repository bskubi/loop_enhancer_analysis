#%%
from pathlib import Path
import polars as pl
import duckdb
import subprocess

#%%
flank = 30000
chrom = "chr11"
region_start = 61718683
region_end = 61901683
print(f"chr11:{region_start}-{region_end}")
gene_list = ["FEN1", "FADS1", "FADS2", "FADS3"]
threshold_bp = 1000
chromsizes = "raw/Reference/hg38.chrom.sizes"
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
        pl.col.transcript_type.is_in(["protein_coding"])
    )
)
print(transcripts)

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
    transcripts.tss,
    transcripts.strand,
    transcripts.gencode_transcript_id,
    transcripts.gencode_transcript_id_full,
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
    .sort("chrom", "start", "end")
    .write_csv(
        mass_screen_enhancers_path, 
        include_header=False, 
        separator="\t",
    )
)
print("Updating enhancer bb")
subprocess.run(f"bedToBigBed {mass_screen_enhancers_path} {chromsizes} output/browser_tracks/mass_screen_enhancers.bb", shell=True)
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
    pa.gene_name,
    ea.start AS enh_start,
    ea.end AS enh_end,
    pa.start AS pro_start,
    pa.end AS pro_end,
    pa.strand AS pro_strand,
    pa.gencode_transcript_id AS gencode_transcript_id,
    pa.gencode_transcript_id_full AS gencode_transcript_id_full,

FROM promoter_anchors AS pa
JOIN enhancer_anchors AS ea
ON
    pa.loop_id = ea.loop_id
    AND pa.anchor_id != ea.anchor_id
"""
).pl()
#%%
transcript_ids = ",".join(bridging_pe_loops["gencode_transcript_id_full"].unique().to_list())
with open("output/browser_tracks/transcript_ids.txt", "w") as file:
    file.write(transcript_ids)

#%%

#228, 208, 10

path = f"output/browser_tracks/bridging_pe_loops.interact"
drop_columns = [col for col in bridging_pe_loops.columns if col not in ["chrom"]]
(
    bridging_pe_loops
    .with_columns(
        chrom = pl.col.chrom1,
        chromStart = pl.col.start1,
        chromEnd = pl.col.end1,
        name = pl.lit("."),
        score = 1000,
        value = 1000,
        exp = pl.lit("."),
        color = 0,
        sourceChrom = pl.col.chrom1,
        sourceStart = pl.col.enh_start,
        sourceEnd = pl.col.enh_end,
        sourceName = pl.lit("."),
        sourceStrand = pl.lit("."),
        targetChrom = pl.col.chrom1,
        targetStart = pl.col.pro_start,
        targetEnd = pl.col.pro_end,
        targetName = pl.lit("."),
        targetStrand = pl.col.pro_strand
    )
    .sort("chrom", "chromStart", "chromEnd")
    .drop(*drop_columns)
    .write_csv(path, separator="\t", include_header=False)
)
chromsizes = "raw/Reference/hg38.chrom.sizes"
output = f"output/browser_tracks/bridging_pe_loops.bb"
subprocess.run(f"bedToBigBed -as=output/browser_tracks/interact.as -type=bed5+13 {path} {chromsizes} {output}", shell=True)
import sys
sys.exit()
# %%
# Extract the needed subset of the bigwigs

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
import glob
from pathlib import Path

bed_files = glob.glob("output/browser_tracks/*.bed")
chromsizes = "raw/Reference/hg38.chrom.sizes"
for file in bed_files:
    bed_df = (
        pl.read_csv(file, separator="\t", new_columns=["chrom", "start", "end"]).select("chrom", "start", "end")
        .sort("chrom", "start", "end")
    )
    print(bed_df)
    temp_file = Path("/tmp") / Path(file).name
    bed_df.write_csv(temp_file, separator="\t", include_header=False)
    output = Path("output/browser_tracks/") / Path(file).name.replace(".bed", ".bb")
    command = f"bedToBigBed {temp_file} {chromsizes} {output}"
    subprocess.run(command, shell=True)
# %%
