import os
import subprocess
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from Bio.SeqUtils import GC


# Step 1: Adapter trimming using cutadapt
def trim_adapters(input_file, output_file, adapter_seq):
    cmd = ["cutadapt", "-a", adapter_seq, "-o", output_file, input_file]
    subprocess.run(cmd, check=True)

# Step 2: Align trimmed reads to reference database using bowtie2
def align_reads(input_file, output_file, reference_db):
    cmd = ["bowtie2", "-x", reference_db, "-U", input_file, "-S", output_file, "--very-sensitive"]
    subprocess.run(cmd, check=True)

# Step 3: Parse the SAM file and calculate statistics
def parse_sam_file(sam_file):
    alignments = SeqIO.parse(sam_file, "sam")
    alignment_stats = []
    for record in alignments:
        query_id = record.id
        ref_id = record.reference_id
        alignment_length = len(record.seq)
        gc_content = GC(record.seq)
        alignment_stats.append([query_id, ref_id, alignment_length, gc_content])
    return pd.DataFrame(alignment_stats, columns=["Query_ID", "Ref_ID", "Alignment_Length", "GC_Content"])

# Step 4: Visualize the alignment statistics
def visualize_alignment_stats(df):
    sns.set(style="whitegrid")
    sns.boxplot(data=df, x="Ref_ID", y="Alignment_Length")
    sns.boxplot(data=df, x="Ref_ID", y="GC_Content")
    sns.pairplot(data=df, hue="Ref_ID", diag_kind="kde", markers="+")
    sns.clustermap(data=df.corr(), annot=True, cmap="coolwarm")

# Run the pipeline
def run_pipeline(input_fastq, adapter_seq, reference_db):
    trimmed_fastq = "trimmed_reads.fastq"
    aligned_sam = "aligned_reads.sam"

    trim_adapters(input_fastq, trimmed_fastq, adapter_seq)
    align_reads(trimmed_fastq, aligned_sam, reference_db)

    alignment_stats_df = parse_sam_file(aligned_sam)
    visualize_alignment_stats(alignment_stats_df)

if __name__ == "__main__":
    input_fastq = "8656_S56_L001_R1_001.fastq"
    adapter_seq = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  # Replace with your adapter sequence
    reference_db = "/root/code/data_science/Gut_data/assembly_summary_refseq.txt"  # Replace with your reference database path

    run_pipeline(input_fastq, adapter_seq, reference_db)
