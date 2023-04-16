from Bio import SeqIO

input_file = "trimmed_reads.fastq"

# Function to split sequence and quality scores into chunks
def split_by_n(seq, n):
    return [seq[i:i + n] for i in range(0, len(seq), n)]

# Read the FASTQ file
with open(input_file, "r") as file:
    for record in SeqIO.parse(file, "fastq"):
        read_id = record.id
        sequence = record.seq
        quality_scores = record.letter_annotations["phred_quality"]

        # Split sequence and quality scores into chunks of 20
        sequence_chunks = split_by_n(sequence, 20)
        quality_chunks = split_by_n(quality_scores, 20)

        # Print read ID
        print(f"Read ID: {read_id}\n")

        # Print sequence and quality scores in a formatted table
        print("Sequence\t\tQuality scores")
        print("--------\t\t-------------")
        for seq_chunk, qual_chunk in zip(sequence_chunks, quality_chunks):
            seq_str = ''.join(seq_chunk)
            qual_str = ' '.join(map(str, qual_chunk))
            print(f"{seq_str}\t{qual_str}")

        print("\n" + "-" * 60 + "\n")
