#!/usr/bin/env bash
sudo apt install -y libncurses5 libncurses5-dev
conda create -y -n samtools
conda activate samtools
conda install -y -c bioconda samtools

#Convert SAM to BAM for better performance and storage.
samtools view -@ 128 -Sb SRR20458000_dedup_trimmed_paired_aligned_hg38.sam >SRR20458000.bam
# Sort the BAM file by reference and position.
samtools sort -@ 128 SRR20458000.bam -o SRR20458000_sorted.bam
# Extract reads where both mates are unmapped.
samtools view -@ 128 -b -f 12 SRR20458000_sorted.bam > SRR20458000_unmapped_read_mate.bam
# Create indexes for quick random access
samtools index -@ 128 SRR20458000_sorted.bam
samtools index -@ 128 SRR20458000_unmapped_read_mate.bam
# Generate a summary of alignments
samtools flagstat SRR20458000_sorted.bam > SRR20458000_sorted_flagstat.txt
samtools flagstat SRR20458000_unmapped_read_mate.bam >SRR20458000_unmapped_read_mate_flagstat.txt
# Generate alignment statistics
samtools idxstats SRR20458000_sorted.bam > SRR20458000_sorted_idxstats.txt
samtools idxstats SRR20458000_unmapped_read_mate.bam >SRR20458000_unmapped_read_mate_idxstats.txt
# Count the number of aligned reads
samtools view -@ 128 -c SRR20458000_sorted.bam >SRR20458000_sorted_count.txt
samtools view -@ 128 -c SRR20458000_unmapped_read_mate.bam >SRR20458000_unmapped_read_mate_count.txt
# Compute coverage depth for each position
samtools depth -@ 128 SRR20458000_sorted.bam >SRR20458000_sorted_depth.txt
samtools depth -@ 128 SRR20458000_unmapped_read_mate.bam >SRR20458000_unmapped_read_mate_depth.txt
# Convert unmapped read pairs back to FASTQ format.
samtools bam2fq -@ 128 -1 SRR20458000_unmapped_R1.fastq -2 SRR20458000_unmapped_R2.fastq SRR20458000_unmapped_read_mate.bam

# copy all of the files to the s3 bucket
ls *.txt | xargs -I {} aws s3 cp {} s3://peru-mummy/
ls *.bam* | xargs -I {} aws s3 cp {} s3://peru-mummy/

conda create -y -n shovill && conda activate shovill && conda install -c conda-forge -c bioconda -c defaults shovill
export TMPDIR=/home/ubuntu/data/tmp

# spades runs into a ulimit issue because each thread requires a bunch of open files. i had to
# limit it to 20 CPUs to get it to run without shitting itself over a ulimit.
# also, pilon in shovill takes FOREVER, so i just don't run it.
mkdir -p tmp && mkdir -p SRR21031366_shovill && shovill --cpus 20 --RAM 480 --force --tmpdir /home/ubuntu/data/tmp  --outdir SRR21031366_shovill --R1 SRR21031366_unmapped_R1.fastq.gz --R2 SRR21031366_unmapped_R2.fastq.gz

# shut down the node
sudo shutdown -h now
