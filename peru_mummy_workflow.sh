# Reports
# This is the NGS report for "Victoria"
#   https://www.the-alien-project.com/wp-content/uploads/2019/04/15-04-2019-Rapport-danalyses-ADN-élargie-dAbraxas-GB.pdf
# Here's the report for all three samples on SRA:
#   https://www.the-alien-project.com/wp-content/uploads/2018/12/ABRAXAS-EN.pdf
# This is the extraction/amplification of "Maria", who is not "Victoria"
# https://www.the-alien-project.com/wp-content/uploads/2018/11/2018-02-06-PALEO-DNA-MARIA-COMPARAISON-ADN.pdf
#
# system notes
# /u/VerbalCant and /u/Big_Tree_Fall_Hard are running this on os x apple silicon, /u/flynnston is on Linux.

# OSX notes
# - Use conda. If you're on Apple Silicon, you'll want Rosetta. There's a lot of stuff
# in bioconda that isn't built for ARM.


#############
# I got this by running AdapterRemoval --identify-adapters on the first million lines of each
# file, which i got with `gunzip -c filename.gz | head -n 100000 >filename_head.fastq`, e.g.:
#      AdapterRemoval --identify-adapters --file1 SRR21031366_1_head.fastq --file2 SRR21031366_2_head.fastq
# Adapter notes:
# SRR21031366 AdapterRemoval notes:
#   --adapter1:  AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG
#   --adapter2:  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
# SRR20458000 AdapterRemoval notes
#   --adapter1:  AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG
#   --adapter2:  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
# SRR20755928 AdapterRemoval notes:
#   --adapter1:  AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG
#   --adapter2:  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
# TruSeq™ single index (previously LT) and TruSeq CD index (previously HT)-based kits
# Matching adapter FASTA: trimmomatic TruSeq3-PE-2.fa

# sra_toolkit prefetch to locally cache the run results
# makes working much faster. these are paired-end runs
bin/prefetch     SRR20458000 --max-size UNLIMITED # ancient0004,  Victoria,  Sample label: Momia 5 -DNA
bin/prefetch     SRR21031366 --max-size UNLIMITED #ancient0002, Victoria, Sample label: Neck Bone Med Seated 00-12 Victoria 4
bin/prefetch     SRR20755928 --max-size UNLIMITED #ancient0003

# fasterq-dump parses the reads into fastq files, used for
# analysis
bin/fasterq-dump SRR20458000 --threads 128
bin/fasterq-dump SRR21031366 --threads 128
bin/fasterq-dump SRR20755928 --threads 128

# Now we use fastqc to check the quality of the data.
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

FastQC

# FastQC checks look good (39%GC) except for high sequence duplication
# Possible causes: PCR relic, degraded sample from ancient DNA, contamination?

# in any case we have to dedup before moving on to the next stage

clumpify.sh in1=SRR20458000_1.fastq in2=SRR20458000_2.fastq out1=SRR20458000_1_dedup.fastq out2=SRR20458000_2_dedup.fastq dedupe
clumpify.sh in1=SRR21031366_1.fastq in2=SRR21031366_2.fastq out1=SRR21031366_1_dedup.fastq out2=SRR21031366_2_dedup.fastq dedupe
clumpify.sh in1=SRR20755928_1.fastq in2=SRR20755928_2.fastq out1=SRR20755928_1_dedup.fastq out2=SRR20755928_2_dedup.fastq dedupe

clumpify.sh in1=SRR20755928_1.fastq.gz in2=SRR20755928_2.fastq.gz out1=SRR20755928_1_dedup.fastq.gz out2=SRR20755928_2_dedup.fastq.gz dedupe
time trimmomatic PE -phred33 -threads 128 SRR20755928_1_dedup.fastq.gz SRR20755928_2_dedup.fastq.gz SRR20755928_R1_dedup_trimmed_paired.fastq.gz SRR20755928_R1_dedup_trimmed_unpaired.fastq.gz  SRR20755928_R2_dedup_trimmed_paired.fastq.gz SRR20755928_R2_dedup_trimmed_unpaired_trimmed.fastq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10   LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
time trimmomatic PE -phred33 -threads 128 SRR21031366_1_dedup.fastq.gz SRR21031366_2_dedup.fastq.gz SRR21031366_R1_dedup_trimmed_paired.fastq.gz SRR21031366_R1_dedup_trimmed_unpaired.fastq.gz  SRR21031366_R2_dedup_trimmed_paired.fastq.gz SRR21031366_R2_dedup_trimmed_unpaired_trimmed.fastq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10   LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
# now i need to get the reference genome to align to
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz && gunzip hg38.analysisSet.fa.gz

# ...and index it for bowtie2
bowtie2-build hg38.analysisSet.fasta human_genome

# align to hg38 w/bowtie2
bowtie2 -p 7 --local --quiet  -S SRR20458000.sam -x human_genome -1 SRR20458000_1_dedup.fastq -2 SRR20458000_2_dedup.fastq human_genome
bowtie2 -p 7 --local --quiet  -S SRR21031366.sam -x human_genome -1 SRR21031366_1_dedup.fastq -2 SRR21031366_2_dedup.fastq human_genome

# (sam to bam and sort)
samtools view -@ 20 -Sb  SRR21031366.sam >SRR21031366.bam
samtools sort -@ 20 SRR21031366.bam -o SRR21031366_sorted.bam
samtools view -@ 20 -b -f 12 SRR21031366_sorted.bam > SRR21031366_unmapped_read_mate.bam
samtools index -@ 20 SRR21031366_sorted.bam
samtools index -@ 20 SRR21031366_unmapped_read_mate.bam
samtools flagstat SRR21031366_sorted.bam > SRR21031366_sorted_flagstat.txt
samtools flagstat SRR21031366_unmapped_read_mate.bam >SRR21031366_unmapped_read_mate_flagstat.txt
samtools idxstats SRR21031366_sorted.bam > SRR21031366_sorted_idxstats.txt
samtools idxstats SRR21031366_unmapped_read_mate.bam >SRR21031366_unmapped_read_mate_idxstats.txt
samtools view -@ 20 -c SRR21031366_sorted.bam >SRR21031366_sorted_count.txt
samtools view -@ 20 -c SRR21031366_unmapped_read_mate.bam >SRR21031366_unmapped_read_mate_count.txt
samtools depth -@ 20 SRR21031366_sorted.bam >SRR21031366_sorted_depth.txt
samtools depth -@ 20 SRR21031366_unmapped_read_mate.bam >SRR21031366_unmapped_read_mate_depth.txt
samtools bam2fq -@ 20 -1 SRR21031366_unmapped_R1.fastq -2 SRR21031366_unmapped_R2.fastq SRR21031366_unmapped_read_mate.bam

samtools view -@ 128 -Sb  SRR20458000_dedup_trimmed_paired_aligned_hg38.sam >SRR21031366.bam
samtools sort -@ 128 SRR21031366.bam -o SRR21031366_sorted.bam
samtools view -@ 128 -b -f 12 SRR21031366_sorted.bam > SRR21031366_unmapped_read_mate.bam
samtools index -@ 128 SRR21031366_sorted.bam
samtools index -@ 128 SRR21031366_unmapped_read_mate.bam
samtools flagstat SRR21031366_sorted.bam > SRR21031366_sorted_flagstat.txt
samtools flagstat SRR21031366_unmapped_read_mate.bam >SRR21031366_unmapped_read_mate_flagstat.txt
samtools idxstats SRR21031366_sorted.bam > SRR21031366_sorted_idxstats.txt
samtools idxstats SRR21031366_unmapped_read_mate.bam >SRR21031366_unmapped_read_mate_idxstats.txt
samtools view -@ 128 -c SRR21031366_sorted.bam >SRR21031366_sorted_count.txt
samtools view -@ 128 -c SRR21031366_unmapped_read_mate.bam >SRR21031366_unmapped_read_mate_count.txt
samtools depth -@ 128 SRR21031366_sorted.bam >SRR21031366_sorted_depth.txt
samtools depth -@ 128 SRR21031366_unmapped_read_mate.bam >SRR21031366_unmapped_read_mate_depth.txt
samtools bam2fq -@ 128 -1 SRR21031366_unmapped_R1.fastq -2 SRR21031366_unmapped_R2.fastq SRR21031366_unmapped_read_mate.bam


#
# Analysis steps
# This is /u/VerbalCant's version, analyzing ancient0002, SRR21031366, "Victoria"
#

# Next step is to run --pa2 for a tax map
#   Build the database first
# kraken2-build --build --standard --db kraken_standard --threads 10
# kraken2-build --db kraken_standard --download-taxonomy --threads 10
# kraken2-build --build --db kraken_standard --threads 10 # ugh i just downloaded it

kraken2 --db k2_standard_20230605  --paired SRR21031366_unmapped_R1.fq.gz SRR21031366_unmapped_R2.fq.gz --use-names  --threads 8 --memory-mapping --gzip-compressed --output SRR21031366_kraken2_output.tab --report SRR21031366_kraken2_report.tab

# First megahit for de novo assembly on hg38-unaligned reads

megahit -1 SRR21031366_unmapped_R1.fq.gz -2 SRR21031366_unmapped_R2.fq.gz --presets meta-large -o SRR21031366_megahit

# then get the final.contigs.fa out of the megahit output directory and start working with that
seqkit stats SRR21031366.megahit_contigs.fa > SRR21031366.megahit_contigs.fa.stats.txt
# file                            format  type  num_seqs      sum_len  min_len  avg_len  max_len
# SRR21031366.megahit_contigs.fa  FASTA   DNA    997,918  614,055,333      200    615.3  279,692

jellyfish count -m 31 -s 100M -t 10 -C -o SRR21031366_megahit_contigs_jellyfish.jf SRR21031366.megahit_contigs.fa
jellyfish histo SRR21031366_megahit_contigs_jellyfish.jf >SRR21031366_megahit_contigs_jellyfish_histogram.txt
# Plot overall histogram - wow, really not a lot of good stuff here
# NOTE: this is from the `bin/` directory of this repo
plot_jellyfish_histo.py SRR21031366_megahit_contigs_jellyfish_histogram.txt -o SRR21031366_megahit_contigs_jellyfish_histogram.png

seqtk comp SRR21031366.megahit_contigs.fa | awk '{sum+=$2; count++} END {print sum/count}'
# This gives an N50 of 615.336

seqkit seq -m 1000 -i SRR21031366.megahit_contigs.fa -o S21031366.megahit_contigs_gt1000.fasta -g

# get a random sample for blastn
seqkit sample -p 0.001 SRR21031366.megahit_contigs.fa > SRR21031366_0.1pct_megahit_sampled_contigs.fa

blastn -query SRR21031366_0.1pct_megahit_sampled_contigs.fa -db /Volumes/Nauvoo/blast/nt -outfmt "6 std staxids" -out SRR21031366_contigs_random_blastn.tsv



# kraken2 the raw reads against nt
kraken2 --db k2_nt --memory-mapping --paired SRR20458000fixed_non_human_R1.fastq SRR20458000fixed_non_human_R2.fastq --output SRR20458000_kraken2_output.txt --report SRR20458000_kraken2_report.txt
kraken2 --db k2_nt --memory-mapping --paired SRR21031366_1.fastq.gz SRR21031366_2.fastq.gz --output SRR21031366_kraken2_output.txt --report SRR21031366_kraken2_report.txt --gzip-compressed

# ancient003
trimmomatic PE -phred33 -threads 128  SRR20755928_1_dedup.fastq.gz SRR20755928_2_dedup.fastq.gz   SRR20755928_1_paired.fastq.gz SRR20755928_1_unpaired.fastq.gz   SRR20755928_2_paired.fastq.gz SRR20755928_2_unpaired.fastq.gz   ILLUMINACLIP:TruSeq3-PE.fa:2:30:10   LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
dedupe.sh  in1=SRR20755928_R1_paired.fasta.gz in2=SRR20755928_R2.paired.fasta.gz out=SRR20755928_paired_deduped.fasta
# clumpify.sh -Xmx24g  in1=SRR20755928_1.fastq.gz in2=SRR20755928_2.fastq.gz out1=SRR20755928_1_dedup.fastq.gz out2=SRR20755928_2_dedup.fastq.gz dedupe threads=128
kraken2 --db /home/ubuntu/data/k2_nt  --paired SRR20755928_1.fastq.gz SRR20755928_2.fastq.gz  --use-names --threads 128 --output SRR20755928_kraken2_raw_reads_output.tab --report SRR20755928_kraken2_raw_reads_report.tab
bowtie2 -p 128 --local --quiet -S SRR20755928_cleaned_deduped_aligned_hg38.sam -x human_genome -1 SRR20755928_1_paired.fastq.gz -2 SRR20755928_2_paired.fastq.gz human_genome \
  | samtools view -@ 128 -bS - \
  | samtools sort  -@ 128 -o SRR20755928_sorted.bam
samtools view -@ 20 -b -f 12 SRR20755928_sorted.bam > SRR20755928_unmapped_read_mate.bam
samtools index -@ 20 SRR20755928_sorted.bam
samtools index -@ 20 SRR20755928_unmapped_read_mate.bam
samtools flagstat SRR20755928_sorted.bam > SRR20755928_sorted_flagstat.txt
samtools flagstat SRR20755928_unmapped_read_mate.bam >SRR20755928_unmapped_read_mate_flagstat.txt
samtools idxstats SRR20755928_sorted.bam > SRR20755928_sorted_idxstats.txt
samtools idxstats SRR20755928_unmapped_read_mate.bam >SRR20755928_unmapped_read_mate_idxstats.txt
samtools view -@ 20 -c SRR20755928_sorted.bam >SRR20755928_sorted_count.txt
samtools view -@ 20 -c SRR20755928_unmapped_read_mate.bam >SRR20755928_unmapped_read_mate_count.txt
samtools depth -@ 20 SRR20755928_sorted.bam >SRR20755928_sorted_depth.txt
samtools depth -@ 20 SRR20755928_unmapped_read_mate.bam >SRR20755928_unmapped_read_mate_depth.txt

## megahit alignment
## spades(via shovill) alignment
mkdir -p tmp && mkdir -p SRR21031366_shovill && shovill --cpus 20 --RAM 480 --force --tmpdir /home/ubuntu/data/tmp  --outdir SRR21031366_shovill --R1 SRR21031366_unmapped_R1.fastq.gz --R2 SRR21031366_unmapped_R2.fastq.gz

## binning megahit
rm -rf workflow && conda activate quast && ./binning_workflow.py   --num_samples 200000 -o workflow SRR21031366_megahit_contigs.fasta.gz --no-checkm  && conda activate checkm && ./binning_workflow.py   --num_samples 200000 -o workflow SRR21031366_megahit_contigs.fasta.gz --no-subset --no-maxbin
## binning spades
rm -rf workflow && conda activate quast && ./binning_workflow.py   --num_samples 200000 -o workflow SRR21031366_spades_contigs.fasta.gz --no-checkm  && conda activate checkm && ./binning_workflow.py   --num_samples 200000 -o workflow SRR21031366_spades_contigs.fasta.gz --no-subset --no-maxbin
## CDS
prodigal -i SRR21031366_spades_contigs.fasta.gz -o SRR21031366_spades_contigs.gbk -a SRR21031366_spades_contigs.faa -p meta
## blastn contig samples
## MEME for motif finding

##
## comparative analysis
##

##
## paper
##




