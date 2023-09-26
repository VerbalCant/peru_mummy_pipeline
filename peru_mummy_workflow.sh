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
# - default gcc compiler is problematic on os x. better bet is to `brew install gcc` and put
#   `/opt/homebrew/bin` in your path. or at the very least, use gcc@13 for these builds.
#   in particular, kraken2 and samtools support multithreading, but default gcc (at least on apple silicon
#   and Ventura?) won't build with it.

# sra_toolkit prefetch to locally cache the run results
# makes working much faster. these are paired-end runs
bin/prefetch     SRR20458000 --max-size UNLIMITED # ancient0004,  Victoria,  Sample label: Momia 5 -DNA
bin/prefetch     SRR21031366 --max-size UNLIMITED #ancient0002, Victoria, Sample label: Neck Bone Med Seated 00-12 Victoria 4
bin/prefetch     SRR20755928 --max-size UNLIMITED #ancient0003

# fasterq-dump parses the reads into fastq files, used for
# analysis
bin/fasterq-dump SRR20458000 --threads 8
bin/fasterq-dump SRR21031366 --threads 8
bin/fasterq-dump SRR20755928 --threads 8

# Now we use fastqc to check the quality of the data.
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

FastQC

# FastQC checks look good (39%GC) except for high sequence duplication
# Possible causes: PCR relic, degraded sample from ancient DNA, contamination?

# in any case we have to dedup before moving on to the next stage

../bbmap/clumpify.sh -Xmx24g  in1=SRR20458000_1.fastq in2=SRR20458000_2.fastq out1=SRR20458000_1_dedup.fastq out2=SRR20458000_2_dedup.fastq dedupe
../bbmap/clumpify.sh -Xmx24g  in1=SRR21031366_1.fastq in2=SRR21031366_2.fastq out1=SRR21031366_1_dedup.fastq out2=SRR21031366_2_dedup.fastq dedupe
../bbmap/clumpify.sh -Xmx24g  in1=SRR20755928_1.fastq in2=SRR20755928_2.fastq out1=SRR20755928_1_dedup.fastq out2=SRR20755928_2_dedup.fastq dedupe

# now i need to get the reference genome to align to
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz && gunzip hg38.analysisSet.fa.gz

# ...and index it for bowtie2
../bowtie2-2.5.1-macos-arm64/bowtie2-build hg38.analysisSet.fasta human_genome

# align to hg38 w/bowtie2
../bowtie2-2.5.1-macos-arm64/bowtie2 -p 7 --local --quiet  -S SRR20458000.sam -x human_genome -1 SRR20458000_1_dedup.fastq -2 SRR20458000_2_dedup.fastq human_genome
../bowtie2-2.5.1-macos-arm64/bowtie2 -p 7 --local --quiet  -S SRR21031366.sam -x human_genome -1 SRR21031366_1_dedup.fastq -2 SRR21031366_2_dedup.fastq human_genome
../bowtie2-2.5.1-macos-arm64/bowtie2 -p 7 --local --quiet  -S SRR20755928.sam -x human_genome -1 SRR20755928_1_dedup.fastq -2 SRR20755928_2_dedup.fastq human_genome

# less permissive bowtie
#../bowtie2-2.5.1-macos-arm64/bowtie2 -p 7 --local --quiet -S SRR21031366_1.sam -x human_genome -1 SRR21031366_1_dedup.fastq -2 SRR21031366_2_dedup.fastq --ma 1 --mp 6,2 --np 1 --rdg 5,3 --rfg 5,3 --met-file SRR21031366_bowtie_alignment_metrics.txt --no-mixed --no-discordant

#
# Analysis steps
# This is /u/VerbalCant's version, analyzing ancient0002, SRR21031366, "Victoria"
#

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

# Next step is to run kraken2 for a tax map
#   Build the database first
# kraken2-build --build --standard --db kraken_standard --threads 10
# kraken2-build --db kraken_standard --download-taxonomy --threads 10
# kraken2-build --build --db kraken_standard --threads 10 # ugh i just downloaded it

kraken2 --db k2_standard_20230605  --paired SRR21031366_unmapped_R1.fq.gz SRR21031366_unmapped_R2.fq.gz --use-names  --threads 8 --memory-mapping --gzip-compressed --output SRR21031366_kraken2_report.tab --unclassified-out SR21031366_kraken2_unclassified_R#.fastq --classified-out SRR21031366_kraken2_classified_R#.fastq

# First megahit for de novo assembly on hg38-unaligned reads

megahit -1 SRR21031366_unmapped_R1.fq.gz -2 SRR21031366_unmapped_R2.fq.gz --presets meta-large -o SRR21031366_megahit

# will also report on the kraken2 tax classifications.
#
# NOTE: see below about running through spades for error correction. gonna try without it
# first.
# use assembled contigs for BLAST and further classification


###########################
# NOTES
###########################
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5826002/#:~:text=De%20Bruijn%20graph–based%20genome,as%20the%20best%20genome%20assemblers.
# https://astrobiomike.github.io/genomics/de_novo_assembly#assembly
# Try both spades and megahit??
# >> Relevant quotes
# "I've had really good results with SPAdes for isolate or enrichment cultures when I’m trying to reconstruct just one or a few genomes.
# But when working with high diversity metagenomic samples, sometimes SPAdes can’t handle it and MEGAHIT is pretty awesome with how
# well it does with such a small memory footprint – and it’s insanely fast."
# ...
# "I will note that I’ve consistently found that incorporating an error-correction step tends improve assembly results, and the one
# I happen to use is available through the SPAdes program. So even when I end up using the assembly from another program, I typically
# run error-correction on the reads with SPAdes first, and then put the output of that into whatever other assembler I’m using. A
# default SPAdes run with the current version (noted at the top of this page) will run the error-correction step and save the reads
# from it so we can then use them elsewhere. If we don’t want to do the assembly with SPAdes, but run the error-correction step, we
# can set the --only-error-correction flag like we do first here."


# quality analysis of de novo contigs

# BLAST a random sample of the de novo contigs

# annotate de novo genome and hg38 alignment

# what next?
