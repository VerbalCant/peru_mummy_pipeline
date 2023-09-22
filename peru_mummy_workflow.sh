# Reports
# https://www.the-alien-project.com/wp-content/uploads/2019/04/15-04-2019-Rapport-danalyses-ADN-élargie-dAbraxas-GB.pdf
# https://www.the-alien-project.com/wp-content/uploads/2018/11/2018-02-06-PALEO-DNA-MARIA-COMPARAISON-ADN.pdf

# system notes
# /u/VerbalCant and /u/Big_Tree_Fall_Hard are running this on os x apple silicon, /u/flynnston is on Linux.
# OSX notes
# - default gcc compiler is problematic on os x. better bet is to `brew install gcc` and put
#   `/opt/homebrew/bin` in your path. or at the very least, use gcc@13 for these builds.
#   in particular, kraken2 and samtools support multithreading, but default gcc won't
#   build it. if you

# sra_toolkit prefetch to locally cache the run results
# makes working much faster. these are paired-end runs
bin/prefetch     SRR20458000 --max-size UNLIMITED # ancient0004
bin/prefetch     SRR21031366 --max-size UNLIMITED #ancient0002
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
# these are still running
../bowtie2-2.5.1-macos-arm64/bowtie2 -p 7 --local --quiet -S SRR20458000_1.sam -x human_genome -1 SRR20458000_1_dedup.fastq -2 SRR20458000_2_dedup.fastq --ma 1 --mp 6,2 --np 1 --rdg 5,3 --rfg 5,3 --met-file SRR20458000_bowtie_alignment_metrics.txt --no-mixed --no-discordant
../bowtie2-2.5.1-macos-arm64/bowtie2 -p 7 --local --quiet -S SRR21031366_1.sam -x human_genome -1 SRR21031366_1_dedup.fastq -2 SRR21031366_2_dedup.fastq --ma 1 --mp 6,2 --np 1 --rdg 5,3 --rfg 5,3 --met-file SRR21031366_bowtie_alignment_metrics.txt --no-mixed --no-discordant
../bowtie2-2.5.1-macos-arm64/bowtie2 -p 7 --local --quiet -S SRR20755928_1.sam -x human_genome -1 SRR20755928_1_dedup.fastq -2 SRR20755928_2_dedup.fastq --ma 1 --mp 6,2 --np 1 --rdg 5,3 --rfg 5,3 --met-file SRR20755928_bowtie_alignment_metrics.txt --no-mixed --no-discordant

#
# Analysis steps
# This is /u/VerbalCant's version, analyzing ancient0002, SRR21031366,
#

# (sam to bam and sort)
## samtools view -Sb SRRxxx.sam > SRRxxx.bam
## samtools sort SRRxxx.bam -o SRRxxx_sorted.bam
samtools view -@ 20 SRR21031366.sam >SRR21031366.bam
##samtools sort -@ 20 SRR21031366.bam -o SRR21031366_sorted.bam
# de novo assembly of unaligned reads.


# get both read and mate unmapped and convert back to paired end fastq for assembly
# 4 is read unmapped, 8 is mate unmapped, so let's max our chances by 4 & 8
#samtools view -@ 20 -b -f 12 SRR21031366_sorted.bam > SRR21031366_unmapped_read_mate.bam

# index the bams for great justice
##samtools index -@ 20 SRR21031366_sorted.bam
##samtools index -@ 20 SRR21031366_unmapped_read_mate.bam

# let's compare
##samtools flagstat SRR21031366_sorted.bam > SRR21031366_sorted_flagstat.txt
##samtools flagstat SRR21031366_unmapped_read_mate.bam >SRR21031366_unmapped_read_mate_flagstat.txt

#samtools view -@ 20 -c SRR20458000_sorted.bam # 732469872
#samtools view -@ 20 -c SRR20458000_unmapped_read_mate.bam  # 635844786

../samtools/bin/samtools bam2fq -1 SRR20458000_unmapped_R1.fastq -2 SRR20458000_unmapped_R2.fastq SRR20458000_unmapped_read_mate.bam

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

# de novo assembly of all reads

# quality analysis of de novo contigs

# taxmap with kraken2
# kraken2 --db kraken_db --output SRRxxx_kraken_output.txt --report SRRxxx_kraken_report.txt contigs.fasta

# BLAST a random sample of the de novo contigs

# annotate de novo genome and hg38 alignment

# what next?
