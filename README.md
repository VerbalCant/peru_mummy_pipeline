# peru_mummy_pipeline

```
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
../bowtie2-2.5.1-macos-arm64/bowtie2 -p 7 --local --quiet -S SRR20458000_1.sam -x human_genome -1 SRR20458000_1_dedup.fastq -2 SRR20458000_2_dedup.fastq --ma 1 --mp 6,2 --np 1 --rdg 5,3 --rfg 5,3 --met-file SRR20458000_bowtie_alignment_metrics.txt --no-mixed --no-discordant
../bowtie2-2.5.1-macos-arm64/bowtie2 -p 7 --local --quiet -S SRR21031366_1.sam -x human_genome -1 SRR21031366_1_dedup.fastq -2 SRR21031366_2_dedup.fastq --ma 1 --mp 6,2 --np 1 --rdg 5,3 --rfg 5,3 --met-file SRR21031366_bowtie_alignment_metrics.txt --no-mixed --no-discordant
../bowtie2-2.5.1-macos-arm64/bowtie2 -p 7 --local --quiet -S SRR20755928_1.sam -x human_genome -1 SRR20755928_1_dedup.fastq -2 SRR20755928_2_dedup.fastq --ma 1 --mp 6,2 --np 1 --rdg 5,3 --rfg 5,3 --met-file SRR20755928_bowtie_alignment_metrics.txt --no-mixed --no-discordant
```
