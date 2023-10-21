# variant calling... for Ancient003 this produces a >300GB VCF
# from a 60GB BAM.
samtools mpileup -uf hg38.analysisSet.fasta SRR20755928_sorted.bam | bcftools call -c --threads 12 -o SRR20755928_variants.vcf

# check qual score distribution
cat SRR20755928_variants.vcf | awk 'BEGIN {OFS="\t"} !/^#/ {print $6}' > SRR20755928_qual_scores.txt
# get stats
awk '{sum+=$1; sumsq+=$1*$1} END {print "Mean = " sum/NR; print "Std Dev = " sqrt(sumsq/NR - (sum/NR)*(sum/NR)); print "Num Rows = " NR}' SRR20755928_qual_scores.txt
# Mean = 89.5416
# Std Dev = 61.7857
# Num Rows = 2.82063e+09

# Randomly sample 100 million scores from the file and sort it
shuf -n 100000000 SRR20755928_qual_scores.txt | sort -n > SRR20755928_qual_scores_sample.txt

# Calculate line numbers for 25th, 50th, 75th, and 90th percentiles
awk 'BEGIN{
    num_lines=100000000;
    print int(num_lines * 0.25);
    print int(num_lines * 0.5);
    print int(num_lines * 0.75)
    print int(num_lines * 0.9)
}' | xargs -I {} awk 'NR=={}' SRR20755928_qual_scores_sample.txt

# 25th percentile: 44.9942
# 50th percentile: 68.9807
# 75th percentile: 110.995
# 90th percentile: 177.996
# Mean = 89.5416
# Std Dev = 61.7857
# Num Rows = 2.82063e+09
# Median is about QUAL 69. a stddev would be ~151.
# Let's try median.

# filter for quality > median
bcftools filter --threads 12 -i'%QUAL>69' SRR20755928_variants.vcf  -Oz -o SRR20755928_filtered_variants.vcf.gz
# filter to SNPs only
bcftools view --threads 12 -v snps -O z -o SRR20755928_filtered_SNPs_only.vcf.gz SRR20755928_filtered_variants.vcf.gz
bcftools index SRR20755928_filtered_SNPs_only.vcf.gz

# whoops
samtools view -h SRR20755928_sorted.bam | \
awk -v OFS='\t' '{ if ($0 ~ /^@/) { print; } else { $12="RG:Z:Group1"; print; } }' | \
samtools view -Sb - > SRR20755928_final.bam

echo "SRR20755928" > sample_name.txt
bcftools reheader -s sample_name.txt -o SRR20755928_filtered_SNPs_only_renamed.vcf.gz SRR20755928_filtered_SNPs_only.vcf.gz

# phase it
whatshap phase \
  --indels \
  --output SRR20755928_filtered_SNPs_only_renamed_phased.vcf.gz \
  --reference ~/data/hg38.analysisSet.fasta \
  --nthreads 12 \
  --memory 60G \
  SRR20755928_filtered_SNPs_only_renamed.vcf.gz \
  SRR20755928_sorted_withRG.bam

bcftools view -r chrM -Oz -o SRR20755928_mtDNA_variants.vcf.gz SRR20755928_filtered_SNPs_only_renamed_phased.vcf.gz
bcftools view -r chrY -Oz -o SRR20755928_Y_variants.vcf.gz SRR20755928_filtered_SNPs_only_renamed_phased.vcf.gz
