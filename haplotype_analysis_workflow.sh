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

bcftools filter -i'%QUAL>69' SRR20755928_variants.vcf  -Oz -o SRR20755928_filtered_variants.vcf.gz
