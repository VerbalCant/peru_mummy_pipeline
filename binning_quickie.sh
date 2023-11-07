# Don't just blindly run this! You'll remove a directory called "workflow"
# My workaround script for dependency hell. Thanks/go to hell, conda
rm -rf workflow && conda activate quast && ./binning_workflow.py   --num_samples 1000000 -o workflow SRR20458000_megahit_contigs.fasta.gz SRR20458000_spades_contigs.fasta.gz --no-checkm && conda activate checkm && ./binning_workflow.py   --num_samples 1000000 -o workflow SRR20458000_megahit_contigs.fasta.gz SRR20458000_spades_contigs.fasta.gz --no-subset --no-maxbin

rm -rf workflow && conda activate quast && ./binning_workflow.py   --num_samples 1000000 -o workflow SRR20755928_megahit_contigs.fasta.gz SRR20755928_spades_contigs.fasta.gz --no-checkm && conda activate checkm && ./binning_workflow.py   --num_samples 1000000 -o workflow SRR20755928_megahit_contigs.fasta.gz SRR20755928_spades_contigs.fasta.gz --no-subset --no-maxbin


xstreme --p SRR20458000_megahit_contigs2.fasta --oc SRR20458000_megahit_xstreme --meme-p 6 --meme-searchsize 200000
