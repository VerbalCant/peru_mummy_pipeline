# This is how I provision a server to run whatever tool I
# have. I use aws ec2.
# I was having a problem with the most recent Anaconda so I just
# use this one.
curl -O https://repo.anaconda.com/archive/Anaconda3-2023.07-2-Linux-x86_64.sh
bash Anaconda3-2023.07-2-Linux-x86_64.sh 
conda update -y -n base -c defaults conda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
sudo apt update
sudo apt install -y awscli # used for connecting to S3 or IDK whatever you want do do on the command line

sudo apt install -y pigz # multithreaded gzipper.

# Now run this command to configure AWS.
aws configure
### Stop coying and pasting now!
###
### FILL IN YOUR INFO WHEN `aws configure` asks for it
###

# Configure some S3 params.
aws configure set default.s3.max_concurrent_requests 20
aws configure set default.s3.multipart_chunksize 32MB
aws configure set default.s3.max_bandwidth 1TB/s

# For example, you can use it to copy files to/from s3:
# aws s3 cp somefile.gz s3://my-bucket/somefile.gz

###
### I set up a different environment for every tool
### because I hate dealing with conflicts. It's kind
### of a pain.

# megahit example
conda create -y -n megahit
conda activate megahit
conda install -y -c bioconda megahit 
# megahit whatever blah blah

# kraken2
conda create -y -n kraken2
conda activate kraken2
conda install -y -c bioconda kraken2

