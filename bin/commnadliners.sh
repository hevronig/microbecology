## Chronological record of command lines used in reasearch

# tblastx
time tblastx -db <path to database> -query <path to fasta file> -out <path to fileout> -outfmt 6 -num_threads 8 -evalue 0.1

# List lengths of sequences in a fasta file (brew install homebrew/science/bioawk)
bioawk -c fastx '{ print $name, length($seq) }' < scaffolds.fasta

# rsync transfer files
rsync -av --numeric-ids --progress -e "ssh -T -c arcfour -o Compression=no -x" <source path> <destination path>

# Start interactive qsub channel
qsub -l vmem=8gb,mem=1gb,nodes=1:ppn=1 -q P -I

# Start screen session
screen -S <session-name>

# List active screen sessions
srcreen -ls

# Resume screen session
screen -r <session-name>

# Kill screen session
screen -X -S <session-name> quit

# SPAdes assembly single-cell (MDA)
time spades.py --careful --sc -o <output dir> -1 <path to R1.fastq> -2 <path to R2.fastq>

# SPAdes hybrid assembly single-cell, extend scaffolds from 1st assembly round
time spades.py --careful --sc -o <output dir> -1 <path to R1.fastq> -2 <path to R2.fastq> --untrusted-contigs contigs.fasta

# Check assemble contigs coverage (http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/)
sh bbmap.sh in=<path to R1.fastq> in2=<path to R2.fastq> ref=<path to reference scaffolds> covstats=<path to output stats txt file>

# Create salmon index
salmon index -t <path to fasta file> -i <name for index>

# Map reads to reference genome with salmon
salmon quant -1 <path to R1.fq.gz> -2 <path to R2.fq.gz> -o <path to output dir> -i <path to index dir> -l IU


