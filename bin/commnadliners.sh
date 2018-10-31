## Chronological record of command lines used in reasearch

# Rename fasta headers with sed
time for file in *.fa; do cut -d ' ' -f2 $file | sed -e 's/^scaffold/>scaffold/' > trim."$file"; done

# Add a string to fasta headers
sed 's/^\(>.*\)$/\1@RS2012/' input.fa > output.fna

# fastq to fasta sed [stakoverflow](https://stackoverflow.com/questions/1542306/converting-fastq-to-fasta-with-sed-awk)
sed -n '1~4s/^@/>/p;2~4p' file.fastq > file.fasta

# Split file
split -l 200000 --numeric-suffixes=0001 --suffix-length=4 --additional-suffix=.fa <files to split> <prefix to use>

# Prodigal (with parallel)
time find . -name '*.fa' | parallel -j4 'prodigal -i {} -d {}.fna -a {}.faa -p meta -o {}.gbk'

# Find and exce built-in operations
find <path dir> -name "<string to find>" -exec cat {} \; > <output file>

# Find and execute a scripts passing I/O variables
find . -name "file.txt" -exec ~/script.sh {} {}.out.tsv \;

# Trim string with awk
for file in <path to files>; do base=$(dirname "$file"); awk -v OFS='\t' -F '\t' '{ $1=substr($1,1,12); print}' "$file" > "$base"/new_file.tsv ; done

# tblastx
time tblastx -db <path to database> -query <path to fasta file> -out <path to fileout> -outfmt 6 -num_threads 8 -evalue 0.1

# List lengths of sequences in a fasta file using bioawk (brew install homebrew/science/bioawk)
bioawk -c fastx '{ print $name, length($seq) }' < scaffolds.fasta

# List lengths of sequences in a fasta file using regular awk
# from https://stackoverflow.com/questions/23992646/sequence-length-of-fasta-file
awk '/^>/ {if (seqlen){print seqlen}; printf $0" " ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' file.fasta > lengths.tsv

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

# Use ATLAS qsub to submit a script
qsub -l vmem=30gb,mem=12gb,nodes=1:ppn=4 -q P <path to script> -o <path to standard output dir> -e <path to standard error dir> -N <job name> -V -v "outdir=<path to results dir>, fq1=<path to R1.fq.gz>, fq2=<path to R2.fq.gz>, contigs=<path to contigs.fasta>"

# IDBA merge convert fastq R1 and R2 to fasta to use with idba assembler. Must unzip the fastq first.
/bin/fq2fa --merge <path to R1.fq> <path to R2.fq> <path to out merged fasta file>

# IDBA-UD. Basic operation
bin/idba_ud -r <path to merged fasta file> --num_threads 4 -o <path to results dir>

# Find all files with common name (e.g. generic name output of an analysis), perform a process on each file and rename the output with a meaningful name (e.g. the name of the sample used in the analysis, in this case the name of the dir)
for file in $(find . -name "some_name.tsv"); do base=$(dirname "$file"); f=$(basename $base); head -n 100 $file > $f".head_100.tsv"; done

# Use abacas to map contigs to a reference genome
for file in <path to contigs dir>/*.fasta; do perl ./abacas.1.3.1.pl -r <path to reference genome>/ref.fasta -q $file -p 'nucmer'; done

# Use quast to generate assembly statistics
python ~/quast-4.5/quast.py -o <outdir> -1 <R1.fastq.gz> -2 <R2.fastq.gz> contigs.fasta

# Linearize a fasta file using awk (from [Biostars](https://www.biostars.org/p/9262/))
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < file.fa > file.linear.fa

# Sum a column of numbers in a file using awk (https://stackoverflow.com/questions/450799/shell-command-to-sum-integers-one-per-line)
cut -d, -f2 file.csv | awk '{s+=$1} END {print s}'

# csvkit (pip install csvkit) to edit csv files in the command line
csvsql -t --query "select Name,NumReads from salmon_counts where NumReads >= 5;" salmon_counts.tsv | csvformat -T > salmon_counts_larger_than_5.tsv

csvjoin -t -c TXNAME,Name annotation_file.tsv salmon_counts_larger_than_5.tsv | csvformat -T > annotation_file_larger_than_5.tsv

csvcut -c TXNAME,GENEID -t annotation_file_larger_than_5.tsv | csvformat -T > tx2gene.txt

# parse string with awk and print/paste as columns (delimiter '_' e.g. vDNA_PO2015_JAR_663_CAAAAG_L007)
awk '{ split($0, a, "_"); printf "%s\t%s\t%s\t%s\t%s\n", $0, a[1], a[2], a[3], a[4] }' sample_name.txt > sample_data.txt
