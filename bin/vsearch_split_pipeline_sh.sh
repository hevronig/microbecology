fasta_extension=".fna"
lines_per_splitted_file=250000
for file in *."$fasta_extension"; do bi=$(basename "$file" "$fasta_extension"); time vsearch --derep_fulllength "$file" --fasta_width 0 --relabel_keep --relabel_md5 --maxseqlength 50000000 --output "$file".derep_fl.fa; done > runVsearch_derep_full.log 2>&1

for file in *.derep_fl.fa; do bi=$(basename "$file" .derep_fl.fa); vsearch --derep_prefix "$file" --fasta_width 0 --relabel_keep --relabel_md5 --maxseqlength 50000000 --minseqlength 500 --output "$file".derep_prefix.fa; done > runVsearch_derep_prefix.log 2>&1

for file in *.derep_prefix.fa; do bi=$(basename $file .derep_prefix.fa); vsearch --shuffle $i --output "$bi".shuffled.fasta --randseed 13 --maxseqlength 50000000 --fasta_width 0; done > runVsearchShuffle.log 2>&1

for file in *.shuffled.fasta; do bi=$(basename $file .shuffled.fasta); mkdir "$bi"; split -l "$lines_per_splitted_file" --numeric-suffixes=1 --additional-suffix=.split.fa $file "$bi".; mv "$bi".*.split.fa "$bi"; done > runFastaSplit.log 2>&1