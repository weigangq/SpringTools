#! /bin/bash
usage="\nBefore running please specify path to pipeline, glimmer and blast\n\****ALL OPTIONS ARE REQUIRED****\nUSAGE:\n-f [input] <FASTA FILE> \n-i [icm] \n-o [name] \n-p [GID] \n-r [type] \n-h [get predicted orfs] <0 or 1> \n \n"

while getopts ":i:o:f:p:r:c:n:" opt; do
    case $opt in 
	i) i=$OPTARG ;; #icm
	o) o=$OPTARG ;; #name eg gid_4
	f) f=$OPTARG ;; #fasta
	p) p=$OPTARG ;; #gid eg 1,2,3,4
	r) r=$OPTARG ;; #type
	c) c=$OPTARG ;;
	n) n=$OPTARG ;; #name eg DK
	esac 
done 


glimmer long-orfs -n -t 1.15 -l $f $o.longorfs 
glimmer extract -t $f $o.longorfs > $o.train
glimmer build-icm -r $o.icm < $o.train
glimmer glimmer3 -o 10 -t 30 -g 300 -l $f $o.icm $o
glimmer extract -w $f $o.predict > $o.fas
cat $o.fas | sed -E 's/^>(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/>\1|'$n'_'$p'|gid_'$p'|\2|\3/' > $o.fas2
head $o.fas2 
~/dnatwizzer-code/biocode -t $o.fas2 > $o.pep

makeblastdb -in $o.pep -dbtype 'prot' -out $o 

printf "Forward Blasting \n" ;
blastp -query $o.pep -db /home/shared/rrahman/feb_2015/springdb/gid_1/gid_1 -outfmt "6 qseqid sseqid pident length qlen evalue" -evalue 1e-5 -out $o.fwd

printf "Reverse Blasting \n" ;

blastp -query /home/shared/rrahman/feb_2015/springdb/gid_1/gid_1.orth_id_for_blasting -db $o -outfmt "6 qseqid sseqid pident length qlen evalue" -evalue 1e-5 -out $o.rev

perl /home/shared/rrahman/check-reciprocal-with-options.pl $o.fwd $o.rev > $o.orth 2> $o.no_orth

perl /home/shared/rliang/scripts/Parsers/update_orth_with_reciprocol_top_hits.pl $o.orth  > $o.tophits 
perl /home/shared/rliang/scripts/Parsers/parse_predict_orf.pl $o.fas2 > $o.txt 

