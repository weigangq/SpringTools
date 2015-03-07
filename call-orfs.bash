#!/bin/bash

usage="SPRING-Tools Usage:\n****ALL OPTIONS ARE REQUIRED****\n-f [input] <FASTA FILE> \n-n [strain name] \n-i [genome ID] \n" 

while getopts ":i:f:n:h" opt; do
  case $opt in 
	h) printf "$usage" ; exit ;;
	f) f=$OPTARG ;; #fasta
	p) i=$OPTARG ;; #gid eg 1,2,3,4
	n) n=$OPTARG ;; #strain name eg DK
	\?) echo "Invalid option: -$OPTARG" >&2 ; printf "$usage" ; exit ;;
	esac 
done 

if [[ $f && ${f} ]] || [[ $i && ${i} ]]|| [[ $n && ${n} ]]
then
    printf "Running SPRINGDB pipeline for $f with gid $i and strain name $n\n"
else
    printf "$usage"
    exit
fi

glimmer long-orfs -n -t 1.15 -l $f gid_$i.longorfs 
glimmer extract -t $f gid_$i.longorfs > gid_$i.train
glimmer build-icm -r gid_$i.icm < gid_$i.train
glimmer glimmer3 -o 10 -t 30 -g 300 -l $f gid_$i.icm gid_$i
glimmer extract -w $f gid_$i.predict > gid_$i.fas
cat gid_$i.fas | sed -E 's/^>(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/>\1|'$n'_'$i'|gid_'$i'|\2|\3/' > gid_$i.fas2
head gid_$i.fas2 
~/dnatwizzer-code/biocode -t gid_$i.fas2 > gid_$i.pep

makeblastdb -in gid_$i.pep -dbtype 'prot' -out gid_$i 

printf "Forward Blasting \n" ;
blastp -query gid_$i.pep -db /home/shared/rrahman/feb_2015/springdb/gid_1/gid_1 -outfmt "6 qseqid sseqid pident length qlen evalue" -evalue 1e-5 -out gid_$i.fwd

printf "Reverse Blasting \n" ;

blastp -query /home/shared/rrahman/feb_2015/springdb/gid_1/gid_1.orth_id_for_blasting -db gid_$i -outfmt "6 qseqid sseqid pident length qlen evalue" -evalue 1e-5 -out gid_$i.rev

perl /scripts/check-reciprocal-with-options.pl gid_$i.fwd gid_$i.rev > gid_$i.orth 2> gid_$i.no_orth

perl /scripts/update_orth_with_reciprocol_top_hits.pl gid_$i.orth  > gid_$i.tophits 
perl /scripts/parse_predict_orf.pl gid_$i.fas2 > gid_$i.txt 

