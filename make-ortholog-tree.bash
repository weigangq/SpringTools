#!/bin/bash
mkdir results
perl /home/shared/rrahman/get_random_core.pl -o 100
for file in ./*.fas
do 
	fas=$( basename "$file" ".fas")
	bioseq -t1 $fas.fas | muscle -out $fas.aln -clwstrict
	/home/shared/rrahman/align_coding_seq.pl $fas.aln $fas.fas > $fas.cds
done
/home/shared/rrahman/AlignConcat.V2.pl ./*.cds > alignedconcat.txt

leng=$(bioaln -i 'fasta' -l alignedconcat.txt)
echo "DNA, gene1codon1 = 1-$leng\3" >> part.txt
echo "DNA, gene1codon2 = 2-$leng\3" >> part.txt 
echo "DNA, gene1codon3 = 3-$leng\3" >> part.txt 

raxmlHPC-SSE3 -f a -x 12345 -p 12345 -# 100 -s alignedconcat.txt -m GTRCAT -n pa_tree -q part.txt

mv RAxML* results/
