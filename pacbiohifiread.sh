#!/bin/bash
# Universitat Potsdam
# Author Gaurav Sablok
# date: 2024-2-22
# a read checker for the pacbiohifi reads making seqtk working for the hifi
# it is completely based on the regular expression and it is absolutely fast
# started writing it before lunch and stuck at a point and while having lunch solved it
read -r -p "pacbiohifi reads:" reads
if [[ "${reads}" == "" ]]; then
          echo "cant process with the fasta manipulations"
fi

if [[ -f "${reads}" ]]; then
              sed "s/@/>@/g" "${reads}" >"${reads%.*}".head.fastq
fi

          echo "storing ids on the run"
grep "^>" "${reads%.*}".head.fastq | 
                   sed "s/ /|/g" | cut -f 1 -d "|" | sed "s/>@//g" >fastaids.txt
cat fastaids.txt | while read line; do 
                    grep -A 1 "${line}" "${reads%.*}".head.fastq >"$line".pacbiohifi.fasta
done

echo "all the pacbiohifi reads are stored as independent fasta files"
cat *.pacbiohifi.fasta >assimilated.pacbiohifi.fasta
echo "the assimilated pacbiohifi fasta is present in pacbiohifi.fasta"
echo "the fastaids for the all pacbiohifi are present in the fastaids.txt"
echo "thank you for using the pacbiohifi fastq fasta read assimilator"
echo "UNIVERSITAT POTSDAM"
