#https://pubmlst.org/data/dbases.xml
python ../src/mlst_populate.py
cat *-alleles.fa > allAlleles.fa
../bin/cd-hit-est -i allAlleles.fa -o temp90.fa -c 0.90 -n 8 -aS 0.8 -g 1 -p 1 -d 0 >db_90.log
python ../src/paddingFasta.py temp90.fa
../bin/makeblastdb -in allAlleles.fa -dbtype nucl
../bin/bwa index temp90_padded.fa
../bin/samtools faidx temp90_padded.fa
cat *-stfile.tsv|grep -v stid > allSTs.tsv
