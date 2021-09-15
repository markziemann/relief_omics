#!/bin/bash

for FQZ1 in *R1_001.fastq.gz ; do
  FQZ2=$(echo $FQZ1 | sed 's/_R1_/_R2_/' )
  ../sw/skewer/skewer -q 10 -t 8 $FQZ1 $FQZ2

  FQ1=$(echo $FQZ1 | sed 's/.gz$/-trimmed-pair1.fastq/')
  FQ2=$( echo $FQ1 | sed 's/pair1.fastq/pair2.fastq/' )
  ../sw/kallisto/kallisto quant \
  -i ../refgenome/gencode.v33.transcripts.fa.idx \
  -o ${FQ1}_kal -t 16 $FQ1 $FQ2

  rm $FQ1 $FQ2
done


for TSV in */*abundance.tsv ; do
  NAME=$(echo $TSV | cut -d '_' -f1) ; cut -f1,4 $TSV | sed 1d | sed "s/^/${NAME}\t/"
done > 3col.tsv
