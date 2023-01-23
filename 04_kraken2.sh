#!/bin/bash

if [ -d ../output/kraken_20220329 ]; then
	rm -rf ../output/kraken_20220329/
	mkdir ../output/kraken_20220329
else
	mkdir ../output/kraken_20220329
fi

for f in `ls ../output/assembled/*/scaffolding/k41.scafSeq`; do
	
	dn=`dirname ${f}`
	bn=`echo ${dn} | cut -d '/' -f 4`

	mkdir ../output/kraken_20220329/${bn}
	
	echo "[Kraken] Obtaining taxonomic classification of ${bn} assemble."

	kraken2 --db /data/db/kraken2/DB/ \
		--memory-mapping --output ../output/kraken_20220329/${bn}/${bn}.txt \
		--threads 30 ${f} \
		--confidence 0.05 \
		--report ../output/kraken_20220329/${bn}/${bn}.report.txt &> ../output/kraken_20220329/${bn}/${bn}.report.log.txt

	/usr/local/bioinfo/anaconda/bin/perl /usr/local/bioinfo/bioinfoutilities/kraken2ToMPA.pl -i ../output/kraken_20220329/${bn}/${bn}.txt  > ../output/kraken_20220329/${bn}/${bn}_labels.txt

done
