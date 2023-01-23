#!/bin/bash

if [ -d ../output/kallisto ]; then
	rm -rf ../output/kallisto
	mkdir ../output/kallisto
	echo ""
else
	mkdir ../output/kallisto
fi

for f in `ls ../output/assembled/*/scaffolding/k41.scafSeq`; do
	
	dn=`dirname ${f}`
	bn=`echo ${dn} | cut -d '/' -f 4`
	
	if [ -d ../output/kallisto/${bn}/${i} ]; then
#                rm -rf ../output/kallisto/${bn}
#                mkdir -p ../output/kallisto/${bn}
 		echo "Index of ${bn} assemble already done!"
	else
		mkdir -p ../output/kallisto/${bn}
		kallisto index --index=../output/kallisto/${bn}/${bn}.index ${f}
	fi

    kallisto index --index=../output/kallisto/${bn}/${bn}.index ${f}

	for i in `grep ${bn} ./Sample_Details_LibrariesSequenced-01092021.txt | cut -f 2`; do

		lista=""	

		samples_R1=`grep ${i} ./Sample_Details_LibrariesSequenced-01092021.txt | cut -f 2 | sed 's/./..\/output\/processed\/prinseq\/&/' | sed 's/$/.atropos_final.prinseq_1.fastq/' | tr '\n' ' '`
		samples_R2=`grep ${i} ./Sample_Details_LibrariesSequenced-01092021.txt | cut -f 2 | sed 's/./..\/output\/processed\/prinseq\/&/' | sed 's/$/.atropos_final.prinseq_2.fastq/' | tr '\n' ' '`
		samples_SG1=`grep ${i} ./Sample_Details_LibrariesSequenced-01092021.txt | cut -f 2 | sed 's/./..\/output\/processed\/prinseq\/&/' | sed 's/$/.atropos_final.prinseq_1_singletons.fastq/' | tr '\n' ' '`
		samples_SG2=`grep ${i} ./Sample_Details_LibrariesSequenced-01092021.txt | cut -f 2 | sed 's/./..\/output\/processed\/prinseq\/&/' | sed 's/$/.atropos_final.prinseq_2_singletons.fastq/' | tr '\n' ' '`

		lista=${lista}${samples_R1}${samples_R2}${samples_SG1}${samples_SG2}

		if [ -d ../output/kallisto/${bn}/${i} ]; then
			rm -rf ../output/kallisto/${bn}/${i}
			mkdir -p ../output/kallisto/${bn}/${i}
			echo "[Kallisto] Counts of ${bn} assemble and ${i} libraries already done!"
		else
			mkdir -p ../output/kallisto/${bn}/${i}
	
			echo "[Kallisto] Obtaining counts of ${bn} assemble and ${i} libraries."
	
			kallisto quant --index=../output/kallisto/${bn}/${bn}.index --plaintext --threads=30 --verbose --bias ${lista} -o ../output/kallisto/${bn}/${i}/Kb_${bn}_${i} &> ../output/kallisto/${bn}/${i}/Kb_${bn}_${i}_log.txt
		
			rm -rf ../output/kallisto/${bn}/${i}/Kb_${bn}_${i}/bs*
	
		fi
		
	done
	
done
