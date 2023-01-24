#!/bin/bash
#
#              INGLÊS/ENGLISH
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  http://www.gnu.org/copyleft/gpl.html
#
#
#             PORTUGUÊS/PORTUGUESE
#  Este programa é distribuído na expectativa de ser útil aos seus
#  usuários, porém NÃO TEM NENHUMA GARANTIA, EXPLÍCITAS OU IMPLÍCITAS,
#  COMERCIAIS OU DE ATENDIMENTO A UMA DETERMINADA FINALIDADE.  Consulte
#  a Licença Pública Geral GNU para maiores detalhes.
#  http://www.gnu.org/copyleft/gpl.html
#
#  Copyright (C) 2019  Universidade Estadual Paulista "Júlio de Mesquita Filho"
#
#  Universidade Estadual Paulista "Júlio de Mesquita Filho" (UNESP)
#  Faculdade de Ciências Agrárias e Veterinárias (FCAV)
#  Laboratório de Bioinformática (LB)
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com
#  http://www.fcav.unesp.br 
#


# input - diretório contendo os arquivos de entrada no formato .fastq
input=$1

# validação do parâmetro "input"
if [ ! ${input} ]
then   
        echo "[ERROR] Missing input directory." 1>&2
        exit
else   
        if [ ! -d ${input} ]
        then   
                echo "[ERROR] Wrong input directory (${input})." 1>&2
                exit
        fi
fi

# output - diretório para armazenar o resultado do processo de montagem
output=$2

# validação do parâmetro "output"
if [ ! ${output} ]
then   
        echo "[ERROR] Missing output directory." 1>&2
        exit
else   
        if [ ! -d ${output} ]
        then   
                echo "[ERROR] Wrong output directory (${output})." 1>&2
                exit
        fi
fi

# Número de CORES para o processamento
# ATENÇÃO: Não exceder o limite da máquina
THREADS=$3

if [ ! ${THREADS} ]; then
	echo "[ERROR] Missing number of threads." 1>&2
	exit
fi

libfile=$4
if [ ! ${libfile} ]; then
	echo "[ERROR] Missing library file with three columns (SAMPLE_ID, SAMPLE_NAME, ASSEMBLY_NAME)" 1>&2
	exit
else
        if [ ! -e ${libfile} ]
        then   
                echo "[ERROR] Wrong library file (${libfile})." 1>&2
                exit
        fi
		
fi

# Quantidade de memória para o processamento com Jellyfish
# ATENÇÃO: Não exceder o limite da máquina
MEM=$5
MEMPARAM=""
if [ ! ${MEM} ]; then
	echo "[WARNING] Using maximum available memory." 1>&2
else	
	MEMPARAM=" --memory ${MEM} "
fi

# Arquivos e diretórios de saída (output) 

for asmname in $(cut -f 3 ${libfile} | sort -u); do
	
	echo "Collecting libraries for the assembly (${asmname}) ..." 1>&2

	asmdir_out="${output}/${asmname}"

	left=()
	right=()
	singleton=()
	singleton_left=()
	singleton_right=()
	
	for libname in $(grep -P "\t${asmname}" ${libfile} | cut -f 2); do
		echo -e "\tLibrary (${libname}) ..." 1>&2
		
		leftlib=()
		rightlib=()
		leftlib_singleton=()
		rightlib_singleton=()

		shopt -s nullglob
		leftlib=(${input}/${libname}.*_1.fastq)
		leftlib_singleton=(${input}/${libname}.*_1_singletons.fastq)
		rightlib_singleton=(${input}/${libname}.*_2_singletons.fastq)
		shopt -u nullglob # Turn off nullglob to make sure it doesn't interfere with anything later
		
		if [ "${#leftlib[@]}" -gt 0 ]; then
			echo -e "\t\tFound left PE reads ${#leftlib[@]}" 1>&2
			
			for lr in "${leftlib[@]}"; do
				#echo ">${lr}<"
				if [ ! -e ${lr}.renamed ]; then
					reallr=`readlink -f ${lr}`
					echo -e "\t\tRenaming reads on ${reallr} ..."
					#awk '{ if (NR%4==1) { if ($1!~/\/1$/) { print $1"/1" } else { print $1 } } else if (NR%4==3) { print "+" } else { print $1 } }' ${reallr} > ${reallr}.temp
					#mv ${reallr}.temp ${reallr}
					touch ${lr}.renamed
				fi
				rr=`echo ${lr} | sed 's/_1.fastq/_2.fastq/'`
				#echo ">${rr}<"
				if [ -e ${rr} ]; then
					if [ ! -e ${rr}.renamed ]; then
						realrr=`readlink -f ${rr}`
						echo -e "\t\tRenaming reads on ${realrr} ..."
						#awk '{ if (NR%4==1) { if ($1!~/\/2$/) { print $1"/2" } else { print $1 } } else if (NR%4==3) { print "+" } else { print $1 } }' ${realrr} > ${realrr}.temp
						#mv ${realrr}.temp ${realrr}
						touch ${rr}.renamed
					fi						
					rightlib+=("${rr}")
				fi
			done
			left+=("${leftlib[@]}")
		else
			echo -e "\t\tNot Found left PE reads" 1>&2
		fi
		
		if [ "${#rightlib[@]}" -gt 0 ]; then
			echo -e "\t\tFound right PE reads ${#rightlib[@]}" 1>&2
			right+=("${rightlib[@]}")
		else
			echo -e "\t\tNot Found right PE reads" 1>&2
		fi
		
		if [ "${#leftlib_singleton[@]}" -gt 0 ]; then
			echo -e "\t\tFound left singleton reads ${#leftlib_singleton[@]}" 1>&2
			singleton_left+=("${leftlib_singleton[@]}")
			singleton+=("${leftlib_singleton[@]}")
		else
			echo -e "\t\tNot Found left singleton reads" 1>&2
		fi

		if [ "${#rightlib_singleton[@]}" -gt 0 ]; then
			echo -e "\t\tFound right singleton reads ${#rightlib_singleton[@]}" 1>&2
			singleton_right+=("${rightlib_singleton[@]}")
			singleton+=("${rightlib_singleton[@]}")
		else
			echo -e "\t\tNot Found right singleton reads" 1>&2
		fi
		
	done

	rm -fr ${asmdir_out}
	
	echo "LEFT.....: " $(IFS=, ; echo "${left[*]}") 1>&2
	echo "RIGHT....: " $(IFS=, ; echo "${right[*]}") 1>&2
	
	UPARAM=""
	
	if [ "${#singleton[@]}" -gt 0 ]; then
		#echo "UNPAIRED.: " $(IFS=, ; echo "${singleton[*]}") 1>&2
		UPARAM=" -r $(IFS=, ; echo "${singleton[*]}")"
	fi
		
	echo "DATASET: ${asmname}"
	echo "Performing metagenomic assembly ..."

	megahit --num-cpu-threads ${THREADS} --presets meta-sensitive ${MEMPARAM} \
		-1 $(IFS=, ; echo "${left[*]}") \
		-2 $(IFS=, ; echo "${right[*]}") \
		${UPARAM} -o ${asmdir_out} &> ${output}/${asmname}.log.txt
		
	rm -fr ${asmdir_out}/intermediate_contigs/

	echo "LEFT singletons..: $(IFS=, ; echo "${singleton_left[*]}")"
	echo "RIGHT singletons.: $(IFS=, ; echo "${singleton_right[*]}")"
	
	echo "Performing SOAP scaffolding using SOAPdenovo-fusion ..."

	mkdir -p ${asmdir_out}/scaffolding
	
	./mkSOAPconfig.pl -pe1 $(IFS=, ; echo "${left[*]}") \
			  -pe2 $(IFS=, ; echo "${right[*]}") \
			  -se1 $(IFS=, ; echo "${singleton_left[*]}") \
			  -se2 $(IFS=, ; echo "${singleton_right[*]}") > ${asmdir_out}/scaffolding/config
	
	cur_dir=`pwd`
	abs_asmdir_out=`readlink -f ${asmdir_out}`
	cd ${asmdir_out}/scaffolding/
	
	SOAPdenovo-fusion -D -s config -p ${THREADS} -K 41 -g k41 -c ${abs_asmdir_out}/final.contigs.fa &> ${abs_asmdir_out}/scaffolding/SOAPdenovo-fusion.log.txt
	SOAPdenovo-127mer map -s config -p ${THREADS} -g k41 &> ${abs_asmdir_out}/scaffolding/SOAPdenovo-127mer-map.log.txt
	SOAPdenovo-127mer scaff -p 40 -g k41 -F &> ${abs_asmdir_out}/scaffolding/SOAPdenovo-127mer-scaff.log.txt
	
	cd ${cur_dir}

	#echo "Performing BWA alignment reads X assembly ..."

	mkdir -p ${output}/aligned
	bwa index ${asmdir_out}/final.contigs.fa &> ${asmdir_out}/bwa_index.log.txt

	for i in $(seq 0 $[${#left[@]}-1]); do 
		bn=`basename ${left[${i}]} | sed 's/\..*$//'`
		bwa mem -M -t ${THREADS} ${asmdir_out}/final.contigs.fa \
		${left[${i}]} ${right[${i}]} | samtools view --threads ${THREADS} -b -h -F 260 - | samtools sort --threads ${THREADS} -n -o  ${output}/aligned/${bn}.bam
	done
	
	echo "Performing metaquast evaluation (contigs) ..."
	
	evaldir_out="${output}/${asmname}/metaquast"
	
	mkdir -p ${evaldir_out}
	
       metaquast.py    ${asmdir_out}/final.contigs.fa \
                       --output-dir  ${evaldir_out} \
       		--pe1 $(IFS=, ; echo "${left[*]}") \
       		--pe2 $(IFS=, ; echo "${right[*]}") \
       		--single $(IFS=, ; echo "${singleton[*]}") \
                       --min-contig 1000 \
                       --threads ${THREADS} \
                       --space-efficient \
                       -L \
                       --gene-finding \
                       --rna-finding \
			--conserved-genes-finding \
			--no-icarus \
			--reuse-combined-alignments \
			--ambiguity-usage all \
			--ambiguity-score 1 \
			--max-ref-number 0 \
                       --plots-format svg &> ${evaldir_out}.log.txt

	echo "Performing metaquast evaluation (scaffolds) ..."
	
	scaff_evaldir_out="${output}/${asmname}/scaffolding/metaquast"
	
	mkdir -p ${evaldir_out}
	
       metaquast.py    ${asmdir_out}/scaffolding/k41.scafSeq \
                       --output-dir  ${scaff_evaldir_out} \
       		--pe1 $(IFS=, ; echo "${left[*]}") \
       		--pe2 $(IFS=, ; echo "${right[*]}") \
       		--single $(IFS=, ; echo "${singleton[*]}") \
                       --min-contig 1000 \
                       --threads ${THREADS} \
                       --space-efficient \
                       -L \
                       --gene-finding \
                       --rna-finding \
			--conserved-genes-finding \
			--no-icarus \
			--reuse-combined-alignments \
			--ambiguity-usage all \
			--ambiguity-score 1 \
			--max-ref-number 0 \
                       --plots-format svg &> ${evaldir_out}.log.txt

done
