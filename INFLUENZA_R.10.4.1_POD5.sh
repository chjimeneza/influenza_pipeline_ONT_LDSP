#!/bin/bash

echo ""
echo -e "\e[93m##############################################################################################\e[0m"
echo -e "\e[93m##\e[0m             \e[92mLABORATORIO DEPARTAMENTAL DE SALUD PUBLICA DE ANTIOQUIA                      \e[93m##\e[0m"
echo -e "\e[93m##\e[0m                       \e[92m....................................                               \e[93m##\e[0m"
echo -e "\e[93m##\e[0m                       \e[92m....................................                               \e[93m##\e[0m"
echo -e "\e[93m##\e[0m    \e[92m                     IDENTIFICACION DE CLADOS DE MONKEYPOX                            \e[93m##\e[0m"
echo -e "\e[93m##\e[0m                       \e[92m....................................                               \e[93m##\e[0m"
echo -e "\e[93m##\e[0m                       \e[92m....................................                               \e[93m##\e[0m"
echo -e "\e[93m##\e[0m                       \e[92m....................................                               \e[93m##\e[0m"
echo -e "\e[93m##\e[0m                       \e[92m....................................                               \e[93m##\e[0m"
echo -e "\e[93m##\e[0m             \e[92mDirección: Carrera 72 A 78 B – 141, Robledo - Medellín                       \e[93m##\e[0m"
echo -e "\e[93m##\e[0m                            \e[92mTeléfono: (+574) 3835400                                      \e[93m##\e[0m"
echo -e "\e[93m##\e[0m             \e[92mCorreo institucional: laboratoriodep@antioquia.gpv.co                        \e[93m##\e[0m"
echo -e "\e[93m##\e[0m                                \e[92mNIT: 890900286-0                                          \e[93m##\e[0m"
echo -e "\e[93m##############################################################################################\e[0m"

echo ""
echo ""
echo ""
echo -e "\e[93m############################\e[0m \e[92mLas variables a ingresar son las siguientes\e[0m  \e[93m################################\e[0m"
echo ""

#echo -e "\e[92m---> Ingrese la ruta donde estan los archivos POD5:\e[0m "
#read POD5

echo -e "\e[92m---> Ingrese el numero del primer barcode con el que se llevo a cabo la secuenciación:\e[0m "
read bar_start

echo -e "\e[92m---> Ingrese el numero del ultimo barcode con el que se llevo a cabo la secuenciación:\e[0m "
read bar_end

echo -e "\e[92m---> Ingresa el nombre del archivo de salida (Ejemplo: 08102025 ó 08102025_MPXV):\e[0m "
read NOMBRE

echo -e "\e[92m---> Ingresa la ruta (Ejemplo: /mnt/HPC-DISK-210/MonkeyPox_Results/) donde quiere que se guarden los archivos:\e[0m "
read RUTA_FINAL

prefijo="_Influenza"
#reference="/root/artic-ncov2019/primer_schemes/mpxv/V2/mpxv.reference.fasta"
RUN_DIR="${RUTA_FINAL}/${NOMBRE}"
#DirScheme="/root/artic-ncov2019/primer_schemes"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate base

cd "${RUTA_FINAL}"
mkdir -p "${NOMBRE}/fastq"
cd "${RUTA_FINAL}/${NOMBRE}/fastq"
OUT_FASTQ="$(pwd)"

echo ""
echo "#############################################################################"
echo "###   Basecaller con Dorado de POD5 $date ==> $(date)  ##"
echo "#############################################################################"

/mnt/HPC-DISK-210/programas/dorado-1.1.1-linux-x64/bin/dorado basecaller /mnt/HPC-DISK-210/programas/dorado-1.1.1-linux-x64/bin/dna_r10.4.1_e8.2_400bps_hac@v4.2.0/ "${POD5}" > "${RUTA_FINAL}/${NOMBRE}/fastq/${NOMBRE}${prefijo}_raw.bam"

echo ""
echo "#############################################################################################"
echo "###   Demultiplexación del archivo concatenado .BAM $date ==> $(date)  ##"
echo "#############################################################################################"

/mnt/HPC-DISK-210/programas/dorado-1.1.1-linux-x64/bin/dorado demux \
        -o "${RUTA_FINAL}/${NOMBRE}/fastq/" \
        --kit-name SQK-NBD114-96 \
        --emit-fastq "${RUTA_FINAL}/${NOMBRE}/fastq/${NOMBRE}${prefijo}_raw.bam"

echo ""
echo "####################################################################"
echo "###   Renombrar archivos fastq $date ==> $(date)  ##"
echo "####################################################################"

for n in $(seq -w "${bar_start}" "${bar_end}"); do
        archivo_origen=$(ls "${OUT_FASTQ}"/*_SQK-NBD114-96_barcode"${n}".fastq 2>/dev/null)

        if [ -f "$archivo_origen" ]; then
                nuevo_nombre="${NOMBRE}${prefijo}_barcode${n}.fastq"
                echo -e "\e[94mRenombrando:\e[0m $archivo_origen → $nuevo_nombre"
                mv "$archivo_origen" "$nuevo_nombre"
        else
                echo -e "\e[91m[ERROR]\e[0m No se encontró archivo para barcode${n}"
        fi

done

echo -e "\e[92mProceso completado. Archivos renombrados guardados en:\e[0m ${OUT_FASTQ}"

echo "-------------------> RUTA ACTUAL: ${OUT_FASTQ}"

echo "-------------------> RUTA ACTUAL: ${RUTA_FINAL}/${NOMBRE}"

echo ""

echo ""
echo "####################################################################"
echo "###   Mapeo de lecturas a referencia $date ==> $(date)  ##"
echo "####################################################################"

	mini_align \
		-i reads.fastq.gz \
		-r reference.fasta \
		-p align_tmp \
		-t 2 \
		-m

echo ""
echo "####################################################################"
echo "###   Mapeo de lecturas a referencia $date ==> $(date)  ##"
echo "####################################################################"

    	samtools view --write-index -F 4 align_tmp.bam -o align.bam##idx##align.bam.bai

    	stats_from_bam -o align.bamstats -s align.bam.summary -t 2 align.bam

	samtools depth -aa align.bam -Q 20 -q 1 > depth.txt

	samtools faidx reference.fasta

	cut -f1-2 reference.fasta.fai > regions.txt

	upper=`echo \$((\${length}+(\${length}*70/100)))`
        lower=`echo \$((\${length}-(\${length}*70/100)))`

        samtools view -bh align.bam \${region} > \${region}.bam;

      # ignore regions with no reads
        count=`samtools view -c \${region}.bam`

        if [ "\${count}" -eq "0" ];
        then
          echo "no reads in \${region} so continuing"
          cp \${region}.bam \${region}_all.bam
          continue;
        fi

        lines=( 5000 / 2 )


























