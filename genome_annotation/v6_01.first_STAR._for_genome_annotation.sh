#!/bin/bash
#BSUB -q long
#BSUB -W 23:59
#BSUB -n 24
#BSUB -R "rusage[mem=2000]"
#BSUB -R "span[hosts=1]" 
#BSUB -J star1
#BSUB -o out/01_star1.out
#BSUB -e out/01_star1.err
#BSUB -B -u ye.duan@umassmed.edu
#BSUB -N -u ye.duan@umassmed.edu

cd /home/ye.duan-umw/pacbio/v6_annotation

# 1 First STAR index
ml star/2.7.10a

mkdir -p v6_index

STAR \
--runMode genomeGenerate \
--runThreadN 24 \
--genomeDir v6_index/ \
--genomeSAindexNbases 13 \
--genomeFastaFiles hofPB_v6.FINAL.fa

# loop started

inputdir="RNAseq/PE"

mkdir -p RNAseq/temp_fq
tempdir="RNAseq/temp_fq"

for file in $inputdir/*_R1.fq; do
  base=$(basename "$file" | sed -E 's/_R1\.fq//')

  read_raw_1="${inputdir}/${base}_R1.fq"
  read_raw_2="${inputdir}/${base}_R2.fq"  

  read_processed_1="${tempdir}/${base}_processed_R1.fq"
  read_processed_2="${tempdir}/${base}_processed_R2.fq"

 # 2 quality filter and size filter
 ml cutadapt
 cutadapt -m 15 -q 10 $read_raw_1 -o $read_processed_1
 cutadapt -m 15 -q 10 $read_raw_2 -o $read_processed_2

 # 3 first STAR mapping
 mkdir -p first_STAR_map/${base}

 ml star/2.7.10a
 STAR \
  --runThreadN 24 \
  --genomeDir v6_index \
  --readFilesIn $read_processed_1 $read_processed_2 \
  --outSAMattributes NH HI AS nM MD jM jI \
  --outFileNamePrefix first_STAR_map/${base}/${base}. \
  --outSAMtype BAM SortedByCoordinate
 
 # 4 move SJ_tab to download
 mkdir -p first_STAR_SJ
 cp first_STAR_map/${base}/${base}.SJ.out.tab first_STAR_SJ

 # 5 index first_BAM
 ml samtools/1.16.1
 samtools index \
  first_STAR_map/${base}/${base}.Aligned.sortedByCoord.out.bam

 # 6 report_1
 echo "$(date): ${base} first_STAR is done" >> v6_annotation.log

done

## Download the first.map.SJ dir to process SJ files in R




