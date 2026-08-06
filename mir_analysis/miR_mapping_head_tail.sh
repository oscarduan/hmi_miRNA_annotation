#!/bin/bash
#BSUB -q long
#BSUB -W 24:00
#BSUB -n 16
#BSUB -R "rusage[mem=2000]"
#BSUB -J reg
#BSUB -o out/regPB.out
#BSUB -e out/regPB.err
#BSUB -B -u ye.duan@umassmed.edu
#BSUB -N -u ye.duan@umassmed.edu

# Name every FASTQ input in format of A1_H00.full.fastq by mv

set -euo pipefail
shopt -s nullglob
cd /home/ye.duan-umw/sRNA.seq

GENOME_FA=/home/ye.duan-umw/RNAseq/ref/hofPB_v6.FINAL.fa
STAR_INDEX=/home/ye.duan-umw/pacbio/v6_index_second
GTF_PATH=/home/ye.duan-umw/sRNA.seq/refs/miRNA_final.gtf

STAR_OUTPUT=/home/ye.duan-umw/sRNA.seq/Star.alignment/regen_half_PB
COUNT_OUTPUT=/home/ye.duan-umw/sRNA.seq/subread/regen_half_PB

TEMP_DIR=/home/ye.duan-umw/sRNA.seq/temp/regPB
mkdir -p "$TEMP_DIR"

mkdir -p log
LOGFILE=log/log.260605.regeneration

mkdir -p "$STAR_OUTPUT" "$COUNT_OUTPUT" "$TEMP_DIR" out log

#1. Extract file name
for files in /home/ye.duan-umw/sRNA.seq/input/240128_regeneration/filtered/*.filtered.fastq

do

file=$(basename "$files" .filtered.fastq)

cd /home/ye.duan-umw/sRNA.seq

# #2 Extract UMI
# ml umi_tools/1.1.2
# umi_tools extract \
# -p NNNNNNNNNNNN \
# --3prime \
# -I input/240128_regeneration/$file\.full.fastq \
# -S temp/$file\.umi.fastq

# echo "$(date): $file\_umi_extract is done." >> log/log.260605.regeneration

# #3 trim off the adaptor
# ml cutadapt
# cutadapt -a AACTGTAGGCACC temp/$file\.umi.fastq \
# -m 18 -M 26  \
# -q 10 \
# -e 0.2 \
# -o temp/$file\.cutadapt.fastq  # Here -o cannot link to input/, so I made next line 

# #4 filter ncRNA contamination
# module load bowtie2/2.4.1
# bowtie2 --end-to-end  -N 0 -D 15 -R 2 --no-1mm-upfront  \
# -S temp/#file\.nc.sam -L 15 -i S,1,1.15 -x refs/nc.v22 \
# -q temp/$file\.cutadapt.fastq \
# --un input/240128_regeneration/filtered/$file\.filtered.fastq \
# 2> log/bowtie2.log/$file\.filter.log

# echo "$(date): $file\_filtering is done." >> log/log.260605.regeneration
# rm temp/#file\.nc.sam

#5 STAR mapping
ml star/2.7.10a

mkdir -p "$STAR_OUTPUT/${file}"

STAR_TMP="$TEMP_DIR/${file}.STAR_tmp"
rm -rf "$STAR_TMP"        # remove leftovers ONLY

STAR \
  --runThreadN 16 \
  --genomeDir "$STAR_INDEX" \
  --readFilesIn "$files" \
  --outFilterMultimapNmax 100 \
  --outSAMattributes NH HI AS nM \
  --alignSJoverhangMin 30 \
  --outTmpDir "$STAR_TMP" \
  --outFileNamePrefix "$STAR_OUTPUT/${file}/${file}." \
  --outSAMtype BAM SortedByCoordinate

rm -rf "$STAR_TMP"

echo "$(date): $file star is done" >> "$LOGFILE"

#6 dedup the BAM file
ml samtools
  samtools sort -n \
    -o "$TEMP_DIR/${file}.STAR.name_sorted.bam" \
    "$STAR_OUTPUT/${file}/${file}.Aligned.sortedByCoord.out.bam"

  samtools view -h \
    -o "$TEMP_DIR/${file}.STAR.name_sorted.sam" \
    "$TEMP_DIR/${file}.STAR.name_sorted.bam"

  echo "$(date '+%F %T')  [$file] dedup-prep finished" >> "$LOGFILE"

  ml umi_tools/1.1.2
  umi_tools dedup \
  --method unique \
  --extract-umi-method=read_id \
  -I "$TEMP_DIR/${file}.STAR.name_sorted.sam" \
  -S "$TEMP_DIR/${file}.STAR.dedup.bam"

samtools sort \
  -o "$STAR_OUTPUT/${file}/${file}.dedup.sorted.bam" \
  "$TEMP_DIR/${file}.STAR.dedup.bam"

samtools index "$STAR_OUTPUT/${file}/${file}.dedup.sorted.bam"

  rm "$TEMP_DIR/${file}.STAR.name_sorted.sam"
  rm "$TEMP_DIR/${file}.STAR.name_sorted.bam"
  rm "$TEMP_DIR/${file}.STAR.dedup.bam"
  
  echo "$(date '+%F %T')  [$file] UMI dedup finished" >> "$LOGFILE"


#8 featureCounts
  ml subread

  featureCounts \
    -t miRNA \
    -M -s 1 -O --fraction \
    -g miRNA \
    -a "$GTF_PATH" \
    --extraAttributes gene_id,pre_miRNA,strand \
    -o "$COUNT_OUTPUT/${file}.miRNA.txt" \
    "$STAR_OUTPUT/${file}/${file}.dedup.sorted.bam"


  echo "$(date '+%F %T')  [$file] featureCounts finished" >> "$LOGFILE"
  echo "$(date '+%F %T')  [$file] DONE" >> "$LOGFILE"
  echo "----------------------------------------" >> "$LOGFILE"

done

echo "===== RUN END $(date '+%F %T') =====" >> "$LOGFILE"