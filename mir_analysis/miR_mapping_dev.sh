#!/bin/bash
#BSUB -q long
#BSUB -W 47:59
#BSUB -n 16
#BSUB -R "rusage[mem=8000]"
#BSUB -J star
#BSUB -o out/star.out
#BSUB -e out/star.err
#BSUB -B -u ye.duan@umassmed.edu
#BSUB -N -u ye.duan@umassmed.edu

set -euo pipefail
shopt -s nullglob
cd /home/ye.duan-umw/sRNA.seq

GENOME_FA=/home/ye.duan-umw/RNAseq/ref/hofPB_v6.FINAL.fa
STAR_INDEX=/home/ye.duan-umw/pacbio/v6_index_second
BOWTIE_INDEX=/home/ye.duan-umw/sRNA.seq/refs/bowtie
GTF_PATH=/home/ye.duan-umw/sRNA.seq/refs/ZIV_full_relcoords.gtf

TEMP_DIR=/home/ye.duan-umw/sRNA.seq/PB_temp
INPUT=/home/ye.duan-umw/sRNA.seq/input/filtered

STAR_OUTPUT=/home/ye.duan-umw/sRNA.seq/Star.alignment/PB_map
COUNT_OUTPUT=/home/ye.duan-umw/sRNA.seq/subread/PB_STAR_map

mkdir -p "$TEMP_DIR" "$COUNT_OUTPUT" "$STAR_OUTPUT" out

LOGFILE="out/star_pipeline.log"
echo "===== RUN START $(date '+%F %T') =====" >> "$LOGFILE"

# Step 1: Loop over all filtered fastq files
for files in "$INPUT"/*.filtered.fastq; do
  filename=$(basename "$files")
  file=${filename%%.filtered.fastq}

  echo "Processing sample: $file"
  echo "$(date '+%F %T')  [$file] START" >> "$LOGFILE"

  # Step 2: STAR mapping
  mkdir -p "$STAR_OUTPUT/${file}"

 STAR_TMP="$TEMP_DIR/${file}.STAR_tmp"
rm -rf "$STAR_TMP"        # remove leftovers ONLY

ml star/2.7.10a
STAR \
  --runThreadN 16 \
  --genomeDir "$STAR_INDEX" \
  --readFilesIn "$files" \
  --outFilterMultimapNmax 100 \
  --outSAMattributes NH HI AS nM \
  --alignSJoverhangMin 30 \
  --outTmpDir "$STAR_TMP" \
  --outFileNamePrefix "$STAR_OUTPUT/${file}/${file}.STAR." \
  --outSAMtype BAM SortedByCoordinate


  rm -rf "$STAR_TMP"
  echo "$(date '+%F %T')  [$file] STAR mapping finished" >> "$LOGFILE"

  # Step 3: Dedup
  ml samtools
  samtools sort -n \
    -o "$TEMP_DIR/${file}.STAR.name_sorted.bam" \
    "$STAR_OUTPUT/${file}/${file}.STAR.Aligned.sortedByCoord.out.bam"

  samtools view -h \
    -o "$TEMP_DIR/${file}.STAR.name_sorted.sam" \
    "$TEMP_DIR/${file}.STAR.name_sorted.bam"

  echo "$(date '+%F %T')  [$file] dedup-prep finished" >> "$LOGFILE"

  ml umi_tools/1.1.2
  umi_tools dedup \
    --method unique \
    --extract-umi-method=read_id \
    -I "$TEMP_DIR/${file}.STAR.name_sorted.sam" \
    -S "$STAR_OUTPUT/${file}/${file}.STAR.dedup.bam"

  rm "$TEMP_DIR/${file}.STAR.name_sorted.sam"
  rm "$TEMP_DIR/${file}.STAR.name_sorted.bam"

  echo "$(date '+%F %T')  [$file] UMI dedup finished" >> "$LOGFILE"

  ml samtools
  samtools index "$STAR_OUTPUT/${file}/${file}.STAR.Aligned.sortedByCoord.out.bam"
  samtools index "$STAR_OUTPUT/${file}/${file}.STAR.dedup.bam"

  echo "$(date '+%F %T')  [$file] BAM indexing finished" >> "$LOGFILE"

  # Step 4: featureCounts
  ml subread
  featureCounts \
    -t hairpin \
    -M -s 1 -O --fraction \
    -g gene_id \
    -a "$GTF_PATH" \
    -o "$COUNT_OUTPUT/${file}.hairpin.txt" \
    "$STAR_OUTPUT/${file}/${file}.STAR.Aligned.sortedByCoord.out.bam" \
    "$STAR_OUTPUT/${file}/${file}.STAR.dedup.bam"

  featureCounts \
    -t mature \
    -M -s 1 -O --fraction \
    -g gene_id \
    -a "$GTF_PATH" \
    -o "$COUNT_OUTPUT/${file}.mature.txt" \
    "$STAR_OUTPUT/${file}/${file}.STAR.Aligned.sortedByCoord.out.bam" \
    "$STAR_OUTPUT/${file}/${file}.STAR.dedup.bam"

  featureCounts \
    -t star \
    -M -s 1 -O --fraction \
    -g gene_id \
    -a "$GTF_PATH" \
    -o "$COUNT_OUTPUT/${file}.star.txt" \
    "$STAR_OUTPUT/${file}/${file}.STAR.Aligned.sortedByCoord.out.bam" \
    "$STAR_OUTPUT/${file}/${file}.STAR.dedup.bam"

  echo "$(date '+%F %T')  [$file] featureCounts finished" >> "$LOGFILE"
  echo "$(date '+%F %T')  [$file] DONE" >> "$LOGFILE"
  echo "----------------------------------------" >> "$LOGFILE"

done

echo "===== RUN END $(date '+%F %T') =====" >> "$LOGFILE"
