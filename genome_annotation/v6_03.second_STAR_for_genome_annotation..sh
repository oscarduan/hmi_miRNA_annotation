#!/bin/bash
#BSUB -q long
#BSUB -W 47:59
#BSUB -n 8
#BSUB -R "rusage[mem=4000]"
#BSUB -R "span[hosts=1]" 
#BSUB -J star2.stringtie.new
#BSUB -o out/03_star2.stringtie.new.out
#BSUB -e out/03_star2.stringtie.new.err
#BSUB -B -u ye.duan@umassmed.edu
#BSUB -N -u ye.duan@umassmed.edu

cd /home/ye.duan-umw/pacbio/

# 1 index again
# ml star/2.7.10a

# mkdir -p v6_index_second

# STAR \
# --runMode genomeGenerate \
# --sjdbFileChrStartEnd first_STAR_SJ/processed_txt/$(ls first_STAR_SJ/processed_txt/*.txt) \
# --runThreadN 8 \
# --genomeDir v6_index_second/ \
# --genomeSAindexNbases 13 \
# --genomeFastaFiles hofPB_v6.FINAL.fa

# loop starts
inputdir="RNAseq/PE"

mkdir -p RNAseq/temp_fq
tempdir="RNAseq/temp_fq"

for file in $inputdir/*_R1.fq; do
  base=$(basename "$file" | sed -E 's/_R1\.fq//')

  read_raw_1="${inputdir}/${base}_R1.fq"
  read_raw_2="${inputdir}/${base}_R2.fq"  

  read_processed_1="${tempdir}/${base}_processed_R1.fq"
  read_processed_2="${tempdir}/${base}_processed_R2.fq"

# 1.4 Clean raw FASTQ using seqtk-compatible awk
mkdir -p RNAseq/cleaned_fq
cleaned_fq_dir="RNAseq/cleaned_fq"
cleaned_1="${cleaned_fq_dir}/${base}_cleaned_R1.fq"
cleaned_2="${cleaned_fq_dir}/${base}_cleaned_R2.fq"

# Remove bad reads with mismatched sequence/quality lengths
awk 'NR%4==1{h=$0} NR%4==2{s=$0} NR%4==3{p=$0} NR%4==0{q=$0; if(length(s)==length(q)) print h"\n"s"\n"p"\n"q}' $read_raw_1 > $cleaned_1
awk 'NR%4==1{h=$0} NR%4==2{s=$0} NR%4==3{p=$0} NR%4==0{q=$0; if(length(s)==length(q)) print h"\n"s"\n"p"\n"q}' $read_raw_2 > $cleaned_2

# 1.5 Quality and length filter with cutadapt
ml cutadapt
cutadapt -m 15 -q 10 $cleaned_1 -o $read_processed_1
cutadapt -m 15 -q 10 $cleaned_2 -o $read_processed_2


# 2 second_STAR
 mkdir -p second_STAR_map/${base}

 ml star/2.7.10a
 STAR \
  --outSAMstrandField intronMotif \
  --runThreadN 8 \
  --genomeDir v6_index_second \
  --readFilesIn $read_processed_1 $read_processed_2 \
  --outSAMattributes NH HI AS nM MD jM jI \
  --outFileNamePrefix second_STAR_map/${base}/${base}. \
  --outSAMtype BAM SortedByCoordinate
 
#  # 3 move SJ_tab to download
#  mkdir -p second_STAR_SJ
#  cp second_STAR_map/${base}/${base}.SJ.out.tab second_STAR_SJ

 # 4 index second_BAM
 ml samtools/1.16.1
 samtools index \
  second_STAR_map/${base}/${base}.Aligned.sortedByCoord.out.bam

 # 5 remove processed file in #2 
 rm $read_processed_1 
 rm $read_processed_2

 # 6 report_2
  echo "$(date): ${base} second_STAR is done" >> v6.241224.log

done

# 7 StringTie individual
ml stringtie/2.2.1

cd /home/ye.duan-umw/pacbio/
mkdir -p stringtie

ST_dir="second_STAR_map"

for dir in $ST_dir/*; do
   base=$(basename $dir)

   inputfile1="${ST_dir}/${base}/${base}.Aligned.sortedByCoord.out.bam" 

 stringtie \
  -p 12 \
  -o stringtie/${base}.stringtie.gtf \
  $inputfile1

 echo "$(date): stringtie ${base} is done." >> v6.241224.log

done

# 8 StringTie merge (need to modify for final)
stringtie --merge -o stringtie.merged.gtf stringtie/*.stringtie.gtf

echo "$(date): stringtie merge is done."  >> v6.241224.log



