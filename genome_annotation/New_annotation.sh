#!/bin/bash

# Set base path
path=/home/duan/pacbio/new_annotation
cd $path

### Part 1: Adding false negative from AQ to new strT 
## 1.1 Duplicate unstranded features with both + and - strands to augment the final GTF with stranded information
cd $path

# Extract unstranded and stranded features separately
awk '$7 == "."' stringtie_merged_s3.75.gtf > unstranded_features_s375.gtf
awk '$7 != "."' stringtie_merged_s3.75.gtf > stranded_features_s375.gtf

# Add strand information (+ and -) to the unstranded features
# NOTE: Using sed to safely replace the 7th field without altering formatting of attributes (V9)
sed 's/^\(\([^\t]*\t\)\{6\}\)[^\t]*/\1+/' unstranded_features_s375.gtf > unstranded_plus_s375.gtf
sed 's/^\(\([^\t]*\t\)\{6\}\)[^\t]*/\1-/' unstranded_features_s375.gtf > unstranded_minus_s375.gtf

# Concatenate stranded, plus-strand-augmented, and minus-strand-augmented GTFs
cat stranded_features_s375.gtf unstranded_plus_s375.gtf unstranded_minus_s375.gtf > \
stringtie.merged_s375_strand_augmented.gtf

## 1.2 GFFcmp of AQ referencing new strT_aug
mkdir -p Gffcmp.out

gffcompare \
    ../AQ_s2g_updated.gtf \
    -r ./stringtie.merged_s375_strand_augmented.gtf\
    -o ./Gffcmp.out/gffcmp_s375.aq


## 1.3 GFFcmp of strT_s375 referencing AQ for class_o
gffcompare \
    ./stringtie_merged_s3.75.gtf \
    -r ../AQ_s2g_updated.gtf\
    -o ./Gffcmp.out/gffcmp__0_s375.aq


## 1.4 Generate AQ gtf subset for strT-false negativbes

# This can also be done Rstudio S24A, which has scripts for graphs and update DEMO
# Below is only for WSL:
cd $path
nano new_annotation_add.R

## Add below to the R script:
#########################################################################################
#!/usr/bin/env Rscript
setwd("/home/duan/pacbio/new_annotation")

library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(ggplot2)
library(rtracklayer)

`%notin%` <- Negate(`%in%`)

# Load and preprocess GTF
gtf <- rtracklayer::import("./stringtie_merged_s3.75.gtf") %>% as.data.frame()
gtf.short <- gtf %>% select(10:11)

# Load featureCounts output
AQ.fc <- read.table("./AQ_merged_featureCounts.txt", header = TRUE) %>%
  mutate(
    Length_kb = Length / 1000,
    RPK = SUM_read / Length_kb,
    TPM = RPK / (sum(RPK, na.rm = TRUE) / 1e6)
  ) %>%
  select(-Length_kb, -RPK)

# Parse AQ tracking file
aq.tracking.df <- read.table("Gffcmp.out/gffcmp_s375.aq.tracking",
                             header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")

aq.tracking_parsed <- aq.tracking.df %>%
  mutate(
    strT_gene_id = if_else(V3 == "-", "-", str_extract(V3, "^[^|]+")),
    strT_t_id    = if_else(V3 == "-", "-", str_extract(V3, "(?<=\\|)[^|]+")),
    AQ_id        = if_else(V5 == "-", "-", str_extract(V5, "(?<=q1:)[^|]+\\|[^|]+"))
  ) %>%
  select(6:8, V4) %>%
  group_by(AQ_id) %>%
  summarize(
    strT_gene_id = strT_gene_id,
    strT_t_id = strT_t_id,
    class_code = paste(unique(V4), collapse = ","),
    class_equal = any(V4 == "="),
    class_match = any(V4 %in% c("=", "j", "k", "c", "m", "n"))
  )

# Class o processing
class_o.df0 <- aq.tracking_parsed %>% filter(class_code == "o")
class_o.count <- class_o.df0 %>% pull(AQ_id) %>% unique() %>% length()

strT_0.tracking.df <- read.table("Gffcmp.out/gffcmp__0_s375.aq.tracking",
                                 header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")

strT_0.tracking_parse1 <- strT_0.tracking.df %>%
  mutate(
    strT_gene_id = if_else(V5 == "-", "-", str_extract(V5, "(?<=:)[^|]+")),
    strT_transcript_id = if_else(V5 == "-", "-", str_extract(V5, "(?<=:)[^|]+\\|[^|]+") %>% str_extract("[^|]+$")),
    AQ_id = if_else(V3 == "-", "-", str_extract(V3, "^[^|]+\\|[^|]+"))
  ) %>%
  select(6:8, 4) %>%
  dplyr::rename(strT_t_id = strT_transcript_id, AQ_id_0 = AQ_id)

class_o.df <- class_o.df0 %>%
  left_join(strT_0.tracking_parse1 %>% select(2:4), by = "strT_t_id") %>%
  mutate(
    only_match = AQ_id_0 == AQ_id,
    other_match_good = V4 %in% c("=", "m", "n", "j", "k", "c", "i", "y", "p")
  )

class_o_not_add <- class_o.df %>% filter(!only_match & other_match_good)
class_o_add <- class_o.df %>% filter(AQ_id %notin% class_o_not_add$AQ_id)

class_match1.df <- aq.tracking_parsed %>%
  filter(class_code %in% c("i", "y", "p", "u"))

new_add.df <- bind_rows(class_o_add %>% select(1:6), class_match1.df) %>%
  left_join(AQ.fc %>% select(Geneid, TPM) %>% dplyr::rename(AQ_id = Geneid), by = "AQ_id")

aq_refed.list <- strT_0.tracking_parse1 %>%
  filter(V4 %in% c("=", "m", "n", "j", "k", "c")) %>%
  pull(AQ_id_0)

new_add.filtered <- new_add.df %>%
  filter(TPM >= 0.5) %>%
  filter(AQ_id %notin% aq_refed.list)

new_add.summary <- new_add.filtered %>%
  group_by(class_code) %>%
  summarise(
    gene_number = n(),
    mean_TPM = mean(TPM)
  ) %>%
  ungroup()

write.csv(new_add.filtered, "./new_add_AQ.csv", row.names = FALSE)

# subset AQ GTF
new_add.filtered <- read.csv("./new_add_AQ.csv", header = T)

aq.input.gtf <- read.table("../AQ_s2g_updated.gtf",
                           header = FALSE,
                           sep = "\t",
                           stringsAsFactors = FALSE,
                           fill = TRUE,
                           quote = "")

aq_parsed <- aq.input.gtf %>%
  mutate(
    transcript_id = str_extract(V9, 'transcript_id "[^"]+"') %>% str_remove_all('transcript_id "|"'),
    read_name = str_extract(V9, 'read_name "[^"]+"') %>% str_remove_all('read_name "|"'),
    query_alignment_start = str_extract(V9, 'query_alignment_start "[^"]+"') %>% str_remove_all('query_alignment_start "|"'),
    query_alignment_end = str_extract(V9, 'query_alignment_end "[^"]+"') %>% str_remove_all('query_alignment_end "|"'),
    cigar = str_extract(V9, 'cigar "[^"]+"') %>% str_remove_all('cigar "|"')
  )

add_gene.df <- aq_parsed %>%
  filter(transcript_id %in% as.vector(new_add.filtered$AQ_id) | str_detect(V9, "98015634|98037435") ) %>%
  mutate(gene_id = transcript_id) %>%
  mutate(transcript_id = paste(transcript_id, ".1", sep = ""),
         gene_name = read_name,
         gene_biotype = "protein_coding",
         class_code = "-",
         inherit_ref = FALSE,
         ref_name = read_name,
         Info = "")  %>%
  mutate(
    new_V9 = paste0(
      'gene_id "', gene_id, '"; ',
      'transcript_id "', transcript_id, '"; '
    )
  ) %>%
  select(V1, V2, V3, V4, V5, V6, V7, V8, new_V9)

write.table(add_gene.df, file = "./AQ_add.gtf", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#######################################################################################

## 1.5 # Merge two GTFs and sort by genomic coordinate using bedtools

# 1. Concatenate AQ_add.gtf and stringtie_merged_s3.75.gtf (excluding header lines)
# 2. Pipe into bedtools sort for proper chromosome and position sorting
# 3. Output to merged_sorted.gtf

(
  cat AQ_add.gtf
  grep -v "^#" stringtie_merged_s3.75.gtf
) | bedtools sort -i - > stringtie_merged_1.gtf



# Part 2 TransDecoder for ORF prediction and homology:
# Local WSL
  
# 2.1 get Stringtie FA
cd $path

/home/duan/pacbio/transdecoder/TransDecoder/util/gtf_genome_to_cdna_fasta.pl \
 /home/duan/pacbio/new_annotation/stringtie_merged_1.gtf \
 /home/duan/pacbio/hofPB_v6.FINAL.fa \
 > /home/duan/pacbio/new_annotation/YD_transcripts_transdecoder_1.fasta 

# 2.2 GTF->GFF3
# 2.2.1 YD: direct convert
/home/duan/pacbio/transdecoder/TransDecoder/util/gtf_to_alignment_gff3.pl \
 /home/duan/pacbio/new_annotation/stringtie_merged_1.gtf \
 > /home/duan/pacbio/new_annotation/stringtie.merged.gff3

# 2.3 Generate gene-to-transcript file (this is not included in the wiki)
awk '$3 == "transcript" {
  match($0, /gene_id "([^"]+)"/, gene);
  match($0, /transcript_id "([^"]+)"/, transcript);
  if (gene[1] != "" && transcript[1] != "")
    print gene[1] "\t" transcript[1];
}' /home/duan/pacbio/new_annotation/stringtie_merged_1.gtf > \
/home/duan/pacbio/new_annotation/gene_trans_map.txt

#  2.4 Generate the longest ORF prediction
mkdir -p /home/duan/pacbio/new_annotation/transdecoder
cd /home/duan/pacbio/new_annotation/transdecoder

/home/duan/pacbio/transdecoder/TransDecoder/TransDecoder.LongOrfs \
 -t /home/duan/pacbio/new_annotation/YD_transcripts_transdecoder_1.fasta \
 --gene_trans_map /home/duan/pacbio/new_annotation/gene_trans_map.txt

# 2.5 Homology search_blastp 
## download https://www.uniprot.org/proteomes/UP000005640
## human proteomes is UP000005640
##  2.5.1 generate db
cd /home/duan/pacbio/transdecoder
makeblastdb -in /home/duan/pacbio/transdecoder/UP000005640_9606.fasta -dbtype prot -out ref/uniprot_sprot_db

## 2.5.2 BLASTP
# The following blastp can be run at SCI, about 3X faster. see v6_11_1
blastp \
  -query /home/duan/pacbio/new_annotation/transdecoder/YD_transcripts_transdecoder_1.fasta.transdecoder_dir/longest_orfs.pep  \
  -db /home/duan/pacbio/transdecoder/ref/uniprot_sprot_db \
  -max_target_seqs 1 \
  -outfmt 6 \
  -evalue 1e-3 \
  -num_threads 16 \
  > /home/duan/pacbio/new_annotation/transdecoder/YD.blastp.outfmt6

# 2.6 Generate best candidate ORF prediction:
cd /home/duan/pacbio/new_annotation/transdecoder

  /home/duan/pacbio/transdecoder/TransDecoder/TransDecoder.Predict \
   -t /home/duan/pacbio/new_annotation/YD_transcripts_transdecoder_1.fasta \
   --retain_blastp_hits /home/duan/pacbio/new_annotation/transdecoder/YD.blastp.outfmt6

# 2.7 Propagate ORFs to genome
cd $path

## 2.7.1 Propagate:
/home/duan/pacbio/transdecoder/TransDecoder/util/cdna_alignment_orf_to_genome_orf.pl \
     ./transdecoder/YD_transcripts_transdecoder_1.fasta.transdecoder.gff3 \
     ./stringtie.merged.sorted.gff3 \
     ./YD_transcripts_transdecoder_1.fasta \
     > YD.transdecoder.genome.gff3 \
  2> ./transdecoder/YD.genome_mapping_250530.log

## 2.7.2  Sort the output gff3 by genome coordinate
(grep '^#' YD.transdecoder.genome.gff3; grep -v '^#' YD.transdecoder.genome.gff3 | grep -v '^[[:space:]]*$' | sort -k1,1 -k4,4n) > YD.transdecoder.genome.sorted.gff3

## 2.7.3 Fix the .ID= issue known by TransDecoder
sed -i 's/\.ID=/\tID=/g' YD.transdecoder.genome.sorted.gff3

## 2.7.4 Fix the ^chr1^+ suffix
# sed -E '
#   s/(ID=MSTRG\.[0-9]+)\^.*?\^[-+]/\1/g;
#   s/(Parent=MSTRG\.[0-9]+)\^.*?\^[-+]/\1/g
# ' YD.transdecoder.genome.sorted.gff3 > YD.transdecoder.genome.cleaned.gff3

sed -E '
  s/([; \t](ID|Parent)=[^; \t\^]+)\^[^\^]+\^[-+]/\1/g
' YD.transdecoder.genome.sorted.gff3 > YD.transdecoder.genome.cleaned.gff3

## 2.7.5 remove the .pN of mRNA features and corresponding Parant=
# sed -E '
#   s/(ID=MSTRG\.[0-9]+\.[0-9]+)\.p[0-9]+/\1/g;
#   s/(Parent=MSTRG\.[0-9]+\.[0-9]+)\.p[0-9]+/\1/g
# ' YD.transdecoder.genome.cleaned.gff3 > YD.transdecoder.genome.final.gff3

sed -E '
  s/([; \t](ID|Parent)=[^; \t]+)\.p[0-9]+/\1/g
' YD.transdecoder.genome.cleaned.gff3 > YD.transdecoder.genome.final.gff3


## 2.8 Process to GTF 
### 2.8.1 change Name= to Info= for tagging
sed -E 's/;Name="/;Info="/g' YD.transdecoder.genome.final.gff3 > YD.transdecoder.genome.gff3.tagged

## 2.8.2 Change score= to score: so AGAT won't mis-interpret
cp YD.transdecoder.genome.gff3.tagged YD.transdecoder.genome.gff3.tagged.bak
sed -i 's/score=/score_/g' YD.transdecoder.genome.gff3.tagged

## 2.8.3 AGAT to convert gff3 to gtf
# only run once for the 1st time
apptainer pull agat_1.1.0.sif docker://quay.io/biocontainers/agat:1.1.0--pl5321hdfd78af_0 # only run once

apptainer exec agat_1.1.0.sif \
  agat_convert_sp_gff2gtf.pl \
  --gff YD.transdecoder.genome.gff3.tagged \
  -o YD.transdecoder.genome.tagged_0.gtf

########### 2.8.2.0 (alternative) Change score= to score: so AGAT won't mis-interpret
# cd /n/boslfs02/LABS/srivastava_lab/Lab/oscar/

# cp YD.transdecoder.genome.gff3.tagged YD.transdecoder.genome.gff3.tagged.bak
# sed -i 's/score=/score_/g' YD.transdecoder.genome.gff3.tagged

# salloc -p test -c 1 -t 00-07:00 --mem=32000

# singularity exec /cvmfs/singularity.galaxyproject.org/a/g/agat:1.1.0--pl5321hdfd78af_0 \
#  agat_convert_sp_gff2gtf.pl \
#   --gff YD.transdecoder.genome.gff3.tagged \
#   -o YD.transdecoder.genome.tagged.gtf 

# # Then download
######################################################################

### 2.8.4 remove original biotype flag:
sed -E 's/[[:space:]]*original_biotype[^;]*;[[:space:]]*//g' YD.transdecoder.genome.tagged_0.gtf > YD.transdecoder.genome.tagged.gtf

### 2.8.5  Remove Info from AGAT output to avoid the "" parsing issue
sed -E '
  s/;[[:space:]]*Info[[:space:]]*".*""//g
' YD.transdecoder.genome.tagged.gtf > YD.transdecoder.genome.cleaned.gtf

### 2.8.6.1 Extract gene_id and Name/Info flags from the original GFF3
awk '
BEGIN { FS=OFS="\t" }

$3 == "gene" {
  # Extract attributes from $9 (keep spaces inside quoted strings)
  attr = $9
  id = ""; name = ""

  split(attr, fields, ";")
  for (i in fields) {
    gsub(/^[ \t]+|[ \t]+$/, "", fields[i])  # trim field
    if (fields[i] ~ /^ID=/) {
      split(fields[i], a, "=")
      id = a[2]
    } else if (fields[i] ~ /^Name=/) {
      match(fields[i], /Name="(.*)"/, b)
      name = b[1]
    }
  }

  if (id != "" && name != "") {
    print id, name
  }
}' YD.transdecoder.genome.final.gff3 > gene_info_map.txt

### 2.8.6.2 Dedup the map
awk -F'\t' '{
  split($2, a, ",");
  seen = ""; out = "";
  for (i in a) {
    if (!match(seen, "(^|,)" a[i] "($|,)")) {
      seen = seen ? seen "," a[i] : a[i];
      out = out ? out "," a[i] : a[i];
    }
  }
  print $1 "\t" out
}' gene_info_map.txt > gene_info_map.dedup.txt

### 2.8.6.3 Inject Info attribute into GTF
awk '
BEGIN {
  FS=OFS="\t"
  while ((getline < "gene_info_map.dedup.txt") > 0) {
    gene2info[$1] = $2
  }
}

{
  if ($0 ~ /^#/) { print; next }

  n = split($9, a, /;[ \t]*/)
  attr = ""
  gene_id = ""
  for (i = 1; i <= n; i++) {
    if (a[i] ~ /^[ \t]*$/) continue
    match(a[i], /^[ \t]*([^ ]+) "(.*)"/, kv)
    if (kv[1] != "" && kv[2] != "") {
      attr = attr kv[1] " \"" kv[2] "\"; "
      if (kv[1] == "gene_id") gene_id = kv[2]
    }
  }

  if (gene_id in gene2info) {
    attr = attr "Info \"" gene2info[gene_id] "\";"
  }

  $9 = attr
  print
}' YD.transdecoder.genome.cleaned.gtf > YD.transdecoder.genome.with_info.gtf


# 2.9 Process not propagated ORF & no ORF genes
cd $path

## 2.9.1 Isolate transcript IDs
### 2.9.1a Get all transcript IDs from FASTA
grep "^>" YD_transcripts_transdecoder_1.fasta | sed 's/>//' | awk '{print $1}' > all_transcript_ids.txt

### 2.9.1b Get predicted ORF transcript IDs from BED file
cut -f1 transdecoder/YD_transcripts_transdecoder_1.fasta.transdecoder.bed | grep -v '^track' | sort -u > predicted_ids.txt

### 2.9.1c Get genome-mapped ORF transcript IDs from final genome GFF3
awk -F'\t' '$3 == "CDS" && $9 ~ /transcript_id/ {
  match($9, /transcript_id "([^"]+)"/, a);
  if (a[1] != "") print a[1];
}' YD.transdecoder.genome.with_info.gtf | sort -u > mapped_ids.txt

### 2.9.1d Transcripts that do NOT have an ORF prediction
comm -23 <(sort all_transcript_ids.txt) predicted_ids.txt > YD_no_orf_transcripts.txt
### 2.9.1e Transcripts that HAVE an ORF prediction but were NOT mapped to genome
comm -23 predicted_ids.txt mapped_ids.txt > YD_orf_not_mapped.txt

## 2.9.2 Subset GTFs of the un-propagated transcripts
# Extract lines for transcripts WITHOUT any ORF prediction
awk 'FNR==NR {ids[$1]; next}
     {
       match($0, /transcript_id "([^"]+)"/, a);
       if (a[1] in ids) print $0;
     }' YD_no_orf_transcripts.txt ./stringtie_merged_1.gtf > stringtie.no_orf.gtf

# Extract lines for transcripts WITH predicted ORFs that FAILED to map to genome
awk 'FNR==NR {ids[$1]; next}
     {
       match($0, /transcript_id "([^"]+)"/, a);
       if (a[1] in ids) print $0;
     }' YD_orf_not_mapped.txt ./stringtie_merged_1.gtf > stringtie.orf_not_mapped.gtf

## 2.9.3 Add gene_biotype
### 2.9.3a "protein_coding"
awk '
BEGIN { FS=OFS="\t" }
{
  if ($0 ~ /^#/) { print; next }

  # Convert attributes field into key-value pairs
  n = split($9, a, /;[ \t]*/)
  attr = ""
  for (i = 1; i <= n; i++) {
    if (a[i] ~ /^[ \t]*$/) continue
    match(a[i], /^[ \t]*([^ ]+) "(.*)"/, kv)
    if (kv[1] != "" && kv[2] != "") {
      attr = attr kv[1] " \"" kv[2] "\"; "
    }
  }

  # Add gene_biotype safely
  attr = attr "gene_biotype \"protein_coding\";"
  $9 = attr
  print
}'  YD.transdecoder.genome.with_info.gtf > YD.transdecoder.genome.with_biotype.gtf

### 2.9.3b "non-coding"
awk '
BEGIN { FS=OFS="\t" }
{
  if ($0 ~ /^#/) { print; next }
  n = split($9, a, /;[ \t]*/)
  attr = ""
  for (i = 1; i <= n; i++) {
    if (a[i] ~ /^[ \t]*$/) continue
    match(a[i], /^[ \t]*([^ ]+) "(.*)"/, kv)
    if (kv[1] != "" && kv[2] != "") {
      attr = attr kv[1] " \"" kv[2] "\"; "
    }
  }
  attr = attr "gene_biotype \"non-coding\";"
  $9 = attr
  print
}' stringtie.no_orf.gtf > stringtie.no_orf.with_biotype.gtf

### 2.9.3c "protein_coding_antisense"
awk '
BEGIN { FS=OFS="\t" }
{
  if ($0 ~ /^#/) { print; next }
  n = split($9, a, /;[ \t]*/)
  attr = ""
  for (i = 1; i <= n; i++) {
    if (a[i] ~ /^[ \t]*$/) continue
    match(a[i], /^[ \t]*([^ ]+) "(.*)"/, kv)
    if (kv[1] != "" && kv[2] != "") {
      attr = attr kv[1] " \"" kv[2] "\"; "
    }
  }
  attr = attr "gene_biotype \"protein_coding_antisense\";"
  $9 = attr
  print
}' stringtie.orf_not_mapped.gtf > stringtie.orf_not_mapped.with_biotype.gtf

# 2.9.4 Concatenate all three GTFs and sort by chromosome and start position
cat \
  YD.transdecoder.genome.with_biotype.gtf \
  stringtie.no_orf.with_biotype.gtf \
  stringtie.orf_not_mapped.with_biotype.gtf \
  > merged.biotype.gtf

sort -k1,1V -k4,4n merged.biotype.gtf > merged.biotype.sorted.gtf

java -jar IGVTools/igvtools.jar index merged.biotype.sorted.gtf

# 2.10 backup 
cp merged.biotype.sorted.gtf ./stringtie_s3.75.transdecoder.final.gtf

java -jar IGVTools/igvtools.jar index stringtie_s3.75.transdecoder.final.gtf



# Part 3 Inherit AQ's gene_name
## 3.1 GFFcmp of strT transdecoder referencing AQ for name inheritance
cd $path

gffcompare \
    ./stringtie_s3.75.transdecoder.final.gtf \
    -r ../AQ_s2g_updated.gtf\
    -o ./Gffcmp.out/gffcmp__0_s375.aq

## 3.2 Inherit gene name by R
nano inherit.AQ.R

###
####
#####
######
#!/usr/bin/env Rscript

# Set working directory to absolute path (WSL-compatible)
setwd("/home/duan/pacbio/new_annotation")

library(dplyr)
library(Biostrings)
library(stringr)
library(tibble)
library(readr)
library(GenomicRanges)
library(rtracklayer)

mkdir -p DEMO

`%notin%` <- Negate(`%in%`)

# 1 load initial input ----
yd.input.gtf <- read.table("stringtie_s3.75.transdecoder.final.gtf",
                           header = FALSE,
                           sep = "\t",
                           stringsAsFactors = FALSE,
                           fill = TRUE,
                           quote = "")
aq.input.gtf <- read.table("AQ_s2g_updated.gtf",
                           header = FALSE,
                           sep = "\t",
                           stringsAsFactors = FALSE,
                           fill = TRUE,
                           quote = "")

strT_parsed <- yd.input.gtf %>%
  mutate(
    transcript_id = str_extract(V9, 'transcript_id "[^"]+"') %>% str_remove_all('transcript_id "|"'),
    gene_id = str_extract(V9, 'gene_id "[^"]+"') %>% str_remove_all('gene_id "|"'),
    gene_biotype = str_extract(V9, 'gene_biotype "[^"]+"') %>% str_remove_all('gene_biotype "|"'),
    Info = str_extract(V9, 'Info "[^"]+"') %>% str_remove_all('Info "|"')
  )

aq_parsed <- aq.input.gtf %>%
  mutate(
    transcript_id = str_extract(V9, 'transcript_id "[^"]+"') %>% str_remove_all('transcript_id "|"'),
    read_name = str_extract(V9, 'read_name "[^"]+"') %>% str_remove_all('read_name "|"'),
    query_alignment_start = str_extract(V9, 'query_alignment_start "[^"]+"') %>% str_remove_all('query_alignment_start "|"'),
    query_alignment_end = str_extract(V9, 'query_alignment_end "[^"]+"') %>% str_remove_all('query_alignment_end "|"'),
    cigar = str_extract(V9, 'cigar "[^"]+"') %>% str_remove_all('cigar "|"')
  )

# 2 process gff_compare output gtf ----
tracking.df <- read.table("Gffcmp.out/gffcmp__0_s375.aq.tracking",
                          header = FALSE,
                          sep = "\t",
                          stringsAsFactors = FALSE,
                          fill = TRUE,
                          quote = "") 

tracking_parse <- tracking.df %>%
  mutate(
    strT_gene_id = case_when(
      V5 == "-" ~ "-",
      str_detect(V5, "MSTRG") ~ str_extract(V5, "(?<=:)[^|]+"),
      TRUE ~ str_extract(V5, "(?<=:)([^|]+\\|[^|]+)")
    ),
    
    strT_transcript_id = case_when(
      V5 == "-" ~ "-",
      str_detect(V5, "MSTRG") ~ str_extract(V5, "(?<=:)[^|]+\\|[^|]+") %>% str_extract("[^|]+$"),
      TRUE ~ str_extract(V5, "(?<=:)([^|]+\\|){2}[^|]+\\|[^|]+") %>% 
        str_extract("([^|]+\\|[^|]+$)")
    ),
    
    AQ_id = if_else(
      V3 == "-", "-",
      str_extract(V3, "^[^|]+\\|[^|]+")
    )
  )

# 3 process by class_code ----
class_by_loci <- tracking_parse %>%
  filter(strT_gene_id != "-") %>%
  group_by(strT_gene_id) %>%
  summarize(
    class_code = paste(unique(V4), collapse = ","),
    class_equal = any(V4 == "="),
    t_number = n_distinct(V1),
    AQ_t_number = n_distinct(AQ_id),
    AQ_orf = any(AQ_id != "-"),
    YD_orf = any(strT_gene_id != "-"),
    class_match = any(V4 %in% c("=", "j", "k", "c", "m", "n")),
    AQ_id = AQ_id,
    YD_id = strT_gene_id
  )

# 4 Translate gene id ----
# One AQ matches multi strT----
multi.strT.count <- class_by_loci %>%
  filter(class_match == TRUE) %>%
  filter(AQ_id != "-") %>%
  group_by(AQ_id) %>%
  summarise(strT_per_AQ = n_distinct(strT_gene_id)) %>%
  ungroup() %>%
  filter(strT_per_AQ > 1)

multi.strT.list <- multi.strT.df %>% pull(strT_gene_id) %>% unique() 

multi.strT.df <- class_by_loci %>% 
  filter(class_match == TRUE) %>%
  filter(AQ_id %in% as.vector(multi.strT.count$AQ_id)) %>%
  distinct()

# extract the AQ_id to search the TransDecoder.gff3
multi.AQ.count <- class_by_loci %>%
  filter(strT_gene_id %notin% as.vector(multi.strT.df$strT_gene_id)) %>% # not one AQ gene shared by two strT genes
  filter(class_match == TRUE) %>%
  filter(AQ_id != "-") %>%
  group_by(strT_gene_id) %>%
  summarise(AQ_per_strT = n_distinct(AQ_id)) %>%
  ungroup() %>%
  filter(AQ_per_strT > 1)

AQ_id_pool <- class_by_loci %>% 
  filter(strT_gene_id %in% as.vector(multi.AQ.count$strT_gene_id)) %>%
  distinct() %>%
  pull(AQ_id)

# Extract from the modified AQ_transdecoder.gtf (see v6_07.sh)
AQ.transdecoder.gtf <- read_tsv("AQ.modified.transdecoder.gtf",
                                comment = "##",
                                col_names = FALSE,
                                col_types = cols(.default = "c"),
                                quote = "") # avoid issues with extra quotes in column 9

AQ.transdecoder_parsed <- AQ.transdecoder.gtf %>%
  filter(X3 == "gene") %>%
  transmute(
    gene_id = str_extract(X9, 'gene_id \\"([^\\"]+)\\"') %>% 
      str_remove('gene_id \\"') %>% 
      str_remove('\"$'),
    trD.score = str_extract(X9, 'score_([^\\"]+)') %>% 
      str_remove('score_'),
    ORF_type = str_extract(X9, 'ORF type:[^ ]+') %>% 
      str_remove('ORF type:') %>% 
      str_trim(),
    length = str_extract(X9, 'len:([0-9]+)') %>% 
      str_remove('len:') %>% 
      as.numeric(),
    homolog = str_extract(X9, 'sp\\|([^\\|]+\\|[^\\"]+)'),
    accession = if_else(!is.na(homolog), str_extract(homolog, "^sp\\|[^\\|]+"), NA_character_),
    homolog_gene = if_else(!is.na(homolog), str_extract(homolog, "(?<=\\|)[^\\|]+(?=_HUMAN)"), NA_character_),
    perc_identity = if_else(!is.na(homolog),
                            str_extract(homolog, "_HUMAN\\|([0-9.]+)") %>% 
                              str_remove("_HUMAN\\|") %>% 
                              as.numeric(),
                            NA_real_),
    E_value = if_else(
      !is.na(homolog),
      str_extract(homolog, "(0|[0-9.]+e[-+][0-9]+)") %>% as.numeric(),
      NA_real_
    )
  )           

AQ_ref.df <- AQ.transdecoder_parsed %>% 
  dplyr::rename(AQ_id = gene_id) %>%
  left_join(class_by_loci %>%
              filter(class_match == TRUE),
            by = "AQ_id") %>%
  select(1:4,6:12,17) %>%
  filter(!is.na(strT_gene_id))

multi.strT.ref <- AQ_ref.df %>% filter(AQ_id %in% as.vector(multi.strT.df$AQ_id)) %>% distinct()

# 4.1.1 One strT matches multi AQ ----
multi.AQ.count <- class_by_loci %>%
  filter(strT_gene_id %notin% as.vector(multi.strT.df$strT_gene_id)) %>% # not one AQ gene shared by two strT genes
  filter(class_match == TRUE) %>%
  filter(AQ_id != "-") %>%
  group_by(strT_gene_id) %>%
  summarise(AQ_per_strT = n_distinct(AQ_id)) %>%
  ungroup() %>%
  filter(AQ_per_strT > 1)

# match name with this policy:
# Perform a TransDecoder and BLASTP for AQ transcriptome against PacBio genome
# a) if AQ genes are not propagated to genome, use the propagated AQ gene.
AQ_best_ref0 <- AQ_ref.df %>% 
  filter(AQ_id %in% AQ_id_pool) %>%
  filter(strT_gene_id %in% as.vector(multi.AQ.count$strT_gene_id)) %>%
  distinct() 

  # b) if AQ genes have complete ORF versus other types, choose complete ORF
AQ_best_ref1 <- AQ_best_ref0 %>%
  group_by(strT_gene_id) %>%
  filter(if (any(ORF_type == "complete")) ORF_type == "complete" else TRUE) %>%
  ungroup() 

  # c) if still multiple, choose the AQ gene with lower BLASTP E-value
AQ_best_ref2 <- AQ_best_ref1 %>%
  group_by(strT_gene_id) %>%
  filter(
    if (n() == 1) TRUE
    else if (all(is.na(E_value))) TRUE
    else if (all(E_value == E_value[1], na.rm = TRUE)) TRUE
    else E_value == min(E_value, na.rm = TRUE)
  ) %>%
  ungroup() 

# d) if still multiple, choose the AQ gene with higher TransDecoder Propagation score
AQ_best_ref3 <- AQ_best_ref2 %>%
  group_by(strT_gene_id) %>%
  mutate(trD.score.num = as.numeric(trimws(trD.score))) %>%
  filter(
    if (n() == 1) TRUE
    else trD.score.num == max(trD.score.num, na.rm = TRUE)
  ) %>%
  select(-trD.score.num) %>%
  ungroup()

# e) if still multiple, choose 9820| or later version
AQ_best_ref4 <- AQ_best_ref3 %>%
  group_by(strT_gene_id) %>%
  filter(
    if (n() == 1) TRUE
    else {
      aq_prefix <- as.numeric(substr(AQ_id, 1, 8))
      aq_prefix == max(aq_prefix, na.rm = TRUE)
    }
  ) %>%
  ungroup()

# 4.1.3 one-to-one match ----
class_equal.df <- class_by_loci %>% 
  filter(class_equal == TRUE) %>%
  filter(YD_orf == TRUE & AQ_orf == TRUE) %>%
  filter(AQ_id != "-") %>%
  filter(strT_gene_id %notin% as.vector(multi.AQ.count$strT_gene_id)) %>%
  filter(strT_gene_id %notin% as.vector(multi.strT.df$strT_gene_id)) %>%
  distinct()

inherit_name_1.df <- class_equal.df %>% mutate(gene_id = AQ_id)

class_non_equal_match.df <- class_by_loci %>%
  filter(class_match == TRUE) %>%
  filter(class_equal == FALSE) %>%
  filter(AQ_id != "-") %>%
  filter(YD_orf == TRUE & AQ_orf == TRUE) %>%
  filter(strT_gene_id %notin% as.vector(multi.AQ.count$strT_gene_id)) %>%
  filter(strT_gene_id %notin% as.vector(multi.strT.df$strT_gene_id)) %>%
  distinct()

inherit_name_2.df <- class_non_equal_match.df %>% mutate(gene_id = AQ_id)

inherit_name_3.df <- AQ_best_ref4 %>%
  filter(strT_gene_id %notin% as.vector(multi.strT.df$strT_gene_id))  # not one AQ gene shared by two strT genes


inherit_name.merge <- bind_rows(inherit_name_1.df, inherit_name_2.df, inherit_name_3.df) %>%
  select(strT_gene_id, AQ_id, class_code) %>%
  mutate(inherit_ref = TRUE) %>%
  mutate(ref_name = AQ_id) %>%
  mutate(gene_name = AQ_id) %>%
  select(-AQ_id) %>%
  mutate(
    homolog_gene = str_extract(ref_name, "\\|(.*?)(?:-|$)") %>%
      str_remove("^\\|") %>%
      str_remove("-$"),
    homolog_gene = if_else(str_detect(homolog_gene, "^unknown"), toupper(homolog_gene), homolog_gene),
    homolog_gene = na_if(homolog_gene, "NA"),
    homolog_gene = na_if(homolog_gene, "UNKNOWN") # Add this line
  )

# 4.2 NEW NAMES not match ----
class_non_match.df <- class_by_loci %>%
  filter(YD_id %notin% as.vector(inherit_name.merge$strT_gene_id)) %>%
  distinct() %>%
  select(strT_gene_id, AQ_id, class_code) %>%
  mutate(inherit_ref = FALSE) %>%
  mutate(ref_name = AQ_id) %>%
  select(-AQ_id)

# get the homology
strT.ref <- strT_parsed %>%
  transmute(
    gene_id = gene_id,
    transcript_id = transcript_id,
    trD.score = str_extract(Info, 'score=([^,]+)') %>%
      str_remove('score=') %>%
      as.numeric(),
    ORF_type = str_extract(Info, 'ORF type:[^ ]+') %>%
      str_remove('ORF type:') %>%
      str_trim(),
    length = str_extract(Info, 'len:([0-9]+)') %>%
      str_remove('len:') %>%
      as.numeric(),
    homolog = str_extract(Info, 'sp\\|([^\\|]+\\|[^,]+)'),
    accession = if_else(!is.na(homolog),
                        str_extract(homolog, "^sp\\|[^\\|]+"),
                        NA_character_),
    homolog_gene = if_else(!is.na(homolog),
                           str_extract(homolog, "(?<=\\|)[^\\|]+(?=_HUMAN)"),
                           NA_character_),
    perc_identity = if_else(!is.na(homolog),
                            str_extract(homolog, "_HUMAN\\|([0-9.]+)") %>%
                              str_remove("_HUMAN\\|") %>%
                              as.numeric(),
                            NA_real_),
    E_value = if_else(!is.na(homolog),
                      str_extract(Info, "(0|[0-9.]+e[-+][0-9]+)") %>%
                        as.numeric(),
                      NA_real_)
  )

# get AQ's name pool
AQ_id_full <- names(Biostrings::readDNAStringSet("AQ_transcriptome.nt"))

AQ_id_name <- sub("^[^|]+\\|([^ ]+).*", "\\1", names(Biostrings::readDNAStringSet("AQ_transcriptome.nt")))

AQ_id.df <- data.frame(
  AQ_id = AQ_id_name,
  AQ_root = sub("-.*$", "", AQ_id_name),
  stringsAsFactors = FALSE
)

# process strT_ref
strT.ref.process0 <- strT.ref %>%
  select(gene_id, homolog_gene) %>%
  mutate(strT_root = if_else(is.na(homolog_gene), "unknown", tolower(homolog_gene))) %>%
  filter(gene_id %notin% as.vector(inherit_name.merge$strT_gene_id)) %>%
  distinct() 

AQ_id_set <- unique(AQ_id_name)  # Make a set for faster matching

strT.ref.process <- strT.ref.process0 %>%
  mutate(strT_num = sprintf("%08d", 98300000 + row_number())) %>%
  
  group_by(strT_root) %>%
  mutate(
    strT_name = {
      candidate_names <- c(strT_root, paste0(strT_root, "-", 2:1000))
      candidate_names[!candidate_names %in% AQ_id_set][seq_len(n())]
    }
  ) %>%
  ungroup() 

# Find which "unknown", and how many
unknown_idx <- which(strT.ref.process$strT_root == "unknown")

if (length(unknown_idx) > 0) {
  # Find suffixes already used by AQ_id_name for "unknown"
  unknown_used_suffix <- AQ_id_name[startsWith(AQ_id_name, "unknown")]
  
  # Extract used numbers
  used_numbers <- c(1, as.integer(str_remove(unknown_used_suffix, "unknown-")))
  used_numbers <- used_numbers[!is.na(used_numbers)]  # remove NAs if any
  
  # Find available numbers
  available_numbers <- setdiff(1:(length(unknown_idx) + max(used_numbers, 1) + 100), used_numbers)
  
  # Assign new strT_name for unknowns
  strT.ref.process$strT_name[unknown_idx] <- ifelse(
    seq_along(unknown_idx) == 1,
    "unknown",
    paste0("unknown-", available_numbers[2:length(unknown_idx)])
  )
}

new_name.df0 <- strT.ref.process %>%
  mutate(gene_name = paste(strT_num, strT_name, sep = "|")) %>%
  dplyr::rename(strT_gene_id = gene_id) %>%
  select(strT_gene_id, homolog_gene, gene_name) %>%
  mutate(inherit_ref = FALSE) %>% 
  left_join(class_by_loci %>% select(strT_gene_id, AQ_id, class_code), by = "strT_gene_id") %>%
  dplyr::rename(ref_name = AQ_id) %>%
  distinct()

# merge if multi AQ ref to one strT_gene_id
new_name.df <- new_name.df0 %>%
  group_by(strT_gene_id) %>%
  summarise(
    homolog_gene = dplyr::first(homolog_gene),
    gene_name = dplyr::first(gene_name),
    inherit_ref = dplyr::first(inherit_ref),
    ref_name = paste(ref_name[ref_name != "-"], collapse = ";"),
    class_code = paste(class_code[ref_name != "-"], collapse = ";")
  ) %>%
  ungroup()

## Generate GTF
info.all <- rbind(new_name.df, inherit_name.merge) %>%
  dplyr::rename(gene_id = strT_gene_id) %>%
  mutate(gene_name = ifelse(gene_name == "98300002|unknown", "98300002|unknown-1", gene_name))

gtf.df00 <- strT_parsed %>% 
  left_join(info.all, by = "gene_id")


#### Revise published genes
## Get published genes
pub_csv <- read.csv("Published Genes.csv", header = T) %>% select(-Old.ID)

pub_all.df <- pub_csv %>%
  dplyr::rename(gene_name = New.ID) %>%
  left_join(info.all %>% select(gene_id, gene_name, ref_name, class_code, inherit_ref), by = "gene_name")

pub_no_match.df <- pub_all.df %>%
  filter(gene_name %notin% as.vector(info.all$gene_name)) %>%
  select(gene_name, Paper) %>%
  dplyr::rename(published_name = gene_name) %>%
  mutate(Current_gene_name = NA_character_)

# write.csv(pub_no_match.df, "published_id_revise.csv",row.names = F)

revised_ref_name <- read.csv("published_id_revise_mod.csv", header = T)  %>%
  select(Current_gene_name, Revised_gene_name) 

gtf.df0 <- gtf.df00 %>%
  left_join(revised_ref_name, by = c("gene_name" = "Current_gene_name")) %>%
  mutate(gene_name = ifelse(!is.na(Revised_gene_name), Revised_gene_name, gene_name)) %>%
  select(-Revised_gene_name) %>%
  left_join(revised_ref_name, by = c("ref_name" = "Current_gene_name")) %>%
  mutate(ref_name = ifelse(!is.na(Revised_gene_name), Revised_gene_name, ref_name)) %>%
  select(-Revised_gene_name)

gtf.df <- gtf.df0 %>% 
  mutate(
    ref_name = ifelse(is.na(ref_name) | ref_name == "", "-", ref_name),
    new_V9 = paste0(
      'gene_name "', gene_name, '"; ',
      'gene_id "', gene_id, '"; ',
      'transcript_id "', transcript_id, '"; ',
      'gene_biotype "', gene_biotype, '"; ',
      'Info "', Info, '"; ',
      'ref_name "', ref_name, '"; ',
      'class_code "', class_code, '"; ',
      'inherit_ref "', inherit_ref, '";'
    )
  ) %>%
  select(V1, V2, V3, V4, V5, V6, V7, V8, new_V9)


gtf.df1 <- rbind(gtf.df) %>%
  arrange(V1, V4, V5, factor(V3, levels = c("gene", "transcript", "exon", "CDS", "five_prime_UTR", "three_prime_UTR")))

# Save to GTF-like format
write.table(gtf.df1, file = "PBv6_annotation_v0.2.gtf", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(gtf.df1, file = "../DEMO/PBv6_annotation_v0.2.gtf", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

library(rtracklayer)
# Import GTF
gtf <- import("PBv6_annotation_v0.2.gtf", format = "gtf") 

# statistics
gtf.df <- gtf %>% as.data.frame()

gtf.df %>% pull(gene_id) %>% unique() %>% length()   #29836
gtf.df %>% pull(transcript_id) %>% unique() %>% length()   #58581

# Export to GFF3
export(gtf, "PBv6_annotation_v0.2.gff3", format = "gff3")
export(gtf, "../DEMO/PBv6_annotation_v0.2.gff3", format = "gff3")

#######
#####
###
#

Rscript inherit.AQ.R

# Part 4 BUSCO
cd $path

## 4.1.1 StrT all transcript export 
gffread PBv6_annotation_v0.3.gtf -g ../hofPB_v6.FINAL.fa -w PBv6_transcriptome_all_dirty_header.fa -W

##4.1.2 Make txid_to_genename and longest_mRNA_id
nano longest_mRNA_and_txid.R

#!/usr/bin/env Rscript

# ---- Load libraries ----
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(rtracklayer)
  library(Biostrings)
})

# ---- Load genome ----
genome <- Biostrings::readDNAStringSet("../hofPB_v6.FINAL.fa")
names(genome) <- str_replace(names(genome), " .*", "")  # Clean FASTA headers

# ---- Load GTF and convert to data.frame ----
gtf <- rtracklayer::import("../PBv6_annotation_v0.3.gtf") %>% as.data.frame()

# Optional: write parsed GTF for inspection
write.csv(gtf, "PBv6_gtf_parsed_v0.3.csv", row.names = FALSE)

# ---- Extract transcript lengths ----
transcripts_df <- gtf %>%
  filter(type == "transcript") %>%
  select(gene_id, transcript_id, start, end) %>%
  mutate(length = abs(end - start) + 1)

# ---- Get longest transcript per gene ----
longest_tx_df <- transcripts_df %>%
  group_by(gene_id) %>%
  slice_max(order_by = length, n = 1, with_ties = FALSE) %>%
  ungroup()

# Save longest transcript IDs to file
write_lines(longest_tx_df$transcript_id, "longest_transcripts_v0.3.txt")

# ---- Extract transcript-to-gene name mapping ----
map_df <- gtf %>%
  filter(type == "transcript") %>%
  select(transcript_id, gene_name) %>%
  distinct()

# Save transcript ID to gene name mapping (for use in AWK, etc.)
write_tsv(map_df, "txid_to_genename.tsv", col_names = FALSE)

###################
Rscript longest_mRNA_and_txid.R

## 4.1.3 clean header and add gene_name to header
awk '
BEGIN {
  FS = "\t"
  # Read transcript-to-gene mapping
  while ((getline < "txid_to_genename.tsv") > 0) {
    tx2gene[$1] = $2
  }
}
/^>/ {
  # Save previous record if longer
  if (tid && length(seq) > len[tid]) {
    seqs[tid] = seq
    len[tid] = length(seq)
  }

  # Extract transcript ID
  match($0, /^>([^ ]+)/, m)
  tid = m[1]
  seq = ""
  next
}
{
  seq = seq $0
}
END {
  # Save last one
  if (tid && length(seq) > len[tid]) {
    seqs[tid] = seq
  }

  # Print results
  for (id in seqs) {
    gene = (id in tx2gene) ? tx2gene[id] : "NA"
    print ">" id " " gene
    print seqs[id]
  }
}
' PBv6_transcriptome_all_dirty_header.fa > PBv6_transcriptome_all.fa

## 4.1.4 export longest mRNA
seqkit grep -f longest_transcripts_v0.3.txt PBv6_transcriptome_all.fa > PBv6_transcriptome_longest_v0.3.fa

# 4.2 Harvard RC 
nano /n/srivastava_lab/Lab/oscar/busco_transcriptome.sh

#!/usr/bin/bash
#SBATCH -J BUSCO2
#SBATCH -N 1                      # Ensure that all cores are on one machine
#SBATCH -n 1                # Use n cores for one job
#SBATCH -t 1-23:59                # Runtime in D-HH:MM
#SBATCH -p shared              # Partition to submit to
#SBATCH --mem=96000            # Memory pool for all cores
#SBATCH -o out/busco_t.out    # File to which STDOUT will be written
#SBATCH -e out/busco_t.err    # File to which STDERR will be written
#SBATCH --mail-type=ALL           # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=yeduan@fas.harvard.edu # Email to which notifications will be se

cd /n/srivastava_lab/Lab/oscar/busco
mkdir -p busco_t

singularity exec /cvmfs/singularity.galaxyproject.org/b/u/busco:5.8.2--pyhdfd78af_0 \
  busco -i ../PBv6_transcriptome_all.fa \
        -l eukaryota_odb12 \
        -o busco_t/busco_output_strT_eukaryota_12_t_all \
        -m transcriptome \
        -f \
        --offline


singularity exec /cvmfs/singularity.galaxyproject.org/b/u/busco:5.8.2--pyhdfd78af_0 \
  busco -i ../PBv6_transcriptome_longest.fa  \
        -l eukaryota_odb12 \
        -o busco_t/busco_output_strT_eukaryota_12_t_longest \
        -m transcriptome \
        -f \
        --offline

singularity exec /cvmfs/singularity.galaxyproject.org/b/u/busco:5.8.2--pyhdfd78af_0 \
  busco -i ../YL_transcriptome.fa \
        -l eukaryota_odb12 \
        -o busco_t/busco_output_YL_eukaryota_12_t_all \
        -m transcriptome \
        -f \
        --offline


singularity exec /cvmfs/singularity.galaxyproject.org/b/u/busco:5.8.2--pyhdfd78af_0 \
  busco -i AQ_transcriptome_clean.fa \
        -l eukaryota_odb12 \
        -o busco_t/busco_output_AQ_eukaryota_12_t_all \
        -m transcriptome \
        -f \
        --offline


# Transciptome BUSCO Metazoa_odb12
singularity exec /cvmfs/singularity.galaxyproject.org/b/u/busco:5.8.2--pyhdfd78af_0 \
  busco -i ../PBv6_transcriptome_all.fa \
        -l metazoa_odb12 \
        -o busco_t/busco_output_strT_metazoa_12_t_all \
        -m transcriptome \
        -f \
        --offline


singularity exec /cvmfs/singularity.galaxyproject.org/b/u/busco:5.8.2--pyhdfd78af_0 \
  busco -i ../PBv6_transcriptome_longest.fa  \
        -l metazoa_odb12 \
        -o busco_t/busco_output_strT_metazoa_12_t_longest \
        -m transcriptome \
        -f \
        --offline

singularity exec /cvmfs/singularity.galaxyproject.org/b/u/busco:5.8.2--pyhdfd78af_0 \
  busco -i ../YL_transcriptome.fa \
        -l metazoa_odb12 \
        -o busco_t/busco_output_YL_metazoa_12_t_all \
        -m transcriptome \
        -f \
        --offline


singularity exec /cvmfs/singularity.galaxyproject.org/b/u/busco:5.8.2--pyhdfd78af_0 \
  busco -i AQ_transcriptome_clean.fa \
        -l metazoa_odb12 \
        -o busco_t/busco_output_AQ_metazoa_12_t_all \
        -m transcriptome \
        -f \
        --offline
