# 1.1 Genome BUSCO Eukaryota_odb12
# Harvard RC 
nano /n/srivastava_lab/Lab/oscar/busco_eukaryota12.sh

#!/usr/bin/bash
#SBATCH -J BUSCO1
#SBATCH -N 1                      # Ensure that all cores are on one machine
#SBATCH -n 1                # Use n cores for one job
#SBATCH -t 1-23:59                # Runtime in D-HH:MM
#SBATCH -p shared              # Partition to submit to
#SBATCH --mem=96000            # Memory pool for all cores
#SBATCH -o out/busco_euk1.out    # File to which STDOUT will be written
#SBATCH -e out/busco_euk1.err    # File to which STDERR will be written
#SBATCH --mail-type=ALL           # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=yeduan@fas.harvard.edu # Email to which notifications will be se

cd /n/srivastava_lab/Lab/oscar/busco

singularity exec /cvmfs/singularity.galaxyproject.org/b/u/busco:5.8.2--pyhdfd78af_0 \
  busco -i ../hofPB_v6.FINAL.fa \
        -l eukaryota_odb12 \
        -o busco_output_v6_eukaryota_12 \
        -m genome \
        -f \
        --offline


singularity exec /cvmfs/singularity.galaxyproject.org/b/u/busco:5.8.2--pyhdfd78af_0 \
  busco -i ../Hmia.030120.fasta  \
        -l eukaryota_odb12 \
        -o busco_output_NP_eukaryota_12 \
        -m genome \
        -f \
        --offline


# 1.2 Genome BUSCO Metazoa_odb12
nano /n/srivastava_lab/Lab/oscar/busco_metazoa12.sh

#!/usr/bin/bash
#SBATCH -J BUSCO2
#SBATCH -N 1                      # Ensure that all cores are on one machine
#SBATCH -n 1                # Use n cores for one job
#SBATCH -t 1-11:59                # Runtime in D-HH:MM
#SBATCH -p shared              # Partition to submit to
#SBATCH --mem=96000            # Memory pool for all cores
#SBATCH -o out/busco_m1.out    # File to which STDOUT will be written
#SBATCH -e out/busco_m1.err    # File to which STDERR will be written
#SBATCH --mail-type=ALL           # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=yeduan@fas.harvard.edu # Email to which notifications will be se

cd /n/srivastava_lab/Lab/oscar/busco

singularity exec /cvmfs/singularity.galaxyproject.org/b/u/busco:5.8.2--pyhdfd78af_0 \
  busco -i ../hofPB_v6.FINAL.fa \
        -l metazoa_odb12 \
        -o busco_output_v6_metazoa_12 \
        -m genome \
        -f \
        --offline


singularity exec /cvmfs/singularity.galaxyproject.org/b/u/busco:5.8.2--pyhdfd78af_0 \
  busco -i ../Hmia.030120.fasta  \
        -l metazoa_odb12 \
        -o busco_output_NP_metazoa_12 \
        -m genome \
        -f \
        --offline


# 2.0 Transcriptome BUSCO export

#Local

# 2.0.1 YL & AQ
gffread final_gffread.gff -g Hmia.030120.fasta -w YL_transcriptome.fa -W
awk '/^>/ {print $1; next} {print}' YL_transcriptome.fa > YL_cleaned_transcriptome.fa

sed 's|/|_|g; s|,|_|g; s|"||g; s| |=|g' AQ_transcriptome.nt > AQ_transcriptome_clean.fa


# 2.0.2 strT all transcripts
gffread PBv6_annotation_v0.2.gtf -g ../hofPB_v6.FINAL.fa -w PBv6_transcriptome_all_dirty_header.fa -W

# clean header and add gene_name to header
awk '
BEGIN {
  FS = "\t"
  while ((getline < "txid_to_genename.tsv") > 0) {
    tx2gene[$1] = $2
  }
}
/^>/ {
  split($0, a, " ")
  gsub(/^>/, "", a[1])
  tid = a[1]
  gene = (tid in tx2gene) ? tx2gene[tid] : "NA"
  print ">" tid " " gene
  next
}
{ print }
' PBv6_transcriptome_all_dirty_header.fa > PBv6_transcriptome_all.fa


## 2.0.3 strT longest transcripts
cd /home/duan/pacbio/new_annotation
Rscript v6_longest_transcripts.R # see the script in the file

seqkit grep -f longest_transcripts_v0.2.txt PBv6_transcriptome_all.fa > PBv6_transcriptome_longest.fa


# 2.1 Transcriptome BUSCO Eukaryota_odb12
# Harvard RC 
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


# 2.2 Transciptome BUSCO Metazoa_odb12
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
