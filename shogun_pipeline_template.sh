#!/bin/bash -l

#PBS -m ae
#PBS -M swandro@ucsd.edu
#PBS -S /bin/bash
#PBS -e ~/logs
#PBS -o ~/logs
#PBS -l walltime=8:00:00
#PBS -l nodes=1:ppn=16
#PBS -l mem=60gb
#PBS -N shogun_pipeline

#Change to whatever environment shogun is installed in
##Software required in shogun environment: shogun, pigz
conda activate shogun

#In directory of fastq files
in_dir=~/FASTQ_DIR

#Set output directory
out_root=~/ANALYSIS_DIR

#Load modules (Should already exist on barnacle)
module load bedtools_2.26.0
module load samtools_1.3.1

#Make sure these files/folders exist.
index_file=/projects/cmi_proj/host_depletion/three_studies/TruSeq3-PE-2_G.fa
filter_db=/databases/bowtie/Human_phiX174/Human_phix174
flash_path=/projects/cmi_proj/host_depletion/three_studies/flash/FLASH-1.2.11-Linux-x86_64/flash
db=/projects/wol/20170307/release/databases/shogun
#Alternatively, use rep82
#db=/databases/genome/rep82/shogun

#LAST THING TO CHECK IS THE INPUT FASTQ FILE SUFFIX BELOW.
#Everything else below should be fine.

# Calculate the number of processors allocated to this run.
cpus=$PBS_NUM_PPN

#Make remp directory in panscratch to work in
export TMPDIR=/panfs/panfs1.ucsd.edu/panscratch/$USER
tmp=$(mktemp -d --tmpdir)
export TMPDIR=$tmp
trap "rm -r $tmp; unset TMPDIR" EXIT

#Move to temp directory
cd $tmp

#Create output directories
out_dir_qual=$out_root/qf_fastq
log_dir_qual=$out_dir_qual/log
out_dir_filter=$out_root/human_filtered_fastq
log_dir_filter=$out_dir_filter/log
out_dir_merge=$out_root/merged_fastq
log_dir_merge=$out_dir_merge/log
out_dir_shogun_alignments=$out_root/shogun_wol_alignments
out_dir_shogun_profiles=$out_root/shogun_wol_profiles

mkdir -p $out_dir_qual $log_dir_qual $out_dir_filter $log_dir_filter $out_dir_merge $log_dir_merge $out_dir_shogun_alignments $out_dir_shogun_profiles

for f_read in $in_dir/*R1_paired.fastq.gz
  do
    ######Get filenames######
    no_root=${f_read##*/}
    ###########################Edit this base name if the suffix dos not match##############################################
    base_name=${no_root%%_L001_001.trimmed_R1_paired.fastq.gz}
    r_read=${f_read/R1/R2}

    ####Check if file has been done####
    if [ ! -f $out_dir_shogun_alignments/"$base_name"_bowtie2_wol_alignment.sam ]
      then
      ##########################
      ## RUN QUALITY TRIMMING ##
      ##########################
      echo "Running trimmomatic on $base_name"
      trimmomatic PE \
      -phred33 \
      $f_read $r_read \
      $tmp/"$base_name"_R1_paired.fastq \
      $tmp/"$base_name"_R1_single.fastq \
      $tmp/"$base_name"_R2_paired.fastq \
      $tmp/"$base_name"_R2_single.fastq \
      ILLUMINACLIP:"$index_file":2:30:7 \
      MINLEN:100 \
      TRAILING:20 \
      AVGQUAL:20 \
      SLIDINGWINDOW:20:20 \
      2> $log_dir_qual/"$base_name"_trimlog.txt
      echo "Finished trimmomatic on $base_name"

      ##########################
      ## RUN HUMAN FILTERING ###
      ##########################
      echo "Running human filter on $base_name"
      #Filter paired human reads
      bowtie2 -p $cpus -x $filter_db -1 $tmp/"$base_name"_R1_paired.fastq -2 $tmp/"$base_name"_R2_paired.fastq --fast-local 2> $log_dir_filter/"$base_name"_human_filter_log.txt |\
      samtools view -f 12 -F 256 |\
      samtools sort -@ $cpus -n |\
      samtools view -bS |\
      bedtools bamtofastq -i - -fq $tmp/"$base_name"_R1_paired_filt.fastq -fq2 $tmp/"$base_name"_R2_paired_filt.fastq
      #Filter unpaired human reads
      bowtie2 -p $cpus -x $filter_db -U $tmp/"$base_name"_R1_single.fastq,$tmp/"$base_name"_R2_single.fastq --fast-local |\
      samtools view -f 4 -F 256 |\
      samtools sort -@ $cpus -n |\
      samtools view -bS |\
      bedtools bamtofastq -i - -fq $tmp/"$base_name"_R1_single_filt.fastq -fq2 $tmp/"$base_name"_R2_single_filt.fastq
      echo "Finished running human filter on $base_name"

      ##########################
      #### Run Flash Merge #####
      ##########################
      echo "Running Flash on $base_name"
      $flash_path \
      -M 150 \
      -m 20 \
      -x 0.01 \
      -o $base_name \
      -d $tmp \
      $tmp/"$base_name"_R1_paired_filt.fastq $tmp/"$base_name"_R2_paired_filt.fastq \
      2>&1 | tee $log_dir_merge/"$base_name"_mergelog.txt
      rm $tmp/"$base_name".hist
      rm $tmp/"$base_name".histogram
      echo "Finished Flash on $base_name"

      ##########################
      ###### Run SHOGUN ########
      ##########################
      echo "Running Shogun on $base_name"
      #Combine merged and unmerged seqs and the singletons
      cat $tmp/"$base_name".extendedFrags.fastq \
      $tmp/"$base_name".notCombined_1.fastq \
      $tmp/"$base_name".notCombined_2.fastq \
      $tmp/"$base_name"_R1_single_filt.fastq \
      $tmp/"$base_name"_R2_single_filt.fastq | seqtk seq -A > $tmp/"$base_name".fasta
      #Run Shogun
      shogun align -t $cpus -d $db -a bowtie2 -i $tmp/"$base_name".fasta -o $tmp/
      #Rename alignment
      mv $tmp/alignment.bowtie2.sam $tmp/"$base_name"_bowtie2_wol_alignment.sam
      #Assign taxonomy
      shogun assign_taxonomy -d $db -a bowtie2 \
      -i $tmp/"$base_name"_bowtie2_wol_alignment.sam \
      -o $tmp/"$base_name"_wol_profile.tsv
      #Cleanup
      rm $tmp/"$base_name".fasta
      echo "Finished Shogun on $base_name"

      ###################################
      ####### Zip and Move Files#########
      ###################################
      echo "Zipping and moving files for $base_name"
      pigz $tmp/*.fastq
      #Qulity
      mv $tmp/"$base_name"_R1_paired.fastq.gz $out_dir_qual
      mv $tmp/"$base_name"_R1_single.fastq.gz $out_dir_qual
      mv $tmp/"$base_name"_R2_paired.fastq.gz $out_dir_qual
      mv $tmp/"$base_name"_R2_single.fastq.gz $out_dir_qual
      #Filtered
      mv $tmp/"$base_name"_R1_paired_filt.fastq.gz $out_dir_filter
      mv $tmp/"$base_name"_R2_paired_filt.fastq.gz $out_dir_filter
      mv $tmp/"$base_name"_R1_single_filt.fastq.gz $out_dir_filter
      mv $tmp/"$base_name"_R2_single_filt.fastq.gz $out_dir_filter
      #Merged
      mv $tmp/"$base_name".extendedFrags.fastq.gz $out_dir_merge
      mv $tmp/"$base_name".notCombined_1.fastq.gz $out_dir_merge
      mv $tmp/"$base_name".notCombined_2.fastq.gz $out_dir_merge
      #Shogun
      mv $tmp/"$base_name"_bowtie2_wol_alignment.sam $out_dir_shogun_alignments
      mv $tmp/"$base_name"_wol_profile.tsv $out_dir_shogun_profiles
      echo "Finished zipping and moving files for $base_name"
    fi
  done
