#Stephen Wandro
#Updated 5/19/2020

#Load modules
module load samtools_1.3.1
module load bedtools_2.29.2
module load bowtie_2.2.3

#Options
human_database=/databases/bowtie/Human_phiX174/Human_phix174
cpus=1

##################
###Paired Reads###
##################
#Input
input_fastq_paired_1=path/to/file_R1.fastq
input_fastq_paired_2=path/to/file_R2.fastq
#Output
filtered_fastq_paired_1=path/to/filtered_file_R1.fastq
filtered_fastq_paired_2=path/to/filtered_file_R2.fastq
log_file=path/to/bt2_log.txt

#Run command to filter paired reads
bowtie2 -p $cpus -x $human_database \
  -1 $input_fastq_paired_1 \
  -2 $input_fastq_paired_2 \
  --fast-local 2> $log_file |\
  samtools view -f 12 -F 256 |\
  samtools sort -@ $cpus -n |\
  samtools view -bS |\
  bedtools bamtofastq -i - \
  -fq $filtered_fastq_paired_1 \
  -fq2 $filtered_fastq_paired_2


####################
###Unpaired Reads###
####################
#Input
input_fastq_unpaired=path/to/file.fastq
#Output
filtered_fastq_unpaired=path/to/filtered_file.fastq
log_file=path/to/bt2_log.txt

#Filter unpaired human reads
bowtie2 -p $cpus -x $human_database \
-U $input_fastq_unpaired \
--fast-local |\
samtools view -f 4 -F 256 |\
samtools sort -@ $cpus -n |\
samtools view -bS |\
bedtools bamtofastq -i - \
-fq $filtered_fastq_unpaired
