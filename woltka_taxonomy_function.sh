#!/bin/bash -l
#PBS -m ae
#PBS -M swandro@ucsd.edu
#PBS -S /bin/bash
#PBS -e ~/logs
#PBS -o ~/logs
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=2
#PBS -l mem=64gb
#PBS -N woltka_ext_metag


###Edit the input and output folder and conda environment below######

source ~/.bash_profile
#Activate conda environment containing woltka
conda activate woltka

#Set root directory for all output folders
out_root=~/ANALYSIS_DIR

#Set input directories
alignment_out=$out_root/alignments
tax_db=/projects/wol/20170307/release/taxonomy
annotation_db=/projects/wol/20170307/release/annotation

########Should not need to edit below this point###################

############
#Run Woltka#
############

#Create temp directory
export TMPDIR=/panfs/panfs1.ucsd.edu/panscratch/$USER
tmp=$(mktemp -d --tmpdir)
export TMPDIR=$tmp
trap "rm -r $tmp; unset TMPDIR" EXIT
cd $tmp

#Count cpus
cpus=$PBS_NUM_PPN

#Create output directories
map_dir=$out_root/maps
profile_dir=$out_root/profile
taxfunc_dir=$out_root/taxfunc
woltka_time_log=$out_root/woltka_time_log.txt
mkdir -p $map_dir $function_dir $profile_dir $taxfunc_dir

#####################
#Create gOTU profile#
#####################
#This output can be used alongside woltka phylogeny for phylogenetic analyses.
woltka gotu -i $alignment_out -o $profile_dir/gotu_profile.biom

################
#Taxa profiling#
################
#Create taxonomy at phylum, genus, species level
#Also create map for taxa stratified functional profiling

start=`date +%s`
woltka classify \
  -i $alignment_out \
  --add-lineage \
  --map $tax_db/taxmap.txt \
  --nodes $tax_db/nodes.dmp \
  --names $tax_db/names.dmp \
  --rank phylum,genus,species \
  --name-as-id \
  --outmap $map_dir \
  -o $profile_dir/
end=`date +%s`
runtime=$((end-start))
echo "Woltka finished making taxa profiles in $runtime seconds." >> $woltka_time_log

##############################
####Functional profiling######
##############################
#Function stratified by genus

start=`date +%s`
woltka classify \
  -i $alignment_out \
  --coords $annotation_db/coords.txt.xz \
  --map $annotation_db/uniref.map.xz \
  --map $annotation_db/metacyc/protein.map \
  --map $annotation_db/metacyc/protein2enzrxn.map \
  --map $annotation_db/metacyc/enzrxn2reaction.map \
  --map $annotation_db/metacyc/reaction2pathway.map \
  --map $annotation_db/metacyc/pathway2class.map \
  --map-as-rank \
  --rank uniref,protein,protein2enzrxn,enzrxn2reaction,reaction2pathway,pathway2class \
  --stratify $map_dir/genus \
  -o $taxfunc_dir/
end=`date +%s`
runtime=$((end-start))
echo "Woltka finished functional profiling in $runtime seconds." >> $woltka_time_log
