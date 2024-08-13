#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l mem=32g
#PBS -l walltime=10:00:00
#PBS -M k.okonechnikov@dkfz-heidelberg.de
#PBS -j oe
#PBS -o /omics/odcf/analysis/OE0290_projects/bioskryp/ResolveOme/Medulloblastoma/rna/Logs/STAR

#TOOL=/home/okonechn/tools/STAR-STAR_2.4.1d/bin/Linux_x86_64/STAR

echo "##########" 
echo "qsub_runSTAR.sh" 
echo `date`
echo `hostname`

echo "Tumor case: $TUMOR"
echo "Cell: $INPUT"
echo "Input reads: $READS"

module load STAR/2.7.6a-foss-2017a
TOOL=STAR

GENOME="/omics/odcf/reference_data/by-species/Homo_sapiens/GRCh38/GRCh38_decoy_ebv_phiX/STAR/2.7.10a/gencode_v39_chr_patch_hapl_scaff"

RESDIR="/omics/odcf/analysis/OE0290_projects/bioskryp/ResolveOme/Medulloblastoma/rna/Result/$TUMOR/STAR"


echo "Entering $RESDIR"
cd $RESDIR

SID=$INPUT

cmd="$TOOL --genomeDir $GENOME --readFilesCommand zcat --runThreadN 8  --outSAMtype BAM SortedByCoordinate --clip5pNbases 26 --quantMode GeneCounts   --outFileNamePrefix "${SID}_" --readFilesIn $READS"
echo $cmd
$cmd

echo "Finished!"

