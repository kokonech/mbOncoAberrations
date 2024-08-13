#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l mem=10g
#PBS -l walltime=5:00:00
#PBS -M k.okonechnikov@dkfz-heidelberg.de
#PBS -j oe
#PBS -o /omics/odcf/analysis/OE0290_projects/bioskryp/ResolveOme/Medulloblastoma/rna/Logs/fCounts

TOOL=$HOME/tools/subread-1.6.4-Linux-x86_64/bin/featureCounts

GENCODE=/omics/odcf/analysis/OE0290_projects/Ependymoma/annotations/star_index_hg38/gencode.v32.annotation.gtf

MAINDIR=/omics/odcf/analysis/OE0290_projects/bioskryp/ResolveOme/Medulloblastoma/rna/Result


echo "##########"
echo "Run featureCounts"
echo `date`
echo `hostname`
echo "INPUT: ${INPUT}"

SID=$INPUT
PID=$TUMOR

DATADIR=$MAINDIR/$PID/STAR
BAMFILE=$DATADIR/${SID}_Aligned.sortedByCoord.out.bam

RESFILE=$MAINDIR/$PID/counts/${SID}.counts



echo "Processing $SID"
 
cd $RESDIR   

# option -t gene allows to use full gene annotation with introns


# default settings
cmd="$TOOL -p -T 8 -t gene --tmpDir $DATADIR -a $GENCODE -o $RESFILE $BAMFILE"

# option -M  mulit-mapped
#cmd="$TOOL -p -T 8 -t gene -M --tmpDir $DATADIR -a $GENCODE -o $RESFILE $BAMFILE"

echo $cmd
$cmd

echo "Finished!"


