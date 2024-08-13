#!/bin/bash

#PBS -M k.okonechnikov@dkfz-heidelberg.de
#PBS -l mem=32g
#PBS -l nodes=1:ppn=8
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -o /omics/odcf/analysis/OE0290_projects/bioskryp/ResolveOme/Medulloblastoma/dna/Logs/dnaQC

# This code runs BaseJumper BJ-DNA-QC pipeline per cell
# https://docs.basejumper.bioskryb.com/pipelines/secondary/bj-dna-qc/1.9.1/docs/


echo "Tumor case: $TUMOR"
READSPATH="${READS}_R{1,2}.fastq.gz"
echo "Input reads: $READSPATH"

OUTDIR="/omics/odcf/analysis/OE0290_projects/bioskryp/ResolveOme/Medulloblastoma/dna/Result/BJ-DNA-QC/${TUMOR}_100Kbp"

mkdir -p $OUTDIR

OUTPUT=$OUTDIR/$INPUT
echo "Sample cell output folder: $OUTPUT"
#mkdir -p $OUTPUT

INST_DIR="/b06x-isi/share/csw/pipelines/bioskryb/basejumper/bj-dna-qc"

export NXF_VER=22.10.7
export NXF_SINGULARITY_CACHEDIR=/b06x-isi/share/csw/pipelines/bioskryb/basejumper/singularity_images
export NXF_OFFLINE='TRUE'

export SENTIEON_LICENSE=b06x-pbs01.inet.dkfz-heidelberg.de:4000

WORK="/omics/odcf/analysis/OE0290_projects/bioskryp/ResolveOme/Medulloblastoma/dna/Result/BJ-DNA-QC/Temp"
GENOME_BASE="/b06x-isi/share/csw/pipelines/bioskryb/basejumper/data"


module load Nextflow/22.10.7

cd $WORK && nextflow run ${INST_DIR}/main.nf -profile singularity,test --max_cpus 8 --max_memory 32.GB --genomes_base ${GENOME_BASE} \
--bin_size 100000 --publish_dir ${OUTPUT} \
--reads "${READSPATH}"

find ${OUTPUT}* -type d -print0 | xargs -0 chmod g+rwx
find ${OUTPUT}* -type f -print0 | xargs -0 chmod g+rw


echo "Done!"

