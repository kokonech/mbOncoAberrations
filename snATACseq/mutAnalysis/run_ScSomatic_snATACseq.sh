#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l mem=48g
#PBS -l walltime=72:00:00
#PBS -M k.okonechnikov@dkfz-heidelberg.de
#PBS -j oe
#PBS -o /b06x-isilon/b06x-m/mbCSF/results/humanTumor/mbSpatial/logs/scomatic


module load Python/3.8.2-foss-2020a
module load R/3.6.3-foss-2017a
module load SAMtools/1.4.1
module load BEDTools


SCOMATIC=/home/ad.dkfz-heidelberg.de/okonechn/tools/SComatic
SCO_DATA=/omics/odcf/analysis/OE0290_projects/Ependymoma/annotations/SComatic_data


echo "INPUT: $INPUT"
echo "BAM: $BAM"
echo "ANN: $ANN"

sample=$INPUT

RESDIR=/b06x-isilon/b06x-m/mbCSF/results/humanTumor/mbSpatial/atacMut


output_dir=$RESDIR/$sample
mkdir -p $output_dir


output_dir1=$output_dir/Step1_BamCellTypes
mkdir -p $output_dir1

python $SCOMATIC/scripts/SplitBam/SplitBamCellTypes.py --bam $BAM \
        --meta $ANN\
        --id ${sample} \
        --n_trim 5 \
        --min_MQ 30 \
     --outdir $output_dir1


REF=/omics/odcf/analysis/OE0290_projects/Ependymoma/annotations/hg_GRCh38/GRCh38_v32.primary_assembly.genome.fa

output_dir2=$output_dir/Step2_BaseCellCounts
mkdir -p $output_dir2

for bam in $(ls -d $output_dir1/*bam);do
  
  # Cell type
  cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')

  res2=$output_dir2/${sample}.${cell_type}.tsv

  if [ -f $res2 ]; then
    echo "Result for $res2 exists, skipping..."
    continue
  fi
       

  # Temp folder
  temp=$output_dir2/temp_${cell_type}
  mkdir -p $temp

  python $SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py --bam $bam \
    --ref $REF \
    --chrom all \
    --out_folder $output_dir2 \
    --min_bq 30 \
    --min_mq 30 \
    --tmp_dir $temp \
    --nprocs 8

  rm -rf $temp
done


output_dir3=$output_dir/Step3_BaseCellCountsMerged
mkdir -p $output_dir3

python $SCOMATIC/scripts/MergeCounts/MergeBaseCellCounts.py --tsv_folder ${output_dir2} \
  --outfile ${output_dir3}/${sample}.BaseCellCounts.AllCellTypes.tsv



output_dir4=$output_dir/Step4_VariantCalling
mkdir -p $output_dir4

# step1
python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step1.py \
          --infile ${output_dir3}/${sample}.BaseCellCounts.AllCellTypes.tsv \
          --outfile ${output_dir4}/${sample} \
          --ref $REF

# step 2

# this could be custom
PON=$SCO_DATA/PoNs/PoN.scATACseq.hg38.tsv

python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step2.py \
          --infile ${output_dir4}/${sample}.calling.step1.tsv \
          --outfile ${output_dir4}/${sample} \
          --pon $PON




