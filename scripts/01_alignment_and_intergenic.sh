################################################################################################
#Here we align to the refernce. This is done with cellranger which is a module in on our cluster. The version used here is 7.1.0. 
cd ./results/
#####Alignement 
samples=(S22_105-2-Striatum_STR S23_105-2-Thalamus_THA S24_105-3-Striatum_STR S25_105-3-Thalamus_THA S26_105-10-Striatum_STR S27_105-10-Thalamus_THA)
#Loop to submit each
for s in "${samples[@]}"; do
  raw_in="./path" #input by sample
  outdir="./Alignement/Align/${s}" #Output directory one per sample
  export TMPDIR=./Alignement/Align/${s}_TMP
  mkdir -p "$outdir" logs #Make the output direcotyr

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${s}
#SBATCH --output=${s}.out
#SBATCH --error=${s}.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --account=general
#SBATCH --partition=general

module load cellranger/v7.1.0
cd $outdir
cellranger count --id=${s} \
  --fastqs=${raw_in} \
  --sample=105-2-Striatum \
  --transcriptome="./Alignement/GENOME/ensembl_sus_scrofa_11.1/" \
  --localcores=32 --localmem=200
EOF
done
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+


################################################################################################
#Here we clean the reads. This is done with cellbender 
# sample list and chromosome list
mkdir -p Align
cd ./Alignement/Align
samples=(S22_105-2-Striatum_STR S23_105-2-Thalamus_THA S24_105-3-Striatum_STR S25_105-3-Thalamus_THA S26_105-10-Striatum_STR S27_105-10-Thalamus_THA)
#Loop to submit each
for s in "${samples[@]}"; do
  raw_in="./path" #input by sample
  outdir="./Alignement/Align/${s}" #Output directory one per sample
  export TMPDIR=./Alignement/Align/${s}_TMP
  mkdir -p "$outdir" logs #Make the output direcotyr

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${s}
#SBATCH --output=${s}.out
#SBATCH --error=${s}.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --account=general
#SBATCH --partition=general

module load miniconda3/4.10.3_gcc_12.3.0
export CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ./CellBender
cd $outdir
cellbender remove-background --input $outdir/outs/raw_feature_bc_matrix.h5 --output $outdir/outs/${s}.h5 --learning-rate 2.5e-5 --fpr 0.01
ptrepack --complevel 5 $outdir/outs/${s}_filtered.h5:/matrix $outdir/outs/${s}_seurat.h5:/matrix
perl "./scripts/1_mt_gene_extract.pl" --GTF ./Alignement/GENOME/Sus_scrofa.112.GFP.gtf --prefix ${s}
EOF
done


#Here we get counts of reads mapping to intergenic regions. This is DNA spilliage and for nuclei is an additional meausre of quality.
#While stats as a whole are provided for these values their per cell values are not. By measuring them directly with the steps below we can include them in QC decisions.
#####Intergenic Alignement
# sample list and chromosome list
samples=(S22_105-2-Striatum_STR S23_105-2-Thalamus_THA S24_105-3-Striatum_STR S25_105-3-Thalamus_THA S26_105-10-Striatum_STR S27_105-10-Thalamus_THA)
#Chromosome list for Pig including the transgene "chromosome"
chroms=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 X Y MT GFP_Transgene Selectable_Marker)
#Loop to submit each BAM to samtools to break it out by chromosome as the interginc counting is to computationally slow to do per genome.
for s in "${samples[@]}"; do
  bam_in="./Alignement/Align/${s}/outs/possorted_genome_bam.bam" #input bam by sample
  outdir="./Alignement/Align/${s}/outs/" #Output directory one per sample
  mkdir -p "$outdir" logs #Make the output direcotyr

  for chr in "${chroms[@]}"; do #loop through each chromosome and resquest the below resourses. load the sam tools module and extract the chromosome with samtools view
    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=split_${s}_${chr}
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/split_${s}_${chr}_%j.out
#SBATCH --error=logs/split_${s}_${chr}_%j.err
#SBATCH -p general
#SBATCH --account=general

module load samtools #module available on our cluster
samtools view -@ 8 -b ${bam_in} ${chr} > ${outdir}/${chr}.bam
EOF
  done
done
#GTF generated with the alignement reference build
GTF="./Alignement/GENOME/ensembl_sus_scrofa_11.1/genes/genes.gtf.gz"

#Loop to run per_barcode_read_distribution.py per sampel per chromosome. 
for s in "${samples[@]}"; do
  outdir="./Alignement/Align/${s}/outs/Contamintion/"
  mkdir -p "$outdir" logs2

  for chr in "${chroms[@]}"; do
    bam="${outdir}/${chr}.bam"
    out="${outdir}/${chr}.tsv"

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=perbc_${s}_${chr}
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=logs2/perbc_${s}_${chr}_%j.out
#SBATCH --error=logs2/perbc_${s}_${chr}_%j.err
#SBATCH -p general
#SBATCH --account=general


module load miniconda3/4.10.3_gcc_12.3.0
export CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ./barcode_reads
python "../scripts/support/per_barcode_read_distribution.py" --bam ${bam} --gtf ${GTF} --out ${out}
EOF
  done
done

#Loop to get final counts.
conda activate ../environments/barcode_reads
for sample in "${samples[@]}"; do
  Rscript ../scripts/support/collapse_per_barcode.R ./Alignement/Align/${sample}/outs/Contamintion/${sample} ./Alignement/Align/${sample}/outs/Contamintion/${sample}_collapsed.tsv
done
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+
