#This serves as the masterworkflow to run each step of the analyses performed
#We build the enviroenemnts first and then directory structure to run each step
#This workflow was generated to recreate results and is not generalized

#if the project directory pig_gfp_comparative is copied it will inlcude a scripts and environement directory
#Set up the rest of the directories from pig_gfp_comparative
mkdir -p results
mkdir -p ./results/Alignement
mkdir -p ./results/Alignement/GENOME
mkdir -p ./results/Alignement/Align
mkdir -p ./results/Processing
mkdir -p ./results/Processing/AnnotationRef
mkdir -p ./results/Processing
mkdir -p ./results/GFP
mkdir -p ./results/comparison


##Environement Built##
#This workflow requires multiple conda environments that were built in a specific way. 
#most can be generated using the yml file that copied the environments builds
#Also included is a recipe shell file that can be run interactivly to allow building from scratch
cd ./environment
conda env create -f Barcode_processing.yml --prefix ./Barcode_processing
conda env create -f CellBender.yml --prefix ./CellBender
conda env create -f enrichment.yml --prefix ./enrichment
conda env create -f Single_Cell.yml --prefix ./Single_Cell #most likely to need to be built with the recipe
conda env create -f Stats.yml --prefix ./Stats



###############################################################################################
#To do anything we have to align to a reference. This chunk is the steps used for building the reference.
#####Alignement Reference Build
#from project directory working directory
sh ./scripts/00_build_reference_genome.sh
#This changes to ./results/Alignement/GENOME, downloads the appropriate reference and GTF
#Then submits using a slurm sbatch
#NOTE: the submission variables should be updated to fit to the HPC the script is running in
#NOTE: also requires cellranger which was available in our cluster as a module so no environement provided
#NOTE: the GFP sequence must be in the directory availabe from data availability statement
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+



###############################################################################################
#Next perform alignement and non-seurat QC
#from project directory working directory
sh ./scripts/01_alignment_and_intergenic.sh
#Then submits using a slurm sbatch
#NOTE: the submission variables should be updated to fit to the HPC the script is running in
#NOTE: also requires cellranger which was available in our cluster as a module so no environement provided
#NOTE: the SNRNA sequence must be in the directory availabe from data availability statement
#NOTE: This likely works best done interactivly
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+


###############################################################################################
#Next perform processing and QC
#from project directory working directory
conda activate ./environment/Seurat
Rscript ./scripts/02_qc_filtering.R
conda deactivate
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+


###############################################################################################
#Next perform integration and clustering
#from project directory working directory
conda activate ./environment/Seurat
Rscript ./scripts/03_integration_clustering.R
conda deactivate
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+



###############################################################################################
#Next build annotation reference and do annotation.
#from project directory working directory
cd ./results/Processing/AnnotationRef

wget https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/expression_matrices/Zeng-Aging-Mouse-10Xv3/20241130/Zeng-Aging-Mouse-10Xv3-raw.h5ad
wget https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/metadata/Zeng-Aging-Mouse-10Xv3/20241130/cell_cluster_annotations.csv

wget https://datasets.cellxgene.cziscience.com/c2f66cd5-4ff4-4578-876c-55783a57cf8f.h5ad #Nerons
wget https://datasets.cellxgene.cziscience.com/99f27be8-9fac-451e-9723-9e4c7191589e.h5ad #non neurons

conda activate ../../../environment/Seurat
Rscript ../../../scripts/04_annotation_strategy.R
conda deactivate
cd ../../../
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+


###############################################################################################
#Next gfp quantification
#from project directory working directory

conda activate ./environment/Seurat
Rscript ./scripts/05_gfp_quantification.R
conda deactivate

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+


###############################################################################################
#Next gfp enrichment
#from project directory working directory

conda activate ./environment/enrichment
Rscript ./scripts/06_deg_enrichment_analysis.R
conda deactivate

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+


###############################################################################################
#Next ortholog expression comparison
#from project directory working directory

conda activate ./environment/Stats
Rscript ./scripts/07_comparative_orthologs.R
conda deactivate

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+


###############################################################################################
#Next pc gene programs
#from project directory working directory

conda activate ./environment/Stats
Rscript ./scripts/08_pc_gene_programs.R
conda deactivate

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+


###############################################################################################
#Next deg gene programs
#from project directory working directory

conda activate ./environment/Seurat
Rscript ./scripts/09_deg_gene_programs.R
conda deactivate

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+