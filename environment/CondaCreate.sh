#Conda recipes and YML files to perform the analyses
##NOTE: Some tools were availalbe as modules on our local cluster and are identified. These will need to be setup for some steps.

##cellbender
#Used for Ambiant and Background read cleaning
if [ ! -d ./CellBender ] ; then
     conda create --yes --prefix ./CellBender python=3.7
     conda config --set env_prompt '({name})'
     # ^ This changes the env name to /conda and not the full path     
     conda config --add channels defaults
     conda config --add channels conda-forge
     conda config --set channel_priority flexible
     # Put the channels in a specific order
     # Download packages flexibly
fi

conda activate CellBender
conda install conda-forge::pytables -p ./CellBender
conda install conda-forge::torchvision -p ./CellBender
conda install pytorch torchvision -c pytorch
git clone https://github.com/broadinstitute/CellBender.git
pip install CellBender
#or
conda env create -f CellBender.yml --prefix ./CellBender



##BARCODE READS
#Used for intergenic and intragenic read assignement
if [ ! -d ./barcode_reads ] ; then
     conda create --yes --prefix ./barcode_reads python=3.10 
     conda config --set env_prompt '({name})'
     # ^ This changes the env name to /conda and not the full path     
fi


conda install bioconda::pysam -p ./barcode_reads
conda install conda-forge::intervaltree -p ./barcode_reads
conda install conda-forge::zstandard
#or
conda env create -f Barcode_processing.yml --prefix ./barcode_reads


##Seurat
#Used for most steps of SC processing and annotation
if [ ! -d ./Seurat ] ; then
     conda create --yes --prefix ./Seurat
     conda config --set env_prompt '({name})'
     # ^ This changes the env name to /conda and not the full path     
     conda config --add channels defaults
     conda config --add channels conda-forge
	 conda config --add channels bioconda
     conda config --set channel_priority flexible
     # Put the channels in a specific order
     # Download packages flexibly
fi

conda install conda-forge::r-essentials -p ./Seurat
conda install conda-forge::pytables -p ./Seurat
conda install conda-forge::r-openssl -p ./Seurat
conda install conda-forge::r-igraph -p ./Seurat
conda install conda-forge::r-httr -p ./Seurat
conda install conda-forge::r-leiden -p ./Seurat
conda install conda-forge::r-plotly -p ./Seurat
conda install conda-forge::r-seurat -p ./Seurat
conda install conda-forge::r-cairo -p ./Seurat
conda install conda-forge::r-ggrastr -p ./Seurat
conda install conda-forge::r-qqconf -p ./Seurat
conda install conda-forge::r-devtools -p ./Seurat
conda install bioconda::bioconductor-limma -p ./Seurat
conda install conda-forge::r-clustree -p ./Seurat
conda install conda-forge::r-hdf5r -p ./Seurat
conda install conda-forge::r-biocmanager -p ./Seurat
conda install bioconda::bioconductor-singler -p ./Seurat
conda install bioconda::bioconductor-scdblfinder -p ./Seurat
conda install bioconda::bioconductor-scater -p ./Seurat
conda install bioconda::bioconductor-multtest -p ./Seurat
conda install bioconda::bioconductor-biomart -p ./Seurat
conda install bioconda::bioconductor-glmgampoi  -p ./Seurat
conda install bioconda::presto -p ./Seurat
conda install conda-forge::r-metap -p ./Seurat
conda install bioconda::r-harmony
conda install -p ./Seurat bioconda::bioconductor-dirichletmultinomial
conda install -p ./Seurat bioconda::bioconductor-tfbstools
conda install -p ./Seurat conda-forge::r-glmnet
conda install -p ./Seurat conda-forge::r-proc
conda install -p ./Seurat bioconda::bioconductor-mast
#For non conda tools. 
#Not exactly best practice to instal within a conda but this setup has worked for me
R
devtools::install_github('immunogenomics/presto') #no conda 
remotes::install_github("bnprks/BPCells/r")
remotes::install_github('satijalab/azimuth', ref = 'master')

#or
conda env create -f Single_Cell.yml --prefix ./Seurat


##CellTypist
#used for annotation
if [ ! -d ./CellTypist ] ; then
     conda create --yes --prefix ./CellTypist
     conda config --set env_prompt '({name})'
     # ^ This changes the env name to /conda and not the full path     
     conda config --add channels defaults
     conda config --add channels conda-forge
     conda config --add channels bioconda
     conda config --set channel_priority flexible
     # Put the channels in a specific order
     # Download packages flexibly
fi

conda install -p ./CellTypist bioconda::celltypist
#or
conda env create -f CellTypist.yml --prefix ./CellTypist

##Enrichment
#for cell type enrichment analyses
if [ ! -d ./Enrichment ] ; then
     conda create --yes --prefix ./Enrichment  -c conda-forge -c bioconda r-base=4.3 bioconductor-clusterprofiler bioconductor-org.ss.eg.db
     conda config --set env_prompt '({name})'
     # ^ This changes the env name to /conda and not the full path     
fi
conda install bioconda::bioconductor-gseabase -p ./Enrichment
conda install bioconda::bioconductor-fgsea -p ./Enrichment
conda install conda-forge::r-gprofiler2 -p ./Enrichment
conda install conda-forge::r-upsetr -p ./Enrichment
conda install -p ./Enrichment/ conda-forge::r-tidyverse
conda install -p ./Enrichment/ bioconda::bioconductor-org.gg.eg.db
#or
conda env create -f enrichment.yml --prefix ./Enrichment


##Stats
#for comparative analyses
if [ ! -d ./Stats ] ; then
     conda create --yes --prefix ./Stats  
     conda config --set env_prompt '({name})'
     # ^ This changes the env name to /conda and not the full path     
fi

conda install --yes conda-forge::r-essentials -p ./Stats
conda install --yes r-lme4 -p ./Stats
conda install --yes r-matchit -p ./Stats
conda install --yes r-caret -p ./Stats
conda install --yes r-mediation -p ./Stats
conda install --yes r-survey -p ./Stats
conda install --yes r-proc -p ./Stats
#or
conda env create -f Stats.yml --prefix ./Stats

