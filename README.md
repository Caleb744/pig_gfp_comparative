# pig_gfp_comparative
The scripts and environments to rerun analyses associated with "Transgene Expression Characterization and Cross-Species Striatal Gene Expression in eGFP-edited Pigs"
## Overview
Scripts are numbered, with a master bash file to run each script. These correspond to the manuscript's Methods sections and their script titles. While each script should be runnable, we suggest a mindful approach for your system and desired directory structure. 
## Requirements
An environments directory includes conda recipes for the environments these analyses were run in. Some software was available as cluster modules, so version numbers are provided. We encourage users to visit the original developers for download instructions. 
### R environments
R: see conda yml files for R versions and supporting software
### Python environments
Python: see conda yml files for Python versions and supporting software

### Major software used not in conda recipes
#### Alignment
CellRanger v7.1.0 (alignment steps)
#### QC and Filtering
samtools v1.20 (Using htslib v1.20)


 
## Data
Raw sequencing: [SRA accession]

Human and mouse datasets are publicly available. Download instructions are available in scripts/04_annotation_strategy.R. 

Human: Siletti, K., Hodge, R., Mossi Albiach, A., Lee, K. W., Ding, S.-L., Hu, L., Lönnerberg, P., Bakken, T., Casper, T., Clark, M., Dee, N., Gloe, J., Hirschstein, D., Shapovalova, N. V., Keene, C. D., Nyhus, J., Tung, H., Yanny, A. M., Arenas, E., … Linnarsson, S. (2023). Transcriptomic diversity of cell types across the adult human brain. _Science, 382_(6667), eadd7046. https://doi.org/10.1126/science.add7046

Mouse: Yao, Z., van Velthoven, C. T. J., Kunst, M., Zhang, M., McMillen, D., Lee, C., Jung, W., Goldy, J., Abdelhak, A., Aitken, M., Baker, K., Baker, P., Barkan, E., Bertagnolli, D., Bhandiwad, A., Bielstein, C., Bishwakarma, P., Campos, J., Carey, D., … Zeng, H. (2023). A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain. _Nature, 624_(7991), 317–332. https://doi.org/10.1038/s41586-023-06812-z
## Running the scripts
From the chosen local directory, update output and input paths as needed, then run the scripts in numeric order (00 through 09). Each script outputs where directed in the script or the master bash file. All can be updated for the user.
