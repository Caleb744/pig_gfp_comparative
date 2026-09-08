# pig_gfp_comparative
The scripts and environments to rerun analyses associated with "Transgene Expression Characterization and Cross-Species Striatal Gene Expression in eGFP-edited Pigs"
## Overview
Scripts are numbered, with a master bash file to run each script. These correspond to the manuscript's Methods sections, with script titles. While each script should be runnable, we suggest a mindful approach for your system and desired directory structure. 
## Requirements
An environments directory is included with conda recipes for the environments these analyses were run in. There are some software that our cluster had as modules so version numbers are available. We encourage users to visit the original developers for download instructions. 

R: see environment/title_sessionInfo.txt that corresponds to the appropriate R script for full package versions.

Python: see conda recipes for python versions and supporting software

CellRanger v7.1.0 (alignement steps)
## Data
Raw sequencing: [SRA accession]

Human and mouse datasets are publicly avaiable. Download instrcutions are avaianable in scripts/04_annotation_strategy.R. 

Human: Siletti, K., Hodge, R., Mossi Albiach, A., Lee, K. W., Ding, S.-L., Hu, L., Lönnerberg, P., Bakken, T., Casper, T., Clark, M., Dee, N., Gloe, J., Hirschstein, D., Shapovalova, N. V., Keene, C. D., Nyhus, J., Tung, H., Yanny, A. M., Arenas, E., … Linnarsson, S. (2023). Transcriptomic diversity of cell types across the adult human brain. _Science, 382_(6667), eadd7046. https://doi.org/10.1126/science.add7046

Mouse: Yao, Z., van Velthoven, C. T. J., Kunst, M., Zhang, M., McMillen, D., Lee, C., Jung, W., Goldy, J., Abdelhak, A., Aitken, M., Baker, K., Baker, P., Barkan, E., Bertagnolli, D., Bhandiwad, A., Bielstein, C., Bishwakarma, P., Campos, J., Carey, D., … Zeng, H. (2023). A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain. _Nature, 624_(7991), 317–332. https://doi.org/10.1038/s41586-023-06812-z
## Running the scripts
From chosen local directory update output and input paths as needed, then run scripts in numeric order (00 through 09). Each script is designed to output where directed in the script or the master bash file. All can be updated for the user.
