###############################################################################################
#To do anything we have to align to a reference. This chunk is the steps used for building the reference.
#####Alignement Reference Build
#In working directory
cd ./results/Alignement/GENOME
#download refrence sequnce
wget https://ftp.ensembl.org/pub/release-112/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
#unzip the reference
gunzip Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
#add GFP to reference
cat ./Sus_scrofa.Sscrofa11.1.dna.toplevel.fa ./NSRRC_RRID_0016_Transgene.txt > Sscrofal_GFP_ref.fa
#regzip new reference
gzip -v -9 ./Sscrofal_GFP_ref.fa

#download refrence annotation
wget https://ftp.ensembl.org/pub/release-112/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.112.gtf.gz
#Make new chromosome and gene
echo -e 'GFP_Transgene\tunknown\texon\t1\t2071\t.\t+\t.\tgene_id "GFP"; transcript_id "GFP"; gene_name "GFP"; gene_biotype "protein_coding";' > GFP.gtf
echo -e 'Selectable_Marker\tunknown\texon\t1\t1108\t.\t+\t.\tgene_id "sel_mark"; transcript_id "sel_mark"; gene_name "sel_mark"; gene_biotype "protein_coding";' >> GFP.gtf
#Add to reference
cat Sus_scrofa.Sscrofa11.1.112.gtf GFP.gtf > Sus_scrofa.112.GFP.gtf

#Make alignement reference with cellranger.

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=Reference_Build
#SBATCH --output=Reference_Build.out
#SBATCH --error=Reference_Build.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --account=general
#SBATCH --partition=general

module load cellranger/v7.1.0
cd $PWD
cellranger mkref --genome=ensembl_sus_scrofa_11.1 \
  --fasta="./Sscrofal_GFP_ref.fa" \
  --genes="./Sus_scrofa.112.GFP.gtf" \
  --memgb 240 \
  --nthreads 64
EOF

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+