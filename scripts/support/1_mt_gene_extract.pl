#!/usr/bin/perl -w
#XQTL_1_1_MATRIX.pl: script for generating Gene Count Matrix. 
use strict;
use autodie;
use Getopt::Long;

my $help = ();				#initialized Help.

my $GTF = ();				#"/storage/hpc/group/UMAG_test/WORKING/cs744/SC_Analysis/GENOME/ensembl_ARS_1.2_110/genes/genes.gtf.gz"
my $out = ".";
my $prefix = ();

my $options = GetOptions (
	'help|?'			=> \$help,
	'GTF=s'				=> \$GTF,
	'prefix=s'			=> \$prefix,
	'out=s'				=> \$out
	);

	
if ($help) {
print <<HELP;
*******************
USAGE
time perl "./2_extract_mt.pl" --GTF "/GTF_file.txt" --OUT_DIR ../OUT --PREFIX TISSSUE_GE
*******************
OPTIONS
	--help				print this message
	--GTF				GTF path
	--OUT_DIR			Full path to output directory (DEFAULT: ./OUT
	--PREFIX			name of cohort
*******************
HELP
exit;
}

my $genes_gtf = "${GTF}/genes/genes.gtf.gz";
my $genes_out = "${out}/${prefix}_mt_genes.txt";
my @mt;

open(genes_gtf, "gunzip -c $genes_gtf |");
open(OUT, ">$genes_out");
while (<genes_gtf>) {
	if (m/gene_id/) {
		my $input = $_;
		chomp $input;
		my @gtf_array = split("\t", $input);
		if ($gtf_array[2] eq 'gene') {
			my $chrm = $gtf_array[0];
			if ($chrm eq 'MT') {
				my $gNAME = $gtf_array[8];
				if ($gNAME =~ /gene_name \"(\S+)"/) {
					$gNAME =$1;
					my $array_push = "${gNAME}";
					push(@mt, $array_push);
					#print OUT "$gNAME,";
				}
			}
			else {next;}
		} else {next;}
	} else {next;}
}
my $mt_out = join(",", @mt);
print OUT "$mt_out\n";
close OUT;
close genes_gtf;