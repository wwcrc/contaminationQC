#!/usr/bin/perl

# contamination_detector.pl - A Perl script to look for contamination in targeted sequencing experiments.
# Copyright (C) 2015 14MG, 2017 University of Glasgow
# Author: Susie Cooke
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Requires two BED files as input, one of SNP positions on X and Y and one of SNP positions
# for the rest of the genome. Y coordinates are ranges, X coordinates are single SNP
# positions and REF and ALT alleles are given in the fourth column as a comma separated list.
# Multi-allelic SNPs are not permitted. Requires a bam/cram/sam file as input.
# Requires a config file. Requires the gender of the sample. Result is appended to a specified json file.

use strict;
use warnings;

use Getopt::Long;
use JSON;
use Scalar::Util qw(looks_like_number);

my $version = '1.0';

# Defaults

my $xy_bed = 'temp_xy.bed'; # temp file created by script
my $aut_bed = 'temp_aut.bed'; # temp file created by script
my $total_depth_cutoff = 30; # minimum cumulative depth required on chrX
my $MQ = 60; # default value can be over-written by json config file
my $BAQ = 13; # default value can be over-written by json config file
my $depth_cutoff = 5; # default value can be over-written by json config file
my $Xchrom = 'chrX';
my $Ychrom = 'chrY';

# Output

my $contamination_level;


### Command line and options ###

my $config_file;
my $gender;
my $bed_file;
my $output;
my $verbose = 0;

GetOptions (	"j=s" => \$output,
		"g=s" => \$gender,
		"c=s" => \$config_file,
		"b=s" => \$bed_file,
		"verbose" => \$verbose,
		"version" => \&version,
		"help" => \&help )
or die "Error in command line arguments\n";


my $bamfile = $ARGV[0];

unless (defined $output && defined $gender && defined $config_file && defined $bed_file && defined $bamfile) {
	print STDERR "Missing one or more parameters. Try --help\n";
	exit 1;
}

if (scalar @ARGV != 1){
	print STDERR "Expected single cram/bam file as input, got @ARGV\n";
	exit 1;
}


### Set parameters from config ###

open (CONFIG, '<', $config_file)
	or die "Cannot open file: $config_file: $!\n";

my $json_content;

while (<CONFIG>) {
	$json_content .= $_;
}
close CONFIG;

my $thresholds = decode_json $json_content;

my $X_mean = $thresholds->{'X_mean'};
my $X_sd = $thresholds->{'X_sd'};
my $Y_mean = $thresholds->{'Y_mean'};
my $Y_sd = $thresholds->{'Y_sd'};
$MQ = $thresholds->{'MQ'}; # Override default
$BAQ = $thresholds->{'BAQ'}; # Override default
$depth_cutoff = $thresholds->{'depth_cutoff'}; # Override default


### Get chromosome nomenclature from bam/cram files ###

open HEADER, "samtools view -H $bamfile|"
    or die "Cannot execute samtools: $!\n";

while (<HEADER>) {
    chomp;
    next unless $_ =~ /^\@SQ/;
    my @data = split ('\t', $_);
    foreach my $i (0 .. $#data) {
        if ($data[$i] =~ /^SN:(([cC][hH][rR])?X$)/) {
            $Xchrom = $1;
        }
        elsif ($data[$i] =~ /^SN:(([cC][hH][rR])?Y$)/) {
            $Ychrom = $1;
        }
    }
}


### Split bed file into temporary XY and autosome bed files ###

open (BED, '<', $bed_file)
	or die "Cannot open file for reading: $bed_file $!\n";
open (XY, '>', $xy_bed)
	or die "Cannot open file for writing: $xy_bed $!\n";
open (AUT, '>', $aut_bed)
	or die "Cannot open file for writing: $aut_bed $!\n";

while (<BED>) {
	if ($_ =~ /^([XxYy]|chr[XxYy])/) {
		print XY $_;
	}
	else {
		next if $_ =~ /^#/;
		print AUT $_;
	}
}
close BED;
close XY;
close AUT;


### Select appropriate flows to run ###

my $json_field;

if (lc($gender) eq 'female') {
	print "Sample is: $gender, running Y_based routine\n" if $verbose; #verbose only
	$contamination_level = Y_based($xy_bed, $bamfile, $Ychrom, $aut_bed, $MQ);
	$json_field->{'CONTAMINATION_IDX'} = $contamination_level;
}

elsif (lc($gender) eq 'male') {
	print "Sample is: $gender, running X_based routine\n" if $verbose; #verbose only
	$contamination_level = X_based($xy_bed, $bamfile, $Xchrom, $BAQ);
	$json_field->{'CONTAMINATION_IDX'} = $contamination_level;
}
else {
	print STDERR "Sample is: $gender, running both metrics.\n";
    my $cont_if_male = X_based($xy_bed, $bamfile, $Xchrom, $BAQ);
    my $cont_if_female = Y_based($xy_bed, $bamfile, $Ychrom, $aut_bed, $MQ);
	$json_field->{'CONTAMINATION_IDX_M'} = $cont_if_male;
	$json_field->{'CONTAMINATION_IDX_F'} = $cont_if_female;
}

# Output result to a json file

open (OUTPUT, '>', $output)
	or die "Cannot open file: $output: $!\n";

print OUTPUT encode_json($json_field);
close OUTPUT;

unlink $xy_bed
	or die "Can't remove $xy_bed: $!\n";
unlink $aut_bed
	or die "Can't remove $aut_bed: $!\n";


##################### SUBROUTINES ########################

### Y-based flow (for female samples) ###

sub Y_based {
	my ($xy_bed, $file, $Ychrom, $aut_bed, $MQ) = ($_[0], $_[1], $_[2], $_[3], $_[4]);
    # calculate counts at low and high MQ thresholds
	my $all_Y_reads = count($xy_bed, $file, 0, $Ychrom);
	my $all_wg_reads = count($aut_bed, $file, 0);
    my $Y_thresholded = count($xy_bed, $file, $MQ, $Ychrom);
    my $wg_thresholded = count($aut_bed, $file, $MQ);

	print "allY:$all_Y_reads, allWG:$all_wg_reads, goodY:$Y_thresholded, goodWG:$wg_thresholded\n" if $verbose; #verbose only

	# Check for divide by zero errors
	if ($all_wg_reads == 0 || $all_Y_reads < 1 || $wg_thresholded < 11) {
		return 'Unknown'; # can't evaluate if insufficient depth
	}

	# Calculate normalised fraction of chrY reads with high MQ
	print "Comparing to normal panel...\n" if $verbose; #verbose only
	my $cont_value = ($Y_thresholded/$wg_thresholded)/($all_Y_reads/$all_wg_reads);
	print "contamination value:$cont_value\n" if $verbose; #verbose only

	# Check normalised depth against unmatched normal panel
	if ($Y_sd == 0) {
		die "Problem with json config file\n";
	}
	if ($cont_value > $Y_mean) {
		return ($cont_value-$Y_mean)/$Y_sd;
	}
	else { return 0; }
}


### X-based flow (for male samples) ###

sub X_based {
	my ($xy_bed, $file, $Xchrom, $BAQ) = ($_[0], $_[1], $_[2], $_[3]);
	my %SNPs;
	open (BEDFILE, '<', $xy_bed)
		or die "Cannot open file: $xy_bed: $!\n";

	# Make a hash of SNP positions and ref/alt alleles
	while (<BEDFILE>) {
		chomp;
		my @data = split('\t', $_);
		if ($data[0] eq $Xchrom) {
			my $location = join ('_', $data[0], $data[2]);
			$SNPs{$location} = $data[3];
		}
	}
	close BEDFILE;

	if ($verbose) {	#verbose only
		foreach (sort keys %SNPs) {
			print "$_\t$SNPs{$_}\n";
		}
	}

	# Generate pileup data
	my $sum_non_genotype = 0;
	my $total_depth = 0;
	my $valid_positions = 0;

	print "Running pileup on chrX...\n" if $verbose; #verbose only
	open PILEUP, "samtools mpileup -l $xy_bed -r $Xchrom -Q $BAQ $file|";

	while (<PILEUP>) {
		chomp;
		my @data = split('\t', $_);
		$data[4] =~ s/\^.//g; #remove read starts and their MQs (represented as ASCII character)
		my @bases = split('', uc($data[4]));
		my %base_counts = ('A' => 0, 'C' => 0, 'G' => 0, 'T' => 0);
		foreach my $i (0 .. $#bases) {
			if ($bases[$i] =~ /[ACGT]/) {
				$base_counts{$bases[$i]}++;
			}
		}
		# Correct for non-ref bases occurring in the representation of indels
		my @indels = ();
		while (defined $data[4]) {
			if ($data[4] =~ /([\+-][0-9]+[ACGTNacgtn]+)/) { #find indel
				push (@indels, $1); #collect indel
				$data[4] =~ s/[\+-][0-9]+[ACGTNacgtn]+//; #remove indel
			}
			else {
				undef $data[4]; #if position contains no indels empty it to exit the loop
			}
		}
		if (scalar @indels > 0) { #if any indels were present in the string subtract them from the base counts
			foreach (@indels) {
				$_ =~ /[\+-]([0-9]+)([ACGTNacgtn]+)/;
				my $indel_length = $1;
				my @indel_seq = split ('', $2);
				foreach my $i (0..($indel_length-1)) {
					next unless $indel_seq[$i] =~ /[ACGT]/;
					$base_counts{$indel_seq[$i]}--;
				}
			}
		}

		# Calculate level of non-genotype allele present
		my $non_genotype = 0;
		my $depth = 0;
		my $lookup_location = join ('_', $data[0], $data[1]);
		my @alleles = split (',', $SNPs{$lookup_location});
		my $first_allele_depth = $base_counts{uc($alleles[0])};
		my $second_allele_depth = $base_counts{uc($alleles[1])};
		$depth = $first_allele_depth + $second_allele_depth;
		print "$alleles[0]:$first_allele_depth, $alleles[1]:$second_allele_depth, depth:$depth\n" if $verbose; #verbose only

		next if $depth < $depth_cutoff; # ignore position if depth too low

		$valid_positions++;
		if ($first_allele_depth >= $second_allele_depth) {
			$non_genotype = $second_allele_depth;
		}
		else { $non_genotype = $first_allele_depth; }
		$sum_non_genotype = $sum_non_genotype + $non_genotype;
		$total_depth = $total_depth + $depth;
	}

	close PILEUP;

	# Check average non-genotype proportion against unmatched normal panel
	print "Comparing to normal panel...\n" if $verbose; #verbose only
	if ($valid_positions == 0 || $total_depth < $total_depth_cutoff) {
		return 'Unknown';
	}
	print "valid:$valid_positions\n" if $verbose; #verbose only
	my $ave_non_genotype = $sum_non_genotype/$total_depth;
	print "contamination:$ave_non_genotype\n" if $verbose; #verbose only
	if ($X_sd == 0) {
		die "Problem with json config file\n";
	}
	if ($ave_non_genotype > $X_mean) {
		return ($ave_non_genotype-$X_mean)/$X_sd;
	}
	else { return 0; }
}


### Count reads ###
sub count {
	print "Counting reads...\n";
	my ($bed, $bam, $mq) = ($_[0], $_[1], $_[2]);
	my $region = '';
	if (defined $_[3]) {
		$region = $_[3];
	}
	open COUNT, "samtools view -c -L $bed -q $mq $bam $region|"
		or die "Cannot execute samtools:$!\n";

	my $count = <COUNT>;
	close COUNT;
	chomp $count;
	return $count;
}

### Version information ###
sub version {
	print "Version $version\n";
	exit 0;
}

### Help information ###
sub help {
	print "Usage: contaminationQC.pl [options] in.cram
	Options:\t-b BED\tBed file of targeted regions on Y and targeted SNP positions on X
	\t\t-g STR\tgender
	\t\t-c FILE\tconfiguration file
	\t\t-j FILE\toutput file name (output file will be json format)
	\t\t--version
	\t\t--verbose
	\t\t--help\n";
	exit 0;
}

