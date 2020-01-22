#!/usr/bin/perl

# Copyright (C) 2015 14MG, 2017 University of Glasgow
# Susie Cooke

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

# This script generates the configuration json file needed for running contaminationQC


use strict;
use warnings;
use Getopt::Long;
use JSON;

my $version = '1.0';
my $xy_bed = 'temp_xy.bed'; # temp file created by script
my $aut_bed = 'temp_aut.bed'; # temp file created by script
my $total_depth_cutoff = 100; # minimum cumulative depth required on chrX
my $MQ = 60; # default value can be over-written by json config file
my $BAQ = 13; # default value can be over-written by json config file
my $depth_cutoff = 5; # default value can be over-written by json config file

### Command line and options ###

my $manifest; # details of the samples to be used
my $config; # a json file containing thresholds for MQ BAQ and depth
my $bed; # a 'bedplus' file of positions to use
my $output; # name of output json file


GetOptions (	"m=s" => \$manifest,
		"c=s" => \$config,
		"b=s" => \$bed,
		"o=s" => \$output,
		"version" => \&version,
		"help" => \&help )
or die "Error in command line arguments. Try --help\n";

my @files = @ARGV; # can be bam/cram format

if (scalar @ARGV <= 1) {
	die "Expected multiple bam/cram files as input, got @ARGV\n";
}

unless (defined $manifest && defined $config && defined $bed && defined $output) {
	die "Missing one or more parameters. Try --help\n";
}


### Set parameters from config ###

open (CONFIG, '<', $config)
	or die "Cannot open $config: $!\n";

my $content;
while (<CONFIG>) {
	$content .= $_;
}
close CONFIG;

my $thresholds = decode_json $content;

$MQ = $thresholds->{'MQ'}; # cutoff to use for mapping quality
$BAQ = $thresholds->{'BAQ'}; # cutoff to use for base quality
$depth_cutoff = $thresholds->{'depth_cutoff'}; # cutoff to use for depth of sequencing on chrX


### Get chromosome nomenclature from bam/cram files ### FIXME, also should check all files are same assembly

my $Xchrom;
my $Ychrom;

open HEADER, "samtools view -H $files[0]|"
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


### Build hash of samples and genders ###

my @X_values;
my @Y_values;
my $female_count = 0;
my $male_count = 0;

open (GENDERS, '<', $manifest)
    or die "Cannot open file: $manifest $!\n";

my %samples;

while (<GENDERS>) {
    chomp;
    my @data = split ();
    $samples{$data[0]} = $data[1];
}
close GENDERS;


### Split bed file into temporary XY and autosome bed files ###

open (BED, '<', $bed)
	or die "Cannot open file for reading: $bed: $!\n";
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

foreach my $file (@files) {
    if (lc($samples{$file}) eq 'female') {
        print "$file is $samples{$file}, running Y_based routine\n";
        $female_count++;
        my $Y_level = Y_based($xy_bed, $file, $Ychrom, $aut_bed, $MQ);
        push (@Y_values, $Y_level);
    }
    elsif (lc($samples{$file}) eq 'male') {
        print "$file is $samples{$file}, running X_based routine\n";
        $male_count++;
        my $X_level = X_based($xy_bed, $file, $Xchrom, $BAQ);
        push (@X_values, $X_level);
    }
    else {
        die "$file had an unrecognised gender: $samples{$file}\n";
    }
}


### Output results ###

print "Normal panel consists of $female_count female samples, $male_count male samples\n";
if ($female_count < 30 || $male_count < 30) {
    print STDERR "Warning: at least 30 male and 30 female samples should be used to approximate normal distributions\n"
}
print "For QC:\n";
print "X_values <- c(", join (',', @X_values), ")\n";
print "Y_values <- c(", join (',', @Y_values), ")\n";

$thresholds->{'X_mean'} = mean(@X_values);
$thresholds->{'Y_mean'} = mean(@Y_values);
$thresholds->{'X_sd'} = sd(@X_values);
$thresholds->{'Y_sd'} = sd(@Y_values);

open (OUTPUT, '>', $output)
    or die "Cannot open file for writing: $output $!\n";
print OUTPUT encode_json $thresholds;
close OUTPUT;


##################### SUBROUTINES ########################

### Y-based flow (for female samples) ###

sub Y_based {
	my ($xy_bed, $file, $Ychrom, $aut_bed, $MQ) = ($_[0], $_[1], $_[2], $_[3], $_[4]);
    # calculate counts at low and high MQ thresholds
    my $all_Y_reads = count($xy_bed, $file, 0, $Ychrom);
    my $all_wg_reads = count($aut_bed, $file, 0);
    my $Y_thresholded = count($xy_bed, $file, $MQ, $Ychrom);
    my $wg_thresholded = count($aut_bed, $file, $MQ);

    # check for divide by zero errors
    if ($all_wg_reads == 0 || $all_Y_reads == 0 || $wg_thresholded == 0) {
        return 'NA'; # can't evaluate if insufficient depth
    }

    # calculate normalised fraction of chrY reads with high MQ
    else { return (($Y_thresholded/$wg_thresholded)/($all_Y_reads/$all_wg_reads)); }
}


### X-based flow (for male samples) ###

sub X_based {
	my ($xy_bed, $file, $Xchrom, $BAQ) = ($_[0], $_[1], $_[2], $_[3]);
	my %SNPs;
    open (BEDFILE, '<', $xy_bed)
        or die "Cannot open bed file of X/Y positions: $xy_bed: $!\n";
    # make a hash of SNP positions and ref/alt alleles
    while (<BEDFILE>) {
        chomp;
        my @data = split ('\t', $_);
        if ($data[0] eq $Xchrom) {
            my $location = join (':', $data[0], $data[2]);
            $SNPs{$location} = $data[3];
        }
    }
    close BEDFILE;

    # generate pileup data
    my $sum_minor_genotype = 0;
    my $total_depth = 0;
    my $valid_positions = 0;

    open PILEUP, "samtools mpileup -l $xy_bed -r $Xchrom -Q $BAQ $file|"
        or die "Cannot execute samtools:$!\n";

    while (<PILEUP>) {
        chomp;
        my @data = split (/\t/, $_);
        $data[4] =~ s/\^.//g; #remove read starts and their MQs (represented as ASCII character)
        my @bases = split ('', uc($data[4]));
        my %base_counts = ('A' => 0, 'C' => 0, 'G' => 0, 'T' => 0);
        foreach (@bases) {
			next unless $_ =~ m/[ACGT]/;
			$base_counts{$_}++;
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
					next unless $indel_seq[$i] =~ m/[ACGTacgt]/;
					$indel_seq[$i] =~ tr/acgt/ACGT/;
					$base_counts{$indel_seq[$i]}--;
				}
			}
		}
        # calculate level of minor allele
        my $minor_level = 0;
        my $depth = 0;
        my $lookup_location = join (':', $data[0], $data[1]);
        my @alleles = split (',', $SNPs{$lookup_location});
        my $first_allele_depth = $base_counts{uc($alleles[0])};
        my $second_allele_depth = $base_counts{uc($alleles[1])};
        $depth = $first_allele_depth + $second_allele_depth;
        next if $depth < $depth_cutoff; # ignore position if depth too low
        $valid_positions++;

        if ($first_allele_depth >= $second_allele_depth) { # assume most prevalent genotype is correct sample
            $minor_level = $second_allele_depth;
        }
        else { $minor_level = $first_allele_depth; }

        $sum_minor_genotype = $sum_minor_genotype + $minor_level; # calculate cumulative level of contamination
        $total_depth = $total_depth + $depth;
    }
    close PILEUP;

    if ($total_depth < $total_depth_cutoff) { # can't evaluate if insufficient depth across all positions
        return 'NA';
    }
    else { return ($sum_minor_genotype/$total_depth); }
}


### Count reads ###

sub count {
    my ($bed, $cram, $mq) = ($_[0], $_[1] ,$_[2]);
    my $region = '';
    if (defined $_[3]) {
        $region = $_[3];
    }

    open COUNT, "samtools view -c -L $bed -q $mq $cram $region|"
        or die "Cannot execute samtools:$!\n";
    my $count = <COUNT>;
    close COUNT;
    chomp $count;
    return $count;
}


### Calculate mean ###

sub mean {
    my $total = 0;
    my $counter = 0;
    foreach (@_) {
        next if $_ =~ /NA/;
        $total += $_;
        $counter++;
    }
    if ($counter > 0) {
        return ($total/$counter);
    }
    else { return 0; }
}


### Calculate standard deviation ###

sub sd {
    my $mean_value = mean(@_);
    my @sqrs;
    foreach (@_) {
        next if $_ =~ /NA/;
        my $value = ($_-$mean_value)*($_-$mean_value);
        push (@sqrs, $value);
    }
    my $total = 0;
    my $counter = 0;
    foreach (@sqrs) {
        $total += $_;
        $counter++;
    }
    if ($counter > 1) {
        my $variance = $total/($counter-1); # Bessel's correction
        return (sqrt($variance));
    }
    else { return 0; }
}

### Version information ###

sub version {
    print "Version $version\n";
    exit 0;
}

### Help information ###

sub help {
    print <<EOF;
    Usage: reference_panel_build.pl [options] in.cram1 in.cram2 ...
    -c JSON confguration file
    -o STR output file name
    -m TXT sample file list
    -b BED regions of interest
EOF
exit 0;
}
