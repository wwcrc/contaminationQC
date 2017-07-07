#! /bin/bash

echo 'Testing help menu...'
perl reference_panel_build.pl --help
echo

echo 'Checking version...'
perl reference_panel_build.pl --version
echo

echo 'Testing missing parameter...'
perl reference_panel_build.pl -m test.manifest -c test.json -o test01.output test_sample01.bam test_sample02.bam
echo

echo 'Testing missing input files...'
perl reference_panel_build.pl -m test.manifest -c test.json -b test.bed -o test01.output
echo

echo 'Testing typo in gender...'
perl reference_panel_build.pl -m typo_gender.manifest -c test.json -b test.bed -o test01.output test_sample01.bam test_sample02.bam
echo

echo 'Testing bed file split...'
perl reference_panel_build.pl -m test.manifest -c test.json -b test.bed -o test01.output test_sample01.bam test_sample02.bam 1>tmp1 2>tmp2
echo 'Expect autosomal and XY in separate temp files'
echo 'temp file 1 (XY SNPs):'
cat temp_xy.bed
echo 'temp file 2 (autosomal SNPs):'
cat temp_aut.bed
echo

echo 'Testing output...'
echo 'values should be x mean 0, y mean 1, x sd 0, y sd 0'
echo 'output file is:'
python -m json.tool test01.output
echo

echo 'Testing warnings...'
echo 'should indicate insufficient number of samples'
echo 'warnings are:'
cat tmp2
echo

echo 'Testing info...'
echo 'should indicate that sample 01 is male and sample 02 is female, Xvalues and Yvalues are 0 and 1'
echo 'info is:'
cat tmp1
echo

echo 'Testing low depth on X...'
perl reference_panel_build.pl -m lowdepth.manifest -c test.json -b test.bed -o test01.output Xdepth99.bam test_sample02.bam 1>tmp1 2>tmp2
echo 'Xvalues should be NA'
cat tmp1
echo

echo 'Testing minor allele present on X...'
perl reference_panel_build.pl -m sample03.manifest -c test.json -b test.bed -o test01.output test_sample03.bam test_sample02.bam 1>tmp1 2>tmp2
echo 'expected X mean is 0.01'
echo 'output file is:'
python -m json.tool test01.output
echo

rm temp_xy.bed temp_aut.bed test01.output tmp1 tmp2

