import os, sys
import re

# create bins from coassembled metagenome

threads = 120
bins_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/metagenomes/assemblies/bins/'
reads_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/metagenomes/reads/'
intermediate_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/metagenomes/bwa_intermediate_files/'
assembly_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/metagenomes/assemblies/'

os.system('cat ' + bins_folder + '*.fa > ' + bins_folder + 'all_bins.fasta')
ref = bins_folder + 'all_bins.fasta'
ref_check = bins_folder + 'bwa_index/all_bins.fasta.bwt'

if not os.path.exists(ref_check):
	command0 = 'bwa index ' + ref
	print(command0 + '\n')
	os.system(command0)
	os.system('mv ' + bins_folder + 'all_bins.fasta.* ' + assembly_folder + 'bwa_index/')

fastq_files = [f for f in os.listdir(reads_folder) if f.endswith(('_R1.fastq', '_R2.fastq'))]
sample_names = list(set([re.sub('_R\d.fastq', '', s) for s in fastq_files]))

database = assembly_folder  + 'bwa_index/all_bins.fasta'

for sample in sample_names:
	
	forward_read_trim = reads_folder + sample + '_R1.fastq'
	reverse_read_trim = reads_folder + sample + '_R2.fastq'
	
	sample_clean = re.sub('_clean', '', sample)
	output_file_loc = intermediate_folder + sample_clean + '.bins.sam'
	output_file_loc_bam = intermediate_folder + sample_clean + '.bins.bam'
	output_file_loc_bam_sort = intermediate_folder + sample_clean + '.bins.sort.bam'
	
	if not os.path.exists(output_file_loc):
		command1 = 'bwa mem -t ' + str(threads) + ' ' + database + ' ' + forward_read_trim + ' ' + reverse_read_trim + ' > ' + output_file_loc
		print(command1 + '\n')
		os.system(command1)
		
		command2 = 'samtools view -S -b ' + output_file_loc + ' > ' + output_file_loc_bam
		print(command2 + '\n')
		os.system(command2)

		command3 = 'samtools sort -@ ' + str(threads) + ' ' + output_file_loc + ' -o ' + output_file_loc_bam_sort
		print(command3 + '\n')
		os.system(command3)

		command4 = 'rm ' + output_file_loc_bam
		print(command4 + '\n')
		os.system(command4)

		command5 = 'samtools idxstats ' + output_file_loc_bam_sort + ' > output/mapping/' + sample_clean + '_bins_readcounts.txt'
		print(command5 + '\n')
		os.system(command5)

bin_files = [f for f in os.listdir(bins_folder) if f.endswith('.fa')]
output_file = open('output/bin_contig_map.csv', 'w')

for bin_file in bin_files:
	file = open(bins_folder + bin_file, 'r')
	for line in file:
		if line[0] == '>':
			contig_id = re.sub('>', '', line)
			bin_id = bin_file.split('.f')[0]
			output_file.write(bin_id + ',' + contig_id)

output_file.close()
