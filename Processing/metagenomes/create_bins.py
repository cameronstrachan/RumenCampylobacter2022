import os, sys
import re

# create bins from coassembled metagenome

threads = 115
raw_data_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/metagenomes/reads/'
output_folder  = '/data/Unit_LMM/selberherr-group/strachan/campy/metagenomes/assemblies/'
intermediate_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/metagenomes/bwa_intermediate_files/'

fastq_files = [f for f in os.listdir(raw_data_folder) if f.endswith(('_R1.fastq', '_R2.fastq'))]
sample_names = list(set([re.sub('_R\d.fastq', '', s) for s in fastq_files]))

ref = output_folder + 'coassembly_min1500.fasta'

command0 = 'bwa index ' + ref
print(command0 + '\n')
os.system(command0)
os.system('mv ' + output_folder + 'coassembly_min1500.fasta.* ' + output_folder + 'bwa_index/')

database = output_folder + 'bwa_index/coassembly_min1500.fasta'

for sample in sample_names:
	
	forward_read_trim = raw_data_folder + sample + '_R1.fastq'
	reverse_read_trim = raw_data_folder + sample + '_R2.fastq'
	
	sample_clean = re.sub('_clean', '', sample)
	output_file_loc = intermediate_folder + sample_clean + '.sam'
	output_file_loc_bam = intermediate_folder + sample_clean + '.bam'
	output_file_loc_bam_sort = intermediate_folder + sample_clean + '.sort.bam'

	
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


output_file_loc = '/home/strachan/master/campy/metagenomes/output/coassembly_min1000_depth.txt'
bin_location = output_folder + 'bins/bin'


if not os.path.exists(output_file_loc):
	
	command5 = 'jgi_summarize_bam_contig_depths --outputDepth ' + output_file_loc + ' ' + intermediate_folder + '*.sort.bam' 
	print(command5 + '\n')
	os.system(command5)	

	command6 = 'metabat2 -i ' + ref + ' -a ' + output_file_loc + ' -o ' + bin_location + ' -m 1500 -t ' + str(threads) + ' -v'
	print(command6 + '\n')
	os.system(command6) 
