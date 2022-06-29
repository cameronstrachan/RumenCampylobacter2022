import os, sys
import re

# co-assemble reads from metagenomes that have already undergone trimming and host read removal

threads = 115
memory = 970
raw_data_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/metagenomes/reads/'
output_folder  = '/data/Unit_LMM/selberherr-group/strachan/campy/metagenomes/assemblies/'

compiled_fastq_R1 = raw_data_folder + 'forward_read_concat.fastq'
compiled_fastq_R1_zip = raw_data_folder + 'forward_read_concat.fastq.gz'

compiled_fastq_R2 = raw_data_folder + 'reverse_read_concat.fastq'
compiled_fastq_R2_zip = raw_data_folder + 'reverse_read_concat.fastq.gz'

if not os.path.exists(compiled_fastq_R1_zip):
	
	command0 = 'cat ' + raw_data_folder + '*_R1.fastq' + ' > ' +  compiled_fastq_R1
	command1 = 'gzip ' + compiled_fastq_R1
	
	print(command0 + '\n')
	os.system(command0)
	print(command1 + '\n')
	os.system(command1)

if not os.path.exists(compiled_fastq_R2_zip):

        command2 = 'cat ' + raw_data_folder + '*_R2.fastq' + ' > ' +  compiled_fastq_R2
        command3 = 'gzip ' + compiled_fastq_R2

        print(command2 + '\n')
        os.system(command2)
        print(command3 + '\n')
        os.system(command3)

spades_output  = output_folder  + 'coassembly'
coassembly_contigs = output_folder + 'coassembly.fasta'

if not os.path.exists(coassembly_contigs):
	command4 = '/home/strachan/bin/SPAdes-3.15.2/metaspades.py -m ' + str(memory) + ' -t ' + str(threads) + ' -1 ' + compiled_fastq_R1_zip + ' -2 ' + compiled_fastq_R2_zip + ' -o ' + spades_output
	print(command4 + '\n') 
	os.system(command4)
	command5 = 'cp ' + spades_output + '/contigs.fasta ' + coassembly_contigs
	print(command5 + '\n')
	os.system(command5)

file_prefix = 'coassembly'
output_file_name_prinseq = output_folder + file_prefix + '_min1500'
output_file_loc = output_folder + file_prefix + '_min1500.fasta'
input_file = output_folder + file_prefix + '.fasta'

if not os.path.exists(output_file_loc):
	command6 = 'perl /home/strachan/miniconda3/envs/qc/bin/prinseq-lite.pl -fasta ' + input_file + ' -min_len 1500 -out_good ' + output_file_name_prinseq + ' -out_format 1 -out_bad null -verbose'
	print(command6 + '\n')
	os.system(command6)



