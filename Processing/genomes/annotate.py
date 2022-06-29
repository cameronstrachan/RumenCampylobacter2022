import os, sys
import re

threads = 80
gtdb_genomes_folder = '/home/strachan/miniconda3/envs/gtdbtk/share/gtdbtk-1.4.1/db/fastani/database/'
aln_file =  open('/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/gdtb_intermediate_files/align/gtdbtk.bac120.msa.fasta', 'r')
assembly_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/assemblies/'
annotations_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/annotations/'
prokka_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/annotations/prokka/'
prokka_ref_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/annotations/prokka_ref/'
prokka_ref2_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/annotations/prokka_ref2/'

for line in aln_file:
	if line[0] == '>':
		header = re.sub('>', '', line)		
		header_split = header.split()

		if len(header_split) > 1:		

			genome_id = re.sub("GB_", "", re.sub('RS_', '',  header_split[0]))
			classification = re.sub('d__Bacteria;p__Campylobacterota;c__Campylobacteria;o__Campylobacterales;f__Campylobacteraceae;', '', header_split[1]) +  " " + header_split[2]

			genome_file_name = genome_id + '_genomic.fna.gz'
			genome_file_loc = gtdb_genomes_folder + genome_file_name
			genome_file_name_unzip = genome_id + '_genomic.fna'
			genome_file_loc_unzip = assembly_folder + 'gdtb_genomes/' + genome_file_name_unzip

			if not os.path.exists(genome_file_loc_unzip):

				command0 = 'cp ' + genome_file_loc + ' ' + assembly_folder + 'gdtb_genomes/'
				print(command0 + '\n')
				os.system(command0)
			
				command1 = 'gunzip ' + assembly_folder + 'gdtb_genomes/' + genome_file_name 
				print(command1 + '\n')
				os.system(command1)


ref = '/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/annotations/ref/Campylobacter_jejuni_NCTC11168.gb'
ref2 = '/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/annotations/ref/Campylobacter_jejuni_ATCC35925.gb'

assembled_files = [f for f in os.listdir(assembly_folder) if f.endswith('.fasta')]

for assembly_file in assembled_files:
	assembly_loc = assembly_folder + assembly_file 
	prefix = assembly_file.split('.f')[0]
	assembly_annotation_output = prokka_folder + prefix
	if not os.path.exists(assembly_annotation_output):
		if 'bin' not in prefix:
			command3a = 'prokka --cpus ' + str(threads) + ' --outdir ' + assembly_annotation_output + ' --prefix ' + prefix  + " " + assembly_loc 
		else:
			command3a = 'prokka --cpus ' + str(threads) + ' --metagenome --outdir ' + assembly_annotation_output + ' --prefix ' + prefix  + " " + assembly_loc 
		print(command3a + '\n')
		os.system(command3a)

	assembly_annotation_ref_output = prokka_ref_folder + prefix
	if not os.path.exists(assembly_annotation_ref_output):
		if 'bin' not in prefix:
			command3b = 'prokka --proteins ' + ref + ' --cpus ' + str(threads) + ' --outdir ' + assembly_annotation_ref_output + ' --prefix ' + prefix  + " " + assembly_loc
		else:
			command3b = 'prokka --proteins ' + ref + ' --cpus ' + str(threads) + ' --metagenome --outdir ' + assembly_annotation_ref_output + ' --prefix ' + prefix  + " " + assembly_loc
		print(command3b + '\n')
		os.system(command3b)

	assembly_annotation_ref2_output = prokka_ref2_folder + prefix
	if not os.path.exists(assembly_annotation_ref2_output):
		if 'bin' not in prefix:
			command3c = 'prokka --proteins ' + ref2 + ' --cpus ' + str(threads) + ' --outdir ' + assembly_annotation_ref2_output + ' --prefix ' + prefix  + " " + assembly_loc
		else:
			command3c = 'prokka --proteins ' + ref2 + ' --cpus ' + str(threads) + ' --metagenome --outdir ' + assembly_annotation_ref2_output + ' --prefix ' + prefix  + " " + assembly_loc
		print(command3c + '\n')
		os.system(command3c)


selected_stinkeris_genomes = ['131980_min1000.fasta', 'JMF_2102_8_0018.fasta', 'JMF18_min1000.fasta']

for assembly_file in selected_stinkeris_genomes:
        assembly_loc = assembly_folder + assembly_file 
        prefix = assembly_file.split('.f')[0]
        assembly_annotation_output = prokka_folder + prefix + '_meta'
        prefix = prefix + '_meta'
        if not os.path.exists(assembly_annotation_output):
                command4a = 'prokka --cpus ' + str(threads) + ' --metagenome --outdir ' + assembly_annotation_output + ' --prefix ' + prefix  + " " + assembly_loc
                print(command4a + '\n')
                os.system(command4a)


gtdb_files = [f for f in os.listdir(assembly_folder + 'gdtb_genomes/') if f.endswith('.fna')]


for gtdb_file in gtdb_files:
	gtdb_loc = assembly_folder + 'gdtb_genomes/' + gtdb_file 
	prefix = gtdb_file.split('.f')[0]
	gtdb_annotation_output = prokka_folder + prefix
	if not os.path.exists(gtdb_annotation_output):
		command5a = 'prokka --cpus ' + str(threads) + ' --outdir ' + gtdb_annotation_output + ' --prefix ' + prefix  + " " + gtdb_loc
		print(command5a + '\n')
		os.system(command5a)

	gtdb_annotation_ref_output = prokka_ref_folder + prefix
	if not os.path.exists(gtdb_annotation_ref_output):
		command5b = 'prokka --proteins ' + ref + ' --cpus ' + str(threads) + ' --outdir ' + gtdb_annotation_ref_output + ' --prefix ' + prefix  + " " + gtdb_loc
		print(command5b + '\n')
		os.system(command5b)

	gtdb_annotation_ref2_output = prokka_ref2_folder + prefix
	if not os.path.exists(gtdb_annotation_ref2_output):
		command5c = 'prokka --proteins ' + ref2 + ' --cpus ' + str(threads) + ' --outdir ' + gtdb_annotation_ref2_output + ' --prefix ' + prefix  + " " + gtdb_loc
		print(command5c + '\n')
		os.system(command5c)


command6 = 'cp ' + prokka_folder + '*/*.faa ' + annotations_folder + 'prots/'
print(command6 + '\n')
os.system(command6)

command7 = 'cp ' + prokka_folder + '*/*.ffn ' + annotations_folder + 'genes/' 
print(command7 + '\n')
os.system(command7)

command8 = 'cp ' + prokka_folder + '*/*.tsv ' + annotations_folder + 'tables/'  
print(command8 + '\n')
os.system(command8)

command9 = 'cp ' + prokka_ref_folder + '*/*.tsv ' + annotations_folder + 'tables_ref/'  
print(command9 + '\n')
os.system(command9)

command10 = 'cp ' + prokka_ref2_folder + '*/*.tsv ' + annotations_folder + 'tables_ref2/'  
print(command10 + '\n')
os.system(command10)

