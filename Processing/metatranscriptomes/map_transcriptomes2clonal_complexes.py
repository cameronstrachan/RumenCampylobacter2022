import os, sys
import re

# map transcriptomes to selected genomes competitively

threads = 80
genomes_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/assemblies/'
transcriptome_reads_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/transcriptomes/reads/'
intermediate_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/mapping/bwa_intermediate_files/'
assembly_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/mapping/assemblies/'
tables_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/mapping/count_tables/'

ref_genomes =  ['131980_spades_pomoxis_polished_min2000', 'JMF18_spades_pomoxis_polished_min2000', 'JMF_2102_8_0001', '131975_min1000', '131974_min1000', '131978_min1000', 'JMF_2102_8_0010', 'JMF_2102_8_0015_filtered', 'JMF_2102_8_0008', 'JMF_2102_8_0006', 'JMF_2102_8_0004_filtered', '131977_min1000']
string_2_cat = ' '

for ref_genome in ref_genomes:
        string_2_cat = string_2_cat + genomes_folder + ref_genome + '.fasta' +  ' '

os.system('cat' + string_2_cat + ' > ' + assembly_folder +  'clonal_contigs.fasta')

ref = assembly_folder +  'clonal_contigs.fasta'
ref_check = assembly_folder +  'bwa_index/clonal_contigs.fasta.bwt'

if not os.path.exists(ref_check):
        command0 = 'bwa index ' + ref
        print(command0 + '\n')
        os.system(command0)
        os.system('mv ' + assembly_folder + 'clonal_contigs.fasta.* ' + assembly_folder + 'bwa_index/')

fastq_files = [f for f in os.listdir(transcriptome_reads_folder) if f.endswith(('_R1.fastq', '_R2.fastq'))]
sample_names = list(set([re.sub('_R\d.fastq', '', s) for s in fastq_files]))

database = assembly_folder  + 'bwa_index/clonal_contigs.fasta'

for sample in sample_names:
        
        forward_read_trim = transcriptome_reads_folder + sample + '_R1.fastq'

        sample_clean = re.sub('_clean', '', sample)
        output_file_loc = transcriptome_reads_folder + sample_clean + '.R1.transcriptome.clones.sam'
        output_file_loc_sam_sort = intermediate_folder + sample_clean + '.R1.transcriptome.clones.sort.sam'

        if not os.path.exists(output_file_loc_sam_sort):
                command1 = 'bwa mem -t ' + str(threads) + ' ' + database + ' ' + forward_read_trim  + ' > ' + output_file_loc
                print(command1 + '\n')
                os.system(command1)
                
                command2 = 'samtools sort -O sam -@ ' + str(threads) + ' ' + output_file_loc + ' -o ' + output_file_loc_sam_sort
                print(command2 + '\n')
                os.system(command2)

prokka_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/annotations/prokka_ref/'

string_2_cat2 = ' '

for ref_genome in ref_genomes:
        string_2_cat2 = string_2_cat2 + prokka_folder + ref_genome + '/' + ref_genome + '.gff' + ' '

os.system('cat' + string_2_cat2  + ' > ' + assembly_folder + 'gff/gff_clones_concat.gff')

input_file = open(assembly_folder + 'gff/gff_clones_concat.gff', 'r')
output_file = open(assembly_folder + 'gff/gff_clones_concat_clean.gff', 'w')

for line in input_file:
        if line.split('_')[0] == 'NODE':
                output_file.write(line)

output_file.close() 


sam_files = [f for f in os.listdir(intermediate_folder) if f.endswith('.R1.transcriptome.clones.sort.sam')]

for sam_file in sam_files:
        output_table = tables_folder + sam_file.split('.sa')[0] + '.txt'
        
        if not os.path.exists(output_table):
                command3 = 'htseq-count -s no -t CDS -i ID --additional-attr=gene --additional-attr=product ' + intermediate_folder + sam_file + ' ' + assembly_folder + 'gff/gff_clones_concat_clean.gff' + ' > ' + output_table 
                print(command3 + '\n')
                os.system(command3)


