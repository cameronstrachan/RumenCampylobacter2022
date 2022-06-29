import os, sys

threads = 60
max_target_seqs = 5
evalue = 0.05

prots_select_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/annotations/prots_select/'

genes_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/annotations/genes/'
genomes_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/assemblies/'
 
genomes_db_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/blast_intermediate_files/genomes_db/'

output_folder =  '/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/blast_intermediate_files/genes_against_genomes/'

files = [f for f in os.listdir(prots_select_folder) if f.endswith('.faa')]


for prot_file in files:
        genome_file = prot_file.split('.fa')[0] + '.fasta'
        input_db_loc = genomes_folder + genome_file
        output_db_loc = genomes_db_folder + genome_file
        output_db_loc_check = genomes_db_folder + genome_file + '.ndb'
        if not os.path.exists(output_db_loc_check):
                command2 = 'makeblastdb -in ' + input_db_loc + ' -dbtype nucl -out ' + output_db_loc
                print(command2 + '\n')
                os.system(command2)

for file1 in files:
        for file2 in files:
        		
                query_prefix = file1.split('.fa')[0]
                db_prefix = file2.split('.fa')[0]
        
                query_file_loc = genes_folder + query_prefix + '.ffn'
                db_file_loc = genomes_db_folder + db_prefix + '.fasta'
                output_loc = output_folder + query_prefix + '_V_' + db_prefix + '.txt'
                
                if not os.path.exists(output_loc):
                        command3 = 'blastn -query ' + query_file_loc + ' -db ' + db_file_loc +  ' -max_target_seqs ' + str(max_target_seqs) + " -evalue " + str(evalue)  + ' -num_threads ' + str(threads) + " -outfmt '6 qseqid sseqid pident sstart send qstart qend evalue bitscore score qlen length sseq'" + ' -out ' + output_loc
                        print(command3 + '\n')
                        os.system(command3)
