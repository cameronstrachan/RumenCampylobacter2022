import os, sys

threads = 60
max_target_seqs = 5
evalue = 0.05

prots_select_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/annotations/prots_select/'

prots_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/annotations/prots/'
genes_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/annotations/genes/'

prots_db_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/blast_intermediate_files/prots_db/' 
genes_db_folder = '/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/blast_intermediate_files/genes_db/'

output_prots_folder =  '/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/blast_intermediate_files/prots_table_select/'

#output_prots_folder =  '/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/blast_intermediate_files/prots_table/'
output_genes_folder =  '/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/blast_intermediate_files/genes_table/'

prots_files = [f for f in os.listdir(prots_select_folder) if f.endswith('.faa')]


for prots_file in prots_files:
	
	input_db_loc = prots_folder + prots_file
	output_db_loc = prots_db_folder + prots_file
	output_db_loc_check = prots_db_folder + prots_file + '.pdb'
	if not os.path.exists(output_db_loc_check):
		command0 = 'makeblastdb -in ' + input_db_loc + ' -dbtype prot -out ' + output_db_loc
		print(command0 + '\n')
		#os.system(command0)

for query_file in prots_files:
	for db_file in prots_files:
	
		query_prefix = query_file.split('.fa')[0]
		db_prefix = db_file.split('.fa')[0]
	
		query_file_loc = prots_folder + query_file
		db_file_loc = prots_db_folder + db_file
		output_loc = output_prots_folder + query_prefix + '_V_' + db_prefix + '.txt'
		
		if not os.path.exists(output_loc):
			command1 = 'blastp -query ' + query_file_loc + ' -db ' + db_file_loc +  ' -max_target_seqs ' + str(max_target_seqs) + " -evalue " + str(evalue)  + ' -num_threads ' + str(threads) + " -outfmt '6 qseqid sseqid pident sstart send qstart qend evalue bitscore score qlen length sseq'" + ' -out ' + output_loc
			print(command1 + '\n')
			os.system(command1)

    	
genes_files = [f for f in os.listdir(genes_folder) if f.endswith('.ffn')]


for genes_file in genes_files:
        
        input_db_loc = genes_folder + genes_file
        output_db_loc = genes_db_folder + genes_file
        output_db_loc_check = genes_db_folder + genes_file + '.ndb'
        if not os.path.exists(output_db_loc_check):
                command2 = 'makeblastdb -in ' + input_db_loc + ' -dbtype nucl -out ' + output_db_loc
                print(command2 + '\n')
                #os.system(command2)

for query_file in genes_files:
        for db_file in genes_files:
        
                query_prefix = query_file.split('.f')[0]
                db_prefix = db_file.split('.fa')[0]
        
                query_file_loc = genes_folder + query_file
                db_file_loc = genes_db_folder + db_file
                output_loc = output_genes_folder + query_prefix + '_V_' + db_prefix + '.txt'
                
                if not os.path.exists(output_loc):
                        command3 = 'blastn -query ' + query_file_loc + ' -db ' + db_file_loc +  ' -max_target_seqs ' + str(max_target_seqs) + " -evalue " + str(evalue)  + ' -num_threads ' + str(threads) + " -outfmt '6 qseqid sseqid pident sstart send qstart qend evalue bitscore score qlen length sseq'" + ' -out ' + output_loc
                        print(command3 + '\n')
                        #os.system(command3)
