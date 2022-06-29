threads=100

# run classification of bins

gtdbtk classify_wf --genome_dir /data/Unit_LMM/selberherr-group/strachan/campy/metagenomes/assemblies/bins/ --out_dir /data/Unit_LMM/selberherr-group/strachan/campy/metagenomes/gdtb_intermediate_files/bins --extension fa --cpus $threads
cp /data/Unit_LMM/selberherr-group/strachan/campy/metagenomes/gdtb_intermediate_files/bins/gtdbtk.bac120.summary.tsv output/


