qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path /data/Unit_LMM/selberherr-group/strachan/campy/amplicons/raw/neubauer2018/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path /data/Unit_LMM/selberherr-group/strachan/campy/amplicons/qiime_int_files/demux-single-end-neubauer2018.qza
  
qiime dada2 denoise-single \
  --i-demultiplexed-seqs /data/Unit_LMM/selberherr-group/strachan/campy/amplicons/qiime_int_files/demux-single-end-neubauer2018.qza \
  --p-trim-left 25 \
  --p-trunc-len 225 \
  --o-representative-sequences /data/Unit_LMM/selberherr-group/strachan/campy/amplicons/qiime_int_files/rep-seqs-dada2-neubauer2018.qza \
  --o-table /data/Unit_LMM/selberherr-group/strachan/campy/amplicons/qiime_int_files/table-dada2-neubauer2018.qza \
  --o-denoising-stats /data/Unit_LMM/selberherr-group/strachan/campy/amplicons/qiime_int_files/stats-dada2-neubauer2018.qza \
  --p-n-threads 100

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path /data/Unit_LMM/selberherr-group/strachan/campy/amplicons/raw/wetzels2017/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path /data/Unit_LMM/selberherr-group/strachan/campy/amplicons/qiime_int_files/demux-single-end-wetzels2017.qza
  
qiime dada2 denoise-single \
  --i-demultiplexed-seqs /data/Unit_LMM/selberherr-group/strachan/campy/amplicons/qiime_int_files/demux-single-end-wetzels2017.qza \
  --p-trim-left 42 \
  --p-trunc-len 242 \
  --o-representative-sequences /data/Unit_LMM/selberherr-group/strachan/campy/amplicons/qiime_int_files/rep-seqs-dada2-wetzels2017.qza \
  --o-table /data/Unit_LMM/selberherr-group/strachan/campy/amplicons/qiime_int_files/table-dada2-wetzels2017.qza \
  --o-denoising-stats /data/Unit_LMM/selberherr-group/strachan/campy/amplicons/qiime_int_files/stats-dada2-wetzels2017.qza \
  --p-n-threads 100


qiime feature-table merge \
  --i-tables /data/Unit_LMM/selberherr-group/strachan/campy/amplicons/qiime_int_files/table-dada2-neubauer2018.qza \
  --i-tables /data/Unit_LMM/selberherr-group/strachan/campy/amplicons/qiime_int_files/table-dada2-wetzels2017.qza \
  --o-merged-table /data/Unit_LMM/selberherr-group/strachan/campy/amplicons/qiime_int_files/table-merged.qza
  
qiime feature-table merge-seqs \
  --i-data /data/Unit_LMM/selberherr-group/strachan/campy/amplicons/qiime_int_files/rep-seqs-dada2-neubauer2018.qza \
  --i-data /data/Unit_LMM/selberherr-group/strachan/campy/amplicons/qiime_int_files/rep-seqs-dada2-wetzels2017.qza \
  --o-merged-data /data/Unit_LMM/selberherr-group/strachan/campy/amplicons/qiime_int_files/rep-seqs-merged.qza

qiime vsearch cluster-features-de-novo \
  --i-table /data/Unit_LMM/selberherr-group/strachan/campy/amplicons/qiime_int_files/table-merged.qza \
  --i-sequences /data/Unit_LMM/selberherr-group/strachan/campy/amplicons/qiime_int_files/rep-seqs-merged.qza \
  --p-perc-identity 1 \
  --o-clustered-table /data/Unit_LMM/selberherr-group/strachan/campy/amplicons/qiime_int_files/table-merged-clustered.qza \
  --o-clustered-sequences /data/Unit_LMM/selberherr-group/strachan/campy/amplicons/qiime_int_files/rep-seqs-merged-clustered.qza

qiime tools export \
  --input-path qiime_int_files/table-merged-clustered.qza \
  --output-path output
  
biom convert -i output/feature-table.biom \
-o output/asv-counts-merged.txt --to-tsv

rm output/feature-table.biom

qiime tools export \
  --input-path qiime_int_files/rep-seqs-merged-clustered.qza \
  --output-path output
  
mv output/dna-sequences.fasta output/asv-seqs-merged.fasta



