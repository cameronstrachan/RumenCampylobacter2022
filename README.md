# RumenCampylobacter2022
This repository is for the analysis code from the Strachan et al. 2023 paper "Differential carbon utilization enables coexistence of recently speciated Campylobacteraceae in the cow rumen epithelial microbiome"

## Statistical analysis and figure generation
The code for the analysis of the various output tables can be found in the folder 'Figures'. The analysis is organized by the figures that it was used to generate. The input data for these analysis are either found in the output folders within 'Processing' or from the raw PCR data in the folder 'PCRdata'. The metadata files used in the analysis are in 'METAdata'. All this code can be run from this repository to reproduce any of the figures.

## Bioinformatic processing
The code used to process sequence data to various output tables (ex. count tables) is found in the folder 'Processing'. This is where the command line settings for various bioinformatic tools can be found. Most time, commands are run within pipelines that are written in python or bash. The code is organized by the data types used in the processing (ex. metagenomes). Due to the file sizes, the raw input data needs to be downloaded and then the locations and names of the input data would need to be changed.

## Analysis not included in this repository
There are a a subset of analysis and plots (ex. alignments and trees) that were implemented in the software Geneious (https://www.geneious.com/). This is described in the Materials and Methods of the paper.

