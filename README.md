# PanHumanIVT
## Summary
Here we present an IVT dataset for 5 human cell lines; HeLa, HepG2, SH-SY5Y, NTERA, and A549. In this GitHub repository we've included all of our analysis code. Due to the size of the component data we're hosting files such bam files, fast5, fastq, and Eventalign files on [Dropbox](https://www.dropbox.com/sh/8r76vesjz01mkyq/AACweGZTB1c_tPCu_qxtB1Hha?dl=0). Additionally, fastq and fast5 files can be found on NIH NCBI SRA under the BioProject accession [PRJNA947135](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA947135). Many of the notebooks and code will not run without these files. If you would like to reproduce our analysis we recommend downloading the entire directory from dropbox and executing the code from there.

## Biological Replicate Analysis
In FIG1github.R script we run a comparative performance analysis of three independent biological IVT replicates from the HeLa cell line. Generation of the input data is derived from the sam/bam files. The specific terminal commands are included in the comments of the R script. We used [NanoPlot](https://github.com/wdecoster/NanoPlot) v1.40.2 and the --raw option to extract percent identity for each read and all other alignment identity statistics. All plots populated for this analysis are included in the R script. 

## SNV Analysis
In the IVT_SNV_analysis.ipynb we identify SNVs in each of cellular datasets from pysamstats pileups. We only considered SNVs that occured at loci with a minimum of 10 reads for the cell line (pooled if there were multiple sequencing runs) and where the SNV in question accounted for 30% of the total reads at a position. 
We repeated this analysis for SNV occurence of 40, 50, 60, 70, 80, and 95% all while maintaining the minimum read threshold of 10. The SNVs identified at each level were further binned into known variants according to Ensembl, variants occuring in low confidence kmers, and putative IVT variants.
These SNV datasets are available as .tsv files and include the relevant pysamstats information for each SNV.
The tsv files have the following format:
|Column Number| Contents |
| :-: | :-:|
| 1  | Chromosome  |
| 2  | Position  |
| 3  | Variant Nucleotide  |
| 4  | Reference Nucleotide  |
| 5  | A Count  |
| 6  | C Count  |
| 7  | T Count  |
| 8  | G Count  |
| 9  | Deletion Count  |
| 10  | Insertion Count |

Additionally .bed files were created for each of the binned variant datasets for visualization in a genome browser. For a broader overview a pooled bed file was created including all of the bins of variants across all of the cell lines.
## Multiple Cell Line Analysis
The multiple cell line analysis examined the relationships between our IVT datasets of each cell line. Subsampling experiments and similarity measurements were calculated in a series of jupyter notebooks. The following notebooks contain the following analysis:
1. Event_align_plots.ipynb - Creating ionic current distribution for target section of MCM5 22:35424403-35424411

2. gm12878_genes.ipynb - Subsampling experiment counting the number of GM12878 genes captured across cell line combinations (HeLa & A549, HepG2 & NTERA, All cell lines).

3. gm12878_saturation - Subsampling experiment counting the number of alignment hits for the top 10 genes of GM12878 (established via subsampling)

4. Saturation_Curves.ipynb - Subsampling experiment showing the power of 4 cell lines to capture the genes represented by a 5th cell line. Repeated with all 5 of our cell lines as the target cell line.

5. Subsample_Gene_Binary.ipynb - Using gene representation binaries for each cell line and completing a pairwise dice similarity measurement

6. unique_gene_count.ipynb - Subsampling experiment to determine the number of genes in each of our cell lines from up to 1,000,000 reads.

## Software & Hardware Dependencies:
samtools 1.16.1 (using htslib 1.16)

python 3

R 4.1.1

NanoPlot 1.40.2

Jupyterlab 3.4.4

## Installation
All analysis can be run directly from provided scripts, no installation beyond dependencies required.

## Demo 
To run code, please download repository found in the Dropbox URL above and execute scripts and / or notebooks.
Run time for notebooks ranges from minutes to multiple hours for subsampling experiments.

## Expected Output 
Data analysis figures, identified variant tsv files, and bed files for identified variants.
