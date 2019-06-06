# CRSSANT: Cross-linked RNA Secondary Structure Analysis using Network Techniques
CRSSANT is an analysis pipeline for sequencing data produced using the PARIS assay described in [Lu et al., Cell 2016](https://www.sciencedirect.com/science/article/pii/S0092867416304226). CRSSANT automates the process of grouping sequencing reads into duplex groups (DGs), extracts stem groups (SGs) from the duplex groups, and reports the potential stem structures that result from folding SGs using ViennaRNA's RNAfold software.

Briefly, the CRSSANT pipeline is as follows: First, CRSSANT uses network analysis methods to cluster PARIS sequencing reads into DGs. After DGs have been found, a percentile rule is used to filter out reads that have one or both arms that are longer than the some percentile arm length, aggregated over all reads in a gene, and over both arms in a read. The remaining reads are referred to as stem groups, or SGs. SGs are piped into ViennaRNA's RNAfold software, and each SG is checked to see if it forms a valid stem structure. Predicted stem structures are reported as lists of base pairs.

CRSSANT is written in Python and available as source code that you can download and run on yuor own machine

For more about the CRSSANT pipeline, please see the [bioRxiv preprint by Fischer-Hwang et al.](LINKLINKLINK).

## Contents
* [Download](https://github.com/ihwang/CRSSANT#download)
* [Run](https://github.com/ihwang/CRSSANT#run)
    * [Specify pipeline parameters](https://github.com/ihwang/CRSSANT#specify-pipeline-parameters)
        * [Output folder](https://github.com/ihwang/CRSSANT#output-folder)
        * [Chimeric reads file](https://github.com/ihwang/CRSSANT#chimeric-reads-file)
        * [Genes for analysis](https://github.com/ihwang/CRSSANT#Genes-for-analysis)
        * [Custering method](https://github.com/ihwang/CRSSANT#clustering-parameters)
        * [Number of threads](https://github.com/ihwang/CRSSANT#Number-of-threads)        
* [Outputs](https://github.com/ihwang/CRSSANT#outputs)
* [Misc](https://github.com/ihwang/CRSSANT#misc)
    * [Creating a reference.bed file](https://github.com/ihwang/CRSSANT#creating-a-referencebed-file)
    * [Help](https://github.com/ihwang/CRSSANT#help)
* [Test](https://github.com/ihwang/CRSSANT#test)

## Download
Navigate to the latest [release](https://github.com/ihwang/CRSSANT/releases), right click on the source code, and save it to a known path/location, e.g. `CRSSANT_path`. You will need Python version 3.6+ and the following Python packages. We recommend downloading the latest versions of these packages using the Ananconda/Bioconda package manager (follow instructions in links in parentheses):
* [ViennaRNA v2.4.7+](https://www.tbi.univie.ac.at/RNA/) ([Anaconda Cloud link](https://anaconda.org/bioconda/viennarna))
* [ushuffle v1.2.2+](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-192) ([Anaconda Cloud link](https://anaconda.org/bioconda/ushuffle))
* [NetworkX v2.1+](https://networkx.github.io/) ([Anaconda Cloud link](https://anaconda.org/anaconda/networkx))
* [NumPy](http://www.numpy.org/) ([Anaconda Cloud link](https://anaconda.org/anaconda/numpy))
* [scikit-learn](http://scikit-learn.org/stable/) ([Anaconda Cloud link](https://anaconda.org/anaconda/scikit-learn))
* [SciPy](https://www.scipy.org/) ([Anaconda Cloud link](https://anaconda.org/anaconda/scipy))

## Run
To run the CRSSANT pipeline using the Python source code, open a command line interface on your machine and make sure to prepend all pipeline commands with a call to Python on your platform, e.g.:
```
python CRSSANT_path/CRSSANT reads.sam reference.fa reference.bed
```
where files `reads.sam`, `reference.fa`, and `reference.bed` include paths to the reads, reference sequence and reference gene files, respectively. See subsection [Specifying pipeline parameters](https://github.com/ihwang/CRSSANT#specifying-pipeline-parameters) for how to specify non-default pipeline parameters using optional flags.

## Specify pipeline parameters
The input to CRSSANT is assumed to be a SAM file of aligned sequencing reads produced by the PARIS assay. The reads are further assumed to be mapped to the same genomic region (e.g. chromosome or mini genome). The SAM file may contain reads from different genes, but all genes must reside in only a single genomic region.

By default, the CRSSANT pipeline analyzes all reads in a SAM file. The pipeline uses the spectral clustering method to cluster reads into DGs with overlap threshold parameter of `t_o=0.5` and eigenratio threshold of `t_eig=5`, and uses 8 threads for parallel processing. The following parameters allow the user to run CRSSANT using options different from the default ones.

### Output folder
The CRSSANT pipeline automatically writes all results to the same path where the reads file is found, but an output path may be specified with the `-o` flag:
```
CRSSANT_path/CRSSANT reads.sam reference.fa reference.bed -o output
```

### Chimeric reads file
The CRSSANT pipeline assumes that all reads in the input SAM file map to a single reference genome location. Reads that map to multiple locations within the same genomic region, or chimeric reads, can be specified using the chimeric reads file flag, `-c`. When using the `-c` flag, the `reads.sam` file is assumed to contain only normally-aligned reads that were mapped to a single reference genome location, as indicated with an `XG:i:0` tag. Using the `-c` flag will create and save a new SAM file with filename ending in `_chimeric.sam` in the `reads.sam` file path. All parsed chimeric reads added to this new file will contain a new chiastic group (XG) field designation, an `XG:i:1` tag. The CRSSANT analysis pipeline is then run on the new file.

To specify a chimeric reads file, run with flag `-c`:
```
CRSSANT_path/CRSSANT reads.sam reference.fa reference.bed -c chimeric.sam -o output
```

### Genes for analysis
By default, CRSSANT analyzes all possible pairs of genes present in the SAM file. The user may also specify a particular pair of genes for analysis using the gene flag `-g`, e.g. `-g gene1,gene2` indicates that the CRSSANT pipeline should analyze only reads whose left arms map to gene1 and whose right arms map to gene2.

To run CRSSANT on a particular gene pair of interest, run with flag `-g`:
```
CRSSANT_path/CRSSANT reads.sam reference.fa reference.bed -g g1,g2 -o output
```

### Clustering method
The default spectral clustering method may be operated with different overlap threshold and eigenratio threshold parameters by specifying one or both with the flags `t_o` and `t_eig`, respectively. `t_o` may be any float between 0 and 1, and `t_eig` may be any positive number. Increasing `t_o` tends to result in more DGs that contain fewer reads, and increasing `t_eig` tends to result in fewer DGs containing more reads. For example, the following command runs CRSSANT with spectral clustering using overlap threshold 0.6 and eigenratio threshold 8:
```
CRSSANT_path/CRSSANT reads.sam reference.fa reference.bed -t_o 0.7 -t_eig 8 -o output
```

The user may also specify the cliques-finding method for clustering DGs by specifying the clustering flag `cluster` with the `cliques` option, e.g. `-cluster cliques`. If the cliques-finding method is specified, `t_o` may also be specified, and again may be any float between and 1. The eigenvalue threshold `t_eig` is not used in the cliques-finding method. By default, for the cliques-finding method the overlap threshold is set to 0.1. For example, the following command runs CRSSANT on reads whose arms both map to gene1, and performs DG clustering with the cliques-finding method using overlap threshold 0.3:

```
CRSSANT_path/CRSSANT reads.sam reference.fa reference.bed -g gene1,gene1 -cluster cliques -t_o 0.3 -o output
```

For more on these parameters, see the bioRxiv preprint referenced at the top of this README.

### Number of threads
CRSSANT runs using a default of 8 threads in parallel. The user may specify a different number of threads with the `-n` flag.

## Outputs
CRSSANT produces up to five output files with the following extensions added to the reads (or combined normal reads and chimeric reads) file(s) name. After the DG clustering step, CRSSANT verifies that the DGs do not contain any non-overlapping reads, i.e. any reads where the start position of its left arm is greater than or equal to the stop position of the right arm of any other read in the DG. If the DGs do not contain any non-overlapping reads, then the following output files ending in the following are written:

1. `_CRSSANT.sam`: SAM file containing reads that were successfully assigned to DGs, plus DG and non-overlapping group (NG) annotations
2. `_CRSSANT_dg.bed`: BED file listing all duplex groups. The file header contains the following columns:
```
region    DG start    DG stop    Group_ID_coverage    # reads in DG    -   DG start    DG start    color    2    DG left arm length,DG right arm length    DG left arm start,DG right arm start
```
where `coverage` is defined as c / sqrt(a\*b) and
* c = number of reads in a given DG
* a = number of reads overlapping the left arm of the DG
* b = number of reads overlapping the right arm of the DG

If DG were successfully clustered, the CRSSANT pipeline performs SG assembly for nine percentile threshold cutoffs: 10th percentile, 20th percentile, ..., up to and including the 90th percentile. After SG assembly for each percentile threshold, CRSSANT checks if there is a non-zero number of SGs that were assembled from the DGs. If there is at least one SG, then output files ending in the following are written:

3. `_CRSSANT_sg_bp.bed`: BED file containing basepairs only for SGs that form valid stem structures. The file header contains the following columns:
```
region    SG bp start    SG bp stop     SG
```
4. `_CRSSANT_sg_arc.bed`: BED file containing arcs for all SGs. Arc start and stop positions are the means of the start and stop indices of left and right SG arms, respectively. The file header contains the following columns:
```
region    SG arc start    SG arc stop     SG
```
5. `_CRSSANT_sg.aux`: auxiliary file containing SG crosslinking and stem length information, and arm statistics for all DGs. If no SG was formed from the DG, or if an SG did not result in a valid structure, some of the following information is replaced with null information. The file header contains the following columns:
```
Group_ID    num_reads     UU_cl,UC_cl,num_basepairs    L_start_min,L_start_max,L_start_std     L_stop_min,L_stop_max,L_stop_std        R_start_min,R_start_max,R_start_std     R_stop_min,R_stop_max,R_stop_std
```
where
* SG Group IDs correspond to DG IDs
* `num_reads` is the number of reads in the SG, i.e. the number of reads in the DG that passed percentile filtering
* `UU_cl` and `UC_cl` are the number of staggered uridine and uridine-cytosine base pairs, respectively, in the SG stem structure
* `stem_length` is the length of SG stem structure
* `_min`, `_max`, `_std` are the minimum arm index, maximum arm index, and standard deviation of all start and stop indices of the SG arms

## Misc

### Creating a `reference.bed` file
CRSSANT assumes that the `reference.bed` file contains minimal gene information in the following 6-column format:
```
genomic region    gene start    gene stop    gene name
```
where columns are separated by tabs. Some specifics on format:
* `genomic region` must match the gene naming format in the `reference.fa` file
* `gene start` and `gene stop` must be integers

### Help
To see specifics on arguments for running CRSSANT, run
```
CRSSANT_path/CRSSANT -h
```
The help information may take ~10 seconds to load.

## Test

You can test CRSSANT using a collection of Homo sapiens ribosomal RNA (rRNA) test data that we have compiled:

1. Download the compressed folder of [test data](https://github.com/ihwang/CRSSANT/tree/master/tests.tar.gz) and decompress using the command `tar -zxvf tests.tar.gz` or by double-clicking on the tar.gz file
2. Specify the path/location where results should be written, e.g. `output`

Run CRSSANT on all rRNA genes in region hs54S:
```
CRSSANT_path/CRSSANT tests/hsrRNA_reads.sam tests/hsrRNA.fa tests/hsrRNA_gene.bed -o output
```
or analyze specific genes, e.g. only reads whose left arms map to gene 5.8S and whose right arms map to gene 28S:
```
CRSSANT_path/CRSSANT tests/hsrRNA_reads.sam tests/hsrRNA.fa tests/hsrRNA_gene.bed -g 5.8S,28S -o output
```
