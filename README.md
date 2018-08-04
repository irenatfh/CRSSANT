# CRSSANT: Crosslinked RNA Secondary Structure Analysis using Network Techniques

CRSSANT is an analysis pipeline for sequencing data produced using the PARIS assay described in [Lu et al., Cell 2016](https://www.sciencedirect.com/science/article/pii/S0092867416304226). CRSSANT automates the process of grouping sequencing reads into duplex groups, and reports the potential stem structures that result from folding duplex groups using ViennaRNA's RNAfold software.


## Install

You can install CRSSANT either as a Python executable for Linux machines, or you can directly download all of the source code. Both the executable and source code files are available at the [release](https://github.com/ihwang/CRSSANT/releases).

### Installing the Python executable:
Navigate to the release page, right click on the executable, and save it to a known path/location, e.g. `CRSSANT_path`. You may need to change the software's permissions so that you can actually execute the pipeline via, e.g. `chmod u+x CRSSANT_path/CRSSANT`.

### Downloading the Python source code:
Navigate to the release page, right click on the source code, and save it to the known path/location, e.g. `CRSSANT_path`. You will need Python version 3.6+, as well as the following Python packages:
* [ViennaRNA v2.4.7+](https://www.tbi.univie.ac.at/RNA/)
* [ushuffle v1.2.2+](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-192)
* [NetworkX v2.1+](https://networkx.github.io/)
as well as typical packages like numpy, scikit-learn, and scipy. We recommend downloading the latest versions of these packages using the Ananconda/Bioconda package manager.

## Run

To run the CRSSANT executable, open a command-line interface and run
```
CRSSANT_path/CRSSANT reads.sam reference.fa reference.bed -r regions -g genes output
```
where files `reads.sam`, `reference.fa`, and `reference.bed` include paths to the files reads, reference sequence and reference gene files, respectively, and `output` is the path/location where outputs should be written. See below for how to specify regions and genes using the optional `-r` and `-g` flags. Users may also specify a chimeric reads file using the optional `-c` flag.

To perform the CRSSANT analysis pipeline on all rRNA regions and genes, run:
```
CRSSANT_path/CRSSANT reads.sam reference.fa reference.bed output
```
To run the CRSSANT pipeline using the Python source code, prepend all commands by calling Python on your platform, e.g. `python CRSSANT_path/CRSSANT reads.sam reference.fa reference.bed -r regions -g genes output`.

### Specifying regions and genes for analysis
Genomic regions and genes of interest may be specified with the region flag `-r` and gene flag `-g`. Note that the regions and genes may be specified independently, and that multiple regions and genes may be specified as comma-separated strings, eg. `-r region1,region2,region3` or `-g gene1,gene2`.

* To perform CRSSANT analysis on a particular region of interest, run with flag `-r`:
```
CRSSANT_path/CRSSANT reads.sam reference.fa reference.bed -r region output
```

* To run CRSSANT on a particular gene of interest, run with flag `-g`:
```
CRSSANT_path/CRSSANT reads.sam reference.fa reference.bed -g gene output
```

* To run CRSSANT on a particular genes and regions of interest, run with both region and gene flags:
```
CRSSANT_path/CRSSANT reads.sam reference.fa reference.bed -r region1,region2 -g gene1,gene2,gene3 output
```
### Creating a `reference.bed` file
CRSSANT assumes that the `reference.bed` file contains minimal gene information in the following 6-column format:
```
genomic region    gene start    gene stop    gene name    1000    +
```
where columns are separated by tabs. Some specifics on format:
* `genomic region` must match the gene naming format in the `reference.fa` file
* `gene start` and `gene stop` must be integers
* The last two columns are default values

### Specifying a chimeric reads file
CRSSANT can also parse and include chimeric reads in the analysis pipeline. If you have a chimeric reads file `chimeric.sam`, you can specify it in the command line using the `-c` flag:
```
CRSSANT_path/CRSSANT reads.sam reference.fa reference.bed -c chimeric.sam output
```
When specifying the `-c` flag, the `reads.sam` file is assumed to contain only normally-aligned reads. Using the `-c` flag will create and save a new SAM file with filename ending in `_chimeric.sam` containing an additional chiastic group (XG) field in the `reads.sam` file path. Aligned reads from `reads.sam` are added to the new file and tagged `XG:i:0`, while paired chimeric reads that do not have any reverse complement components are parsed and appended to the file with an `XG:i:1` tag. The CRSSANT analysis pipeline is then run on the new file.

### Outputs
To understand the output files that CRSSANT produces, it is important to understand the analysis that CRSSANT performs. After DGs have been found, a 90th percentile rule is used to filter out reads that have one or both arms that are longer than the 90th percentile arm length, aggregated over all reads in a gene, and over both arms in a read. The remaining reads are referred to as stem groups, or SGs. SGs are piped into ViennaRNA's RNAfold software, and each SG is checked to see if it forms a valid stem structure.

CRSSANT will produce the following 6 output files with these extensions added to the original file name:

1. `_CRSSANT.log`: logfile recording which regions and genes were analyzed
2. `_CRSSANT.sam`: SAM file containing reads that were successfully assigned to DGs, plus DG and non-overlapping group (NG) annotations
3. `_CRSSANT_info.bed`: BED file listing all duplex groups. The file header contains the following columns:
```
region    DG start    DG stop    Group_ID_coverage    # reads in DG    -   DG start    DG start    color    2    DG left arm length,DG right arm length    DG left arm start,DG right arm start
```
where `coverage` is defined as c / sqrt(a\*b) and
* c = number of reads in a given DG
* a = number of reads overlapping the left arm of the DG
* b = number of reads overlapping the right arm of the DG
4. `_sg.aux`: auxiliary file containing crosslinking and stem length information, and arm statistics for SGs. Note that if a DG did not result in an SG due to 90th percentile filtering, or due to the SG not resulting in a valid structure, some of the following information is replaced with null information. The file header contains the following columns:
```
DG_coverage    num_reads     UU_cl,UC_cl,stem_length    L_start_min,L_start_max,L_start_std     L_stop_min,L_stop_max,L_stop_std        R_start_min,R_start_max,R_start_std     R_stop_min,R_stop_max,R_stop_std
```
where
* `num_reads` is the number of reads in the SG, i.e. the number of reads in the DG that passed 90th percentile filtering
* `UU_cl` and `UC_cl` are the number of staggered uridine and uridine-cytosine base pairs, respectively, in the SG structure
* `stem_length` is the length of SG structure
* `_min`, `_max`, `_std` are the minimum arm index, maximum arm index, and standard deviation of all start and stop indices of the SG arms
5. `_sg_arc.bed`: BED file containing arcs for all SGs. Arcs are drawn between The file header contains the following columns:
```
region    SG start    SG stop     SG    1    +    SG start    SG stop     0,0,0
```
6. `_sg_bp.bed`: BED file containing basepairs only for SGs that form valid stem structures. The file header contains the following columns:
```
region    SG start    SG stop     DG    1    +    SG start    SG stop     0,0,0
```

### Help
To see specifics on arguments for running CRSSANT, run
```
CRSSANT_path/CRSSANT -h
```
The help information may take a little bit to load (~10 seconds).

## Test

You can test CRSSANT using a collection of Homo sapiens ribosomal RNA (rRNA) test data that we have compiled:

1. Download the compressed folder of [test data](https://github.com/ihwang/CRSSANT/tree/master/tests.tar.gz) and decompress using the command `tar -zxvf tests.tar.gz` or by double-clicking on the tar.gz file
2. Specify the path/location where results should be written, e.g. `output`

Run CRSSANT on all rRNA regions and genes:
```
CRSSANT_path/CRSSANT tests/hsrRNA_reads.sam tests/hsrRNA.fa tests/hsrRNA_gene.bed output
```

or analyze specific regions and genes, e.g. all genes in region hs12S and only genes 5.8S and 28S in region hs45S:
```
CRSSANT_path/CRSSANT tests/hsrRNA_reads.sam tests/hsrRNA.fa tests/hsrRNA_gene.bed -r hs12S -g 5.8S,28S output
```
