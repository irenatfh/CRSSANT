# CRSSANT: Computational RNA Secondary Structure Analysis using Network Techniques

CRSSANT is an analysis pipeline for sequencing data produced using the PARIS assay described in [Lu et al., Cell 2016](https://www.sciencedirect.com/science/article/pii/S0092867416304226). CRSSANT automates the process of grouping sequencing reads into duplex groups (DGs), and tests the DGs to find potential secondary structures.


## Install

CRSSANT is packaged as a Python executable, so no prerequisites are needed. Download the appropriate executable from the [release](https://github.com/ihwang/CRSSANT/releases) page and save it to a known path/location, e.g. `CRSSANT_path`.

## Run

To run CRSSANT, open a command-line interface and run
```
CRSSANT_path/CRSSANT.exe reads.sam reference.fa reference.bed -r regions -g genes output
```
where files `reads.sam`, `reference.fa`, and `reference.bed` include paths to the files reads, reference sequence and reference gene files, respectively, and `output` is the path/location where outputs should be written.

To perform the CRSSANT analysis pipeline on all rRNA regions and genes, run:
```
CRSSANT_path/CRSSANT.exe reads.sam reference.fa reference.bed output
```

### Specifying regions and genes for analysis
Genomic regions and genes of interest may be specified with the region flag `-r` and gene flag `-g`. Note that the regions and genes may be specified independently, and that multiple regions and genes may be specified as comma-separated strings, eg. `-r region1,region2,region3` or `-g gene1,gene2`.

* To perform CRSSANT analysis on a particular region of interest, run with flag `-r`:
```
CRSSANT_path/CRSSANT.exe reads.sam reference.fa reference.bed -r region output
```

* To run CRSSANT on a particular gene of interest, run with flag `-g`:
```
CRSSANT_path/CRSSANT.exe reads.sam reference.fa reference.bed -g gene output
```

* To run CRSSANT on a particular genes and regions of interest, run with both region and gene flags:
```
CRSSANT_path/CRSSANT.exe reads.sam reference.fa reference.bed -r region1,region2 -g gene1,gene2,gene3 output
```
### Outputs
CRSSANT will produce the following 5 output files with these extensions added to the original file name:

1. `_CRSSANT.log`: logfile recording which regions and genes were analyzed
2. `_CRSSANT.sam`: SAM file containing reads that were successfully assigned to DGs, plus DG and non-overlapping group (NG) annotations
3. `_CRSSANT_info.bed`: BED file listing all duplex groups. Column information is
```
region    DG start    DG stop    Group_ID_coverage    # reads in DG    -   DG start    DG start    color    2    DG left arm length,DG right arm length    DG left arm start,DG right arm start
```
where `coverage` is defined as c / sqrt(a\*b) and
* c = number of reads in a given DG
* a = number of reads overlapping the left arm of the DG
* b = number of reads overlapping the right arm of the DG
4. `_CRSSANT.aux`: auxiliary file containing crosslinking and stem length information, and arm statistics for each DG--see file header
5. `_CRSSANT_bp.bed`: BED file containing basepairs for only the DGs that pass structure tests

### Test

You can test CRSSANT using a collection of Homo sapiens ribosomal RNA (rRNA) test data that we have compiled:

1. Download the [folder](https://github.com/ihwang/CRSSANT/tree/master/tests) of test data to a known path/location, e.g. `downloads`
2. Specify the path/location where results should be written, e.g. `results`

Run CRSSANT on all rRNA regions and genes:
```
CRSSANT_path/CRSSANT.exe downloads/hsrRNA_reads.sam hsrRNA.fa hsrRNA_gene.bed output
```

or analyze specific regions and genes, e.g. all genes in region hs12S and only genes 5.8S and 28S in region hs45S:
```
CRSSANT_path/CRSSANT.exe downloads/hsrRNA_reads.sam hsrRNA.fa hsrRNA_gene.bed -r hs12S -g 5.8S,28S output
```
