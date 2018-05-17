# CRSSANT: Computational RNA Secondary Structure Analysis using Network Techniques

CRSSANT is an analysis pipeline for sequencing data produced using the PARIS assay described in [Lu et al., Cell 2016](https://www.sciencedirect.com/science/article/pii/S0092867416304226). CRSSANT automates the process of grouping sequencing reads into duplex groups (DGs), and tests the DGs to find potential secondary structures.


## Getting Started

### Install

CRSSANT is packaged as a Python executable, so no prerequisites are needed. Download the appropriate executable below for your operating system:

### Run

To run CRSSANT, open a command-line interface and enter the following:

```
CRSSANT_path/CRSSANT.exe reads.sam reference.seq reference.bed -r regions -g genes output
```
Where files `reads.sam`, `reference.seq`, and `reference.bed` include paths to the files, and `output` is a path to where outputs should be written.

### Test

You can test CRSSANT using a collection of test data that we have compiled:

1. Download the [folder](https://github.com/ihwang/CRSSANT/CRSSANT/tests) of test data
2. 
```
Give an example
```
