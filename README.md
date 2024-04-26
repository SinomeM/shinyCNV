# README

## Description

Small shiny app to perform visual validation of CNV calls.
Plots are generated on the fly. CNVs, samples and SNPs table
are in the usual format for my packages and pipelines,
described in details
[here](https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.621).


## Interface

The plot has three row:

- middle row show the BAF pattern for the CNV and 4-6 lengths on each side
- bottom row shows the LRR patter for this same region plus a dotted
  segment to indicate where the CNV call is.
- top row shows the same as the bottom but even more zoomed out

If a specific locus is selected, locus boundaries are marked with tick marks
in all three rows.


## Usage

To use the app just run
`Rscript app.R /path/to/workingdir /path/to/cnvs.txt /path/to/samplee.txt /path/to/snps.txt`.
Input files can be anywhere, all path must be full (starting from `/` with
no symlinks), the results will be saved in `workingdir`.

It is possible to filter CNVs based on previous visual inspection or
on CNV type (deletions or duplications).

The buttons do what you expect them to do ;)

Save will write the results table in the `workingdir` as `project_name`\_vi\_res.txt.

To filter CNVs follow the instructions in on the left and the then click
"Run Filtering!" at the very bottom.

If the app freezes and the console report and error on these lines:

```
Warning: Error in if: missing value where TRUE/FALSE needed
  3: runApp
  2: print.shiny.appobj
  1: <Anonymous>
```

then it means there are no CNVs in the selected locus with the current
filters set. If the SNPs and length filters results in an empty table,
a message will appear in the middle of the plot.



## Visual inspection codes

The `vo` column is coded as follows:

-  1 : true
-  2 : false
-  3 : unknown
-  4 : true but incorrect boundaries
- \-9 : new
- \-7 : error


## Install

Just clone the repo, `app.R` is all is needed to run the program.

R dependencies:

- shiny
- data.table
- ggplot
- CNValidatron

Linux dependencies:

- tabix
- a browser


## Citation

If you use this tool in your work please cite "Accurate and Effective
Detection of Recurrent Copy Number Variants in Large SNP Genotype Datasets",
DOI: [https://doi.org/10.1002/cpz1.621]( https://doi.org/10.1002/cpz1.6210).

