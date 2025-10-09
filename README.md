# README

## Description

Small shiny app to perform visual validation of CNV calls. It features four different
modes:

1. Normal: the app will go through all CNVs in the provided table
2. Simple Filtering: the CNVs table will be filtered based on the provided
   criteria (vo, GT, length, numsnp).
3. Fixed locus: the app will select putative carriers in the provided
   locus. If any filter is provided, it will also be applied (before selecting).
   **NB**: in this mode the locus is the units of analysis, not the CNVs.
4. Select CNVs in region: the app will allow the user to specify a genomic region
   and will only show CNVs that overlap with this region (above a certain IOU
   threshold). Also in this mode any provided filter will be applied first.


## Inputs

This app takes as input:

- a folder to save results
- a table of CNVs (sample_ID, chr, start, end, GT, vo, CN, length, numsnp)
- a table linking each sample_ID to a tabix-indexed file (sample_ID, file_path_tabix)
- a table of filtered SNP (Name, Chr, Position)

The tabix-indexed files should contain is expected to be a tab-separated file
with the following columns (no header): chr, position, position, LRR, BAF, LRR_adj
(if LRR_adj is not present, the app will use LRR instead).
The input files are in the usual format for my packages and pipelines, described in details
[here](https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.621).   
Example for all input files are provided in `data`.


## Citation

If you use this tool in your work please cite "Accurate and Effective
Detection of Recurrent Copy Number Variants in Large SNP Genotype Datasets",
DOI: [https://doi.org/10.1002/cpz1.621]( https://doi.org/10.1002/cpz1.621).


## Interface

The app interface is divided in two main parts:

- A sidebar on the left with the controls and input for filtering, fixed locus, and select
  region. The sidebar can be hidden by clicking on the small arrow at the top left.
- The main panel on the right shows the CNV line as a table on the top and the
  plot in the middle of the screen. Below the plot there are buttons to
  navigate through CNVs, store the visual validation for the CNV and to refine boundaries.
  The plot has two rows. BAF is on the top, LRR on the bottom. Each dot represent
  a SNP. Semi-transparent rectangles indicate the CNV call(s) in the region, while
  a gray box represents the locus of region in fixed locus or select region mode.
  In all modes the app will plot all CNVs for the sample (before any filtering)
  on the same chromosome as the current CNV. The CNV of interest is highlighted with
  a thin border. In fixed locus and select region modes, the locus/region is also
  highlighted with a gray box and thicker border. The plot can be zoomed in and out
  using the scroll wheel or the `+`/`-` buttons on the top right of the plot.


# Update Boundaries

This new version implements also the ability to refine CNV boundaries.
When the user clicks on the `Refine boundaries` button, the app will ...


## Usage

To use the app just run
`Rscript app.R /path/to/workingdir /path/to/cnvs.txt /path/to/samplee.txt /path/to/snps.txt`.
Input files can be anywhere, all path must be full (starting from `/` with
no symlinks), the results will be saved in `workingdir`.
Every times the app the app moves from one CNV to the next (or previous) the
output table is saved in the provided `workingdir` as `project_name.tsv`. If
`workingdir` the app will attempt to create it. If the also the parent folder does
exist, the app will assume there is a typo in the provided path and will stop.   
**NB**: the app will overwrite any existing `project_name.tsv` file in the provided
`workingdir` without warning.


## Visual inspection codes

The `vo` column is coded as follows:

-  1 : true
-  2 : false
-  3 : unknown/unclear
- \-9 : new
- \-7 : error


## Install

Just clone the repo, `app.R` is all is needed to run the program (assuming the
dependencies are met). 

R dependencies:

- shiny
- data.table
- ggplot2
- bslib
- DT
- plotly

Linux dependencies:

- tabix
- a browser


# OLD versions

An old version of the app (and of this `README`) is available in the `old_version` folder.
This version is simpler and uses a static plot, for this reason it should also
be faster. However, it does not allow to zoom in and out of the plot, and
it does not have the boundary refinement feature. It is also less polished in general.