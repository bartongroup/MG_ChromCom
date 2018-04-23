# Chromatin compaction: a simple 3-state model

This repository contains software used in [Eykelenboom et al. 2018](http://biorxiv.org/cgi/content/short/305391v1).

To run the software and create the markdown document you need to do the following. First, create a conda environment with all required software.

```
conda create --name chromcom --file conda-spec.txt
source activate chromcom
```

Then, run `snakemake`:

```
snakemake -c "qsub -V -cwd -o snakelog -e snakelog -pe smp {threads}" --jobs=50
```

When all the calculations are finished (which make take several days on a cluster) you can "knit" the r-markdown document `doc/mod3.Rmd`. First, you need to modify the `projectDir` in `R/setup.R`, to point to where your results are. Once this is done open `doc/mod3.Rmd` in Rstudio and knit it.
