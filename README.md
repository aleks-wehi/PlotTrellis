# PlotTrellis

Pay the setup once, not per figure. One text file describes your experiment. Every gene
thereafter is a few clicks or lines of R from a clean, fair figure, perfect for exploring your data.
Share the file, and your lab can make their own.

## What you need

- R 4.5 or later, on Linux or macOS (Windows cannot read bigWig files).
- Processed files: bigWig, Salmon output, junctions, Hi-C, peaks. See
  [Input data](https://aleks-wehi.r-universe.dev/articles/PlotTrellis/inputs.html).
- Mouse or human, with a GENCODE or Ensembl GTF.

## Install

Installing from scratch takes 40 to 50 minutes on Linux:

```r
install.packages("BiocManager")
install.packages("PlotTrellis", repos = c("https://aleks-wehi.r-universe.dev",
                                          BiocManager::repositories()))
```

Use this line, not the shorter one on the R-universe page: that one lacks the Bioconductor
packages PlotTrellis needs.

## Start

```r
BiocManager::install("org.Mm.eg.db", update = FALSE)   # mouse gene names, once
library(PlotTrellis)
run_plottrellis(project_dir = "PlotTrellis-example", samples_txt = pt_example())
```

This opens the bundled mouse example in the app. The
[Quick Start](https://aleks-wehi.r-universe.dev/articles/PlotTrellis/PlotTrellis.html) takes
you from there to your own data, in the app or in R.

Found a bug? See [CONTRIBUTING.md](https://github.com/aleks-wehi/PlotTrellis/blob/main/CONTRIBUTING.md).

---

**AI use.** PlotTrellis was built with AI coding assistants (Claude, Anthropic). Most of the
code, tests and documentation were drafted by the assistant and directed, reviewed and tested
by the author.
