# PlotTrellis

Build genomic track figures (coverage, Hi-C, splice junctions, transcript structure,
protein domains) in a point-and-click app, without writing code. Every figure can be
downloaded as a standalone R script that reproduces it.

PlotTrellis needs R 4.5 or later (check with `R.version.string`). On older R the install
stops with `Error in install.packages : unknown input format`.

Install from R-universe:

```r
install.packages("BiocManager")
install.packages("PlotTrellis", repos = c("https://aleks-wehi.r-universe.dev",
                                          BiocManager::repositories()))
```

Then follow the guide, which starts with your first figure in the app:
[PlotTrellis vignette](https://aleks-wehi.r-universe.dev/articles/PlotTrellis/PlotTrellis.html).

Found a bug? See [CONTRIBUTING.md](https://github.com/aleks-wehi/PlotTrellis/blob/main/CONTRIBUTING.md).
