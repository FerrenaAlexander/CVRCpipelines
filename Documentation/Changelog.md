# Version changelog

For all changes, please update changelog and use Year-Month-Day

## 0.1.5
2025.02.19
- add `quickpdf()` function
- add my celltype clustering function `celltype_cluster_map()` (developed for a prior project still unpublished)
- update docs for `save_plotlist_suffix_listnames()` function
- suggest a fix for installing aPEAR from cran archive

## 0.1.4
2024.11.26
- add defaults for `verbose` (F) and `workernum` (1) in the `gsea_apear_pipeline()` function
- update docs of the `gsea_apear_pipeline()` function
- add tutorial vignette for GSEA pipeline

## 0.1.3
2024.11.18
- added option "filename_prefix" to GSEA functions, an extra string appended to all saved file names to make excel integration easier


## 0.1.2
2024.10.23
- added options "title_separator" and "verbose" to `save_plotlist_suffix_listnames()` function

## 0.1.1
2024.10.09
- add pval threshold option
- add pathway categories as an optional argument in gsea script
- general documentation update
- fix general wrapper `...` argument passing to just gsea script

## 0.1
2024.10.09
- initial commit, get gsea pipeline running reproducibly across runs