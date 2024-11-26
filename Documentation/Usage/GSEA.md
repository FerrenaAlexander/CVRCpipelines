# GSEA Vignette

This package contains a pipeline for easily running Gene Set Enrichment Analysis (GSEA). The pipeline is a wrapper around the [FGSEA](https://www.biorxiv.org/content/10.1101/060012v3) package using [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8) results as the input.

Additionally, network analysis of the enriched pathways is also provided via a wrapper around the [aPEAR](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8) package.

The pipeline requires two inputs:

1.  **Differential expression summary statistics** for genes formatted in a style similar to DESeq2 "results", including gene symbols, Log2FoldChanges, and pvalues. GSEA should be used with all detected genes, not just thresholded / significant DEGs.

2.  **Pathway / gene sets, ie from a database like MSIGDB** - a function to easily prepare this is also provided.

<br>
<br>

## Basic Pipeline overview

GSEA is run to find enrichment in your DESeq2 result of pathways / genesets from the MSIGDB.

GSEA first ranks the genes, and then checks where the database-derived gene sets/pathways fall along the ranking. The ranking is thus important. This pipeline ranks genes as follows: **-log10(pvalue) \* sign(Log2FoldChange)**. The ranks are then sorted high to low. The ends of this list correspond the most significantly differentially expressed genes. GSEA works by randomizing the order, and comparing the random "expected" distribution of the pathway genes to the true "observed" ordered distribution. Pathways whose genes are distributed at one end of the list are considered *significantly enriched* and *significantly differentially expressed*.

The Gene sets are drawn from a database called the [Molecular Signatures DataBase (MSIGDB)](https://www.gsea-msigdb.org/gsea/msigdb) which is organized into various "categories" of pathways / genesets, including:

-   the **Hallmarks** pathways (1);

-   Gene Ontology (GO) which includes: the Biological Process **(GOBP)** (2), Molecular Function **(GOMF)** (3), and Cellular Component **(GOCC)** (4) subsets;

-   Kyoto Encyclopedia of Genes and Genomes **(KEGG)** (5);

-   **Reactome** (6);

-   and two gene set databases of Transcription Factor Targets (TFT); including the CHIP-seq based **"TFT GTRD"** (7), and the motif / predicted binding-based **"TFT Legacy"** (8) database.

All eight (8) of these are run by default in the pipeline.

After running GSEA for all eight categories, aPEAR is used to detected enriched networks of "functional modules" after. Many pathways share a large intersect of gene members, and many pathways are truly involved in biological cross-talk. This analysis can identify broad molecular effects driven by specific sets of genes and reveal information about their upstream transcriptional regulation.

Various plotting visualizations are also provided, including dotplots for pathway-level enrichments and a network plot showing the aPEAR clusters.

<br>
<br>

## Quick start - whole pipeline wrapper

The pipeline is composed of separate modules. A single pipeline function to launch all the modules at once is also provided for convenience.

The pipeline is parallelized over "gene set categories", like the eight MSIGDB categories mentioned above, so multi-threading with up to eight CPU cores is supported. A typical run on a DESeq2 input may take about 5 minutes with eight cores and less than 5GB.

### Checking function parameters

Please note, the documentation of each function is extensive, please check the function details like so in R:

`?CVRCpipelines::gsea_apear_pipeline()`

<br>

### Prepare Inputs

##### 1. DESeq2 result input format

Prepare DESeq2 result data.frame. All that is needed is something like this below, the column names are important:

```         
 gene_symbol  log2FoldChange   pvalue
        Apoe       6.293962 1.58e-55
        Ctss       8.995505 8.91e-44
       Abcg1       6.037940 3.76e-39
```

GSEA should be used with all detected genes, not just thresholded / significant DEGs.

<br>


##### 2. Prepare the pathways

Run the following command to download the MSIGDB database and parse it.

Considerations:

-   This function prepares and saves a large data.frame of the MSIGDB pathways. It can be used repeatedly for multiple projects, but the database sometimes gets updated over time. Think of it a bit like a reference genome.

-   Run `msigdbr::msigdbr_species()` in R, for a list of available species.

```         
pathways <- CVRCpipelines::preppathways_pathwayanalysis_crosscondition_module(
  species = 'Mus musculus', # replace with your species
  outdir_int = 'path/to/directory' # replace with your desired path
  )
```

This function will save a .rds file in the outdir_int, which you can read into memory via:

```         
pathways <- readRDS("path/to/msigdb_file.rds")
```

<br>

### Run the whole pipeline at once

To run the whole pipeline, you can use this:

```         
library(CVRCpipelines)

set.seed(54321)

gsea_result_list <- gsea_apear_pipeline(
    deseqres,              # your DESeq2 res data.frame
    pathways,              # data.frame of pathways
    outdir,                # directory path to save outputs to
    verbose=T,             # fun to look at
    workernum=1,           # up to 8 is supported
    filename_prefix = ''   # a string to append to saved output csv files
)
```

Note the default threshold for enriched pathways is **adjusted p value < 0.1**. This can be set via the `pathway_padj_thres` parameter.


The function will return a complex named list with all the main results. They are accessible via the the dollar sign operator, such as `gsea_result_list$apear_plot`. This list is not as important as the saved outputs, detailled below. But the list elements are as follows:

-   gseareslist - a list of of dataframes containing results for each pathway category in table format

-   dp_l - a list of ggplot dotplots, showing the top enriched pathways from each category

-   apearoutlist - a list containing the aPEAR clusters and similarity matrix

-   apear_plot - a ggplot object, the aPEAR network plot

-   clust_dp_l - a list of ggplot dotplots, showing the top enriched pathways from each aPEAR cluster

More important than the output list, however, is the actual saved output.

<br>

### Check the saved outputs

The outdir will look something like this:

```         
outdir/
├── aPEAR/
├── GSEA_Category_Dotplots/
├── GSEA_tables/
├── logfile.txt
└── SignificantPathwaysTable_WithClusters.csv
```

The outputs consist of:

-   aPEAR - a subfolder containing the aPEAR clustering results, including the graph plot, and dotplots with the strongest enriched pathways in each cluster

-   GSEA_Category_Dotplots - a subfolder containing dotplots with the strongest enriched pathways from each pathway category

-   GSEA_tables - a subfolder containing .csv files (openable with excel) with the strongest enriched pathways from each pathway category. Note that excel does not allow you to open multiple files with the same name. If you are running the pipeline on multiple DE tables, this can be a bit annoying, so the `filename_prefix` argument allows you to label each file with an optional string suffix, making them unique.

-   logfile.txt - contains some basic info about software versions, number of enriched pathways, etc

-   SignificantPathwaysTable_WithClusters.csv - a final table containing all the significant pathways from all categories, with an additonal column for the aPEAR clusters. This is the most useful table to look at.

#### The way I check the outputs is typically like this:

1. Look at the file **logfile.txt** (for example in terminal with `less` or `cat`) - this will provide some general information about the pipeline, number of enriched pathways / aPEAR clusters detected, etc.

2.  Open the file **aPEAR/aPEAR_GraphPlots/aPEAR_GraphPlot.date-time.pdf** - This will give you a broad overview of the differentially enriched pathways.

3.  Open the file **GSEA_Category_Dotplots/GSEA_Dotplot_HALLMARK.pdf** - The Hallmarks pathways are high quality annotations, manually curated to optimize uniqueness, so it is always useful to check these for a broad overview of the DE pathways.

4. Open the dotplots in **aPEAR/aPEAR_ClusterPathways_Dotplots** - each dotplot will show the strongest enriched pathways from each aPEAR cluster.

5.  Open the table **SignificantPathwaysTable_WithClusters.csv** - this table has all the main outputs.



Some important terminology:

- "Normalized Enrichment Score" (NES). This is the main "effect size" of enrichment used for interpretation. It is a directional value usually between -3 to 3. The sign of the NES will match the sign of the DESeq2 Log2FC in terms of your experimental design.

- Leading Edge: the genes that drive the enrichment of the pathway.


<br>

## Common errors / issues

#### aPEAR cluster is NA for many pathways

This is expected. Only the top 500 (250 for each signed "direction" of NES) are inputted, by default, otherwise aPEAR often fails. The rest are given NA. Additionally, some pathways may have too few genes to check similarity. Furthermore, some pathways are just too "unique" in terms of their gene set to cluster with any other pathway.



#### aPEAR error: No clusters detected

Sometimes, you may get an error related to aPEAR that says something along the lines of "No clusters found". This can be because few significant pathways were detected, or the detected pathways had little overlap in the genes making up the pathways.

In this case, you can try just running with a try statement. It will run the basic GSEA, but will not save any aPEAR output. Note the **SignificantPathwaysTable_WithClusters.csv** will also not be saved - you will need to analyze each category's output individually.

```         
try(
    gsea_result_list <- gsea_apear_pipeline(
        deseqres,              # your DESeq2 res data.frame
        pathways,              # data.frame of pathways
        outdir,                # directory path to save outputs to
        verbose=T,             # fun to look at
        workernum=1,           # up to 8 is supported
        filename_prefix = ''   # a string to append to saved output csv files
    )
)
```

#### aPEAR visualization looks wierd

The aPEAR graph plot sometimes looks weird. The plot generation is randomized, so you can try running it repeatedly with the code below, to maybe get a better generation. The date-time of the plots is appended to the name by default, so it will not overwrite.

```         
apear_obj_path <- 'outdir/aPEAR/robject_apearlist.rds'
apear_obj <- readRDS(apear_obj_path)

plot_apear_clusters(apear_obj, outdir = 'outdir')
```

The plot and labels may still look weird. Sometimes, big clusters have their labels thrown half-way across the plot. I am not sure why this happens. You can try tweaking the parameters, turning off `repelLabels`, or just editing them with a PDF editor.

For other aPEAR questions, reach out to the authors of aPEAR.

<br>
<br>

## Pipeline modules

The pipeline is organized into the following modules, each with their own extensive documentation and flexible parameters with sensible defaults, run in order as below:

0.  `preppathways_pathwayanalysis_crosscondition_module()` - Download the MSIGDB database and do some pruning and formatting. This function is a wrapper around the [msigdbr](https://cran.r-project.org/web/packages/msigdbr/index.html) package. This function is pulled directly from a package I wrote for single-cell RNAseq analysis called [scDAPP](https://github.com/bioinfoDZ/scDAPP)

1.  `deseq_to_gsea()` - given DESeq2 results and pathways, find the enriched pathways.

2.  `gseares_dotplot_listwrap()` - make dotplots of the strongest enriched pathways from each of the MSIGDB categories

3.  `gsea_to_aPEAR_clusters()` - using the enriched pathways, run aPEAR to check the gene overlap between pathways (Jaccard distance), and find clusters of enriched pathways with similar genes

4.  `plot_apear_clusters()` - plot the aPEAR network graph plot

5.  `plot_apear_cluster_pathway_dotplots()` - make dotplots of the strongest enriched pathways from each of the aPEAR clusters
