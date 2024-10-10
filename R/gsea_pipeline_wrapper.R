#' Full pipeline for DEGs to GSEA and aPEAR pathway clustering.
#' 
#' A wrapper pipeline function for GSEA pathway analsis and aPEAR clustering of enriched pathways, starting from differential expression (DE) summary statistics. Ideally, GSEA should be used with the full result list, not just the subsetted "significant" DEGs. For more details, see the documentation for the pipeline modules: `deseq_to_gsea()`, `gseares_dotplot_listwrap()`, `gsea_to_aPEAR_clusters()`, `plot_apear_clusters()`, `plot_apear_cluster_pathway_dotplots()`
#'
#' @param deseqres data.frame of differential expression, minimally must have columns named "log2FoldChange", "pvalue", and a gene column (with column name denoted by `gene_identifier_type`)
#' @param pathways data.frame of pathways, output of `preppathways_pathwayanalysis_crosscondition_module()`, recommend you run this once and store use the pathways the way you would for a reference genome
#' @param outdir string, path to output directory. Will write result tables and plots to subdirectories inside of here. 
#' @param verbose T/F, default F
#' @param workernum integer, default 1, number of CPUs
#' @param ... other parameters passed on to pipeline modules. For more details, see the documentation for the pipeline modules: `deseq_to_gsea()`, `gseares_dotplot_listwrap()`, `gsea_to_aPEAR_clusters()`, `plot_apear_clusters()`, `plot_apear_cluster_pathway_dotplots()`
#'
#'
#' @return
#' @export
#'
#' @examples
gsea_apear_pipeline <- function(deseqres, 
                                pathways,
                                outdir,
                                verbose,
                                workernum,
                                ...
                                ){
  
  
  
  
  require(CVRCpipelines)
  
  
  
  if(verbose == T){ message('\nA. RUNNING GSEA\n') }
  
  ## run GSEA
  gseareslist <- deseq_to_gsea(deseqres = deseqres,
                               pathways = pathways,
                               
                               
                               outdir = outdir,
                               workernum = workernum,
                               verbose = verbose,
                               ...
                               
  )
  
  
  # get dotplots
  if(verbose == T){ message('\nB. Making Dotplots\n') }
  
  dp_l <- gseares_dotplot_listwrap(gseareslist, 
                                   outdir = outdir,
                                   ...)
  
  
  
  # run apear clustering
  if(verbose == T){ message('\nC. aPEAR Clustering\n') }
  
  apearoutlist <- gsea_to_aPEAR_clusters(gseareslist,
                                         outdir = outdir,
                                         ...)
  
  
  
  #plot apear
  if(verbose == T){ message('\nD. aPEAR Enrichment Network Plot\n') }
  
  apear_plot <- plot_apear_clusters(apearoutlist,
                                    outdir = outdir,
                                    ...)
  
  # apear_plot
  if(verbose == T){ message('\nE. aPEAR clusters pathways dotplot\n') }
  clust_dp_l <- plot_apear_cluster_pathway_dotplots(apearoutlist,
                                                    outdir = outdir,
                                                    ...
  )
  
  
  
  
  
  fulloutlist <- list(gseareslist = gseareslist,
                      dp_l = dp_l,
                      apearoutlist = apearoutlist,
                      apear_plot = apear_plot,
                      clust_dp_l = clust_dp_l
                      )
  
  
  
  
  return(fulloutlist)
  
}