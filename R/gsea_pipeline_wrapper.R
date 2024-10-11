#' Full pipeline for DEGs to GSEA and aPEAR pathway clustering
#' 
#' A wrapper pipeline function for GSEA pathway analysis and aPEAR clustering of enriched pathways, starting from differential expression (DE) summary statistics. Ideally, GSEA should be used with the full result list.
#'
#' @param deseqres data.frame of differential expression, minimally must have columns named "log2FoldChange", "pvalue", and a gene column (with column name denoted by `gene_identifier_type`)
#' @param pathways data.frame of pathways, output of [CVRCpipelines::preppathways_pathwayanalysis_crosscondition_module()], recommend you run this once and store use the pathways the way you would for a reference genome. Minimally, this is a data.frame where each row is a gene and the pathwy it is part of, with the following column names: "gs_name", the pathways; "gs_subcat", the categories/MSIGDB sub-categories of the pathways; and a column matching the parameter `gene_identifier_type`, such as "gene_symbol" for the genes in each pathway
#' @param outdir string, path to output directory. Will write result tables and plots to subdirectories inside of here. 
#' @param aPEAR_cluster_min_size integer, default 3, minimum aPEAR cluster size after clustering the pathways
#' @param min_pathway_gene_size integer, default 3, minimum number of genes in pathway to be considered for clustering
#' @param num_input_sig_pathways_updn integer, default 250, max number of positive and negative NES pathways included, respectively. for example, when set to 250, 500 pathways maximum will be used, 250 positive and 250 negative
#' @param verbose T/F, default F
#' @param workernum integer, default 1, number of CPUs
#' @param ... Other parameters passed on to [CVRCpipelines::deseq_to_gsea]
#'
#' @return
#' @export
#'
#' @examples
gsea_apear_pipeline <- function(deseqres, 
                                pathways,
                                outdir,
                                aPEAR_cluster_min_size = 3,
                                min_pathway_gene_size = 3,
                                num_input_sig_pathways_updn = 250,
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
                                   outdir = outdir
                                   )
  
  
  
  # run apear clustering
  if(verbose == T){ message('\nC. aPEAR Clustering\n') }
  
  apearoutlist <- gsea_to_aPEAR_clusters(gseareslist,
                                         outdir = outdir
                                         )
  
  
  
  #plot apear
  if(verbose == T){ message('\nD. aPEAR Enrichment Network Plot\n') }
  
  apear_plot <- plot_apear_clusters(apearoutlist,
                                    outdir = outdir)
  
  # apear_plot
  if(verbose == T){ message('\nE. aPEAR clusters pathways dotplot\n') }
  clust_dp_l <- plot_apear_cluster_pathway_dotplots(apearoutlist,
                                                    outdir = outdir
  )
  
  
  
  
  
  fulloutlist <- list(gseareslist = gseareslist,
                      dp_l = dp_l,
                      apearoutlist = apearoutlist,
                      apear_plot = apear_plot,
                      clust_dp_l = clust_dp_l
                      )
  
  
  
  
  return(fulloutlist)
  
}