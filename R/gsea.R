#' Prep MSIGDB pathways for pathway analysis
#'
#' This is a modular component of the scRNAseq analysis pipeline. Prep pathways from MSIGDB via the msigdbr package. We include the Hallmarks category; Gene Ontology BP, MF and CC; Reactome; KEGG; transcription factor CHIP-seq targets in the Gene Transcription Regulation Database (TFT_GTRD); inferred transcription factor targets via motif analysis from Xie et al Nature 2005 (TFT_Legacy). Only gene sets with < 500 genes are included.
#'
#' @param species string, species such as "Homo sapeins" or "Mus musculus"
#' @param outdir_int string, directory to save pathways to. Will create a sub-directory called "pathwayanalysis_crosscondition" and save inside of there. We save pathways since the database updates over time.
#'
#' @return a data.frame similar to the output of `msigdbr::msigdbr`, but filtering for some specific categories / subcategories.
#' @export
#'
#' @examples
#' \dontrun{
#' pathways <- preppathways_pathwayanalysis_crosscondition_module(
#' species = 'Mus musculus',
#' outdir_int = 'path/to/directory')
#' }
preppathways_pathwayanalysis_crosscondition_module <- function(species,
                                                               outdir_int)
{
  
  
  require(msigdbr)
  # require(msigdbf)
  
  #prep the pathways
  # make sure to save it. database can update over time
  pwayoutdir <- paste0(outdir_int, '/pathwayanalysis_crosscondition/')
  dir.create(pwayoutdir, recursive = T)
  
  #read if already there
  if( file.exists( paste0(pwayoutdir, '/msigdb_pathways.rds') ) ){
    
    message('Reading cached msigdbr pathways')
    pathways <- readRDS(paste0(pwayoutdir, '/msigdb_pathways.rds') )
    
    
  }
  
  
  
  message('Accessing MSIGDBR database')
  
  #read pathways
  pathways <- msigdbr::msigdbr(species = species)
  
  
  
  ## 2025.03.31: MSIGDB V10 was released in mid march 2025. it changed the format a lot. 
  # I think for now we can just add the old column names; gs_subcat and gs_cat
  #check if the old columns names are in; if not, add the new columns as the old
  msigdbrcolnames <- colnames(pathways)
  if( any(!c('gs_subcat', 'gs_cat') %in% msigdbrcolnames) ){
    
    pathways$gs_cat <- pathways$gs_collection
    pathways$gs_subcat <- pathways$gs_subcollection
    
    
    #for ease, do this here..
    pathways$gs_subcat <- gsub(':', '_', pathways$gs_subcat)
    
    ## replace "TFT_TFT_LEGACY" with "TFT_TFT_Legacy", as per the old name...
    pathways[pathways$gs_subcat == 'TFT_TFT_LEGACY', "gs_subcat"] <- 'TFT_TFT_Legacy'
    
    #they added a new kegg medicus and old kegg is now kegg_legacy; set the kegg_legacy as CP_KEGG as before
    pathways[pathways$gs_subcat == 'CP_KEGG_LEGACY', "gs_subcat"] <- 'CP_KEGG'
    
  }
  
  
  #replace : with _ in actual pathway names:
  pathways$gs_subcat <- gsub(':', '_', pathways$gs_subcat)
  
  
  
  
  
  #### picking default categories
  # because hallmark is a "category" and rest are "subcategories", it is hard to make this automated
  # guess it may be possible if we set missing subcat as cat...
  # table( pathways[pathways$gs_subcat=='',"gs_cat"] )
  # for now hardcode these
  pwaycats <- c("HALLMARK", "GO_BP", "GO_MF", "GO_CC", "CP_REACTOME", "CP_KEGG", "TFT_GTRD", "TFT_TFT_Legacy")
  
  
  ### try to replace ':' with "_"
  # in pwaycats, user provided subcategories:
  pwaycats <- gsub(':', '_', pwaycats)
  names(pwaycats) <- pwaycats
  
  #in actual pathway names:
  pathways$gs_subcat <- gsub(':', '_', pathways$gs_subcat)
  
  #also get hallmarks...
  pathways[pathways$gs_cat == 'H', 'gs_subcat'] <- "HALLMARK"
  
  #prep pathways using categories defined by user
  pathways <- as.data.frame( pathways[pathways$gs_subcat %in% pwaycats,] )
  
  #fgsea recommends no pathways over 500 genes
  pathways <- pathways[table(pathways$gs_name) <= 500,]
  
  #let's also remove pathways with less than 3 genes
  pathways <- pathways[table(pathways$gs_name) >= 3,]
  
  
  #purge mem
  invisible(gc(full = T, reset = F, verbose = F))
  
  #save it
  
  saveRDS(pathways, paste0(pwayoutdir, '/msigdb_pathways.rds') )
  
  
  return(pathways)
  
  
}





#' Fix underflow in a sorted list of -log pvalues
#'
#' Given a list of P values from Seurat Wilcox or EdgeR test, one should multiply P values times -log10 times sign of lfc. With this input, sometimes underflow occurs. This function will fix this underflow (ie for input to GSEA) by adding +1 to all Inf values and -1 to all -Inf values.
#'
#' @param scores named vector of -log or -log10(pvalues) times sign(LFC), sorted hi to low. For the +inf and -inf values, one suggestion is to sort them by logFoldChange. names are gene sybols or IDs, values are -logP values * sign of LFC
#' @param logFC_vec named vector of some effect size, such as log2fc, by which the underflow P value DEGs may be sorted, or else the sorting may be arbitrary. names are gene symbols or IDs, values are L2FC (or some other effect size such as L2FC * pct.diff)
#'
#' @return a vector of P value scores, where underflow -logP values have been ranked appropriated
#' @export
#'
#' @examples
fix_underflow <- function(scores,
                          logFC_vec,
                          verbose = F){
  
  
  if(missing(logFC_vec)){logFC_vec = NULL; warning('it would be optimal to provide a vector of logFC values or some other effect size to logFC_vec, or else sorting of top DEGs may be arbitrary')}
  
  bool = ( scores == Inf | scores == -Inf)
  tbl = table(factor(bool, levels = c(F,T)))
  
  if(verbose == T){
    message(tbl['TRUE'], ' underflow genes detected')
  }
  
  
  #underflow, positive
  if(any(scores == Inf)){
    scores_uf <- scores[scores == 'Inf']
    
    #sort INF values by lfc
    if(!is.null(logFC_vec)){
      logFC_vec_uf <- logFC_vec[names(logFC_vec) %in% names(scores_uf)]
      logFC_vec_uf <- sort(abs(logFC_vec_uf), decreasing = T)
      scores_uf <- scores_uf[match(names(logFC_vec_uf), names(scores_uf))]
    }
    
    last <- scores[length(scores_uf) + 1]
    
    replacement <- c()
    for(i in 1:length(scores_uf)){
      to_replace = ifelse(i == 1, yes = last, no = replacement[i-1])
      to_replace = to_replace + 1 #for pos + 1
      replacement <- c(replacement, to_replace)
    }
    
    # reverse , lo to hi
    replacement <- rev(replacement)
    names(replacement) <- names(scores_uf)
    
    
    scores[names(replacement)] <- replacement
    
  }
  
  
  #underflow, negative
  if(any(scores == -Inf)){
    
    scores <- rev(scores)
    scores_uf <- scores[scores == '-Inf']
    
    #sort INF values by lfc
    if(!is.null(logFC_vec)){
      logFC_vec_uf <- logFC_vec[names(logFC_vec) %in% names(scores_uf)]
      logFC_vec_uf <- sort(abs(logFC_vec_uf), decreasing = T)
      scores_uf <- scores_uf[names(logFC_vec_uf)]
    }
    
    
    last <- scores[length(scores_uf) + 1]
    
    replacement <- c()
    for(i in 1:length(scores_uf)){
      to_replace = ifelse(i == 1, yes = last, no = replacement[i-1])
      to_replace = to_replace - 1 #for negative only, minus 1...
      replacement <- c(replacement, to_replace)
    }
    
    #for negative only
    replacement <- rev(replacement)
    names(replacement) <- names(scores_uf)
    
    
    
    scores[names(replacement)] <- replacement
    
    scores <- rev(scores)
    
  }
  
  scores <- sort(scores, decreasing = T)
  
  
  return(scores)
}










#' Wrapper script for FGSEA applied to MSIGDB subcategories
#' 
#' Will run FGSEA over the following MSIGB subcategories: "HALLMARK" (technical a category but coded as a sub-category here), "GO_BP", "GO_MF", "GO_CC", "CP_REACTOME", "CP_KEGG", "TFT_GTRD", "TFT_TFT_Legacy". Makes extensive use of the fgsea package. ["Fast Gene Set Enrichment Analysis", Korotkevich et al bioRxiv 2021](https://www.biorxiv.org/content/10.1101/060012v3). Ideally, GSEA should be used with the full result list.
#'
#' @param deseqres data.frame, minimally must have columns named "log2FoldChange", "pvalue", and a gene column (with column name denoted by `gene_identifier_type`)
#' @param gene_identifier_type string, default "gene_symbol", column name of gene identifiers. must be one of 'gene_symbol', 'gene_name', 'ensembl_gene', 'entrez_gene'. if 'gene_name' is passed, it will create and use 'gene_symbol'
#' @param pathways data.frame of pathways, output of [CVRCpipelines::preppathways_pathwayanalysis_crosscondition_module()], recommend you run this once and store use the pathways the way you would for a reference genome. Minimally, this is a data.frame where each row is a gene and the pathwy it is part of, with the following column names: "gs_name", the pathways; "gs_subcat", the categories/MSIGDB sub-categories of the pathways; and a column matching the parameter `gene_identifier_type`, such as "gene_symbol" for the genes in each pathway
#' @param outdir string, path to output directory. will create a subdir called "GSEA_tables"
#' @param pathway_padj_thres numeric, default 0.1, alpha for adjusted pvalue threshold for gsea pathways. Set to Inf (without quotes) to ignore the padj filter.
#' @param pathway_pval_thres numeric, default Inf, alpha for nominal pvalue threshold for gsea pathways, normally not used. If set to any value other than Inf, then pathway_padj_thres will be set to Inf. 
#' @param workernum integer, default 1, number of CPUs, parallelization occurs over the eight database categories, max is 8
#' @param verbose T/F, default F, verbosity
#' @param pwaycats string or character vector. List of pathway subcategories to run. it should match the column of `pathways$gs_subcat`. Default is: c("HALLMARK", "GO_BP", "GO_MF", "GO_CC", "CP_REACTOME", "CP_KEGG", "TFT_GTRD", "TFT_TFT_Legacy")
#' @param filename_prefix string, default is empty. An optional string to add to all saved filenames, useful to open multiple tables later in excel.
#'
#' @return
#' @export
#'
#' @examples
deseq_to_gsea <- function(deseqres,
                          gene_identifier_type = 'gene_symbol',
                          pathways,
                          outdir,
                          pathway_padj_thres = 0.1,
                          pathway_pval_thres = Inf,
                          workernum = 1,
                          verbose = F,
                          pwaycats = c("HALLMARK", "GO_BP", "GO_MF", "GO_CC", "CP_REACTOME", "CP_KEGG", "TFT_GTRD", "TFT_TFT_Legacy"),
                          filename_prefix = ''
                          
){
  
  require(fgsea)
  require(foreach)
  require(doParallel)
  require(parallel)
  require(stringr)
  require(ggplot2)
  
  
  
  
  ######################  Input Exceptions   ###################### 
  
  if(verbose == T){
    # message('Running GSEA - ', paste(date()) )
    message('Running GSEA')
    message('- 0 running input checks')
  }
  
  #check input gene type and adjust if nonstandard ones are given
  accepted_gene_idents = c('gene_symbol', 'gene_name',
                           'ensembl_gene',
                           'entrez_gene')
  
  if(!(gene_identifier_type %in% accepted_gene_idents)){
    
    stop('gene_identifier_type was set to "', gene_identifier_type,'", invalid',
         '\n expected one of ', paste(accepted_gene_idents, collapse = ', '))
    
  }
  
  
  
  
  ## for nonstandard ones, do some adjusting
  
  #gene_name; use gene_symbol, and make sure to add this column to deseq2 name...
  if(gene_identifier_type == 'gene_name'){
    warning('gene_identifier_type was set to "gene_name", will use "gene_symbol"',
            '\nalso setting deseqres$gene_symbol <- deseqres$gene_name')
    deseqres$gene_symbol <- deseqres$gene_name
    gene_identifier_type <- 'gene_symbol'
  }
  
  
  
  ## ensembl_gene, make sure version is stripped
  if(gene_identifier_type == 'ensembl_gene'){
    
    if( any(grepl('\\.',deseqres[,gene_identifier_type])) ){
      
      warning('period symbols (".") deteced in deseqres$ensembl_gene, potentially denoting ensembl version identifiers in additon to gene IDs. Will strip to ensembl ID preceeding "." symbols.')
      
      deseqres[,gene_identifier_type] <- sapply(strsplit( deseqres[,gene_identifier_type], '\\.'), function(x){x[1]})
      
    }
    
  }
  
  
  
  #check colnames of deseq res input
  important_colnames <- c(gene_identifier_type, 'log2FoldChange', 'pvalue')
  
  if(any(!(important_colnames %in% colnames(deseqres)))){
    
    missingimportant_colnames <- important_colnames[!(important_colnames %in% colnames(deseqres))]
    
    missingimportant_colnames <- paste(missingimportant_colnames, collapse = ', ')
    
    stop('column(s) "',missingimportant_colnames, '" are missing from the column names of deseqres.\nPlease make sure to format the input deseq2 res with the following:\n', paste(important_colnames, collapse = ', '))
    
  }
  
  
  
  
  #make sure workernum is not more than 8
  if(workernum > 8){
    warning('workernum was set too high to ', workernum, ', reducing to 8 as this function is parallelized over 8 database categories')
    
    workernum <- 8
  }
  
  
  
  #warn about outdir not being set
  if(missing(outdir)){warning('no outdir set; will run but will not write to disc')}
  
  #warn about dir creation
  if(!missing(outdir)){
    
    if(!dir.exists(outdir)){
      
      warning('outdir was set to:\n', outdir, '\nthis folder does not exist; creating')
      
      dir.create(outdir, recursive = T)
      
    }
    
  }
  
  
  
  
  
  
  
  
  # check number of genes from deseqres that are in the database
  
  tab <- table(deseqres[,gene_identifier_type] %in% pathways[,gene_identifier_type])
  perctab <- round(tab/sum(tab)*100, 1)
  
  
  if( all(!(deseqres[,gene_identifier_type] %in% pathways[,gene_identifier_type])) ){
    
    stop('No genes from deseqres$', gene_identifier_type, ' were detected in pathways$', gene_identifier_type,'.',
         '\nPlease make sure you are using the correct gene_identifier_type and check all inputs'
    )
    
  }
  
  if(verbose == T){
    
    message(tab['TRUE'], ' genes from deseqres$', gene_identifier_type, ' (',perctab['TRUE'],'%)', ' were detected in pathways$', gene_identifier_type,'.')
    
  }
  
  if(perctab['TRUE'] < 50){
    warning('Only ',tab['TRUE'], ' genes from deseqres$', gene_identifier_type, ' (',perctab['TRUE'],'%)', ' were detected in pathways$', gene_identifier_type,'.')
    
  }
  
  
  
  ## make sure pwaycats are in pathways
  # this is the pre-defined list of subcategories we will loop over
  # pwaycats <- c("HALLMARK", "GO_BP", "GO_MF", "GO_CC", "CP_REACTOME", "CP_KEGG", "TFT_GTRD", "TFT_TFT_Legacy")
  names(pwaycats) <- pwaycats
  
  unique_subcats <- unique(pathways$gs_subcat)
  
  if(any(!(pwaycats %in% unique_subcats))){
    
    
    missing_pwaycats <- pwaycats[!(pwaycats %in% unique_subcats)]
    
    stop('the following database subcategories "',paste(missing_pwaycats, collapse = ', '), '" are missing from unique(pathways$gs_subcat).\nExpected to find ', paste(pwaycats, collapse = ', '), 
         '\nPlease run the function: preppathways_pathwayanalysis_crosscondition_module()')
    
  }
  
  
  
  #if using pval thres, do not use padj thres
  # update 2025.07.31; set this to inf instead of 1 otherwise pathways get filtered out
  
  if(pathway_pval_thres != Inf){
    
    # update 2025.07.31; set this to inf instead of 1 otherwise pathways get filtered out
    pathway_padj_thres <- Inf
    
    
  }
  
  
  
  ###################### End Input Exceptions   ###################### 
  
  
  
  if(verbose == T){
    message('- 1 formatting inputs')
  }
  
  #subset and keep just important columns of pathways
  pathways <- pathways[,c('gs_subcat', 'gs_name', gene_identifier_type)]
  
  
  #keep just key cols of deseq res
  # gene_name, log2FoldChange, pvalue
  res <- deseqres[,important_colnames]
  
  #remove NAs
  res <- res[complete.cases(res),]
  
  
  
  #make sure they are in the gene univerise of pathways
  res <- res[res[,gene_identifier_type] %in% pathways[,gene_identifier_type],]
  
  
  
  
  
  
  ## prep weighted list ##
  
  #we have to deal with underflow...
  # get -log10 pvalues, sort
  scores <- -log10(res$pvalue)
  scores <- scores * sign(res$log2FoldChange)
  names(scores) <- res[,gene_identifier_type]
  
  #sort by log pval with names
  scores <- sort(scores,decreasing = T)
  
  
  #also get logFC vector; for INf, we will sort them by LFC...
  logFC_vec <- res$log2FoldChange; names(logFC_vec) <- res[,gene_identifier_type]
  
  
  # fix the underflow...
  scores <- CVRCpipelines::fix_underflow(scores, logFC_vec, verbose = verbose)
  
  
  
  #make sure scores is in order of genes
  scores <- scores[match(res[,gene_identifier_type], names(scores))]
  
  
  #add to res
  res$weight <- scores
  
  #order by weight
  res <- res[order(res$weight, decreasing = T),]
  
  
  #just getr scores
  scores <- res$weight
  names(scores) <- res[,gene_identifier_type]
  
  
  #clean env
  Sys.sleep(2)
  rm(res, deseqres, logFC_vec)
  gc(full = T)
  
  
  
  
  
  
  if(verbose == T){
    message('- 2 initiating GSEA\n')
  }
  
  
  #loop over pwaycats and do GSEA
  # loop over pre-defined list of pathway categories, defined above
  cat <- pwaycats[1] #for testing
  
  
  
  ## initiate parallelization
  cl <- parallel::makeCluster(workernum)
  doParallel::registerDoParallel(cl)
  
  
  
  
  gseareslist <- foreach(cat = pwaycats, 
                         .export = c('pathways', 'scores', 'outdir', 'pathway_padj_thres', 'pathway_pval_thres', 'gene_identifier_type'),
                         .packages = c('fgsea', 'stringr', 'ggplot2'),
                         .verbose = verbose
  ) %dopar% {
    
    # gseareslist <- lapply(cats, function(type){
    
    
    set.seed(42069)
    
    message(cat)
    
    
    t2g <- pathways[pathways$gs_subcat == cat,c('gs_name', gene_identifier_type)]
    pwayl = split(t2g[,gene_identifier_type], t2g$gs_name)
    rm(t2g)
    
    
    ## run GSEA ##
    # first run multilevel, then try npermsimple
    gseares <- fgsea::fgsea(pathways=pwayl, stats=scores, nproc = 1)
    
    invisible(gc(full = T, reset = F, verbose = F))
    
    
    #sometimes there are NAs due to "severely unbalanced pathways", try to fix
    if( any(is.na(gseares$NES)) ){
      rm(gseares)
      
      gseares <- fgsea::fgsea(pathways=pwayl, stats=scores, nPermSimple=10000, nproc = 1)
      
      invisible(gc(full = T, reset = F, verbose = F))
      
    }
    
    
    
    #format as data.frame instead of data.table
    gseares <- as.data.frame(gseares)
    
    #ensure no NAs are kept
    gseares <- gseares[complete.cases(gseares[,1:7]),,drop=F]
    
    #order by NES
    gseares <- gseares[order(gseares$NES, decreasing = T),]
    
    #apply cutoff of pathway_pval_thres, then pathway_padj_thres
    gseares <- gseares[gseares$pval < pathway_pval_thres, ,drop=F]
    gseares <- gseares[gseares$padj < pathway_padj_thres, ,drop=F]
    
    #select pathways with more than just 1 gene in the list
    # update 2025.07.31; do not apply this filter; we will filter this in aPEAR anyway
    # gseares <- gseares[gseares$size > 2,,drop=F]
    
    #add the category, if there are any rows at all...
    if( nrow(gseares) > 0 ){
      gseares$category <- cat
      gseares <- gseares[,c('pathway', 'category', 'pval', 'padj', 'log2err', 'ES', 'NES', 'size', 'leadingEdge')]
    }
    
    
    #also, adjust the leading egde
    gseares$leadingEdge <- sapply(gseares$leadingEdge, function(x){ paste(x, collapse = '/') })
    
    #return the result table and the plot
    return(gseares)
    
    
    
    
    
    
    
    # })
    
  }
  
  names(gseareslist) <- pwaycats
  
  
  parallel::stopCluster(cl)
  
  
  #count number of enriched pathways
  numpways <- sapply(gseareslist, nrow)
  numpways_updn <- sapply(gseareslist, function(x){
    table(factor(sign(x$NES), c(1, -1)))
  })
  
  if(verbose == T){
    
    message('Num pathways sig:\n',
            paste0(names(numpways), ' = ', numpways, ', up = ', as.data.frame.matrix(numpways_updn)[1,], ', dn = ',as.data.frame.matrix(numpways_updn)[2,],  '\n')
    )
    
  }
  
  
  
  
  ## remove empty cats
  numpways <- numpways[numpways > 0]
  
  gseareslist <- gseareslist[names(gseareslist) %in% names(numpways)]
  
  
  
  #write out the tables
  if(!missing(outdir)){
    
    
    
    outdir_tables <- paste0(outdir, '/GSEA_tables/')
    dir.create(outdir_tables, recursive = T)
    
    
    if(verbose == T){message(' - 4 writing outputs to ', outdir_tables)}
    
    lapply( names(gseareslist) , function(cat){
      gseares <- gseareslist[[cat]]
      resfile <- paste0(outdir_tables, '/GSEA_results_', cat, '.csv')
      
      if(filename_prefix != ''){
        resfile <- paste0(outdir_tables, '/GSEA_results_', cat, '_', filename_prefix, '.csv')
      }
      
      write.csv(gseares, resfile, row.names = F)
      return(cat)
    })
    
  }
  
  
  
  #write some log information out
  
  logfile <- paste0(outdir, '/logfile.txt')
  
  write_lines_l <- list(date(),
                        paste0('fgsea version - ', packageVersion('fgsea')),
                        paste0('CVRCpipelines version - ', packageVersion('CVRCpipelines')),
                        
                        paste0('\n'),
                        
                        paste0('pathway_pval_thres - ', pathway_pval_thres),
                        paste0('pathway_padj_thres - ', pathway_padj_thres),
                        
                        paste0('\n\n\n'),
                        
                        paste0(tab['TRUE'], ' genes from deseqres$', gene_identifier_type, ' (',perctab['TRUE'],'%)', ' were detected in pathways$', gene_identifier_type,'.'),
                        
                        paste0('\n'),
                        
                        paste0('Num pathways sig:'),
                        paste( paste0(names(numpways), ' = ', numpways, ', up = ', as.data.frame.matrix(numpways_updn)[1,], ', dn = ',as.data.frame.matrix(numpways_updn)[2,]), collapse='\n'),
                        
                        paste0('\n\n\n')
                        
                        
  )
  
  write('FGSEA log\n\n', logfile)
  
  lapply(write_lines_l, function(x){write(x, logfile, append = T)})
  
  
  
  return(gseareslist)
  
  
}













#' Make a dotplot from fgsea result dataframe
#'
#' @param gseares fgsea result data.frame
#' @param n_up integer, default 5, max number of upregulated pathways to plot (positive NES)
#' @param n_dn integer, default 5, max number of dnregulated pathways to plot (positive NES)
#' @param line_char_width integer, default 35, max character width per line of pathway y axis labels, newline added after words, see [stringr::str_wrap()]
#' @param up_color color string name ort code, default 'red'
#' @param midpoint_color color string name ort code, default 'white'
#' @param dn_color color string name ort code, default 'steelblue'
#' @param dotplot_fontsize numeric, size of y axis pathway label font
#' @param include_pways vector of characters; if included, will always look for and include these even if not in n_up or n_dn
#'
#' @return a ggplot dotplot
#' @export
#'
#' @examples
gseares_dotplot <- function(gseares,
                            n_up = 5,
                            n_dn = 5,
                            line_char_width = 35,
                            up_color = 'red',
                            midpoint_color = 'white',
                            dn_color = 'steelblue',
                            dotplot_fontsize = 5,
                            include_pways
                            
){
  
  
  

  
  #update 2025.03.31; do this before cutoff
  # #make sure order is by -log(padj) * NES# update; do this before filtering
  gseares$weight <- (-log(gseares$padj)) * sign(gseares$NES)
  gseares <- gseares[order(gseares$weight, decreasing = T),]
  # gseares_plot$pathway <- factor(gseares_plot$pathway, levels = rev(gseares_plot$pathway)  )
  
  
  #if more than 20 ,select just 20
  # gseares <- gseares_plot #overwrite this; this is important for later...
  gseares_plot <- gseares
  gseares_plot <- rbind( head( gseares_plot[gseares_plot$NES>0,,drop=F], n_up) ,
                         tail( gseares_plot[gseares_plot$NES<0,,drop=F], n_dn) )
  
  
  #also look for and include the include_pways
  if(!missing(include_pways)){
    
    
    
    gseares_plot <- rbind(gseares_plot,
                          gseares[gseares$pathway %in% include_pways,]
    )
    
    #make sure no dups
    gseares_plot <- gseares_plot[!duplicated(gseares_plot$pathway),]
    
  }
  
  
  #make pathway names more readable by using spaces instead of underscores
  gseares_plot$pathway <- gsub(gseares_plot$pathway, pattern = '_', replacement = ' ')
  
  #make pathway names more readable by splitting long ones to multiple lines
  gseares_plot$pathway <- stringr::str_wrap(gseares_plot$pathway, width = line_char_width)
  
  # #make sure order is by -log(padj) * NES# update; do this before filtering now
  ## REPEAT IT HERE TO MAKE SURE FACTOR IS WORKING PROPERLY
  gseares_plot$weight <- -log(gseares_plot$padj) * sign(gseares_plot$NES)
  gseares_plot <- gseares_plot[order(gseares_plot$weight, decreasing = T),]
  gseares_plot$pathway <- factor(gseares_plot$pathway, levels = rev(gseares_plot$pathway)  )
  
  
  
  #plot it
  
  #fix color issue when just 1 obs
  if(nrow(gseares_plot) == 1){
    dp_single_col <- sign(gseares_plot$NES)
    dp_single_col <- ifelse(dp_single_col==1, yes = up_color, no = dn_color)
    
    
    dp <- ggplot(gseares_plot, aes(-log10(padj), pathway, col=NES, size = size))+
      geom_point()+
      theme_linedraw()+
      theme(axis.text=element_text(size=dotplot_fontsize) )+
      scale_color_gradientn(colors = dp_single_col) +
      scale_size(range=c(2,6))
    
  } else{
    
    
    dp <- ggplot(gseares_plot, aes(-log10(padj), pathway, col=NES, size = size))+
      geom_point()+
      theme_linedraw()+
      theme(axis.text=element_text(size=dotplot_fontsize) )+
      scale_color_gradient2(low = dn_color, high = up_color, mid = midpoint_color, midpoint = 0, name = 'Normalized\nEnrichment\nScore')+
      scale_size(range=c(2,6))
    
  }
  
  
  return(dp)
  
}






#' Make a list of dotplots from a list of fgsea dataframes
#' 
#' list wrapper around [CVRCpipelines::gseares_dotplot()]. see that function for details 
#'
#' @param gseareslist output of `deseq_to_gsea()`
#' @param outdir string, path to output directory, will create a subdir called "GSEA_Category_Dotplots"
#' @param plotprefix string, default is 'GSEA_Dotplot', will name each of the category dotplot files with this prefix, ie "GSEA_Dotplot_HALLMARK.pdf". If changed, no need to put a "_" at the end, it is added automatically
#' @param pdfheight numeric, default 7
#' @param pdfwidth numeric, default 7
#' @param filename_prefix string, default is empty. An optional string to add to all saved filenames, useful to open multiple tables later in excel.
#' @param ... other arguments passed to [CVRCpipelines::gseares_dotplot()]
#'
#' @return a list of ggplot dotplots
#' @export
#'
#' @examples
gseares_dotplot_listwrap <- function(gseareslist, 
                                     outdir, 
                                     plotprefix = 'GSEA_Dotplot',
                                     pdfheight=7, 
                                     pdfwidth=7, 
                                     filename_prefix = '',
                                     ...){
  
  
  
  dp_l <- lapply(names(gseareslist), function(cat){
    gseares <- gseareslist[[cat]]
    gseares_dotplot(gseares, ...) + labs(subtitle = cat)
  })
  names(dp_l) <- names(gseareslist)
  
  
  
  if(!missing(outdir)){
    
    outdir_dp_l <- paste0(outdir, '/GSEA_Category_Dotplots/')
    dir.create(outdir_dp_l, recursive = T)
    
    
    if(filename_prefix != ''){
      plotprefix <- paste0(plotprefix, '_', filename_prefix)
    }
    
    
    save_plotlist_suffix_listnames(plotfilename = plotprefix,
                                   plotlist = dp_l,
                                   plotdir = outdir_dp_l,
                                   pdfheight = pdfheight,
                                   pdfwidth = pdfwidth  
    )
    
    
  }
  
  
  return(dp_l)
  
}









#' aPEAR Wrapper
#' 
#' Using the gsea result list from [CVRCpipelines::deseq_to_gsea()], run aPEAR. ["Advanced Pathway Enrichment Analysis Representation", Kerseviciute and Gordevicius Bioinformatics 2023.](https://academic.oup.com/bioinformatics/article/39/11/btad672/7342237). Note many pathways may not be clustered due to uniqueness or being assigned to clusters below minimum cluster size.
#'
#' @param gseareslist output of [CVRCpipelines::deseq_to_gsea()]
#' @param aPEAR_cluster_min_size integer, default 3, minimum aPEAR cluster size after clustering the pathways
#' @param min_pathway_gene_size integer, default 3, minimum number of genes in pathway to be considered for clustering
#' @param num_input_sig_pathways_updn integer, default 250, max number of positive and negative NES pathways included, respectively. for example, when set to 250, 500 pathways maximum will be used, 250 positive and 250 negative
#' @param outdir string, path to output directory, will create a subdir called "aPEAR"
#' @param filename_prefix string, default is empty. An optional string to add to all saved filenames, useful to open multiple tables later in excel.
#'
#' @return
#' @export
#'
#' @examples
gsea_to_aPEAR_clusters <- function(gseareslist,
                                   aPEAR_cluster_min_size = 3,
                                   min_pathway_gene_size = 3,
                                   num_input_sig_pathways_updn = 250,
                                   outdir,
                                   filename_prefix = ''
                                   
                                   
){
  
  
  require(dplyr)
  require(aPEAR)
  
  
  ## get a df of gseares
  gdf <- dplyr::bind_rows(gseareslist)
  
  
  
  
  ## prep input pathways ##
  
  apear_input <- gdf
  
  #get pways above gene size threshold
  apear_input <- apear_input[apear_input$size >= min_pathway_gene_size, ,drop=F]
  
  
  
  #get top/bottom pathways
  apear_input$score <- -log(apear_input$pval) * sign(apear_input$NES)
  apear_input <- apear_input[order(apear_input$score, decreasing = T),]
  
  apear_input <- rbind(head(apear_input[apear_input$NES>0,], num_input_sig_pathways_updn),
                       tail(apear_input[apear_input$NES<0,], num_input_sig_pathways_updn))
  
  
  
  
  
  ## run aPEAR
  
  # reformat input #
  
  # colnames(apear_input)
  
  apear_input_formatted <- data.frame(Description = apear_input$pathway,
                                      pathway = apear_input$pathway,
                                      NES = apear_input$NES,
                                      pathwayGenes = apear_input$leadingEdge,
                                      Size = apear_input$size,
                                      pvalue = apear_input$pval
  )
  
  ## format pathway names a bit ##
  
  #whitespace instead of _
  apear_input_formatted$Description <- gsub('_', ' ', apear_input_formatted$pathway)
  
  #title case
  apear_input_formatted$Description <- stringr::str_to_title(apear_input_formatted$Description)
  
  # keep GOBP, GOCC, GOMF, and KEGG all caps
  apear_input_formatted$Description <- gsub('Gobp', 'GOBP', apear_input_formatted$Description)
  apear_input_formatted$Description <- gsub('Gocc', 'GOCC', apear_input_formatted$Description)
  apear_input_formatted$Description <- gsub('Gomf', 'GOMF', apear_input_formatted$Description)
  apear_input_formatted$Description <- gsub('Kegg', 'KEGG', apear_input_formatted$Description)
  apear_input_formatted$Description <- gsub('rna', 'RNA', apear_input_formatted$Description)
  apear_input_formatted$Description <- gsub('dna', 'DNA', apear_input_formatted$Description)
  
  
  #wrap long names
  # apear_input_formatted$Description <- stringr::str_wrap(apear_input_formatted$Description, width = line_char_width)
  # apear plotting seems to do this already
  
  
  
  ### CLUSTER ###
  apearclusts <- aPEAR::findPathClusters(apear_input_formatted, 
                                         minClusterSize = aPEAR_cluster_min_size,
                                         verbose = T
  )
  
  
  
  
  # # #test plot
  # apear_plot <- aPEAR::plotPathClusters(apear_input_formatted,
  #                                       apearclusts$similarity, apearclusts$clusters,
  #                                       colorBy = 'NES',
  #                                       nodeSize='Size',
  #                                       pCutoff = min(log10(apear_input_formatted$pvalue)),
  #                                       repelLabels = T, drawEllipses = F, max.overlaps = Inf,
  #                                       fontSize = 3
  # )
  # 
  # apear_plot
  
  
  
  
  #add apear clusters to main gdf result
  apearinput_MATCHNAME <- apear_input_formatted[match(apearclusts$clusters$Pathway, apear_input_formatted$Description),]
  apearclusts$clusters$OrigDesc <- apearinput_MATCHNAME$pathway
  
  gdf$aPEAR_Clusters = NA
  gdf[match(apearclusts$clusters$OrigDesc, gdf$pathway),'aPEAR_Clusters'] <- apearclusts$clusters$Cluster
  
  #remove the weird characters...
  gdf$aPEAR_Clusters <- gsub('\n', ' ', gdf$aPEAR_Clusters)
  
  
  #order the output by NES
  gdf = gdf[order(gdf$NES , decreasing = T),]
  
  
  
  apearoutlist <- list(apear_input_formatted=apear_input_formatted,
                       apearclusts=apearclusts,
                       gdf=gdf
  )
  
  
  
  if(!missing(outdir)){
    
    
    outdir_apear <- paste0(outdir, '/aPEAR/')
    dir.create(outdir_apear, recursive = T)
    
    
    #write it in the main outdir
    gdffile <- paste0(outdir, '/SignificantPathwaysTable_WithClusters.csv')
    
    if(filename_prefix != ''){
      gdffile <- paste0(outdir, '/SignificantPathwaysTable_WithClusters_', filename_prefix, '.csv')
    }
    
    
    write.csv(gdf, gdffile, row.names = F)
    
    apearoutlist_file <- paste0(outdir_apear, '/robject_apearlist.rds')
    saveRDS(apearoutlist, apearoutlist_file)
    
    
    
    
    #write some log information out
    numpways <- table(gdf$aPEAR_Clusters)
    numpways['UNCLUSTERED'] = table(factor(is.na(gdf$aPEAR_Clusters), c('FALSE', 'TRUE')))['TRUE']
    
    
    logfile <- paste0(outdir, '/logfile.txt')
    
    write_lines_l <- list(# date(),
                          # paste0('fgsea version - ', packageVersion('fgsea')),
                          paste0('aPEAR version - ', packageVersion('aPEAR')),
                          
                          paste0('\n'),
                          
                          paste0('aPEAR_cluster_min_size - ', aPEAR_cluster_min_size),
                          paste0('min_pathway_gene_size - ', min_pathway_gene_size),
                          paste0('num_input_sig_pathways_updn - ', num_input_sig_pathways_updn),
                          
                          paste0('\n\n\n'),
                          
                          paste0('Clusters and num pathways:'),
                          paste( paste0(names(numpways), ' = ', numpways, ' pathways'), collapse='\n') 
                          
                          
                          
                          
    )
    
    write('aPEAR log\n\n', logfile, append = T)
    
    lapply(write_lines_l, function(x){write(x, logfile, append = T)})
    
    
  }
  
  
  
  return(apearoutlist)
  
  
}





#' Wrapper for aPEAR plotting
#'
#' A wrapper for plotting the aPEAR enrichment network. Sometimes looks weird, so it may be useful to run this a few times, will produce a different plot each time.
#'
#' @param apearoutlist the output of [CVRCpipelines::gsea_to_aPEAR_clusters()]
#' @param outdir string, path to output directory, will create a subdir called "aPEAR/aPEAR_GraphPlots"
#' @param date_plots T/F, default T, whether to append datetime to apear graph plots, useful to remake plots in case they look weird (which is not uncommon)
#' @param repelLabels T/F, default T, adjust labels to avoid overlaps, highly imperfect
#' @param fontSize numeric, default 2.7, cluster label font size
#' @param pdfheight numeric, default 7
#' @param pdfwidth numeric, default 10
#' @param filename_prefix string, default is empty. An optional string to add to all saved filenames, useful to open multiple tables later in excel.
#' @param ... other arguments passed to [aPEAR::plotPathClusters()]
#'
#' @return
#' @export
#'
#' @examples
plot_apear_clusters <- function(apearoutlist,
                                outdir,
                                date_plots = T,
                                repelLabels = T,
                                fontSize = 2.7,
                                pdfheight = 7,
                                pdfwidth = 10,
                                filename_prefix = '',
                                ...
){
  
  
  
  
  
  apear_input_formatted <- apearoutlist$apear_input_formatted
  apearclusts <- apearoutlist$apearclusts
  
  
  
  # make plot
  apear_plot <- aPEAR::plotPathClusters(apear_input_formatted,
                                        apearclusts$similarity, apearclusts$clusters,
                                        colorBy = 'NES',
                                        nodeSize='Size',
                                        pCutoff = min(log10(apear_input_formatted$pvalue)),
                                        repelLabels = repelLabels, 
                                        drawEllipses = F, max.overlaps = Inf,
                                        fontSize = fontSize,
                                        ...
  )
  
  
  
  #add whitespace for the long pway names
  apear_plot <- apear_plot+ scale_x_continuous(expand = expansion(mult = .2))
  
  
  
  if(!missing(outdir)){
    
    
    
    outdir_apear_plot <- paste0(outdir, '/aPEAR/aPEAR_GraphPlots/')
    dir.create(outdir_apear_plot, recursive = T)
    
    
    
    plotname <- paste0(outdir_apear_plot, '/aPEAR_GraphPlot')
    
    
    if(filename_prefix != ''){
      plotname <- paste0(plotname, '_', filename_prefix)
    }
    
    
    
    if(date_plots == T){
      
      plot_date <- gsub(' ', '_', date())
      plot_date <- gsub('__', '.', plot_date)
      plot_date <- gsub(':', '_', plot_date)
      
      
      plotname <- paste0(plotname, '.', plot_date)
      
      
      apear_plot <- apear_plot + labs(subtitle = plot_date)
      
    }
    
    
    
    plotfile <- paste0(plotname, '.pdf')
    
    
    pdf(plotfile, height = pdfheight, width = pdfwidth)
    
    print(apear_plot)
    
    dev.off()
    
    
  }
  
  
  
  
  
  return(apear_plot)
  
  
}








#' aPEAR cluster pathways dotplots
#' 
#' For each aPEAR cluster, make a dotplot of the pathways in that cluster
#'
#' @param apearoutlist the output of [CVRCpipelines::gsea_to_aPEAR_clusters()]
#' @param outdir string, path to output directory, will create a subdir called "aPEAR/aPEAR_ClusterPathways_Dotplots"
#' @param pdfheight integer, default 7
#' @param pdfwidth integer, default 7
#' @param filename_prefix string, default is empty. An optional string to add to all saved filenames, useful to open multiple tables later in excel.
#' @param ... other arguments passed to [CVRCpipelines::gseares_dotplot()]
#'
#' @return
#' @export
#'
#' @examples
plot_apear_cluster_pathway_dotplots <- function(
    apearoutlist,
    outdir,
    pdfheight = 7,
    pdfwidth = 7, 
    filename_prefix = '',
    ...
){
  
  
  
  gdf <- apearoutlist$gdf
  
  
  #get clusters, including the unclustered
  
  clusters <- names(sort(table(gdf$aPEAR_Clusters), decreasing = T))
  names(clusters) <- clusters
  
  
  clusters['UNCLUSTERED'] <- 'UNCLUSTERED'
  
  gdf[is.na(gdf$aPEAR_Clusters),'aPEAR_Clusters'] <- 'UNCLUSTERED'
  
  clust = clusters[1] #test
  
  
  
  clust_dp_l <- lapply(clusters, function(clust){
    
    
    gdfclust <- gdf[gdf$aPEAR_Clusters == clust,]
    
    
    #make sure the actual clustname is in the plot...
    origclustname <- clust
    origclustname <- gsub(' ', '_', clust)
    origclustname <- gdfclust$pathway[grepl(origclustname, gdfclust$pathway, ignore.case = T)]
    
    #add title, but some are long so wrap
    clust_title <- stringr::str_wrap(clust, width = 100)
    
    gseares_dotplot(gdfclust, include_pways = origclustname) + labs(subtitle = clust_title)
    
    
    
  })
  
  
  #for saving, filenames replace whitespaces with '_
  names(clust_dp_l) <- gsub(' ', '_', names(clust_dp_l))
  
  
  
  if(!missing(outdir)){
    
    
    outdir_apear_dotplots <- paste0(outdir, '/aPEAR/aPEAR_ClusterPathways_Dotplots/')
    dir.create(outdir_apear_dotplots, recursive = T)
    
    plotprefix = 'Dotplot_aPEARcluster'
    
    if(filename_prefix != ''){
      plotprefix <- paste0(plotprefix, '_', filename_prefix)
    }
    
    
    
    save_plotlist_suffix_listnames(plotfilename = plotprefix,
                                   plotlist = clust_dp_l,
                                   plotdir = outdir_apear_dotplots,
                                   pdfheight = pdfheight,
                                   pdfwidth = pdfwidth  
    )
    
    
    
    
  }
  
  
  
  
  return(clust_dp_l)
  
  
  
  
}







