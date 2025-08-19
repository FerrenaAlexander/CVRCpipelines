#' Save a named list of plots to individual pdf files, one file per named list element
#' 
#' given a named list of plots, save a pdf file for each list element
#'
#' @param plotfilename string, prefix that will be used to
#' @param plotlist the list to be plotted. must be named
#' @param plotdir the directory to plot to. will *not* create a new dir (for safety)
#' @param pdfheight height of pdf
#' @param pdfwidth width of pdf
#' @param title_separator string, default '_', separator string between `plotfilename` and list element names in pdf file name
#' @param verbose T/F, default F, whether to report which element is being plotted during loop, useful to detect errors
#'
#' @return a list of files in the `plotdir`
#' @export
#'
#' @examples
#' \dontrun{
#' #set plotfilename, a string prefix
#' plotfilename <- 'ScatterPlots'
#' 
#' #set plotlist, a named list of plots, one pdf per list element;
#' # if any element is itself a list, the pdf will have pages for each sublist element
#' plotlist <- list(plot1 = plot(rnorm(100)), plot2 = plot(rnorm(100)))
#' 
#' #set the output folder
#' plotdir <- 'output_plots_will_go_here'
#' 
#' #title_separator; a string character that separates the 
#' # "plotfilename" prefix from the "plotlist" named list names
#' # here, the plots will be named like "ScatterPlots_plot1.pdf" and "ScatterPlots_plot2.pdf"
#' plotfilename = '_'
#' 
#' CVRCpipelines::save_plotlist_suffix_listnames(plotfilename,
#'                                               plotlist=plotlist,
#'                                               plotdir=plotdir)
#' }
save_plotlist_suffix_listnames <- function(
    plotfilename,
    plotlist,
    plotdir,
    pdfheight = NULL,
    pdfwidth = NULL,
    title_separator = '_',
    verbose = F
){
  
  
  #create dir
  # dir.create(plotdir, recursive = T)
  # nah, do this outside just to be super safe
  
  if(is.null(names(plotlist))){
    stop('plotlist lacks names. please use something like names(plotlist) <- c() to provide names to plotlist')
  }
  
  
  lapply(1:length(plotlist), function(i){
    
    #get list element name and append as suffix to general plotfile name
    plotlistelementname <- names(plotlist)[i]
    
    if(verbose == T){message(plotlistelementname)}
    
    
    # plotfilename <- names(plotlist)[i]
    saveplot <- plotlist[[i]]
    # plotfilename <- paste0(plotfilename,      '_',      plotlistelementname, '.pdf')
    plotfilename <- paste0(plotfilename, title_separator, plotlistelementname, '.pdf')
    
    plotfilepathname <- paste0(plotdir, '/', plotfilename)
    pdf(plotfilepathname, height = pdfheight, width = pdfwidth)
    print(saveplot)
    dev.off()
    
    return()
  })
  
  message('\nFiles in: ', plotdir)
  
  
  return(list.files(plotdir))
}







#' Quick one-line wrapper for saving pdf
#'
#' @param plotobject any storable plot object than can be printed to pdf (ie ggplot)
#' @param filepath string, path to save to
#' @param pdfheight numeric, height of pdf
#' @param pdfwidth numeric, width of pdf
#' @param automakedir T/F, default F, whether to detect the folder path from `filepath` and recursively create the dir or not, can be dangerous and messy if not careful
#'
#' @return
#' @export
#'
#' @examples
quickpdf <- function(plotobject, filepath, pdfheight=7, pdfwidth=7, automakedir=F){
  
  if(automakedir == F){
    if(!dir.exists(dirname(filepath))){
      stop('the folder of the input filepath seems to not be found: ', dirname(filepath) ) 
    }
  } else{
    if(!dir.exists(dirname(filepath))){
      warning('the folder of the input filepath seems to not be found: ', dirname(filepath), '\ncreating path and saving' ) 
    }
    dir.create(dirname(filepath), recursive = T)
  }
  
  
  graphics.off()
  
  pdf(filepath, height = pdfheight, width = pdfwidth)
  
  print(plotobject)
  
  graphics.off()
  
  
}



#' Quick one-line wrapper for saving png
#'
#' @param plotobject any storable plot object than can be printed to pdf (ie ggplot)
#' @param filepath string, path to save to
#' @param plotheight numeric, height of plot file in units, default 7
#' @param plotwidth numeric, height of plot file in units, default 7
#' @param units string, default 'in', see `?png`
#' @param automakedir T/F, default F, whether to detect the folder path from `filepath` and recursively create the dir or not, can be dangerous and messy if not careful
#' @param pdfheight numeric, holdover from other plotting functions, will transfer to plotheight
#' @param pdfwidth numeric, holdover from other plotting functions, will transfer to plotwidth
#'
#' @return
#' @export
#'
#' @examples
quickpng <- function(plotobject, filepath, plotheight=7, plotwidth=7, automakedir=F, pdfheight=NULL, pdfwidth=NULL, units='in'){
  
  
  if(!is.null(pdfheight)){
    warning('For PNG use plotheight and plotwidth, not pdfheight and pdfwidth. \nSetting plotheight<-pdfheight')
    plotheight <- pdfheight
  }
  if(!is.null(pdfwidth)){
    warning('For PNG use plotheight and plotwidth, not pdfheight and pdfwidth. \nSetting plotwidth<-pdfwidth')
    plotwidth <- pdfwidth
  }
  
  
  
  if(automakedir == F){
    if(!dir.exists(dirname(filepath))){
      stop('the folder of the input filepath seems to not be found: ', dirname(filepath) ) 
    }
  } else{
    if(!dir.exists(dirname(filepath))){
      warning('the folder of the input filepath seems to not be found: ', dirname(filepath), '\ncreating path and saving' ) 
    }
    dir.create(dirname(filepath), recursive = T)
  }
  
  
  graphics.off()
  
  png(filepath, height = plotheight, width = plotwidth, units = units)
  
  print(plotobject)
  
  graphics.off()
  
  
}






#' Map Celltypes to Clusters using single cell per-cell metadata data.frame
#'
#' @param celltype_cluster data.frame with two columns, celltype (factor) and cluster (factor) 
#' @param cell_min_num integer, default 10, min number of cells for celltype to be considered among the cluster
#' @param singlect_thres numeric 0 to 1, default 0.9. cutoff to assign a single celltype below this assign combined or just "Mixed"
#' @param multict_lothres numeric 0 to 1, default 0.5, cutoff to assign multiple celltypes if none pass singlect_thres, below this just assign "Mixed"
#'
#' @return a vector of length equal to the number of clusters, with the mapping of celltype(s) for each cluster. Can be used with `plyr::mapvalues()` with from = levels(cluster) to assign a new column of metadata.
#' @export
#'
#' @examples
#' \dontrun{
#' #get md from serat obj
#' md <- sobj@meta.data
#' 
#' #prep the input df
#' celltype_cluster <- md[,c('Celltype', 'seurat_clusters')]
#' 
#' #make sure both are factorized
#' celltype_cluster[,1] <- factor(celltype_cluster[,1])
#' celltype_cluster[,2] <- factor(celltype_cluster[,2])
#' 
#' # run it
#' celltype_cluster_mappedlevs <- celltype_cluster_map(celltype_cluster)
#' 
#' }
celltype_cluster_map <- function(celltype_cluster, 
                                 cell_min_num = 10,
                                 singlect_thres = 0.9,
                                 multict_lothres = 0.5
){
  
  
  
  
  #loop thru clusters and get prop of cells
  clusts <- levels(celltype_cluster[,2])
  ctv_l <- lapply(clusts, function(cl){
    
    #get this cluster metadata
    ccc <- celltype_cluster[celltype_cluster[,2] == cl,]
    
    #get table vector
    ctv <- table(ccc[,1])
    
    #sort, remove zeros
    ctv <- ctv[ctv>0]
    ctv <- sort(ctv,decreasing = T)
    
    #cluster must have at least 10 cells?
    # celltype must contribute at lest cell_min_num cells
    ctv_cut <- ctv[ctv>cell_min_num]
    
    #get prop
    ctv <- ctv/sum(ctv)
    
    #filter pfull props by min cell cutoff
    ctv <- ctv[names(ctv) %in% names(ctv_cut)]
    
    #return ctv so we can see how it looks...
    
    return(ctv)
    
    
  })
  names(ctv_l) <- paste0('Cluster_', clusts)
  
  #based on prop of cells, try to assign celltype
  # if first is above 0.9 (singlect_thres), assign 
  # if first is below 0.9, try to assign any with more than 0.1 (multict_lothres)
  ctmapping <- sapply(clusts, function(cl){
    clustlab <- paste0('Cluster_', cl)
    
    ctv <- ctv_l[[clustlab]]
    
    #if above single cutoff thres, just pick it as celltype
    if(ctv[1] >= singlect_thres){
      cm <- names(ctv)[1]
    }
    
    #if below thres, use a combination of celltypes
    
    if(ctv[1] < singlect_thres){
      
      ctv <- ctv[ctv>multict_lothres]
      cm <- names(ctv)
      cm <- paste(cm, collapse = '_')
      
      #sometimes, no cell types pass 50%, just call it mixed
      if( cm == '' ){cm = 'Mixed'}
    }
    
    cm <- paste0(clustlab, '--', cm)
    
    
    
  })
  
  return(ctmapping)
  
}




## code for initiating parallel threads, but maybe this is not really necessary
# initiate_parallel_threads <- function(workernum = 1,
#                                       verbose=T){
#   
#   
#   require(foreach)
#   require(doParallel)
#   require(parallel)
#   
#   
#   
#   ## initiate parallelization
#   cl <- parallel::makeCluster(workernum)
#   doParallel::registerDoParallel(cl)
#   
#   
#   
#   message('\nParallelization intitiated with workernum = ', workernum,
#           '\nTo stop it, please run:\n',
#           '"parallel::stopCluster(cl)"'
#           )
#   
#   
#   return(cl)
#   
# }
