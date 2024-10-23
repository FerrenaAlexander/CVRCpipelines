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

