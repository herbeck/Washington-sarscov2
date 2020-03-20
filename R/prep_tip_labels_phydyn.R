##This is copied and modified from Erik Volz's phydyn SEIJR Coronavirus package
##Public repo here: https://github.com/emvolz/weifang-sarscov2
##I did not write this script.

#' Change the name of sequences to recognize numeric time of sampling and deme 
#'
#' @param path_to_algn A DNAbin alignment *or* a system path (type character) where original alignment can be found
#' @param deme The name of the deme to add to the end of each label 
#' @param regexprs A vector of regular expressions that can be used to match tip labels and categorize in to demes. These must include *all* sequence IDs, and each ID should match at most one regular expression
#' @param invert_regexpr A logical vector specifying if the i'th regular expression should be an inverse match 
#' @param demes A character vector of deme names to be appended to corresponding sequence IDs 
#'
#" @return DNAbin alignment 
#' @export 
#~ prep_tip_labels_phydyn <- function( path_to_align, path_to_save = NULL, deme = 'Il'  ){

prep_tip_labels_phydyn <- function( path_to_align, path_to_save = NULL
                                    , regexprs = c( '.*/WA/.*', '.*/WA/.*' ) 
                                    , invert_regexpr = c( FALSE, TRUE )
                                    , demes = c( 'Il'  , 'exog'  )
){

  library( ape ) 
  library( treedater )
  library( lubridate )
  
  if ( inherits( path_to_align, 'DNAbin' ) )
    d = path_to_align
  else
    d = read.dna( path_to_align, format = 'fasta')
  
  sids = rownames(d) 
  if ( length( regexprs ) != length( demes ))
    stop('Must provide equal numbers of regex and deme names ') 
  
  demegroups = lapply( 1:length(demes), function(k) {
    x = regexprs[k]
    if ( invert_regexpr[k] ) {
      return( sids[ !grepl( pattern = x , sids ) ] )
    }else{
      return( sids[ grepl( pattern = x , sids ) ] )
    }
  })
  int <- do.call( intersect, demegroups )
  uni = do.call( c, demegroups )
  if ( length( int ) > 0 )
    stop( 'Intersection of deme groups is non-empty. Each regex must match a unique set.' )
  if ( length( uni ) < length(sids) ){
    print( setdiff( sids, uni ))
    stop( 'There were some sequence IDs that did not match a regex. ' )
  }
  
  deme <- setNames( rep(demes[1], nrow(d)), sids )
  for ( k in 1:length( demes )){
    deme[ demegroups[[k]] ] <- demes[k]
  }
  
  sts <- sapply( strsplit( rownames(d), '\\|' ) , function(x){
    round(decimal_date( ymd( tail(x,1))), digits = 3)
  })
  rownames(d) <- paste(sep='|', rownames(d), sts, paste0('_', deme) )
  rownames(d) <- gsub( rownames(d), pattern = '\\s' , replacement = '_')
  row.names(d) <- gsub( rownames(d), pattern = "hCoV-19/", replacement = '')
  if ( !is.null( path_to_save ))
    write.dna( d, file = path_to_save, format = 'fasta' )
  d
}
