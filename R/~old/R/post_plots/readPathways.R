#' Parse GMT file and return pathways as list
#'
#' @param fname - (character) path to pathway file in gmt format
#'	pathway score to include pathway in the filter list
#' @param MIN_SIZE, MAX_SIZE: (integer) - min/max num. genes allowed in a pathway.
#'	Pathways with gene counts outside these bounds are excluded
#' @param EXCLUDE_KEGG: (boolean) - if TRUE exclude pathways from KEGG
#' @param IDasName: (boolean) - if FALSE, uses pathway name as key; if TRUE
#'	uses db name and ID as name (e.g.  KEGG:hsa04940)
#' @return List with pathway name as key, vector of genes as value
#' @export
readPathways <- function(fname,MIN_SIZE=10, MAX_SIZE=500, EXCLUDE_KEGG=TRUE,
						 IDasName=FALSE,verbose=TRUE)
{

# change locale to accommodate nonstandard chars in pathway names
oldLocale   <- Sys.getlocale("LC_ALL")
Sys.setlocale("LC_ALL","C")

    out <- list()
    # read list of master pathways
	if (verbose) cat("---------------------------------------\n")
	if (verbose) cat(sprintf("File: %s\n\n", basename(fname)))
    f	<- file(fname,"r")
    # TODO: deal with duplicate pathway names

    #pName <- list()
	ctr <- 0
	options(warn=1)
    repeat {
        s	<- scan(f, what="character",nlines=1,quiet=TRUE,sep="\t")
        if (length(s)==0) break;

        pPos<- gregexpr("%",s[1])[[1]];
		src <- ""
		src_id	<- ""
        if (pPos[1]==-1) {
            #cat("\n\n% symbol not found in pathway name")
			s[1]	<- s[1]
        } else {

			src		<- substr(s[1],pPos[1]+1,pPos[2]-1)
			src_id	<- substr(s[1],pPos[2]+1,nchar(s[1]))
			if (IDasName) {
				s[1]	<- paste(src,src_id,sep=":")
			} else {
        		s[1]	<- substr(s[1],1,pPos[1]-1)
			}
		}
		if (!EXCLUDE_KEGG || (src!="KEGG")) {
			idx <- which(s=="") # remove trailing blank rows.
			if (any(idx)) s <- s[-idx]
        	out[[s[1]]]	<- s[3:length(s)]
        	#pName[[s[1]]] <- s[2] # stores pathway source - prob not needed
		}
		ctr <- ctr+1
    }
    close(f)
	if (verbose) {
		cat(sprintf("Read %i pathways in total, internal list has %i entries\n",

				ctr, length(out)))
    	cat(sprintf("\tFILTER: sets with num genes in [%i, %i]\n",MIN_SIZE,
				MAX_SIZE))
	}
    ln <- unlist(lapply(out, length))
    idx	<- which(ln < MIN_SIZE | ln >= MAX_SIZE)
    out[idx] <- NULL
    #pName[idx] <- NULL
    if (verbose) cat(sprintf("\t  => %i pathways excluded\n\t  => %i left\n",
				length(idx),length(out)))

	# clean pathway names
	# names(out)	<- cleanPathwayName(names(out))

#Sys.setlocale(oldLocale)


    return(out)
}
