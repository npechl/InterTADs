#
#' Title
#'
#' @param entrez.ids
#' @param tech
#'
#' @importFrom AnnotationDbi select
#' @importFrom annotables grch37 grch38
#' 
#' @return
#' @exports
#'
#' @examples


map_entrez_ids <- function(entrez.ids, tech = "hg19"){

    genes <- entrez.ids

    mapping <- select(
        org.Hs.eg.db,
        key = genes,
        columns = c("SYMBOL"),
        keytype = "ENTREZID"
    )

    mapping <- as.data.table(mapping)
    colnames(mapping) <- c("entrez.id", "hgnc.symbol")

    not.found <- mapping[which(is.na(mapping$hgnc.symbol)), ]$entrez.id
    mapping <- mapping[which(!(is.na(mapping$hgnc.symbol))), ]

    if(length(not.found) > 0){
        
        if(tech == "h19"){
            
            x <- grch37
            
        } else {
            
            x <- grch38
            
        }

        x <- as.data.table(x)
        x <- x[,c(2,3)]
        colnames(x) <- c("entrez.id", "hgnc.symbol")

        merge <- x[which(x$entrez.id %in% not.found), ]

        mapping <- rbind(mapping, merge)

        not.found <- not.found[which(!(not.found %in% merge$entrez.id))]
      }

    mapping <- rbind(
        mapping, 
        data.table(
            "entrez.id" = not.found,
            "hgnc.symbol" = NA
        )
    )
    
    mapping <- mapping[order(mapping$hgnc.symbol), ]

    mapping <- unique(mapping)

    return(mapping)
}









