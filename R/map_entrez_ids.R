#' map_entrez_ids
#'
#' @param entrez.ids
#' @param tech Human Genome Reference used
#'
#' @importFrom AnnotationDbi select
#' @import  annotables
#'
#' @return
#' @export
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
    colnames(mapping) <- c("entrez.ids", "hgnc.symbol")

    not.found <- mapping[which(is.na(mapping$hgnc.symbol)), ]$entrez.ids
    mapping <- mapping[which(!(is.na(mapping$hgnc.symbol))), ]

    if(length(not.found) > 0){

        if(tech == "hg19"){

            x <- annotables::grch37

        } else {

            x <- annotables::grch38

        }

        x <- as.data.table(x)
        x <- x[,c(2,3)]
        colnames(x) <- c("entrez.ids", "hgnc.symbol")

        merge <- x[which(x$entrez.ids %in% not.found), ]

        mapping <- rbind(mapping, merge)

        not.found <- not.found[which(!(not.found %in% merge$entrez.ids))]
      }

    mapping <- rbind(
        mapping,
        data.table(
            "entrez.ids" = not.found,
            "hgnc.symbol" = NA
        )
    )

    mapping <- mapping[order(mapping$hgnc.symbol), ]

    mapping <- unique(mapping)

    return(mapping)
}









