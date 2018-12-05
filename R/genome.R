#' Get system path to genome directory
#'
#' @export
genome_dir <- function(genome='', root='~/data/genome') {
    file.path(root, genome)
}

#' Read Genome Configuration
#'
#' @export
read_genome_conf <- function(genome, dirg = '~/data/genome') {
    fi = file.path(dirg, genome, '55.rda')
    if(!file.exists(fi)) {
        stop(sprintf("cannot read %s\n", fi))
    } else {
        readRDS(fi)
    }
}
