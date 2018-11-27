dirr = '~/git/rmaize/R'
source(file.path(dirr, 'header.R'))
source(file.path(dirr, 'plot.R'))

dirg = '~/data/genome'

attach_genome <- function(genome) {
#    if(! genome %in% search()) {
        attach(file.path(dirg, genome, '55.rda'), name = genome)
#    } else {
#        cat(sprintf("already loaded: %s\n", genome))
#    }
}

detach_genome <- function(genome) {
    detach(genome, character.only = T)
}
