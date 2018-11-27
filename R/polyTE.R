dirr = '~/git/rmaize/R'
source(file.path(dirr, 'header.R'))
source(file.path(dirr, 'plot.R'))
source(file.path(dirr, 'genomes.R'))

dirp = '~/projects/polyTE'
dird = file.path(dirp, 'data')

make_tile <- function(cbeg, cend, winsize = 10, winstep = 5) {
    size = cend - cbeg + 1
    nwin = ceil((size - winsize) / winstep) + 1
    tibble(beg = cbeg + winstep * (1:nwin-1),
           end = cbeg + winstep * (1:nwin-1) + winsize - 1) %>%
        mutate(beg = as.integer(beg), end = as.integer(pmin(end, cend)))
}


