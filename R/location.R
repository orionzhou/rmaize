options(scipen=999)

#' make sliding window/tile
#'
#' @export
make_tile <- function(cbeg, cend, winsize = 10, winstep = 5) {
    #{{{
    size = cend - cbeg + 1
    nwin = ceil((size - winsize) / winstep) + 1
    tibble(beg = cbeg + winstep * (1:nwin-1),
           end = cbeg + winstep * (1:nwin-1) + winsize - 1) %>%
        mutate(beg = as.integer(beg), end = as.integer(pmin(end, cend)))
    #}}}
}

#' bedtools intersect
#'
#' @export
intersect_basepair <- function(gr1, gr2) {
    #{{{
    t1 = data.frame(chr = seqnames(gr1), beg = start(gr1)-1, end = end(gr1),
    idx = 1:length(gr1), stringsAsFactors = F)
    gr2 = reduce(gr2)
    t2 = data.frame(chr = seqnames(gr2), beg = start(gr2)-1, end = end(gr2),
    idx = 1:length(gr2), stringsAsFactors = F)

    fbd1 = 'xtest1.bed'
    fbd2 = 'xtest2.bed'
    fres = 'xout.bed'
    options(scipen = 999)
    write.table(t1, fbd1, sep = "\t", row.names = F, col.names = F, quote = F)
    write.table(t2, fbd2, sep = "\t", row.names = F, col.names = F, quote = F)
    options(scipen = 0)
    system(sprintf("intersectBed -wao -a %s -b %s > %s", fbd1, fbd2, fres))

    t3 = read.table(fres, sep = "\t", header = F, as.is = T)
    colnames(t3) = c('chr', 'beg', 'end', 'idx', 'qchr', 'qbeg', 'qend', 'qidx', 'olen')
    system(sprintf("rm %s %s %s", fbd1, fbd2, fres))

    gp = group_by(t3, idx)
    t4 = dplyr::summarise(gp, olen = sum(olen))
    as.numeric(t4$olen)
    #}}}
}

#' bedtools intersect
#'
#' @export
intersect_score <- function(gr1, gr2) {
    #{{{
    t1 = data.frame(chr = seqnames(gr1), beg = start(gr1)-1, end = end(gr1),
    idx = 1:length(gr1), stringsAsFactors = F)
    t2 = data.frame(chr = seqnames(gr2), beg = start(gr2)-1, end = end(gr2),
    idx = 1:length(gr2), score = mcols(gr2)$score, stringsAsFactors = F)

    fbd1 = 'xtest1.bed'
    fbd2 = 'xtest2.bed'
    fres = 'xout.bed'
    options(scipen = 999)
    write.table(t1, fbd1, sep = "\t", row.names = F, col.names = F, quote = F)
    write.table(t2, fbd2, sep = "\t", row.names = F, col.names = F, quote = F)
    options(scipen = 0)
    system(sprintf("intersectBed -wao -a %s -b %s > %s", fbd1, fbd2, fres))

    t3 = read.table(fres, sep = "\t", header = F, as.is = T)
    colnames(t3) = c('chr', 'beg', 'end', 'idx', 'qchr', 'qbeg', 'qend', 'qidx', 'score', 'olen')
    t3$score[t3$olen == 0] = NA
    system(sprintf("rm %s %s %s", fbd1, fbd2, fres))

    gp = group_by(t3, idx)
    t4 = dplyr::summarise(gp, score = sum(score, na.rm = T))
    as.numeric(t4$score)
    #}}}
}

#' bedtools intersect
#'
#' @export
intersect_count <- function(gr1, gr2) {
    #{{{
    t1 = data.frame(chr = seqnames(gr1), beg = start(gr1)-1, end = end(gr1),
    idx = 1:length(gr1), stringsAsFactors = F)
    t2 = data.frame(chr = seqnames(gr2), beg = start(gr2)-1, end = end(gr2),
    idx = 1:length(gr2), stringsAsFactors = F)

    fbd1 = 'xtest1.bed'
    fbd2 = 'xtest2.bed'
    fres = 'xout.bed'
    options(scipen = 999)
    write.table(t1, fbd1, sep = "\t", row.names = F, col.names = F, quote = F)
    write.table(t2, fbd2, sep = "\t", row.names = F, col.names = F, quote = F)
    options(scipen = 0)
    system(sprintf("intersectBed -wao -a %s -b %s > %s", fbd1, fbd2, fres))

    t3 = read.table(fres, sep = "\t", header = F, as.is = T)
    colnames(t3) = c('chr', 'beg', 'end', 'idx', 'qchr', 'qbeg', 'qend', 'qidx', 'olen')
    system(sprintf("rm %s %s %s", fbd1, fbd2, fres))

    gp = group_by(t3, idx)
    t4 = dplyr::summarise(gp, cnt = sum(qidx != '.'))
    as.numeric(t4$cnt)
    #}}}
}

#' bedtools intersect
#'
#' @export
intersect_idx <- function(gr1, gr2) {
    #{{{
    t1 = data.frame(chr = seqnames(gr1), beg = start(gr1)-1, end = end(gr1),
    idx = 1:length(gr1), stringsAsFactors = F)
    t2 = data.frame(chr = seqnames(gr2), beg = start(gr2)-1, end = end(gr2),
    idx = 1:length(gr2), stringsAsFactors = F)

    fbd1 = 'xtest1.bed'
    fbd2 = 'xtest2.bed'
    fres = 'xout.bed'
    options(scipen = 999)
    write.table(t1, fbd1, sep = "\t", row.names = F, col.names = F, quote = F)
    write.table(t2, fbd2, sep = "\t", row.names = F, col.names = F, quote = F)
    options(scipen = 0)
    system(sprintf("intersectBed -wao -a %s -b %s > %s", fbd1, fbd2, fres))

    t3 = read.table(fres, sep = "\t", header = F, as.is = T)
    colnames(t3) = c('chr', 'beg', 'end', 'idx', 'qchr', 'qbeg', 'qend', 'qidx', 'olen')
    system(sprintf("rm %s %s %s", fbd1, fbd2, fres))

    t3[t3$qidx != '.', c('chr','beg','end','idx','qidx','olen')]
    #}}}
}

#' bedtools intersect
#'
#' @export
intersect_idx_maxovlp <- function(gr1, gr2) {
    #{{{
    t1 = data.frame(chr = seqnames(gr1), beg = start(gr1)-1, end = end(gr1),
    idx = 1:length(gr1), stringsAsFactors = F)
    t2 = data.frame(chr = seqnames(gr2), beg = start(gr2)-1, end = end(gr2),
    idx = 1:length(gr2), stringsAsFactors = F)

    fbd1 = 'xtest1.bed'
    fbd2 = 'xtest2.bed'
    fres = 'xout.bed'
    options(scipen = 999)
    write.table(t1, fbd1, sep = "\t", row.names = F, col.names = F, quote = F)
    write.table(t2, fbd2, sep = "\t", row.names = F, col.names = F, quote = F)
    options(scipen = 0)
    system(sprintf("intersectBed -wao -a %s -b %s > %s", fbd1, fbd2, fres))

    t3 = read.table(fres, sep = "\t", header = F, as.is = T)
    colnames(t3) = c('chr', 'beg', 'end', 'idx', 'qchr', 'qbeg', 'qend', 'qidx', 'olen')
    system(sprintf("rm %s %s %s", fbd1, fbd2, fres))

    gp = group_by(t3, idx)
    t4 = dplyr::summarise(gp, qidx = qidx[which(olen == max(olen))[1]])
    as.numeric(t4$qidx)
    #}}}
}

#' bedtools intersect
#'
#' @export
intersect_idx_t <- function(gr1, gr2) {
    #{{{
    t1 = data.frame(chr = seqnames(gr1), beg = start(gr1)-1, end = end(gr1),
    idx = 1:length(gr1), stringsAsFactors = F)
    t2 = data.frame(chr = seqnames(gr2), beg = start(gr2)-1, end = end(gr2),
    idx = 1:length(gr2), stringsAsFactors = F)

    fbd1 = 'xtest1.bed'
    fbd2 = 'xtest2.bed'
    fres = 'xout.bed'
    options(scipen = 999)
    write.table(t1, fbd1, sep = "\t", row.names = F, col.names = F, quote = F)
    write.table(t2, fbd2, sep = "\t", row.names = F, col.names = F, quote = F)
    options(scipen = 0)
    system(sprintf("intersectBed -wo -a %s -b %s > %s", fbd1, fbd2, fres))

    t3 = read.table(fres, sep = "\t", header = F, as.is = T)
    colnames(t3) = c('chr', 'beg', 'end', 'idx', 'qchr', 'qbeg', 'qend', 'qidx', 'olen')
    system(sprintf("rm %s %s %s", fbd1, fbd2, fres))

    as.numeric(t3$idx)
    #}}}
}

#' bedtools intersect
#'
#' @export
intersect_cds_sv <- function(gr1, gr2, cdsidx) {
    #{{{
    t1 = data.frame(chr = seqnames(gr1), beg = start(gr1)-1, end = end(gr1),
    idx = cdsidx, stringsAsFactors = F)
    t2 = data.frame(chr = seqnames(gr2), beg = start(gr2)-1, end = end(gr2),
    idx = 1:length(gr2), stringsAsFactors = F)

    fbd1 = 'xtest1.bed'
    fbd2 = 'xtest2.bed'
    fres = 'xout.bed'
    options(scipen = 999)
    write.table(t1, fbd1, sep = "\t", row.names = F, col.names = F, quote = F)
    write.table(t2, fbd2, sep = "\t", row.names = F, col.names = F, quote = F)
    options(scipen = 0)
    system(sprintf("intersectBed -wao -a %s -b %s > %s", fbd1, fbd2, fres))

    t3 = read.table(fres, sep = "\t", header = F, as.is = T)
    colnames(t3) = c('chr', 'beg', 'end', 'idx', 'qchr', 'qbeg', 'qend', 'qidx', 'olen')
    system(sprintf("rm %s %s %s", fbd1, fbd2, fres))

    gp = group_by(t3, idx, qidx)
    t4 = dplyr::summarise(gp, olen = sum(olen))
    gp2 = group_by(t4, idx)
    t5 = dplyr::summarise(gp2, qidx = qidx[which(olen == max(olen))][1], olen = sum(olen))
    t5
    #}}}
}

#' bedtools intersect
#'
#' @export
intersect_gene <- function(gr1, gr2, cdsidx) {
    #{{{
    t1 = data.frame(chr = seqnames(gr1), beg = start(gr1)-1, end = end(gr1),
    idx = 1:length(gr1), stringsAsFactors = F)
    t2 = data.frame(chr = seqnames(gr2), beg = start(gr2)-1, end = end(gr2),
    idx = cdsidx, stringsAsFactors = F)

    fbd1 = 'xtest1.bed'
    fbd2 = 'xtest2.bed'
    fres = 'xout.bed'
    options(scipen = 999)
    write.table(t1, fbd1, sep = "\t", row.names = F, col.names = F, quote = F)
    write.table(t2, fbd2, sep = "\t", row.names = F, col.names = F, quote = F)
    options(scipen = 0)
    system(sprintf("intersectBed -wo -a %s -b %s > %s", fbd1, fbd2, fres))

    t3 = read.table(fres, sep = "\t", header = F, as.is = T, quote = "")
    colnames(t3) = c('chr', 'beg', 'end', 'idx', 'qchr', 'qbeg', 'qend', 'qidx', 'olen')
    system(sprintf("rm %s %s %s", fbd1, fbd2, fres))

    gp = group_by(t3, idx, qidx)
    t4 = dplyr::summarise(gp, olen = sum(olen))
    t4
    #}}}
}

