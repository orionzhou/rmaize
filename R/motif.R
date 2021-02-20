#' rename a motif
#'
#' @export
rename_mtf <- function(mtf, name) { mtf['name'] = name; mtf }


#' cluster a set of motifs two pass
#'
#' @export
cluster_motifs <- function(ti, min_sites=5, cutHeight1=.1, cutHeight2=.1) {
    #{{{
    #require(dynamicTreeCut)
    # trimming
    cat("input: ", nrow(ti), " motifs\n")
    trim_motifs2 <- function(mtf)
        tryCatch(trim_motifs(mtf), error = function(c) NULL)
    ti2 = ti %>%
        mutate(mtf = map(mtf, trim_motifs2)) %>%
        mutate(j = map_dbl(mtf, length)) %>% filter(j>0) %>% select(-j) %>%
        mutate(conseq=map_chr(mtf, 'consensus')) %>%
        mutate(icscore = map_dbl(mtf, 'icscore')) %>%
        mutate(nsites=map_int(conseq, nchar)) %>%
        filter(nsites >= min_sites) %>% select(-nsites) %>%
        arrange(conseq, desc(icscore))
    cat("after trimming: ", nrow(ti2), " motifs\n")
    tm = ti2 %>% group_by(conseq) %>% slice(1) %>%
        ungroup() %>% select(mid, conseq, icscore, mtf)
    cat("after consensus merging: ", nrow(tm), " clusters\n")
    #
#{{{ # 1-pass clustering
mtfs = tm$mtf
cmp0 = compare_motifs(mtfs, method="PCC", min.mean.ic=.0,
                      min.overlap=5, score.strat="a.mean")
dst0 = as.dist(1 - cmp0)
hcu0 = hclust(dst0, method='average')
#
x = cutree(hcu0, h=cutHeight1)
tx0 = tibble(mid = hcu0$labels, grp = as.integer(x))
tx0$grp[tx0$grp==0] = max(tx0$grp) + 1:sum(tx0$grp==0)
tx = tx0 %>%
    inner_join(tm, by='mid') %>%
    arrange(grp, desc(icscore)) %>%
    select(grp, mid, icscore, mtf)# %>% print(n=30)
#
txs = tx %>% arrange(grp, desc(icscore)) %>%
    group_by(grp) %>%
    summarise(mids = list(mid)) %>% ungroup()
#tx %>% count(grp) %>% arrange(desc(grp))
cat("1st pass clustering: ", nrow(txs), " clusters\n")
#}}}
    # tx txs
    #
#{{{ # 2-pass clustering
tx2 = tx %>% arrange(grp, desc(icscore)) %>%
    group_by(grp) %>%
    summarise(mid=mid[1], mtf=mtf[1]) %>% ungroup()
#
cmp1 = compare_motifs(tx2$mtf, method="PCC", min.mean.ic=.0,
                      min.overlap=6, score.strat="a.mean")
dst = as.dist(1 - cmp1)
hcu = hclust(dst, method='average')
#
x = cutree(hcu, h=cutHeight2)
tx3 = tibble(mid = hcu$labels, grp2 = as.integer(x)) %>%
    inner_join(tx2, by='mid') %>% select(grp, grp2)
#
tr = tx %>% inner_join(tx3, by='grp') %>%
    select(-grp) %>% rename(grp=grp2) %>%
    arrange(grp, desc(icscore)) %>% select(mid, grp)
tr %>% count(grp)
cat("2nd pass clustering: ", length(unique(tr$grp)), " clusters\n")
#}}}
    tr %>% inner_join(tm, by='mid') %>% select(conseq, grp) %>%
        inner_join(ti2, by='conseq') %>%
        select(mid, icscore, grp)
    #}}}
}

