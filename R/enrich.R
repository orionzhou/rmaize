#{{{ old interprosan GO output
if (FALSE) {
    fg = '/home/springer/zhoux379/data/genome/Zmays_v4/61.interpro/15.tsv'
    tg = read.table(fg, sep = "\t", as.is = T, header = F, quote = '')
    colnames(tg) = c("gid", "goid")

    fd = '/home/springer/zhoux379/data/genome/Zmays_v4/61.interpro/16.go.tsv'
    #td = read.table(fd, sep = "\t", header = T, as.is = T, quote = '')

    fi = '/home/springer/zhoux379/data/genome/Zmays_v4/61.interpro/gids.txt'
    ti = read.table(fi, as.is = T)
    gids.all = ti$V1
}
#}}}

#' Read GO annotation for B73 AGP_v4
#'
#' @export
read_go <- function(src = 'uniprot.plants', fgo = '~/projects/genome/data2/Zmays_B73/61_functional/01.go.tsv') {
    #{{{
    ti = read_tsv(fgo)
    if(src == 'all')
        ti
    else if(src %in% ti$ctag)
        ti %>% filter(ctag == src) %>% select(-ctag)
    else
        stop(sprintf("unknown source: %s\n", src))
    #unique(tgo$ctag)
    #gids_all = unique(tgo$gid)
    #tgo_ipr = tgo %>% filter(ctag == 'Interproscan5') %>% select(-ctag)
    #tgo_uni = tgo %>% filter(ctag == 'uniprot.plants') %>% select(-ctag)
    #tgo_ath = tgo %>% filter(ctag == 'arabidopsis') %>% select(-ctag)
    #tgo_arg = tgo %>% filter(ctag == 'argot2.5') %>% select(-ctag)
    #list(go=tgo, ipr=tgo_ipr, uni=tgo_uni, ath=tgo_ath, arg=tgo_arg)
    #}}}
}

#' Perform GO enrichment analysis using hypergeometric test
#'
#' @export
go_enrich <- function(gids, tgo) {
    #{{{
    tgn = tgo %>% distinct(goid, goname, gotype, level)
    tgs = tgo %>% count(goid) %>% select(goid, hitInPop = n)
    #
    gids_all = tgo %>% distinct(gid) %>% pull(gid)
    gids = unique(gids)
    gids = gids[gids %in% gids_all]
    sampleSize = length(gids)
    #
    tz = tgo %>% filter(gid %in% gids) %>%
        count(goid) %>%
        select(goid, hitInSample = n)
    #
    tz %>% inner_join(tgs, by = 'goid') %>%
        mutate(sampleSize = length(gids), popSize = length(gids_all),
            pval.raw = phyper(hitInSample-1, hitInPop, popSize-hitInPop, sampleSize, lower.tail = F)) %>%
            #pval.adj = p.adjust(pval.raw, method = "BH")) %>%
        #filter(pval.raw < 0.05) %>%
        mutate(ratioInSample = sprintf("%d/%d", hitInSample, sampleSize)) %>%
        mutate(ratioInPop = sprintf("%d/%d", hitInPop, popSize)) %>%
        select(goid, ratioInSample, ratioInPop, pval.raw) %>%
        left_join(tgn, by = 'goid') %>% arrange(pval.raw)
    #}}}
}

#' Perform hypergeometic enrichment test on a set of genes
#'
#' @export
hyper_enrich <- function(gids, tgrp) { # grp, gid, note
    #{{{
    tgn = tgrp %>% distinct(grp, note)
    tgs = tgrp %>% count(grp) %>% select(grp, hitInPop = n)
    #
    gids_all = tgrp %>% distinct(gid) %>% pull(gid)
    gids = unique(gids)
    gids = gids[gids %in% gids_all]
    sampleSize = length(gids)
    #
    tz = tgrp %>% filter(gid %in% gids) %>%
        count(grp) %>%
        select(grp, hitInSample = n)
    #
    tz %>% inner_join(tgs, by = 'grp') %>%
        mutate(sampleSize = length(gids), popSize = length(gids_all),
            pval.raw = phyper(hitInSample-1, hitInPop, popSize-hitInPop, sampleSize, lower.tail=F)) %>%
            #pval.adj = p.adjust(pval.raw, method = "BH")) %>%
        #filter(pval.raw < 0.05) %>%
        mutate(ratioInSample = sprintf("%d/%d", hitInSample, sampleSize)) %>%
        mutate(ratioInPop = sprintf("%d/%d", hitInPop, popSize)) %>%
        select(grp, ratioInSample, ratioInPop, pval.raw) %>%
        left_join(tgn, by = 'grp') %>% arrange(pval.raw)
    #}}}
}

#' Perform hypergeometric enrichment test on multiple gene sets
#'
#' @export
hyper_enrich_m <- function(tg, tgrp, method='BH', pval.cutoff=.05) { # group + gid
    #{{{
    tg = tg %>% filter(gid %in% tgrp$gid)
    tgs = tg %>% count(group) %>% filter(n >= 3) %>% select(-n)
    tg = tg %>% filter(group %in% tgs$group)
    tg %>% group_by(group) %>%
        summarise(ng = length(gid), gids = list(gid)) %>% ungroup() %>%
        mutate(data = map(gids, hyper_enrich, tgrp=tgrp)) %>%
        mutate(n_group = nrow(tgs)) %>%
        #filter(!is.na(data)) %>% filter(!is.null(data)) %>%
        select(n_group, group, ng, data) %>% unnest() %>%
        mutate(pval.adj = p.adjust(pval.raw, method=method)) %>%
        #filter(pval.adj < pval.cutoff) %>%
        arrange(pval.adj)
    #}}}
}

#' Generic function to perform hypergeometric/fisher exact test for enrichment
#'
#' @export
test_enrich <- function(hitInSample, sampleSize, hitInPop, popSize,
                        deplete=F, method='phyper') {
    #{{{
    if(method == 'phyper') {
        if(deplete) {
            phyper(hitInSample, hitInPop, popSize-hitInPop, sampleSize, lower.tail=T)
        } else {
            phyper(hitInSample-1, hitInPop, popSize-hitInPop, sampleSize, lower.tail=F)
        }
    } else if (method == 'fet') {
        mat = matrix(c(hitInSample, hitInPop-hitInSample, sampleSize-hitInSample,
                       popSize-hitInPop-sampleSize+hitInSample), 2, 2)
        if(deplete) {
            fisher.test(mat, alternative='less')$p.value
        } else {
            fisher.test(mat, alternative='greater')$p.value
        }
    } else {
        stop("unknown test method:", method)
    }
    #}}}
}



