#' read project meta table
#'
#' @export
read_projects <- function(genome='Zmays_B73', diri = '~/projects/barn') {
    #{{{
    gdic = c('Athaliana' = 'arabidopsis', 'Osativa' = 'rice', 'Zmays_B73' = 'maize')
    stopifnot(genome %in% names(gdic))
    org = gdic[genome]
    read_xlsx(glue("{diri}/{org}.xlsx"))
    #}}}
}

#' get RNA-Seq sample meta table
#'
#' @export
rnaseq_sample_meta <- function(yid='rn18g',dird='~/projects/rnaseq/data') {
    #{{{
    fh1 = sprintf("%s/05_meta/%s.tsv", dird, yid)
    #fh2 = sprintf("%s/06_meta_final/%s.c.tsv", dird, mid)
    #fh = ifelse(file.exists(fh2), fh2, fh1)
    read_tsv(fh1)
    #}}}
}

#' get RNA-Seq mapping stats
#'
#' @export
rnaseq_mapping_stat <- function(yid, dird='~/projects/rnaseq/data') {
    #{{{ mapping stats
    th = rnaseq_sample_meta(yid) %>% select(SampleID,paired)
    diri = file.path(dird, 'raw', yid)
    fi = file.path(diri, 'trimming.tsv')
    tt1 = read_tsv(fi) %>% rename(SampleID = sid) %>%
        separate(SampleID, c('SampleID','suf'), sep='[\\_]', extra='merge', fill='right')
    if('passed_filter_reads' %in% colnames(tt1)) {
        tt1 = tt1 %>% group_by(SampleID) %>%
            summarise(passed = sum(passed_filter_reads),
                      lowQ = sum(low_quality_reads),
                      manyN = sum(too_many_N_reads),
                      short = sum(too_short_reads),
                      long = sum(too_long_reads)) %>% ungroup() %>%
            mutate(failed = lowQ + manyN + short + long) %>%
            inner_join(th, by='SampleID') %>%
            mutate(passed = ifelse(paired, passed/2, passed)) %>%
            mutate(failed = ifelse(paired, failed/2, failed)) %>%
            mutate(total = passed+failed) %>%
            select(SampleID, total, passed, failed)
    } else if('readsIn' %in% colnames(tt1)) {
        tt1 = tt1 %>% rename(total=readsIn, passed=readsOut, failed=readsRemoved) %>%
            group_by(SampleID) %>%
            summarise(total=sum(total), passed=sum(passed), failed=sum(failed)) %>%
            select(SampleID, total, passed, failed)
    } else {
        stop("cannot detect triming.tsv format")
    }
    fi = file.path(diri, 'bamstats.tsv')
    tt2 = read_tsv(fi) %>% rename(SampleID=sid)
    #
    tt = tt1 %>%
        left_join(tt2, by = 'SampleID') %>%
        mutate(mapped=pair_map+pair_orphan+unpair_map,
               mappedu=pair_map_hq+pair_orphan_hq+unpair_map_hq) %>%
        mutate(p.passed=passed/total, p.mapped=mapped/total, p.mappedu=mappedu/total) %>%
        mutate(total=total/1e6, passed=passed/1e6, mapped=mapped/1e6, mappedu=mappedu/1e6) %>%
        select(SampleID,total,passed,mapped,mappedu,p.passed,p.mapped,p.mappedu)
    tt
    #}}}
}

#' read one RNA-Seq study result (raw)
#'
#' @export
rnaseq_cpm_raw <- function(yid, genome='Zmays_B73', diri='~/projects/rnaseq/data/11_qc') {
    #{{{
    fi = sprintf("%s/%s/%s/00.raw.rds", diri, genome, yid)
    cat(fi,'\n')
    stopifnot(file.exists(fi))
    readRDS(fi)
    #}}}
}

#' read one RNA-Seq study result (final)
#'
#' @export
rnaseq_cpm <- function(yid, diri='~/projects/rnaseq/data/11_qc', genome='Zmays_B73') {
    #{{{
    fi = sprintf("%s/%s/%s/01.rds", diri, genome, yid)
    stopifnot(file.exists(fi))
    readRDS(fi)
    #}}}
}

#' filter genes from an RNA-Seq dataset
#'
#' @export
filter_expr <- function(ti, wide=F, min_cpm=1, num_sam_on=0, pct_sam_on=0,
                        min_var_p=0, transform='no') {
    #{{{
    if (wide) ti = ti %>% gather(sid, val, -gid)
    cat(glue("{length(unique(ti$gid))} genes in {length(unique(ti$sid))} samples read\n"), "\n")
    #
    sd2 <- function(x) sd(x[!(is.na(x) | is.infinite(x) | is.nan(x))])
    tis = ti %>% group_by(gid) %>%
        summarise(nsam_on = sum(abs(val) >= min_cpm),
                  psam_on = nsam_on/n(),
                  val_sd = sd2(val)) %>%
        ungroup()
    tis = tis %>% filter(nsam_on >= num_sam_on)
    tis = tis %>% filter(psam_on >= pct_sam_on)
    tis = tis %>% filter(val_sd > 0)
    gids = tis %>% pull(gid)
    #
    tis2 = tis %>% filter(gid %in% gids)
    min_sd = quantile(tis2$val_sd, min_var_p)
    gids = tis2 %>% filter(val_sd >= as.numeric(min_sd)) %>% pull(gid)
    #
    to = ti %>% filter(gid %in% gids)
    #
    if (transform == 'asinh') {
        to = to %>% mutate(val = asinh(val))
    } else if (transform == 'log2') {
        to = to %>% mutate(val = log2(val))
    } else {
        stopifnot(transform == 'no')
    }
    to = to %>% spread(sid, val)
    cat(glue("{nrow(to)} genes passed filtering\n"), "\n")
    to
    #}}}
}

#' normalize RNA-Seq read count matrix
#'
#' @export
readcount_norm <- function(t_rc, t_gs = F) {
    #{{{ normalize
    smMap = t_rc %>% distinct(SampleID) %>% dplyr::rename(oSampleID=SampleID) %>%
        mutate(nSampleID = str_replace_all(oSampleID, '[^a-zA-z0-9_]', '.'))
    tm = t_rc
    if (! identical(smMap$oSampleID, smMap$nSampleID))
        tm = tm %>% inner_join(smMap, by=c('SampleID'='oSampleID')) %>%
            select(-SampleID) %>% select(SampleID=nSampleID, everything())
    tw = tm %>%
        select(SampleID, gid, ReadCount) %>%
        spread(SampleID, ReadCount) %>%
        replace(., is.na(.), 0)
    tm = tw %>% gather(SampleID, ReadCount, -gid)
    gids = tw$gid
    twd = data.frame(tw[,-1])
    rownames(twd) = tw$gid
    # nRC with DESeq2
    require(DESeq2)
    th = tm %>% distinct(SampleID) %>% arrange(SampleID)
    th2 = th %>% mutate(sid = SampleID, SampleID = factor(SampleID))
    thd = column_to_rownames(as.data.frame(th2), var = 'sid')
    stopifnot(identical(rownames(thd), colnames(twd)))
    dds = DESeqDataSetFromMatrix(countData=twd, colData=thd, design = ~ 1)
    dds = estimateSizeFactors(dds)
    sf = sizeFactors(dds)
    t_sf = tibble(SampleID = names(sf), sizeFactor = as.numeric(sf))
    t_nrc = counts(dds, normalized = T) %>% as_tibble() %>%
        mutate(gid = names(dds)) %>% gather(SampleID, nRC, -gid)
    # rCPM and CPM with edgeR
    require(edgeR)
    y = DGEList(counts = twd)
    y = calcNormFactors(y, method = 'TMM')
    t_nf = y$samples %>% as_tibble() %>%
        mutate(SampleID = rownames(y$samples)) %>%
        select(SampleID, libSize=lib.size, normFactor=norm.factors)
    t_cpm1 = cpm(y) %>% as_tibble() %>% mutate(gid = rownames(cpm(y))) %>%
        select(gid, everything()) %>%
        gather(SampleID, CPM, -gid)
    t_cpm2 = cpm(y, normalized = F) %>% as_tibble() %>%
        mutate(gid = rownames(cpm(y))) %>%
        gather(SampleID, rCPM, -gid)
    t_cpm = t_cpm1 %>% inner_join(t_cpm2, by=c('SampleID','gid'))
    # rFPKM & FPKM
    if(is.list(t_gs)) {
        t_tpm = tm %>% left_join(t_gs, by='gid') %>%
            mutate(RPK = ReadCount / (size/1000)) %>%
            group_by(SampleID) %>% mutate(rTPM = RPK / sum(RPK) * 1e6) %>%
            ungroup() %>% inner_join(t_nf, by = 'SampleID') %>%
            mutate(TPM = rTPM / normFactor) %>%
            select(SampleID, gid, TPM, rTPM)
        t_cpm = t_cpm %>% left_join(t_gs, by = 'gid') %>%
            mutate(FPKM = CPM / (size / 1000)) %>%
            mutate(rFPKM = rCPM / (size / 1000)) %>%
            select(-size) %>%
            inner_join(t_tpm, by=c('SampleID','gid'))
    }
    #
    tl = th %>% inner_join(t_sf, by = 'SampleID') %>%
        inner_join(t_nf, by = 'SampleID')
    tm = tm %>%
        left_join(t_nrc, by = c('SampleID','gid')) %>%
        left_join(t_cpm, by = c('SampleID','gid'))
    stopifnot(nrow(tm) == length(gids) * nrow(th))
    if (! identical(smMap$oSampleID, smMap$nSampleID)) {
        tl = tl %>% inner_join(smMap, by=c('SampleID'='nSampleID')) %>%
            dplyr::select(-SampleID) %>% dplyr::select(SampleID=oSampleID, everything())
        tm = tm %>% inner_join(smMap, by=c('SampleID'='nSampleID')) %>%
            dplyr::select(-SampleID) %>% dplyr::select(SampleID=oSampleID, everything())
    }
    list(tl = tl, tm = tm, dds=dds)
    #}}}
}

#' merge replicates
#'
#' @export
merge_reps <- function(th, tm, sid) {
    #{{{
    if("sid" %in% colnames(th)) {
        ths = th %>% distinct(sid, Tissue, Genotype, Treatment) %>%
            mutate(nSampleID = sprintf("%s_%d", !!sid, 1:length(Tissue)))
        th = th %>% inner_join(ths, by = c("sid", "Tissue", "Genotype", "Treatment"))
        t_map = th %>% select(SampleID, nSampleID)
        th = ths %>% select(SampleID=nSampleID, sid, Tissue, Genotype, Treatment)
    } else {
        ths = th %>% distinct(Tissue, Genotype, Treatment) %>%
            mutate(nSampleID = sprintf("%s_%d", !!sid, 1:length(Tissue)))
        th = th %>% inner_join(ths, by = c("Tissue", "Genotype", "Treatment"))
        t_map = th %>% select(SampleID, nSampleID)
        th = ths %>% select(SampleID=nSampleID, Tissue, Genotype, Treatment)
    }
    #
    tm = tm %>% inner_join(t_map, by = 'SampleID') %>%
        mutate(SampleID = nSampleID) %>%
        group_by(gid, SampleID) %>%
        #summarise(ReadCount = sum(ReadCount), nRC = sum(nRC), rCPM = mean(rCPM),
        #          rFPKM = mean(rFPKM), CPM = mean(CPM), FPKM = mean(FPKM)) %>%
        summarise(CPM = mean(CPM), FPKM = mean(FPKM)) %>%
        ungroup()
    list(th = th, tm = tm)
    #}}}
}

#' plot hclust
#'
#' @export
plot_hclust <- function(tm, th, min.value=1, pct.exp=.5,
                        cor.opt='pearson', hc.opt='ward.D',
                        var.lab='lab', var.col='Genotype',
                        lab.size=2.5, expand.x=.2, pal.col='aaas') {
    #{{{ hclust
    tw = tm %>% select(SampleID, gid, value) %>% spread(SampleID, value)
    t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(value>=min.value))
    gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * pct.exp) %>% pull(gid)
    e = tw %>% filter(gid %in% gids) %>% select(-gid)
    dim(e)
    #
    hc_title = sprintf("dist: %s\nhclust: %s", cor.opt, hc.opt)
    edist <- as.dist(1-cor(e, method = cor.opt))
    ehc <- hclust(edist, method = hc.opt)
    tree = as.phylo(ehc)
    lnames = ehc$labels[ehc$order]
    #
    tp = th %>% mutate(taxa = SampleID) %>%
        select(taxa, everything())
    p1 = ggtree(tree, layout = 'rectangular') +
        scale_x_continuous(expand = expansion(mult=c(1e-2,expand.x))) +
        scale_y_discrete(expand = expansion(mult=c(.01,.01)))
    p1 = p1 %<+%
        tp + geom_tiplab(aes(label=get(var.lab), color=get(var.col)), size=lab.size) +
        get(str_c('scale_color', pal.col, sep="_"))() +
        guides(color=F)
    p1
    #}}}
}

#' plot PCA
#'
#' @export
plot_pca <- function(tm, th, min.value=1, pct.exp=.5, pca.center=F, pca.scale=F, ...) {
    #{{{
    tw = tm %>% select(SampleID, gid, value) %>%
        filter(SampleID %in% th$SampleID) %>% spread(SampleID, value)
    t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(value>=min.value))
    gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * pct.exp) %>% pull(gid)
    #
    e = tw %>% filter(gid %in% gids) %>% select(-gid)
    dim(e)
    pca <- prcomp(e, center = pca.center, scale. = pca.scale)
    tx = pca['rotation'][[1]]
    imp = summary(pca)$importance
    #
    xlab = sprintf("PC1 (%.01f%%)", imp[2,1]*100)
    ylab = sprintf("PC2 (%.01f%%)", imp[2,2]*100)
    tp = as_tibble(tx[,1:5]) %>%
        add_column(SampleID = rownames(tx)) %>%
        rename(x=PC1, y=PC2) %>%
        inner_join(th, by = 'SampleID')
    plot_clustering_2d(tp, xtitle=xlab, ytitle=ylab, ...)
    #}}}
}

#' plot tsne
#'
#' @export
plot_tsne <- function(tm, th, min.value=1, pct.exp=.5, perp=6, iter=1200,
                      seed=42, ...) {
    #{{{ tSNE
    require(Rtsne)
    tw = tm %>% select(SampleID, gid, value) %>%
        filter(SampleID %in% th$SampleID) %>% spread(SampleID, value)
    t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(value>=min.value))
    gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * pct.exp) %>% pull(gid)
    tt = tw %>% filter(gid %in% gids)
    dim(tt)
    tsne <- Rtsne(t(as.matrix(tt[-1])), dims=2, verbose=T, perplexity=perp,
                  pca = T, max_iter = iter, check_duplicates = F)
    #
    tp = tsne$Y %>% as.data.frame() %>% as_tibble() %>%
        rename(x=V1, y=V2) %>%
        add_column(SampleID = colnames(tt)[-1]) %>%
        inner_join(th, by = 'SampleID')
    plot_clustering_2d(tp, xtitle='tSNE-1', ytitle='tSNE-2', ...)
    #}}}
}

#' plot umap
#'
#' @export
plot_umap <- function(tm, th, min.value=1, pct.exp=.5,
                      seed=42, ...) {
    #{{{ UMAP
    require(umap)
    tw = tm %>% select(SampleID, gid, value) %>%
        filter(SampleID %in% th$SampleID) %>% spread(SampleID, value)
    t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(value>=min.value))
    gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * pct.exp) %>% pull(gid)
    tt = tw %>% filter(gid %in% gids)
    res <- umap(t(as.matrix(tt[-1])), random_state = seed)
    #
    tp = res$layout %>% as.data.frame() %>% as_tibble() %>%
        rename(x=1, y=2) %>%
        add_column(SampleID = colnames(tt)[-1]) %>%
        inner_join(th, by = 'SampleID')
    plot_clustering_2d(tp, xtitle='UMAP-1', ytitle='UMAP-2', ...)
    #}}}
}


plot_clustering_2d <- function(tp, xtitle='PC1', ytitle='PC2',
    var.lab='lab', var.col='', var.shape='', var.ellipse='',
    shapes = 0:10, leg.col=T, leg.shape=T,
    legend.pos='top.right', legend.dir='v', legend.box='v', legend.title=T,
    point.size=2, lab.size=2, pal.col='aaas') {
    #{{{
    x.max=max(tp$x)
    p = ggplot(tp, aes(x,y))
    # point shape and color aesthetics
    if (var.shape == '' & var.col == '')
        p = p + geom_point(size=point.size)
    else if(var.shape == '' & var.col != '')
        p = p + geom_point(aes(color=get(var.col)), size=point.size) +
            get(str_c('scale_color', pal.col, sep="_"))(name=var.col)
    else if(var.shape != '' & var.col == '')
        p = p + geom_point(aes(shape=get(var.shape)), size=point.size) +
            scale_shape_manual(name=var.shape, values = shapes)
    else
        p = p + geom_point(aes(color=get(var.col), shape=get(var.shape)), size=point.size) +
            get(str_c('scale_color', pal.col, sep="_"))(name=var.col) +
            scale_shape_manual(name=var.shape, values = shapes)
    # ellipse
    if (var.ellipse != '')
        p = p + geom_mark_ellipse(aes(fill=get(var.ellipse)),
            expand=unit(2,'mm'), alpha=0, size = .2,
            con.type='none',label.fontsize=0,label.minwidth=unit(0,'mm'),
            label.buffer=unit(0,'mm'),label.margin = margin(0,0,0,0,"mm"))
    # point label
    if (var.lab != '')
        p = p + geom_text_repel(aes(label=get(var.lab)), size=lab.size)
    p = p +
        scale_x_continuous(name = xtitle) +
        scale_y_continuous(name = ytitle) +
        otheme(legend.pos=legend.pos, legend.dir=legend.dir,
               legend.box=legend.box, legend.title=legend.title,
               xtitle=T, ytitle=T,
               margin = c(.2,.2,.2,.2)) +
        theme(axis.ticks.length = unit(0, 'lines')) +
        guides(fill=F)
    if (!leg.col)
        p = p + guides(color=F)
    else
        p = p# + guides(color=guide_legend(nrow=1, byrow=T))
    if (!leg.shape)
        p = p + guides(shape=F)
    else
        p = p# + guides(shape=guide_legend(nrow=1, byrow=T))
    p
    #}}}
}

#' plot ASE allele-frequency boxplot
#'
#' @export
plot_ase <- function(ta, th, min_rc=20, drc='h',
                     val.col='black', val.fill='white', pal.col='aaas') {
    #{{{
    tp = ta %>% filter(allele1 + allele2 >= min_rc) %>%
        mutate(af = allele1/(allele1 + allele2)) %>%
        inner_join(th, by='SampleID') %>%
        mutate(lab = factor(lab, levels=rev(th$lab)))
    tps = tp %>% group_by(lab) %>%
        summarise(n=n(),af_avg = mean(af), af_med=median(af)) %>%
        mutate(labn=sprintf("%s (%s)", lab, number(n, big.mark = ",", accuracy=1))) %>%
        ungroup() %>% arrange(lab)
    #tp %>% group_by(lab) %>%
        #summarise(q50=median(af), m50=sum(allele1)/sum(allele1+allele2)) %>%
        #ungroup() %>% print(n=70)
    ort = ifelse(drc=='r', 'reverse', ifelse(drc=='v', 'vertical', 'horizontal'))
    p = ggboxplot(tp, x='lab', y='af', color=val.col, fill=val.fill,
                  orientation=ort, palette = pal.col,
                  outlier.shape = NA, bxp.errorbar = T) +
    #p = ggplot(tp) +
        #geom_violin(aes(x=lab, y=af, fill=get(val.fill)), color=val.col,
                    #trim=F) +
        #geom_segment(aes(x=match(lab, levels(lab))-0.1,
                     #xend=match(lab, levels(lab))+0.1,
                     #y=af, yend=af), color='white')+
        geom_hline(yintercept = .5, color='black',linetype='dashed') +
        scale_x_discrete(breaks=tps$lab, labels=tps$labn) +
        scale_y_continuous(name='% allele1',expand=expansion(mult=c(.02,.02))) +
        #coord_flip() +
        otheme(xtext=T, ytext=T, xtick=T, ytick=T, xtitle=T) +
        guides(color=F, fill=F)
    p
    #}}}
}

#' plot RIL haplotype blocks
#'
#' @export
plot_ril_genotype <- function(cp, th, sids_red='', gts=c('B73','Mo17','het')) {
#{{{ RIL genotype
    fw = '~/projects/genome/data/Zmays_B73/15_intervals/20.win11.tsv'
    tw = read_tsv(fw)
    offs = c(0, cumsum(tw$size)[-nrow(tw)]) + 0:10 * 10e6
    tx = tw %>% mutate(off = offs) %>%
        mutate(gstart=start+off, gend=end+off, gpos=(gstart+gend)/2) %>%
        select(rid,chrom,gstart,gend,gpos,off) %>% filter(chrom!='B99')
    #
    tp = cp %>% inner_join(tx, by='rid') %>%
        mutate(start=start+off, end=end+off) %>%
        inner_join(th, by=c('sid'='SampleID'))
    tz = cp %>% mutate(size=end-start) %>% group_by(sid, gt) %>%
        summarise(size = sum(size)) %>%
        mutate(total_size = sum(size)) %>%
        mutate(prop = size / total_size) %>% ungroup() %>%
        select(sid, gt, prop) %>% spread(gt, prop) %>%
        replace_na(list(a=0,b=0,h=0))
    ty = tp %>% distinct(sid, lab) %>% arrange(desc(lab)) %>%
        mutate(y = 1:n()) %>%
        mutate(col = ifelse(sid %in% sids_red, 'red','black')) %>%
        inner_join(tz, by='sid') %>%
        mutate(lab=sprintf("%s (%.1f%%)", lab, a*100))
    tp = tp %>% inner_join(ty, by=c('sid'))
    #
    xmax = max(tp$end)
    p = ggplot(tp) +
        geom_rect(aes(xmin=start,xmax=end,ymin=y-.3,ymax=y+.3, fill=gt)) +
        #geom_text(data=ty, aes(x=xmax+5e6, y=y, label=lab), hjust=0, size=2.5) +
        scale_x_continuous(breaks=tx$gpos, labels=tx$chrom, expand=expansion(mult=c(.01,.01))) +
        scale_y_continuous(breaks=ty$y, labels=ty$lab, expand=expansion(mult=c(.002,.002))) +
        scale_fill_manual(values=pal_simpsons()(8)[c(1,2,5)], labels=gts) +
        otheme(legend.pos='top.center.out', legend.dir='h',
               xtext=T, xtick=T, ytext=T, ytick=T) +
        theme(axis.text.y = element_text(color=ty$col, size=7))
#}}}
}



