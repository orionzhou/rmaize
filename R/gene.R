#' plot kmers
#'
#' @export
plot_kmer <- function(gt, gid, cid, srd, pb, pe) {
    #{{{
    if(is.na(gid)) {
        return(ggplot()+otheme(panel.border=F,xtick=F,xtext=F,panel.spacing=0,margin=c(0,0,0,0)))
    }
    #{{{ prepare
    ggid = glue("{gt}_{gid}")
    ti = tibble(gid=ggid,st=1)
    write_tsv(ti, 'tmp.tsv', col_names=F)
    fm = glue("{dird}/41_ml/00_nf/03_motif_lists/{cid}.tsv")
    cmd = glue("kmer.py prepare_ml tmp.tsv {fm} ",
               "--bin 'TSS:-/+2k,TTS:-/+2k' ",
               "--epi umr --nfea top100 --mod anr ",
               "--fmt long out.tsv")
    system(cmd)
    tp = read_tsv('out.tsv')
    if (nrow(tp)==0) 
        return(ggplot()+otheme(panel.border=F,xtick=F,xtext=F,panel.spacing=0,margin=c(0,0,0,0)))
    tp = tp %>% mutate(pos = round(start+end)) %>%
        select(fid,pos,srd)
    if(srd == '-') tp = tp %>% mutate(pos = pe - pos)
    system("rm tmp.tsv out.tsv")
    #}}}
    ggplot(tp) +
        #geom_dotplot(aes(x=pos), binwidth=(pe-pb)/100, fill="white", stroke=.1) +
        geom_jitter(aes(x=pos,y=0), height=1, width=50, size=.2) +
        scale_x_continuous(limits=c(pb,pe),expand=expansion(mult=c(.01,.01))) +
        scale_y_continuous(expand=expansion(mult=c(.1,.1))) +
        otheme(panel.border=F, xtick=F,xtext=F,margin=c(0,0,0,0))
    #}}}
}

#' obtain txdb gene tibble
#'
#' @export
txdb_gene <- function(txdb, gid, gb=0) {
    #{{{
    if(is.na(gid)) return(NA)
    ti = AnnotationDbi::select(txdb, keys=gid, columns=columns(txdb), keytype='GENEID') %>%
        rename(gid=GENEID, tid=TXNAME, chrom=TXCHROM, tb=TXSTART, te=TXEND,
            srd=TXSTRAND, eb=EXONSTART, ee=EXONEND, cb=CDSSTART, ce=CDSEND)
    ti1 = ti %>% distinct(gid,tid,beg=tb,end=te,srd) %>% mutate(type='rna')
    ti2 = ti %>% filter(!is.na(ee)) %>%
        select(gid,tid,beg=eb,end=ee,srd) %>% mutate(type='exon')
    ti3 = ti %>% filter(!is.na(cb)) %>%
        select(gid,tid,beg=cb,end=ce,srd) %>% mutate(type='cds')
    rbind(ti1,ti2,ti3) %>%
        mutate(beg=beg-gb, end=end-gb)
    #}}}
}

#' plot genes using tibble
#'
#' @export
plot_genes <- function(tg, pb, pe, n_arrows=10, ht.exon=.2, ht.cds=.4, x=T, label.top=T) {
    #{{{
    if(length(tg)==0 || is.na(tg)) {
        p = ggplot() +
        otheme(panel.border=F, xtick=F,xtext=F,panel.spacing=0,margin=c(0,0,0,0))
        return(p)
    }
    #{{{ prepare
    if (! 'y' %in% colnames(tg)) {
        tg0 = tg %>% filter(type=='rna') %>% arrange(gid) %>% mutate(y = 1:n())
        tg = tg %>% inner_join(tg0 %>% dplyr::select(tid,y), by='tid')
    }
    ir = function(xs) tibble(b=xs[-length(xs)], e = xs[-1])
    tgr = tg %>% filter(type=='rna') %>% mutate(pos=(beg+end)/2)
    tgx = tgr %>% mutate(size=end-beg) %>% arrange(desc(size))
    brks = c(tgx$beg[[1]], tgx$end[[1]])
    labs = c("TSS", "TTS")
    if (tgx$srd[[1]] == '-')  labs = c("TTS", "TSS")
    tga = tgr %>% mutate(x=map2(beg,end,seq,length.out=n_arrows)) %>%
        mutate(x1 = map(x, ir)) %>%
        dplyr::select(gid,tid,srd,y,x1) %>% unnest(x1)
    tga1 = tga %>% filter(srd=='+')
    tga2 = tga %>% filter(srd=='-')
    tge = tg %>% filter(type=='exon')
    tgc = tg %>% filter(type=='cds')
    #}}}
    #{{{
    col.exon = 'royalblue'; col.cds = 'royalblue'
    col.syn = 'grey'
    wd.arrow1 = .1; wd.arrow2 = .08
    arrow11 = arrow(length=unit(wd.arrow1,'cm'), angle=30, ends='last',type="open")
    arrow12 = arrow(length=unit(wd.arrow1,'cm'), angle=30, ends='first',type="open")
    arrow21 = arrow(length=unit(wd.arrow2,'cm'), angle=30, ends='last',type="open")
    arrow22 = arrow(length=unit(wd.arrow2,'cm'), angle=30, ends='first',type="open")
    #}}}
    p = ggplot() +
        geom_segment(data=tgr,aes(x=beg,xend=end,y=y,yend=y),size=.5) +
        geom_segment(data=tga1,aes(x=e-5,xend=e,y=y,yend=y),
                     color='black', size=1, arrow=arrow11) +
        geom_segment(data=tga2,aes(x=b,xend=b+5,y=y,yend=y),
                     color='black', size=1, arrow=arrow12) +
        geom_rect(data=tge,aes(xmin=beg,xmax=end,ymin=y-ht.exon,ymax=y+ht.exon),fill=col.exon,color=NA,alpha=1) +
        geom_rect(data=tgc,aes(xmin=beg,xmax=end,ymin=y-ht.cds,ymax=y+ht.cds),fill=col.cds,color=NA,alpha=1) +
        geom_segment(data=tga1,aes(x=e-2,xend=e-1,y=y,yend=y),
                     color='white', size=.3, arrow=arrow21) +
        geom_segment(data=tga2,aes(x=b+1,xend=b+2,y=y,yend=y),
                     color='white', size=.3, arrow=arrow22) +
        scale_x_continuous(label=label_number(accuracy=1),
                           #breaks=brks, labels=labs,
                           limits=c(pb,pe),expand=expansion(mult=c(.01,.01))) +
        #scale_y_continuous(expand=expansion(mult=c(.01,.01))) +
        scale_y_continuous(expand=expansion(add=c(0,.4),mult=c(0,0))) +
        otheme(panel.border=F, xtick=T,xtext=T,margin=c(0,0,0,0)) +
        theme(axis.text.x = element_text(size=6))
    if (label.top) {
        p + geom_text(data=tgr,aes(x=pos,y=y+ht.cds*1.1,label=tid), vjust=0, size=2)
    } else {
        p + geom_text(data=tgr,aes(x=pos,y=y-ht.cds*1.1,label=tid), vjust=1, size=2)
    }
    #}}}
}

#' prepare synteny plots
#'
#' @export
prepare_syn <- function(ta,qgid,tgid,qtss,ttss) {
    #{{{
    if (is.na(qgid))
        return(tibble(srd='*', aln=list(tibble())))
    ta1 = ta %>% filter(id==tgid)
    o1 = ttss - ta1$b1; o2 = qtss - ta1$b2
    srd = ta1$srd[[1]]
    aln = ta1$aln[[1]] %>%
        mutate(rb1=rb1-o1, re1=re1-o1, rb2=rb2-o2, re2=re2-o2)
    tibble(srd=srd, aln=list(aln))
    #}}}
}

#' plot synteny
#'
#' @export
plot_syn <- function(ta, srd, pb, pe, qtop=F,size=.2) {
    #{{{
    if(is.na(ta) || nrow(ta)==0) {
        p = ggplot() +
        otheme(panel.border=F, xtick=F,xtext=F,panel.spacing=0,margin=c(0,0,0,0))
        return(p)
    }
    #{{{
    ta = ta %>% filter(rb1<pe,re1>pb, rb2<pe, re2>pb) %>%
        mutate(ob1 = pmax(0, pb-rb1)) %>%
        mutate(oe1 = pmax(0, re1-pe)) %>%
        mutate(rb1 = rb1 + ob1, re1 = re1 - oe1) %>%
        #mutate(rb2 = rb2 + ob1, re2 = re2 - oe1) %>%
        mutate(ob2 = pmax(0, pb-rb2)) %>%
        mutate(oe2 = pmax(0, re2-pe)) %>%
        #mutate(rb1 = rb1 + ob2, re1 = re1 - oe2) %>%
        mutate(rb2 = rb2 + ob2, re2 = re2 - oe2)
    if(qtop) {
        tay = tibble(xt=c('rb1','re1','rb2','re2'), y=c(-1,-1,1,1))
    } else {
        tay = tibble(xt=c('rb1','re1','rb2','re2'), y=c(1,1,-1,-1))
    }
    if(srd == "-") {
        tay = tay %>% mutate(xo=c(1,2,1,2))
    } else {
        tay = tay %>% mutate(xo=c(1,2,2,1))
    }
    tp = ta %>% mutate(i=1:n()) %>% gather(xt, x, -i) %>%
        inner_join(tay, by='xt') %>%
        arrange(i,y,xo)
    #}}}
    col.syn='coral'
    ggplot() +
        geom_polygon(data=tp, aes(x=x,y=y,group=i), fill=col.syn, alpha=.7,
                     color=ifelse(size>0,'black',NA),size=size) +
        scale_x_continuous(limits=c(pb,pe),expand=expansion(mult=c(.01,.01))) +
        scale_y_continuous(limits=c(-1,1), expand=expansion(mult=c(.0,.0))) +
        otheme(panel.border=F, xtick=F,xtext=F,panel.spacing=0,margin=c(0,0,0,0))
    #}}}
}

#' plot synteny combo plots
#'
#' @export
plot_combo <- function(gid,cid,cond,drc,st,fo, tss, cmps, gts, off=2000) {
    #{{{
    #drc = ifelse(str_detect(st, "\\+"), "+", "-")
    tit = glue("{gid} {cond} {drc} (B|M|W) {st}")
    cat(tit,'\n')
    #tss %>% filter(gid==!!gid) %>% print(width=Inf)
    cfg = tibble(gt=gts) %>% mutate(gid=!!gid) %>%
        mutate(gt=factor(gt,levels=gts)) %>% arrange(gt) %>%
        left_join(tss, by=c('gt','gid')) %>%
        select(gt,gid,gid2,chrom,gb,ge,tss,srd) %>%
        mutate(pb=gb-tss-off, pe=ge-tss+off) %>%
        mutate(pb = min(pb,na.rm=T)) %>%
        mutate(pe = max(pe,na.rm=T)) %>%
        mutate(pk = pmap(list(gt,gid,cid,srd,pb,pe), plot_kmer)) %>%
        inner_join(txdbs, by='gt') %>%
        mutate(tg = pmap(list(txdb, gid2, tss), txdb_gene))
    #
    o = cfg %>%
        mutate(pg = pmap(list(tg,pb,pe), plot_genes, n_arrows=10, ht.exon=.2, ht.cds=.4))
    #
    cmp1 = cmps %>%
        inner_join(cfg %>% select(qry=gt,qgid=gid2,qtss=tss),by='qry') %>%
        inner_join(cfg %>% select(tgt=gt,tgid=gid2,ttss=tss,pb,pe),by='tgt')
    oc = cmp1 %>%
        mutate(x = pmap(list(ta,qgid,tgid,qtss,ttss), prepare_syn)) %>%
        unnest(x) %>% mutate(qtop=c(T,F)) %>%
        mutate(pc = pmap(list(aln,srd,pb,pe,qtop), plot_syn, size=0))
    #
    p = o$pg[[1]] + o$pk[[1]] +
        oc$pc[[1]] +
    o$pg[[2]] + o$pk[[2]] +
        oc$pc[[2]] +
    o$pg[[3]] + o$pk[[3]] +
        plot_layout(ncol=1, heights=c(1,1,1.5, 1,1,1.5, 1,1)) +
        plot_annotation(title=tit,theme=theme(plot.title=element_text(size=8)))
    ggsave(p, filename=fo, width=5, height=3)
    p
    #}}}
}
