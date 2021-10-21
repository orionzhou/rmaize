#' Generate a shape palette (<=15)
#'
#' @export
pal_shapes <- function(n=15) {
    #{{{
    shapes = c(15:18, 0:4, 7:13)
    if(n <= 15) {
        shapes[1:n]
    } else {
        stop(sprintf("error: >15 shapes not allowed\n"))
    }
    #}}}
}

#' Generate a color palette (n=23)
#'
#' @export
pal_23 <- function() c(pal_ucscgb()(23)[c(1:11,13:23)], pal_uchicago()(6))

#' Generate a color palette (n=23)
#'
#' @export
pal_23b <- function() c(pal_igv()(23))

#' Generate a color palette (n=13)
#'
#' @export
pal_13 <- function() c(pal_igv()(13))

#' Generate a color palette (n=12)
#'
#' @export
pal_12 <- function() pal_futurama()(12)

#' Generate a color palette (n=10)
#'
#' @export
pal_10 <- function() pal_aaas()(10)

#' Generate a gradient color palette (n=100)
#'
#' @export
pal_gradient<- function(opt = 1, reverse = F) {
    #{{{
    if(opt == 1) {
        cols_raw = brewer.pal(n = 7, name = "RdYlBu")
    } else {
        stop("unknown opt", opt)
    }
    if(reverse) cols_raw = rev(cols_raw)
    colorRampPalette(cols_raw)(100)
    #}}}
}

#' generate point density
#'
#' @export
get_density <- function(x, y, ...) {
    #{{{
    require(MASS)
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
    #}}}
}

#' use hclust to order matrix rows
#'
#' @export
hc_order_row <- function(ti, cor.opt='euclidean', hc.opt='ward.D') {
    #{{{ gid followed by 2+ value columns
    ti = ti %>% rename(gid = 1)
    e = ti %>% select(-gid) %>% as.data.frame()
    rownames(e) = ti$gid
    if(nrow(ti) <= 1) {
        ti$gid
    } else {
        edist = idist(e, opt='row', method=cor.opt)
        ehc = hclust(edist, method = hc.opt)
        if (cor.opt == 'gower') {
            labs = rownames(e)
        } else {
            labs = ehc$labels
        }
        labs[ehc$order]
    }
    #}}}
}

#' flexible custom ggplot theme
#'
#' @export
otheme <- function(margin = c(.5,.5,.5,.5),
                   strip.size = 8, strip.style = 'white', strip.compact = F,
                   legend.pos='', legend.dir='v', legend.box='v', legend.vjust=NA,
                   legend.title = F, legend.border = F,
                   legend.spacing.x = 0, legend.spacing.y = 0,
                   panel.border = T, panel.spacing = .02,
                   xticks = F, yticks = F, xtitle = F, ytitle = F,
                   xtext = F, ytext = F, xgrid = F, ygrid = F,
                   xsize = 8, ysize = 8, nostrip = T) {
    #{{{ custom theme
    o = theme_bw()
    #{{{ strip style
    strip.color = ifelse(strip.style=='dark', 'white', 'black')
    strip.fill = ifelse(strip.style=='dark', '#5D729D', 'white')
    strip.box.color = ifelse(strip.style=='dark', '#4A618C', '#4A618C')
    if(strip.style == 'dark') {
        o = o +
            theme(strip.background = element_blank(),
                  strip.text = element_textbox(size = strip.size,
                  color = strip.color, fill = strip.fill, box.color = strip.box.color,
                  halign=.5, linetype = 1, r = unit(5, "pt"), width = unit(1, "npc"),
                  padding = margin(2, 0, 1, 0), margin = margin(.05, .05, .05, .05)))
    } else {
        o = o +
            theme(strip.background = element_blank(),
                  strip.text = element_text(size = strip.size))
        if(strip.compact)
            o = o + theme(strip.text=element_text(margin=margin(.05,.05,.05,.05)))
    }
    #}}}
    #{{{ legend
    o = o +
        theme(legend.background = element_blank(),
              legend.key.size = unit(.8, 'lines'),
              legend.text = element_text(size = 8),
              legend.spacing.x = unit(legend.spacing.x, "lines"),
              legend.spacing.y = unit(legend.spacing.y, "lines"),
              legend.margin = margin(.1,.1,.1,.1,'lines'))
    if(legend.title == F) {
        o = o + theme(legend.title = element_blank())
    } else {
        o = o + theme(legend.title = element_text(size=8))
    }
    if(legend.pos == 'none') {
        o = o + theme(legend.position = 'none')
    } else if(legend.pos == 'top.center.out') {
        if(is.na(legend.vjust)) legend.vjust=0
        o = o + theme(legend.position = c(.5,1), legend.justification = c(.5,legend.vjust))
        margin[1] = margin[1] + 1
    } else if(legend.pos == 'top.center') {
        o = o + theme(legend.position = c(.5,1), legend.justification = c(.5,1))
    } else if(legend.pos == 'top.left') {
        o = o + theme(legend.position = c(0,1), legend.justification = c(0,1))
    } else if(legend.pos == 'top.right') {
        o = o + theme(legend.position = c(1,1), legend.justification = c(1,1))
    } else if(legend.pos == 'bottom.right') {
        o = o + theme(legend.position = c(1,0), legend.justification = c(1,0))
    } else if(legend.pos == 'bottom.left') {
        o = o + theme(legend.position = c(0,0), legend.justification = c(0,0))
    } else if(legend.pos == 'left') {
        o = o + theme(legend.position = 'left')
    }
    if(legend.dir == 'h')
        o = o + theme(legend.direction='horizontal')
    if(legend.box == 'h')
        o = o + theme(legend.box='horizontal')
    if(legend.border == T)
        o = o + theme(legend.box.background = element_rect())
    #}}}
    #{{{ axis tick/text/title
    if(xtitle) {
        o = o + theme(axis.title.x = element_text(size = 9))
    } else {
        o = o + theme(axis.title.x = element_blank())
    }
    if(ytitle) {
        o = o + theme(axis.title.y = element_text(size = 9))
    } else {
        o = o + theme(axis.title.y = element_blank())
    }
    if(xtext) {
        o = o + theme(axis.text.x = element_text(size = xsize))
    } else {
        o = o + theme(axis.text.x = element_blank())
    }
    if(ytext) {
        o = o + theme(axis.text.y = element_text(size = ysize))
    } else {
        o = o + theme(axis.text.y = element_blank())
    }
    if(!xticks)
        o = o + theme(axis.ticks.x = element_blank(), axis.ticks.length.x=unit(0,'pt'))
    if(!yticks)
        o = o + theme(axis.ticks.y = element_blank(), axis.ticks.length.y=unit(0,'pt'))
    #}}}
    #{{{ panel & grid
    o = o + theme(panel.grid.minor = element_blank()) +
        theme(panel.spacing = unit(panel.spacing, "lines"))
    if(!xgrid)
        o = o + theme(panel.grid.major.x = element_blank())
    if(!ygrid)
        o = o + theme(panel.grid.major.y = element_blank())
    if(!panel.border)
        o = o + theme(panel.border = element_blank())
    #}}}
    o + theme(plot.margin = unit(margin, "lines")) +
        theme(title=element_text(margin=margin(0,0,0,0,'cm')))
    #}}}
}

#' remove x axis
#'
#' @export
no_x_axis <- function(title=F,text=F,tick=F) {
    #{{{
    p = theme()
    if(!title) p = p + theme(axis.title.x=element_blank())
    if(!text) p = p + theme(axis.text.x=element_blank())
    if(!tick) p = p + theme(axis.ticks.x=element_blank(), axis.ticks.length.x=unit(0, "pt"))
    p
    #}}}
}

#' remove y axis
#'
#' @export
no_y_axis <- function(title=F,text=F,tick=F) {
    #{{{
    p = theme()
    if(!title) p = p + theme(axis.title.y=element_blank())
    if(!text) p = p + theme(axis.text.y=element_blank())
    if(!tick) p = p + theme(axis.ticks.y=element_blank(), axis.ticks.length.y=unit(0, "pt"))
    p
    #}}}
}

#' margin
#'
#' @export
o_margin <- function(t=.2, r=.2, b=.2, l=.2)
  theme(plot.margin = margin(t, r, b, l, 'lines'))

#' quick plot histogram
#
#' @export
plot_hist <- function(x, fo='~/tmp.pdf', xlab='xlab', ylab='count') {
    #{{{
    p1 = ggplot() +
        geom_histogram(aes(x = x)) +
        otheme_bw(legend.pos='right', legend.dir='v',
                  xtitle=T, ytitle=T, xtext=T, ytext=T)
    ggsave(p1, filename = fo, width = 6, height = 6)
    #}}}
}

#' plot proportions of a set of categories
#
#' @export
cmp_cnt1 <- function(ti, xtitle='', ytitle='', ytext=F, xangle=0,
                     margin.l=.1, ypos='left',
    legend.pos='top.center.out', legend.dir='h', legend.title='', barwidth=.8,
    oneline=F, expand.x=c(.02,.02), expand.y=c(.01,.04), fills=pal_npg()(8)) {
    #{{{
    tags1 = levels(ti$tag1)
    tags2 = levels(ti$tag2)
    sep = ifelse(oneline, " ", "\n")
    col1 = 'black'; col2 = 'black'
    tp = ti %>%
        arrange(tag1, desc(tag2)) %>%
        group_by(tag1) %>% mutate(n_tot = sum(n)) %>%
        mutate(prop = n) %>%
        mutate(y = cumsum(prop)) %>%
        mutate(y = y - prop/2) %>%
        mutate(lab = glue("{number(n,accuracy=1)}")) %>%
        ungroup() %>%
        mutate(tag2 = factor(tag2, levels=tags2))
    tpx = tp %>% group_by(tag1) %>%
        summarise(y=sum(n)*1.01, lab=number(sum(n))) %>% ungroup()
    if (ypos == 'right') { tagm = tpx$tag1[nrow(tpx)] }
    else { tagm = tpx$tag1[1] }
    tpy = tp %>% filter(tag1==tagm)
    pp = ggplot(tp) +
        geom_bar(aes(x=tag1,y=prop,fill=tag2), stat='identity', position='stack', width=barwidth, alpha=.8) +
        geom_text(aes(x=tag1,y=y,label=lab),color=col1,size=2.5,lineheight=.8) +
        geom_text(data=tpx, aes(tag1,y,label=lab), color=col2,size=3, vjust=0) +
        scale_x_discrete(name=xtitle, expand=expansion(mult=expand.x)) +
        scale_y_continuous(name=ytitle, breaks=tpy$y, labels=tpy$tag2, expand=expansion(mult=expand.y), position=ypos) +
        scale_fill_manual(name=legend.title, values=fills) +
        otheme(legend.pos=legend.pos,legend.dir=legend.dir,
               legend.title=legend.title, legend.vjust=.2,
            xtick=T, ytick=!!ytext, xtitle=T, xtext=T, ytext=!!ytext, panel.border=F,
            margin = c(.1, .1, .1, margin.l))
    if(xangle != 0) pp = pp + theme(axis.text.x=element_text(angle=xangle, hjust=1,vjust=1))
    pp
    #}}}
}

#' plot proportions of a set of categories
#
#' @export
cmp_proportion1 <- function(ti, xangle=0, ypos='left', alph=.8, acc=1, barwidth=.8,
    xtitle='', ytitle='', legend.title='',
    lab.size=2.5, oneline=F, expand.x=c(.02,.02), expand.y=c(.01,.04),
    fills=pal_npg()(8), ...) {
    #{{{
    tags1 = levels(ti$tag1)
    tags2 = levels(ti$tag2)
    sep = ifelse(oneline, " ", "\n")
    col1 = 'black'; col2 = 'black'
    tp = ti %>%
        arrange(tag1, desc(tag2)) %>%
        group_by(tag1) %>% mutate(n_tot = sum(n)) %>%
        mutate(prop = n/n_tot) %>%
        mutate(y = cumsum(prop)) %>%
        mutate(y = y - prop/2) %>%
        mutate(lab = glue("{number(n,accuracy=acc)}{sep}({percent(prop,accuracy=1)})")) %>%
        ungroup() %>%
        mutate(tag2 = factor(tag2, levels=tags2))
    tpx = tp %>% group_by(tag1) %>% summarise(n=number(sum(n), accuracy=acc)) %>% ungroup()
    if (ypos == 'right') { tagm = tpx$tag1[nrow(tpx)] }
    else { tagm = tpx$tag1[1] }
    tpy = tp %>% filter(tag1==tagm)
    pp = ggplot(tp) +
        geom_bar(aes(x=tag1,y=prop,fill=tag2), alpha=alph, stat='identity', position='fill', width=barwidth) +
        geom_text(aes(x=tag1,y=y,label=lab),color=col1,size=lab.size,lineheight=.8) +
        geom_text(data=tpx, aes(tag1,1.01,label=n), color=col2,size=lab.size+.5, vjust=0) +
        scale_x_discrete(name=xtitle, expand=expansion(mult=expand.x)) +
        scale_y_continuous(name=ytitle, breaks=tpy$y, labels=tpy$tag2, expand=expansion(mult=expand.y), position=ypos) +
        scale_fill_manual(name=legend.title, values=fills) +
        otheme(...)
    if(xangle != 0) pp = pp + theme(axis.text.x=element_text(angle=xangle, hjust=1,vjust=1))
    pp
    #}}}
}

#' plot proportions of a set of categories
#
#' @export
cmp_proportion <- function(ti, xangle=0, barwidth=.8, 
    alph=.8, oneline=F, lab.size=2.5, prop.only=F, strip.compact=F, acc=1,
    xtitle='', ytitle='', legend.title='',
    nc = 5, expand.x=c(.02,.02), expand.y=c('.02,.02'), fills = pal_npg()(8), ...) {
    #{{{
    tags1 = levels(ti$tag1)
    tags2 = levels(ti$tag2)
    sep = ifelse(oneline, " ", "\n")
    col1 = 'black'; col2 = 'black'
    tp = ti %>%
        arrange(tag1, desc(tag2)) %>%
        group_by(pnl, tag1) %>% mutate(n_tot = sum(n)) %>%
        mutate(prop = n/n_tot) %>%
        mutate(y = cumsum(prop)) %>%
        mutate(y = y - prop/2) %>%
        ungroup() %>%
        mutate(tag2 = factor(tag2, levels=tags2))
    if (prop.only)
        tp = tp %>% mutate(lab = glue("({percent(prop,accuracy=1)})"))
    else
        tp = tp %>% mutate(lab = glue("{number(n,accuracy=acc)}{sep}({percent(prop,accuracy=1)})"))
    tpx = tp %>% group_by(pnl, tag1) %>%
        summarise(n=number(sum(n), accuracy=acc)) %>% ungroup()
    pp = ggplot(tp) +
        geom_bar(aes(x=tag1,y=prop,fill=tag2), alpha=alph, stat='identity', position='fill', width=barwidth) +
        geom_text(aes(x=tag1,y=y,label=lab),color=col1,size=lab.size,lineheight=.8) +
        geom_text(data=tpx, aes(tag1,1.01,label=n), color=col2,size=lab.size+.5, vjust=0) +
        scale_x_discrete(name=xtitle, expand=expansion(mult=expand.x)) +
        scale_y_continuous(name=ytitle, expand=expansion(mult=expand.y)) +
        scale_fill_manual(name=legend.title, values=fills) +
        facet_wrap(~pnl, ncol=nc, scale='free') +
        otheme(...)
    if(xangle != 0) pp = pp + theme(axis.text.x=element_text(angle=xangle, hjust=1,vjust=1))
    pp
    #}}}
}

#' single-panel piechart
#
#' @export
pie1 <- function(ti, fills=pal_npg()(8), lab.size=2.5, nudge.x=0, alph=.6) {
    #{{{
    tags = levels(ti$tag)
   n_tot = sum(ti$n)
    sep = " "
    col1 = 'black'
    tp = ti %>%
        arrange(desc(tag)) %>%
        mutate(prop = n/n_tot) %>%
        mutate(y = cumsum(prop) - prop/2) %>%
        mutate(lab = glue("{tag}\n{number(n,accuracy=1)}{sep}({percent(prop,accuracy=1)})"))
    tit = glue("All ({number(n_tot, accurary=1)})")
    ggplot(tp) +
        geom_bar(aes(x='',y=prop,fill=tag), stat='identity',
                 alpha=alph, color='white') +
        coord_polar("y", start=0) +
        geom_text_repel(aes(x='',y=y,label=lab),color=col1,size=lab.size,lineheight=.8,nudge_x=nudge.x) +
        #scale_x_discrete(name=xtitle, expand=expansion(mult=expand.x)) +
        #scale_y_continuous(name=ytitle, expand=expansion(mult=expand.y)) +
        scale_fill_manual(values=fills) +
        ggtitle(tit) +
        otheme(legend.pos='none',
               xtick=F, ytick=F, xtitle=F, xtext=F, ytext=F, panel.border=F) +
        theme(plot.title=element_text(hjust=.5,size=9))
    #}}}
}

#' multi-panel piechart
#
#' @export
multi_pie <- function(ti, legend.pos='top.center.out', legend.dir='h',
                      strip.compact=F, legend.title='', legend.vjust=.5,
                      alph=.6, ncol=5, fills=pal_npg()(8)) {
    #{{{
    tags1 = levels(ti$tag1)
    tags2 = levels(ti$tag2)
    sep = " "
    col1 = 'black'; col2 = 'black'
    tp = ti %>%
        arrange(tag1, desc(tag2)) %>%
        group_by(tag1) %>% mutate(n_tot = sum(n)) %>%
        mutate(prop = n/n_tot) %>%
        mutate(y = cumsum(prop) - prop/2) %>%
        mutate(lab = glue("{number(n,accuracy=1)}{sep}({percent(prop,accuracy=1)})")) %>%
        ungroup() %>%
        mutate(pnl = glue("{tag1} ({n_tot})"))
    ggplot(tp) +
        geom_bar(aes(x='',y=prop,fill=tag2), stat='identity', alpha=alph, color='white') +
        coord_polar("y", start=0) +
        geom_text_repel(aes(x='',y=y,label=lab),color=col1,size=2.5,lineheight=.8,nudge_x=.3) +
        #scale_x_discrete(name=xtitle, expand=expansion(mult=expand.x)) +
        #scale_y_continuous(name=ytitle, expand=expansion(mult=expand.y)) +
        scale_fill_manual(name=legend.title, values=fills) +
        facet_wrap(~pnl, ncol=ncol) +
        otheme(legend.pos=legend.pos,legend.dir=legend.dir,strip.compact=T,
               legend.title=legend.title,
               legend.vjust=legend.vjust,
               strip.compact=!!strip.compact,
            xtick=F, ytick=F, xtitle=F, xtext=F, ytext=F, panel.border=F)
    #}}}
}

#' plot a heatmap with hclust tree shown on left
#
#' @export
heatmap_hc <- function(ti, top=2,bottom=5, text.size=2.5, ratio=3,
                       r.top = .5,
                       leg='pairwise similarity', tiplab=F) {
    #{{{ tibble(aid, bid, dist, algd, blgd, acol, bcol)
    require(cluster)
    require(ape)
    require(ggtree)
    #{{{ make triangular matrix
    tiw = ti %>% mutate(ul = aid <= bid) %>%
        mutate(aid2 = ifelse(ul, aid, bid)) %>%
        mutate(bid2 = ifelse(ul, bid, aid)) %>%
        select(-aid, -bid) %>%
        rename(aid=aid2, bid=bid2) %>%
        arrange(aid, bid, ul) %>%
        group_by(aid, bid) %>%
        summarise(dist=mean(dist), algd=algd[1], acol=acol[1],
                  blgd=blgd[1], bcol=bcol[1]) %>%
        ungroup() %>% select(aid, bid, dist) %>%
        arrange(aid, bid) %>% spread(bid, dist)
    #}}}
    ids = tiw$aid
    dist_mat = as.matrix(tiw[,-1])
    nums = as.numeric(dist_mat[upper.tri(dist_mat)])
    dist_mat = t(dist_mat)
    dist_mat[upper.tri(dist_mat)] = nums
    rownames(dist_mat) = ids
    #
    hc = hclust(as.dist(dist_mat), method = "ward.D")
    oidx = hc$order
    tree = as.phylo(hc)
    ids_sort = ids[oidx]
    #
    tp = ti %>%
        mutate(aid = factor(aid, levels = ids_sort)) %>%
        mutate(bid = factor(bid, levels = ids_sort)) %>%
        mutate(dist = ifelse(aid==bid, NA, dist)) %>%
        mutate(sim = 1-dist) %>%
        mutate(lab = str_remove(sprintf("%.02f", sim), '^0+')) %>%
        mutate(lab = ifelse(lab=='NA', '', lab))
    tpa = tp %>% distinct(aid, algd, acol) %>% arrange(aid)
    tpb = tp %>% distinct(bid, blgd, bcol) %>% arrange(bid)
    #
    tpt = tpa %>% select(taxa = aid, lgd = algd, col = acol)
    p1 = ggtree(tree, ladderize=F) +
        scale_x_continuous(expand=expand_scale(mult=c(.05,0)))
    if (tiplab)
        p1 = p1 %<+% tpt +
        geom_tiplab(aes(label=lgd, col=col), size=text.size) +
        scale_x_continuous(expand=expand_scale(mult=c(0,.5)))
    p1 = p1 +
        theme(plot.margin = margin(top,0,bottom,0, 'lines'))
    #
    cols = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
    p2 = ggplot(tp) +
        geom_tile(aes(x=aid, y=bid, fill=sim)) +
        geom_text(aes(x=aid, y=bid, label=lab), size=text.size) +
        scale_x_discrete(position='bottom', labels=tpa$algd, expand=c(0,0)) +
        scale_y_discrete(labels=tpb$blgd, expand = c(0,0)) +
        scale_fill_gradientn(name=leg, na.value='grey50', colors=cols) +
        #scale_fill_viridis(name=leg) +
        otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
               xtick=T,ytick=T,xtext=T,ytext=T,
               margin = c(r.top,.5,.5,.05)) +
        theme(panel.border = element_blank()) +
        theme(axis.text.x=element_text(color=tpa$acol,angle=45,size=7,hjust=1,vjust=1)) +
        theme(axis.text.y=element_text(color=tpb$bcol))
    ggarrange(p1, p2, nrow = 1, ncol = 2, widths = c(1,ratio))
    #}}}
}

