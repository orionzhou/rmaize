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

#' flexible custom ggplot theme
#'
#' @export
otheme <- function(strip.size = 8, margin = c(.5,.5,.5,.5),
                   legend.pos='', legend.dir='v', legend.box='v', legend.vjust=NA,
                   legend.title = F, legend.border = F,
                   panel.border = T, panel.spacing = .02,
                   xticks = F, yticks = F, xtitle = F, ytitle = F,
                   xtext = F, ytext = F, xgrid = F, ygrid = F,
                   xsize = 8, ysize = 8, nostrip = T) {
    #{{{ custom theme
    o = theme_bw()
    if(nostrip) o = o +
        theme(strip.background = element_blank(),
              strip.text = element_text(size = strip.size,
                                        margin = margin(.05,.05,.05,.05,'lines')))
    #{{{ legend
    o = o +
        theme(legend.background = element_blank(),
              legend.key.size = unit(.8, 'lines'),
              legend.text = element_text(size = 8))
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
        o = o + theme(axis.ticks.x = element_blank())
    if(!yticks)
        o = o + theme(axis.ticks.y = element_blank())
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
    o + theme(plot.margin = unit(margin, "lines"))
    #}}}
}

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
cmp_proportion1 <- function(ti, xtitle='', ytitle='', xangle=0, margin.l=.1,
    legend.pos='top.center.out', legend.dir='h', legend.title='', barwidth=.8,
    oneline=F, expand.x=c(.01,.01), expand.y=c(.01,.04), pal='npg') {
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
        mutate(lab = str_c(number(n,accuracy=1), sep, '(', percent(prop,accuracy=1), ')')) %>%
        ungroup() %>%
        mutate(tag2 = factor(tag2, levels=tags2))
    tpx = tp %>% group_by(tag1) %>% summarise(n=number(sum(n))) %>% ungroup()
    pp = ggplot(tp) +
        geom_bar(aes(x=tag1,y=prop,fill=tag2), stat='identity', position='fill', width=barwidth, alpha=.8) +
        geom_text(aes(x=tag1,y=y,label=lab),color=col1,size=2.5,lineheight=.8) +
        geom_text(data=tpx, aes(tag1,1.01,label=n), color=col2,size=3, vjust=0) +
        scale_x_discrete(name=xtitle, expand=expansion(mult=expand.x)) +
        scale_y_continuous(name=ytitle, expand=expansion(mult=expand.y)) +
        get(str_c('scale_fill', pal, sep="_"))(name=legend.title) +
        otheme(legend.pos=legend.pos,legend.dir=legend.dir,
               legend.title=legend.title, legend.vjust=.2,
            xtick=T, ytick=F, xtitle=T, xtext=T, ytext=F, panel.border=F,
            margin = c(.1, .1, .1, margin.l))
    if(xangle != 0) pp = pp + theme(axis.text.x=element_text(angle=xangle, hjust=1,vjust=1))
    pp
    #}}}
}

#' plot proportions of a set of categories
#
#' @export
cmp_proportion <- function(ti, xtitle='', ytitle='', xangle=0, margin.l=.1,
    legend.pos='top.center.out', legend.dir='h', legend.title='', barwidth=.8,
    oneline=F, expand.x=c(.01,.01), expand.y=c(.01,.04), pal='Pastel2', nc=5) {
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
        mutate(lab = str_c(number(n), sep, '(', percent(prop,accuracy=1), ')')) %>%
        ungroup() %>%
        mutate(tag2 = factor(tag2, levels=tags2))
    tpx = tp %>% group_by(pnl, tag1) %>% summarise(n=number(sum(n))) %>% ungroup()
    pp = ggplot(tp) +
        geom_bar(aes(x=tag1,y=prop,fill=tag2), stat='identity', position='fill', width=barwidth) +
        geom_text(aes(x=tag1,y=y,label=lab),color=col1,size=2.5,lineheight=.8) +
        geom_text(data=tpx, aes(tag1,1.01,label=n), color=col2,size=3, vjust=0) +
        scale_x_discrete(name=xtitle, expand=expansion(mult=expand.x)) +
        scale_y_continuous(name=ytitle, expand=expansion(mult=expand.y)) +
        get(str_c('scale_fill', pal, sep="_"))(name=legend.title) +
        facet_wrap(~pnl, ncol=nc) +
        otheme(legend.pos=legend.pos,legend.dir=legend.dir,
               legend.title=legend.title, legend.vjust=-.3,
            xtick=T, ytick=F, xtitle=T, xtext=T, ytext=F,
            margin = c(.1, .1, .1, margin.l))
    if(xangle != 0) pp = pp + theme(axis.text.x=element_text(angle=xangle, hjust=1,vjust=1))
    pp
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

#' quick linear regression
#'
#' @export
lg <- function(a, b, name1, name2, f_png) {
    #{{{
    png(filename=f_png, width=500, height=500, units='px');
    plot(a, b, type="p", xlab=name1, ylab=name2)
    fit = lm(b~a)
    abline(fit, col="blue")
    fit.sum = summary(fit)
    ann = paste("adjusted Rsquare = ", sprintf("%.04f", fit.sum$adj.r.squared), sep="")
    text(0.8*min(a)+0.2*max(a), 0.8*min(b)+0.2*max(b), ann, col='red')
    dev.off();
    #}}}
}
