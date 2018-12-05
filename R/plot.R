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

#' flexible custom ggplot theme
#'
#' @export
otheme <- function(strip.size = 8, margin = c(.5,.5,.5,.5),
                   legend.pos = '', legend.dir = 'v', legend.border = F, 
                   xticks = F, yticks = F, xtitle = F, ytitle = F, 
                   xtext = F, ytext = F, xgrid = F, ygrid = F,
                   xsize = 8, ysize = 8) {
    #{{{ custom theme
    o = theme_bw() +
        theme(strip.background = element_blank(),
              strip.text = element_text(size = strip.size,
                                        margin = margin(0,0,0,0,'lines')))
    #{{{ legend
    o = o +
        theme(legend.background = element_blank(),
              legend.title = element_blank(),
              legend.key.size = unit(.8, 'lines'), 
              legend.text = element_text(size = 8))
    if(legend.pos == 'none') {
        o = o + theme(legend.position = 'none')
    } else if(legend.pos == 'top.center.out') {
        o = o + theme(legend.position = c(.5,1), legend.justification = c(.5,0))
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
    }
    if(legend.dir == 'h') 
        o = o + theme(legend.direction = 'horizontal')
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
    #{{{ grid
    o = o + theme(panel.grid.minor = element_blank())
    if(!xgrid)
        o = o + theme(panel.grid.major.x = element_blank())
    if(!ygrid)
        o = o + theme(panel.grid.major.y = element_blank())
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
