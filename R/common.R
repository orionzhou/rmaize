#' Get summary statistics for a numeric vector
#'
#' @export
desc_stat <- function(x) {
    #{{{
    if(typeof(x) == 'list') x = x %>% pull(1)
    x = x[!is.na(x)]
    tibble(statK = c('q0','q5','q25','q50','q75','q95','q100',
                     'avg','sd','iqr','mad'),
           statV = c(map_dbl(c(0,.05,.25,.5,.75,.95,1), .f <- function(y) quantile(x,y)),
                     mean(x), sd(x), IQR(x), mad(x)))
    #}}}
}

#' simple t-test to mark significance
#'
#' @export
ttest_signif <- function(x, y) t.test(x, y)$p.value
map_signif <- function(p) ifelse(p<0.001, "***",
    ifelse(p<0.01, "**", ifelse(p<0.05, '*', 'NS')))

#' Modified standard deviation funciton after removing missing values
#'
#' @export
sd2 <- function(x) sd(x[!(is.na(x) | is.infinite(x) | is.nan(x))])

#' Get summary statistics for the 2nd+ columns of a tibble
#'
#' @export
sum_stat_tibble <- function(ti) {
    #{{{
    statks = colnames(ti)[-1]
    to = ti %>% rename(sid = 1) %>%
        gather(stat, v, -sid) %>%
        mutate(stat = factor(stat, levels=statks)) %>%
        group_by(stat) %>%
        summarise(n=n(), q0=quantile(v,0), q5=quantile(v,.05), q25=quantile(v,.25),
            q50=quantile(v,.5), q75=quantile(v,.75), q95=quantile(v,.95),
            q100=quantile(v,1), avg=mean(v), std=sd(v),
            iqr=IQR(v), mad=mad(v)) %>% ungroup()
    to
    #}}}
}

#' get most frequent item in a vector
#'
#' @export
Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

.ls.objects <- function(pos=1, pattern, order.by, decreasing=F, head=F, n=5) {
    #{{{
    napply <- function(names, fn) sapply(names, function(x) fn(get(x, pos=pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.prettysize <- napply(names, function(x) {
                           capture.output(print(object.size(x), units = "auto")) })
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
    names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
    #}}}
}

#' Improved List of Objects
#'
#' @export
lsos <- function(..., n=10) .ls.objects(..., order.by="Size", decreasing=T, head=T, n=n)

#' sort utility
#'
#' @export
sortC <- function(...) {
    #{{{
    a <- Sys.getlocale("LC_COLLATE")
    on.exit(Sys.setlocale("LC_COLLATE", a))
    Sys.setlocale("LC_COLLATE", "C")
    sort(...)
    #}}}
}

#' set bound to a value/vector
#'
#' @export
set_bound <- function(x, minV, maxV) min(max(x,minV),maxV)

#' a unified distance function with custome methods
#'
#' @export
idist <- function(m, opt='row', method='euclidean', ...) {
    #{{{ be default clusters by column
    if( method %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski") ) {
        m2 = if (opt == 'row') m else t(m)
        dist(m2, method = method, ...)
    } else if(method %in% c("pearson",'spearman','kendall')) {
        m2 = if (opt == 'row') t(m) else m
        as.dist(1-cor(m2, method = method, ...))
    } else if (method == 'gower') {
        as.dist(as.matrix(daisy(m, metric=method)))
    } else {
        stop("unsupported dist method: \n", method)
    }
    #}}}
}
