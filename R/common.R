#' Get summary statistics for a numeric vector
#'
#' @export
desc_stat <- function(x) {
    #{{{
    x = x[!is.na(x)]
    tibble(q0 = quantile(x,0), q5=quantile(x,.05), q25=quantile(x,.25),
           q50 = quantile(x,.5), q75=quantile(x,.75), q95=quantile(x,.95),
           q100 = quantile(x,1), mean=mean(x), sd=sd(x),
           iqr = IQR(x), mad=mad(x))
    #}}}
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


