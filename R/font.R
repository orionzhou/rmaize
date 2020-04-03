#' Generate load custom fonts
#'
#' @export
load_fonts <- function(dira = '~/projects/assets') {
    #{{{ add fonts
    font_paths(file.path(dira, 'fonts'))
    #font_files()
    font_add('emoji', "emoji.ttf")
    font_add('fab', "fab.ttf")
    font_add('far', "far.ttf")
    font_add('fas', "fas.ttf")
    font_add('icm', "icomoon.ttf")
    showtext_auto()
    #
    fas = c('\uf4d8', '\uf00d')
    fa1 = fas[1]; fa2 = fas[2]
    #}}}
}

