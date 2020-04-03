#' read Google Spreadsheet "maize_studies"
#'
#' @export
read_gspread_master <- function(book='maize_studies', sheet=1, lib='rnaseq',
                        cred='~/.config/google_account_token.json') {
    #{{{
    drive_auth(path=cred)
    sheets_auth(path=cred)
    gbook = drive_get(book)
    t_cfg = read_sheet(gbook, sheet=sheet, col_names=T) %>%
        rename(ase=ASE,ril=RIL,run=Run) %>%
        mutate(ase=as.logical(ase)) %>%
        mutate(stress=as.logical(stress)) %>%
        mutate(ril=as.logical(ril)) %>%
        replace_na(list(ase=F,stress=F,ril=F)) %>%
        mutate(lgd = sprintf("%s%d %s [%d]", str_to_title(author),year,study,n)) %>%
        filter(libtype == !!lib) %>%
        select(yid,author,year,study,genotype,tissue,n,ase,stress,ril,run)
    t_cfg
    #}}}
}

#' read a generic googlesheet
#'
#' @export
read_gspread <- function(book='maize_studies', sheet=1,
                        cred='~/.config/google_account_token.json') {
    #{{{
    drive_auth(path=cred)
    sheets_auth(path=cred)
    gbook = drive_get(book)
    read_sheet(gbook, sheet=sheet, col_names=T)
    #}}}
}

