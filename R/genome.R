#' Get system path to genome directory
#'
#' @export
genome_dir <- function(genome='', root='~/projects/s3/zhoup-genome') {
    file.path(root, genome)
}

#' Read Chromosomze size
#'
#' @export
read_chrom_size <- function(genome='Zmays_B73v5', dirg = '~/projects/s3/zhoup-genome') {
    #{{{
    fi = file.path(dirg, genome, '15_intervals/01.chrom.sizes')
    if(!file.exists(fi)) {
        stop(sprintf("cannot read %s\n", fi))
    } else {
        read_tsv(fi, col_names=c('chrom','size'))
    }
    #}}}
}

#' Read Genome Configuration
#'
#' @export
read_genome_conf <- function(genome='Zmays_B73v5', dirg = '~/projects/s3/zhoup-genome') {
    #{{{
    fi = glue("{dirg}/{genome}/55.rds")
    if(!file.exists(fi)) {
        stop(sprintf("cannot read %s\n", fi))
    } else {
        readRDS(fi)
    }
    #}}}
}

#' laod txdb
#'
#' @export
load_txdb <- function(org, primary=F) {
    #{{{
    require(GenomicFeatures)
    fn = ifelse(primary, '15', '10')
    fd = glue('~/projects/s3/zhoup-genome/{org}/50_annotation/{fn}.sqlite')
    if (! file.exists(fd) )
        fd = glue('~/projects/s3/zhoup-genome/Zmays_{org}/50_annotation/{fn}.sqlite')
    loadDb(fd)
    #}}}
}

#' laod multiple txdbs
#'
#' @export
load_txdbs <- function(orgs) {
    #{{{
    txdbs = list()
    for (org in orgs) {
        txdbs[[org]] = load_txdb(org)
    }
    txdbs
    #}}}
}

#' Read maize v3 to v4 gene ID mapping
#'
#' @export
v3_to_v4 <- function(fm="~/projects/genome/data2/Zmays_B73/gene_mapping/maize.v3TOv4.geneIDhistory.txt", opt='tibble') {
    #{{{
    fi2 = "~/projects/genome/data2/Zmays_B73/gene_mapping/v3_v4_extra.xlsx"
    ti2 = read_xlsx(fi2, col_names=c('ogid','gid')) %>% mutate(type='extra')
    to = read_tsv(fm, col_names = F) %>%
        filter(!str_detect(X1, '^#')) %>%
        transmute(ogid = X1, gid = X2, change = X3, method = X4, type = X5) %>%
        select(ogid, gid, type) %>%
        bind_rows(ti2)
    if (opt == 'tibble') {
        to
    } else if (opt == 'dict_all') {
        to2 = to %>% filter(type != '1-to-0') %>%
            group_by(ogid) %>%
            slice(1) %>% ungroup()
        gd = to2$gid
        names(gd) = to2$ogid
        gd
    } else if (opt == 'dict_11') {
        to2 = to %>% filter(type %in% c('1-to-1', 'extra')) %>%
            group_by(ogid) %>%
            slice(1) %>% ungroup()
        gd = to2$gid
        names(gd) = to2$ogid
        gd
    } else {
        stop("unknown opt", opt, '\n')
    }
    #}}}
}

#' prepare flattern 2D genome coordinates
#'
#' @export
flattern_gcoord_prepare <- function(t_size, gap=2e7) {
    #{{{
    offsets = c(0, cumsum(t_size$size + gap)[-nrow(t_size)])
    t_size %>% mutate(offset = offsets) %>%
        mutate(start=1+offset, end=size+offset, pos=(start+end)/2)
    #}}}
}

#' flattern 2D genome coordinates tibble[chrom, pos] to 1D coordinates
#'
#' @export
flattern_gcoord <- function(ti, tz) {
#{{{ flattern tibble:[chrom, pos] to 1D coordinates
    ti %>%
        inner_join(tz[,c('chrom','offset')], by = 'chrom') %>%
        mutate(gpos = pos + offset) %>%
        pull(gpos)
#}}}
}

#' read genotype/genome panels
#'
#' @export
read_gt_panels <- function(fp = "~/projects/s3/zhoup-genome2/panels.xlsx") {
    #{{{
    sheets = excel_sheets(fp)
    t_pnl = tibble(panel=sheets) %>%
        mutate(x = pmap(list(fp, sheet=panel), read_xlsx, col_names="gt")) %>%
        unnest(x) %>%
        group_by(panel) %>% summarise(gts = list(gt), ngt = n()) %>% ungroup()
    t_pnl
    #}}}
}

#' get a genotype/genome panel
#'
#' @export
get_gt_panel <- function(pnl, pnls) pnls %>% filter(panel==pnl) %>% pluck('gts', 1)

#' Get maize gene symbols
#'
#' @export
read_loci <- function(org='Zmays_B73', opt='tibble') {
    #{{{
    fi = glue('~/projects/genome/data2/loci/{org}.tsv')
    ti = read_tsv(fi, col_types='cccccccccc')
    if(opt == 'dict') {
        to = ti %>% select(gid, symbol) %>% filter(!is.na(gid))
        sdic = to$symbol
        names(sdic) = to$gid
        sdic
    } else {
        ti
    }
    #}}}
}

#' Get maize syntenic genes
#'
#' @export
read_syn <- function(gcfg, opt=1) {
    #{{{
    if(opt == 2) {
        fi = '~/projects/wgc/data/21_Sbicolor_B73/10.maize.pairs.tsv'
        read_tsv(fi)
    } else if(opt == 3) { # James Schable's table
        fi = '~/projects/genome/data2/Zmays_B73/gene_mapping/synteny_maize.tsv'
        read_tsv(fi)
    } else {
    ftypes = c('non-syntenic',
               'subgenome1-fractionated',
               'subgenome2-fractionated',
               'subgenome1-retained',
               'subgenome2-retained')
    fsyn = '~/projects/wgc/data/21_Sbicolor_B73/10.maize.subgenome.tsv'
    ti = read_tsv(fsyn)
    gcfg$gene %>% filter(ttype == 'mRNA') %>% select(gid) %>%
        left_join(ti, by = 'gid') %>%
        replace_na(list(ftype = 'non-syntenic')) %>%
        mutate(ftype = str_replace(ftype, "_", "-")) %>%
        mutate(ftype = factor(ftype, levels = ftypes))
    }
    #}}}
}

#' Get maize gene tissue specificity scores
#'
#' @export
read_tis_spec <- function(fi='~/projects/rnaseq/data/15_tissue_spec/10.tau.rds') {
    #{{{
    readRDS(fi)
    #}}}
}

#' Read TF IDs
#'
#' @export
read_tf <- function(fi = '~/projects/genome/data/Zmays_B73/61_functional/06.tf.tsv')
    read_tsv(fi)

#' Get maize ACR
#'
#' @export
read_regions <- function(opt='umr', dirw = '~/projects/genome/data2/Zmays_B73/chromatin') {
    #{{{
    if(opt == 'acrL') {
        fi = file.path(dirw, 'acrL.bed')
    } else if(opt == 'acrE') {
        fi = file.path(dirw, 'acrE.bed')
    } else if(opt == 'umr') {
        fi = file.path(dirw, 'umr.bed')
    } else {
        stop("unknown opt\n")
    }
    read_tsv(fi, col_names=F) %>%
        select(chrom=1,start=2,end=3)
    #}}}
}

#' Get maize cross-reference gene table
#'
#' @export
read_xref <- function(opt='v5', dirw = '~/projects/genome/data2/syntelog') {
    #{{{
    if(opt == 'v4') {
        fi = file.path(dirw, 'xref.maize.v4.tsv')
    } else if(opt == 'v5') {
        fi = file.path(dirw, 'xref.maize.v5.tsv')
    } else if(opt == 'maizeGDB') {
        fi = file.path(dirw, 'xref.maizeGDB.tsv')
    } else {
        stop("unknown opt\n")
    }
    read_tsv(fi, col_names=T)
    #}}}
}




