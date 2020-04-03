#' Get system path to genome directory
#'
#' @export
genome_dir <- function(genome='', root='~/projects/genome/data') {
    file.path(root, genome)
}

#' Read Genome Configuration
#'
#' @export
read_genome_conf <- function(genome='Zmays_B73', dirg = '~/projects/genome/data') {
    #{{{
    fi = file.path(dirg, genome, '55.rds')
    if(!file.exists(fi)) {
        stop(sprintf("cannot read %s\n", fi))
    } else {
        readRDS(fi)
    }
    #}}}
}

#' Read maize v3 to v4 gene ID mapping
#'
#' @export
v3_to_v4 <- function(fm="~/data/genome/Zmays_B73/gene_mapping/maize.v3TOv4.geneIDhistory.txt") {
    #{{{
    fi2 = "~/data/genome/Zmays_B73/gene_mapping/v3_v4_extra.xlsx"
    ti2 = read_xlsx(fi2, col_names=c('ogid','gid')) %>% mutate(type='extra')
    read_tsv(fm, col_names = F) %>%
        transmute(ogid = X1, gid = X2, change = X3, method = X4, type = X5) %>%
        select(ogid, gid, type) %>%
        bind_rows(ti2)
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

#' Get maize gene symbols
#'
#' @export
read_symbol <- function(opt=1) {
    #{{{
    fi = '~/projects/genome/data/Zmays_B73/61_functional/00.symbols.tsv'
    ti = read_tsv(fi)
    if(opt == 1) {
        ti %>% select(gid, symbol)
    } else {
        ti %>% select(gid, symbol)
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
        fi = '~/projects/genome/data/Zmays_B73/gene_mapping/synteny_maize.tsv'
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

#' Read TF IDs
#'
#' @export
read_tf <- function(fi = '~/projects/genome/data/Zmays_B73/61_functional/06.tf.tsv')
    read_tsv(fi)

#' Get maize ACR
#'
#' @export
read_regions <- function(opt='umr') {
    #{{{
    dirw = '~/projects/genome/data/Zmays_B73/chromatin'
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



