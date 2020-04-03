extract_read_stats <- function(stat) {
    #{{{
    spots = as.integer(xml_attr(stat, 'nspots'))
    x1 = xml_attrs(xml_find_all(stat, './Read'))
    x2 = as_tibble(data.frame(t(sapply(x1, c)), stringsAsFactors = F)) %>%
        mutate(index=as.integer(index), count=as.integer(count),
               average=as.numeric(average), stdev=as.numeric(stdev))
    x2a = x2 %>% filter(index == 0)
    x2b = x2 %>% filter(index == 1)
    stopifnot(nrow(x2a) == 1)
    stopifnot(nrow(x2b) == 0 | nrow(x2b) == 1)
    nread1 = x2a$count[1]; avgLength1 = x2a$average
    nread2=0; avgLength2=0
    spots_with_mate = 0
    if(nrow(x2b) == 1) {
        nread2 = x2b$count[1]; avgLength2 = x2b$average
        if (spots == nread1 | spots == nread2) {
            spots_with_mate = min(nread1, nread2)
        }
    }
    paired = spots_with_mate / spots >= .5
    avgLength = round((avgLength1*nread1+avgLength2*nread2) / spots)
    tibble(paired=paired, spots=spots, spots_with_mate=spots_with_mate, avgLength=avgLength)
    #}}}
}

extract_lib_stats <- function(stat) {
    #{{{
    ptns = c(
        lib_name = './LIBRARY_NAME',
        lib_strategy = './LIBRARY_STRATEGY',
        lib_source = './LIBRARY_SOURCE',
        lib_selection = './LIBRARY_SELECTION',
        lib_layout = './LIBRARY_LAYOUT'
    )
    tag2str <- function(ptn, xml) {
        r = xml_find_all(xml, ptn)
        ifelse(length(r) > 0, xml_text(r), '')
    }
    tibble(tag = names(ptns), ptns = ptns) %>%
        mutate(value = map_chr(ptns, tag2str, xml=stat)) %>%
        select(-ptns) %>% spread(tag, value)
    #}}}
}

#' Query SRA using BioProject ID(s) and return meta table of Runs
#'
#' @export
get_sra_meta <- function(acc, yid) {
    #{{{ sra_ids = 'PRJNA494874 PRJNA289143'
    require(rentrez)
    require(xml2)
    cat(sprintf("working on %s\n", yid))
    accs = unlist(str_split(acc, ' '))

    n = length(accs)
    batch = ceiling(n/100)
    ids = c()
    for (i in 1:batch) {
        #{{{
        idxb = (i-1) * 100 + 1
        idxe = min(i*100, n)
        accs1 = accs[idxb:idxe]
        cat(sprintf("  searching batch %2d: %d - %d\n", i, idxb, idxe))
        term = str_c(sprintf("%s", accs1), collapse = ' OR ')
        res = entrez_search(db="sra", term=term, retmax=1e5)
        ids = c(ids, res$ids)
        #}}}
    }

    n = length(ids)
    batch = ceiling(n/100)
    for(i in 1:batch) {
        #{{{
        idxb = (i-1) * 100 + 1
        idxe = min(i*100, n)
        ids1 = ids[idxb:idxe]
        cat(sprintf("  retrieving batch %2d: %d - %d\n", i, idxb, idxe))
        recs = entrez_fetch(db="sra", id=ids1, rettype="xml", parsed=F, retmax=1e5)
        xml1 = read_xml(recs)
        if(i == 1) {
            xml = xml1
        } else {
            for (child in xml_children(xml1))
                xml_add_child(xml, child)
        }
        #}}}
    }

    srrs = xml_attr(xml_find_all(xml, '//RUN'), 'accession')
    nr = length(srrs)
    srr_alias = xml_attr(xml_find_all(xml, '//RUN'), 'alias')
    stopifnot(length(srrs)==nr & length(srr_alias)==nr)
    srxs = xml_attr(xml_find_all(xml, '//RUN/EXPERIMENT_REF'), 'accession')
    stats = xml_find_all(xml, "//RUN/Statistics")
    stopifnot(length(srxs)==nr & length(stats)==nr)
    ti = tibble(srr=srrs, srr_alias=srr_alias,
                srx = srxs, stats = c(stats)) %>%
        mutate(data = map(stats, extract_read_stats)) %>%
        select(-stats) %>% unnest(data)

    nx = length(xml_children(xml))
    srxs = xml_attr(xml_find_all(xml, '//EXPERIMENT'), 'accession')
    srx_alias = xml_attr(xml_find_all(xml, '//EXPERIMENT'), 'alias')
    stopifnot(length(srxs)==nx & length(srx_alias)==nx)
    srx_title = xml_text(xml_find_all(xml, '//EXPERIMENT/TITLE'))
    srx_design = xml_text(xml_find_all(xml, '//EXPERIMENT/DESIGN/DESIGN_DESCRIPTION'))
    smps = xml_attr(xml_find_all(xml, '//SAMPLE'), 'alias')
    studies = xml_text(xml_find_all(xml, '//STUDY/DESCRIPTOR/STUDY_TITLE'))
    stats = xml_find_all(xml, '//EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR')
    stopifnot(length(srx_title)==nx & length(stats)==nx)
    tx = tibble(srx=srxs, srx_alias=srx_alias,
                srx_title = srx_title, srx_design = srx_design,
                study_title = studies,
                sample_alias=smps, stats = c(stats)) %>%
        mutate(data = map(stats, extract_lib_stats)) %>%
        select(-stats) %>% unnest(data)

    ti %>% left_join(tx, by = 'srx')
    #}}}
}

#' Read SRA run table, format to tibble
#'
#' @export
read_sra_run <- function(fi, fi2) {
    #{{{
    ti0 = read_csv(fi)
    ti = ti0 %>% select(Run, spots, spots_with_mates, avgLength,
                        LibraryName, LibraryStrategy, LibraryLayout, SampleName,
                        BioSample, Sample, Experiment) %>%
        mutate(paired = ifelse(spots_with_mates/spots >= .5, T, F)) %>%
        print(width = Inf)
    ti %>% count(LibraryLayout, paired) %>% print(n=5)
    if(file.exists(fi2)) {
        ti2 = read_csv(fi2) %>%
            transmute(Experiment = `Experiment Accession`,
                      Title = `Experiment Title`,
                      Title2 = `Study Title`
            )
        ti = ti %>% left_join(ti2, by = 'Experiment') %>% select(-Experiment)
    }
    ti
    #}}}
}

#' Fill in 'Replicate' column of a sample table
#'
#' @export
sra_fill_replicate <- function(th) {
    #{{{
    cmap = c()
    for (i in 1:nrow(th)) {
        tis = th$Tissue[i]; gt = th$Genotype[i]; trt = th$Treatment[i]
        key = sprintf("%s-%s-%s", tis, gt, trt)
        if (key %in% names(cmap))
            cmap[key] = cmap[key] + 1
        else
            cmap[key] = 1
        th$Replicate[i] = cmap[key]
    }
    th = th %>% mutate(Replicate = as.integer(Replicate))
    th %>% dplyr::count(Replicate) %>% print(n=10)
    th
    #}}}
}

#' helper function
generate_sample_id <- function(n, pre='s') {
    #{{{
    ndig = floor(log10(n)) + 1
    fmt = sprintf("%s%%0%dd", pre, ndig)
    sprintf(fmt, 1:n)
    #}}}
}

#' complete sample list
#'
#' @export
complete_sample_list <- function(ti, fillrep=T) {
    #{{{
    if(!'Treatment' %in% colnames(ti)) ti = ti %>% mutate(Treatment='')
    if(!'Replicate' %in% colnames(ti)) ti = ti %>% mutate(Replicate='')
    if( is.na(ti$SampleID[1]) ) ti$SampleID = generate_sample_id(nrow(ti))
    if('directory' %in% colnames(ti)) ti = ti %>% fill(directory)
    if('MergeID' %in% colnames(ti)) ti = ti %>% select(-MergeID)
    to = ti %>% fill(Tissue, Genotype) %>%
        replace_na(list(Treatment='')) %>%
        select(SampleID, Tissue, Genotype, Treatment, Replicate, everything())
    if(fillrep)
        to = to %>% group_by(Tissue,Genotype,Treatment) %>%
            mutate(Replicate = 1:n()) %>%
            ungroup()
    to = to %>% arrange(SampleID, Tissue, Genotype, Treatment, Replicate)
    tos = to %>% distinct(Tissue, Genotype, Treatment)
    tos = tos %>% mutate(MergeID = generate_sample_id(nrow(tos), pre='m'))
    to %>% inner_join(tos, by = c('Tissue','Genotype','Treatment')) %>%
        select(SampleID, Tissue, Genotype, Treatment, Replicate, MergeID, everything())
    #}}}
}

