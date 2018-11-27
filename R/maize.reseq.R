dirr = '~/git/rmaize/R'
source(file.path(dirr, 'header.R'))
source(file.path(dirr, 'plot.R'))
source(file.path(dirr, 'genomes.R'))

dird = '~/projects/maize.reseq/data'
dirp = dird
dirc = '/scratch.global/zhoux379/maize.reseq'
f_cfg = file.path(dird, '01.cfg.xlsx')
t_cfg = read_xlsx(f_cfg, sheet=1, col_names=T) %>% replace_na(list('meta'=F))
f_yml = file.path(dird, '11.cfg.yaml')
Sys.setenv("R_CONFIG_FILE" = f_yml)
#}}}




