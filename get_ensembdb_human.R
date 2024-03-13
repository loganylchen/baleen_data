library(ensembldb)
gtf<-'Homo_sapiens.GRCh38.104.gtf'
DB <- ensDbFromGtf(gtf=gtf)
save.image('./working.RData')