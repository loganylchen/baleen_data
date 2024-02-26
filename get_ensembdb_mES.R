library(ensembldb)
gtf<-'Mus_musculus.GRCm39.104.gtf'
DB <- ensDbFromGtf(gtf=gtf)
save.image('./working_GRCm39.RData')