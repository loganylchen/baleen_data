library(ensembldb)
gtf<-'Homo_sapiens.GRCh38.104.gtf'
DB <- ensDbFromGtf(gtf=gtf)
EDB<-EnsDb(DB)

saveRDS(EDB,'./ensembldb.rds')