library(ensembldb)
gtf<-'Homo_sapiens.GRCh38.104.gtf'
DB <- ensDbFromGtf(gtf=gtf)
EDB<-EnsDb(DB)
df <- read.csv('hek293t_m6ace/Hek293T_m6aceSeq_results.csv')
grange_range <- makeGRangesFromDataFrame(df,seqnames.field=c('Chr'),start.field='Start',end.field='End',strand.field='Strand')
transcript_loc <- genomeToTranscript(grange_range,EDB)
write.csv(as.data.frame(transcript_loc),'hek293t_m6ace/Hek293T_m6aceSeq_results_annotated.csv',quote=FALSE)