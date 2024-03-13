args<- commandArgs(trailingOnly = TRUE)
library(ensembldb)
library(dplyr)
contig <- args[1]
f <- args[2]
outdir <- args[3]
print(paste0('Working on ',contig))
load('working.RData')
EDB<-EnsDb(DB)
df <- read.csv(f)
working_contig_df <- df %>% filter(Chr==contig)
grange_range <- makeGRangesFromDataFrame(working_contig_df,seqnames.field=c('Chr'),start.field='Start',end.field='End',strand.field='Strand')
transcript_loc <- genomeToTranscript(grange_range,EDB)
write.csv(as.data.frame(transcript_loc),paste0(outdir,'/tmp_',contig,'.csv'),quote=FALSE,row.names=FALSE)
