# Readme

## How to run

```bash
wget https://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz
gzip -d Homo_sapiens.GRCh38.104.gtf.gz
docker run --rm -w `pwd` -v `pwd`:`pwd` btrspg/baleen_data:latest bash -c "Rscript convert_genomic2transcriptomic.R"
docker run --rm -w `pwd` -v `pwd`:`pwd` btrspg/baleen_data:latest bash -c "python get_final_benchmark_hek293t.py"
docker run --rm -w `pwd` -v `pwd`:`pwd` btrspg/baleen_data:latest bash -c "bedtools merge -i hek293t_m6ace/Hek293T_m6aceSeq_annotated_m6aceincluded.bed > hek293t_m6ace/Hek293T.bed"
```

