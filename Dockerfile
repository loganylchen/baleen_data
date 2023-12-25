FROM btrspg/vscode-base:0.0.4

ENV DEBIAN_FRONTEND noninteractive
ENV PATH /opt/bin:$PATH

RUN mkdir -p /opt/bin && \
    wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static -O /opt/bin/bedtools && \
    chmod +x /opt/bin/* && \
    R -e 'BiocManager::install("ensembldb");IRkernel::installspec(name = "VSCODER", displayname = "VSCODER",user=FALSE)'
ADD requirements.txt /tmp/requirements.txt
RUN pip3 install -r /tmp/requirements.txt && rm /tmp/requirements.txt
