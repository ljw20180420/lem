#!/bin/bash

download_genome() {
    wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz -O ${DATA_DIR}/chromFa.tar.gz
    genomepy plugin enable blacklist
    genomepy plugin enable bowtie2
    genomepy plugin enable bwa
    genomepy install \
        -p local \
        -l GRCm38 \
        ${DATA_DIR}/chromFa.tar.gz
}

download_hic() {
    wget \
        https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279296/suppl/GSE279296%5Funtagged%5Funtreated%5Fmerge%2Drep1%2D2%2Emcool \
        -O ${DATA_DIR}/wt.mcool
}

download_genome

download_hic
