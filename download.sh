#!/bin/bash

plugin() {
    genomepy plugin enable blacklist
    genomepy plugin enable bowtie2
    genomepy plugin enable bwa
}

download_genome() {
    local url=$1
    local name=$2
    wget "${url}" -O ${DATA_DIR}/${name}.tar.gz
    plugin
    genomepy install \
        -p local \
        -l ${name} \
        ${DATA_DIR}/${name}.tar.gz
}

download_geo() {
    local gse=$1
    geofetch -i ${gse} \
        --processed \
        --data-source all \
        --geo-folder ${DATA_DIR} \
        -m ${DATA_DIR}
}

# download_genome "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz" "GRCm38"
# download_genome "https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz" "GRCm37"

# hic
download_geo GSE279296

# ctcf chip-exp
# download_geo GSE235386