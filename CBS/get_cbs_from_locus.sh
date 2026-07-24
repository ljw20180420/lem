#!/bin/bash

# change to the dir of the script
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

download_motif() {
    mkdir -p jaspar
    wget https://jaspar.elixir.no/api/v1/matrix/MA0139/versions.meme -O jaspar/versions.meme
}

extract_mm10_DNA() {
    local chrom=$1
    local start=$2
    local end=$3
    local strand=$4
    local genome=$5

    bed2fasta -s \
        <(
            printf "%s\t%d\t%d\t.\t.\t%s\n" ${chrom} ${start} ${end} ${strand}
        ) \
        ${genome}
}

find_motif() {
    extract_mm10_DNA chr18 36929998 36930213 + /home/ljw/.local/share/genomes/GRCm38/GRCm38.fa |
    fimo --text \
        --norc \
        --motif MA0139.1 \
        jaspar/versions.meme -
}

# download_motif

find_motif

