#!/bin/bash

extract_DNA() {
    local chrom=$1
    local start=$2
    local end=$3

    bed2fasta -s \
        <(
            printf "%s\t%d\t%d\t.\t.\t+\n" ${chrom} ${start} ${end}
        ) \
        ${genome}
}

find_jaspar_motif() {
    local chrom=$1
    local start=$2
    local end=$3

    if [[ ! -f "versions.meme" ]]
    then
        wget https://jaspar.elixir.no/api/v1/matrix/MA0139/versions.meme
    fi

    extract_DNA ${chrom} ${start} ${end} |
    fimo --text \
        --thresh 0.001 \
        --motif MA0139.1 \
        versions.meme -
}

find_eCBS_motif() {
    local chrom=$1
    local start=$2
    local end=$3

    extract_DNA ${chrom} ${start} ${end} |
    fimo --text \
        --norc \
        ${data_dir}/result/eCBS.meme -
}

find_motif() {
    local meme=$1
    local chrom=$2
    local start=$3
    local end=$4

    if [ "${meme,,}" == "jaspar" ]
    then
        find_jaspar_motif ${chrom} ${start} ${end}
        return
    fi

    extract_DNA ${chrom} ${start} ${end} |
    fimo --text \
        --thresh 0.001 \
        ${meme} -
}

data_dir="/home/ljw/sdc1/cpcdh"
genome="/home/ljw/.local/share/genomes/GRCm38/GRCm38.fa"

find_motif "$@"
