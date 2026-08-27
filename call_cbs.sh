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

data_dir="/home/ljw/sdc1/cpcdh"
genome="/home/ljw/.local/share/genomes/GRCm38/GRCm38.fa"

if [[ "${1,,}" == "jaspar" ]]
then
    cmd="find_jaspar_motif"
elif [[ "${1,,}" == "ecbs" ]]
then
    cmd="find_eCBS_motif"
else
    echo "command is jaspar or ecbs" >&2
    exit 1
fi

shift

${cmd} "$@"
