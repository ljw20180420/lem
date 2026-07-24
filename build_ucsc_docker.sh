#!/bin/bash

config_docker_build_and_run_proxy() {
    local port=$1

    mv ~/.docker/config.json ~/.docker/config.json.bak
    jq -s '.[0] * .[1]' ~/.docker/config.json.bak - \
        > ~/.docker/config.json \
        << EOF
{
    "proxies": {
            "default": {
            "httpProxy": "http://127.0.0.1:${port}",
            "httpsProxy": "http://127.0.0.1:${port}",
            "noProxy": "localhost,127.0.0.1/8"
        }
    }
}
EOF
}

config_docker_daemon_proxy() {
    local port=$1

    mv ~/.config/docker/daemon.json ~/.config/docker/daemon.json.bak
    jq -s '.[0] * .[1]' ~/.config/docker/daemon.json.bak - \
        > ~/.config/docker/daemon.json \
        << EOF
{
    "proxies": {
        "http-proxy": "http://127.0.0.1:${port}",
        "https-proxy": "http://127.0.0.1:${port}",
        "no-proxy": "localhost,127.0.0.0/8"
    }
}
EOF
}

config_docker() {
    local port=$1

    config_docker_build_and_run_proxy ${port}
    config_docker_daemon_proxy ${port}
    systemctl --user daemon-reload
    systemctl --user restart docker
}

config_docker
if ! [ -f Dockerfile ]
then
    wget https://raw.githubusercontent.com/ucscGenomeBrowser/kent/master/src/product/installer/docker/Dockerfile
fi
docker build . -t ljw/ucsc_genomebrowser_image
