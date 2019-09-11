FROM brentp/musl-hts-nim:latest

RUN cd /  &&      \
    git clone -b master --depth 1 git://github.com/brentp/duphold.git && \
    cd duphold && \
    nimble install -y && \
    nim c -d:release -d:danger -o:/usr/bin/duphold --passC:-flto src/duphold && \
    rm -rf /duphold && duphold --help

