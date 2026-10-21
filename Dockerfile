FROM ubuntu:24.04 AS build

RUN apt-get update \
    && apt-get install -y --no-install-recommends g++ make \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY Makefile main.cpp ./
COPY indexer/ indexer/
COPY search/ search/

RUN make -j2 && ./skippy --help

FROM ubuntu:24.04

RUN apt-get update \
    && apt-get install -y --no-install-recommends libstdc++6 libgcc-s1 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /src/skippy /usr/local/bin/skippy

WORKDIR /data
ENTRYPOINT ["/usr/local/bin/skippy"]
CMD ["--help"]
