# syntax=docker/dockerfile:1.6

############################
# Build stage
############################
FROM ubuntu:24.04 AS build

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends                \
    build-essential g++ make pkg-config                                         \
    libgmp-dev libmpfr-dev libmpc-dev                                           \
    git ca-certificates curl                                                    \
 && rm -rf /var/lib/apt/lists/*

# mpfrc++ (mpreal.h).
WORKDIR /opt
RUN git clone --depth=1 https://github.com/advanpix/mpreal.git mpfrcpp
ENV MPREAL_INC=/opt/mpfrcpp

# Copy sources after deps to maximize caching
WORKDIR /src
COPY . .

# Build args: NATIVE=1 turns on -march=native
ARG NATIVE=1
# Common compile-time flags for mpfrc++
# ENV CPPFLAGS="-DOS_USE_MPFR -DMPFR_USE_NO_MACRO -DMPFR_USE_INTMAX_T -I${MPREAL_INC} -include stdint.h"
ENV CPPFLAGS="-DMPFR_USE_NO_MACRO -DMPFR_USE_INTMAX_T -I${MPREAL_INC} -include stdint.h"

# Build
RUN set -eux;                                                                   \
    if [ "$NATIVE" = "1" ]; then                                                \
        EXTRA_CXXFLAGS="-O3 -march=native -mtune=native";                       \
    else                                                                        \
        EXTRA_CXXFLAGS="-O3";                                                   \
    fi;                                                                         \
    make clean || true;                                                         \
    make TARGET=bin/z -j"$(nproc)"                                                           \
        CXXFLAGS="$EXTRA_CXXFLAGS -std=c++20 -MMD -MP -fopenmp $CPPFLAGS"       \
        LDFLAGS="-fopenmp"                                                      \
        MPREAL_INC="$MPREAL_INC"                                                

############################
# Runtime stage
############################
FROM ubuntu:24.04 AS runtime

# update libs
RUN apt-get update && apt-get install -y --no-install-recommends                \
    libgmp10 libmpfr6 libmpc3 libgomp1                                           \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY --from=build /src/bin /app/bin
# If you have config files:
COPY conf /app/conf

# Default OpenMP runtime hints; override at run-time as needed
ENV OMP_PROC_BIND=close \
    OMP_PLACES=cores

# Default entry; you’ll likely pass args explicitly in `docker run`
COPY run.sh /app/run.sh
RUN chmod +x /app/run.sh
ENTRYPOINT ["/bin/bash"]

############################
# Executing build
############################
# docker build -t jor_z:native --build-arg NATIVE=1 .
# docker tag jor_z:native jacobore/jor_z:native
# docker run --rm -e OMP_NUM_THREADS -e OMP_PROC_BIND -e OMP_PLACES jacobore/jor_z:native /app/run.sh z 1122334455667788.758122007832190107859999875437
# docker run --rm -e OMP_NUM_THREADS=4 -e OMP_PROC_BIND=close -e OMP_PLACES jacobore/jor_z:native /app/run.sh z 4400112233445565.831381180135852831954333609817