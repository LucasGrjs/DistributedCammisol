# 1. OS DE BASE
FROM ubuntu:24.04

# Ãviter les questions interactives pendant l'installation
ENV DEBIAN_FRONTEND=noninteractive

# 2. INSTALLATION DES OUTILS SYSTÃME
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    tar \
    ca-certificates \
    build-essential \
    libevent-dev \
    zlib1g-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# INSTALLATION de MAVEN 3.8.6
WORKDIR /opt/maven
RUN curl -L https://archive.apache.org/dist/maven/maven-3/3.8.6/binaries/apache-maven-3.8.6-bin.tar.gz -o maven.tar.gz \
    && tar -xzf maven.tar.gz \
    && rm maven.tar.gz

ENV PATH="/opt/maven/apache-maven-3.8.6/bin:${PATH}"

# 3. INSTALLATION DE JAVA 17
WORKDIR /opt/java
RUN curl -L --retry 5 --retry-delay 3 \
    "https://api.adoptium.net/v3/binary/latest/17/ga/linux/x64/jdk/hotspot/normal/eclipse?project=jdk" \
    -o jdk.tar.gz \
    && tar -xzf jdk.tar.gz \
    && mv jdk-17* jdk-latest \
    && rm jdk.tar.gz

# Configuration des variables d'environnement Java
ENV JAVA_HOME=/opt/java/jdk-latest
ENV PATH=$JAVA_HOME/bin:$PATH

# 4. COMPILATION D'OPEN MPI 4.1.4 (Support Java activÃ©)
WORKDIR /tmp/openmpi
RUN curl -L https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.4.tar.gz -o openmpi.tar.gz \
    && tar -xzf openmpi.tar.gz \
    && cd openmpi-4.1.4 \
    && ./configure --enable-mpi-java --prefix=/usr/local \
    && make -j $(nproc) \
    && make install \
    && ldconfig

# 5. SOURCES DU PROJET
# gama et gama.experimental doivent rester frÃ¨res (gama.parent rÃ©fÃ©rence
# gama.experimental via ../../gama.experimental/...), comme dans le dÃ©pÃ´t.
WORKDIR /app
COPY ./gama ./gama
COPY ./gama.experimental ./gama.experimental
COPY ./Cammisol ./Cammisol

ENV MAVEN_OPTS="-Djdk.xml.maxGeneralEntitySizeLimit=0 -Djdk.xml.totalEntitySizeLimit=0 -Djdk.xml.entityExpansionLimit=0"

# 6. COMPILATION DE GAMA (inclut l'extension MPI via le rÃ©acteur gama.parent)
# sed: les scripts viennent d'un checkout Windows (CRLF), ce qui casse le
# shebang une fois copiÃ© dans le conteneur Linux ("./build.sh: not found").
WORKDIR /app/gama/travis
RUN sed -i 's/\r$//' build.sh && chmod +x build.sh && ./build.sh

# 7. PRÃPARATION DES SCRIPTS CAMMISOL
WORKDIR /app/Cammisol/models
RUN sed -i 's/\r$//' startHeadless startDistributedCammisol \
    && chmod +x startHeadless startDistributedCammisol

# 8. UTILISATEUR ET DROITS
WORKDIR /app
RUN useradd -m mpiuser && chown -R mpiuser:mpiuser /opt /app
USER mpiuser
