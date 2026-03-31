FROM rocker/geospatial:4.5.2

# install packages, some of these are required to make the container behave like ldap client
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        curl \
        dnsutils \
        git \
        jq \
        ldap-utils \
        libnss-ldapd \
        libpam-ldapd \
        less \
        nano \
        nodejs \
        nscd \
        nslcd \
        rsync \
        unzip \
        vim \
        wget \
        zip \
    && rm -rf /var/lib/apt/lists/*

# install kubectl, required for running on the k8s cluster
ARG KUBECTL_VERSION=v1.28.5
RUN curl -LO https://dl.k8s.io/release/$KUBECTL_VERSION/bin/linux/amd64/kubectl \
    && mv kubectl /usr/local/bin/kubectl \
    && chmod +x /usr/local/bin/kubectl

# copy in ood k8s utils
ARG UTILS_HASH=6298fb01f7a7c66a8454e3f0fd74437a32491423
RUN git clone https://github.com/nesi/training-environment-k8s-utils.git /opt/ood-k8s-utils \
    && cd /opt/ood-k8s-utils \
    && git checkout $UTILS_HASH \
    && chmod +x /opt/ood-k8s-utils/files/* \
    && mv /opt/ood-k8s-utils/files/* /bin/ \
    && rm -rf /opt/ood-k8s-utils

# make a dummy module command to avoid warnings from ondemand job_script_content.sh
RUN echo "#!/bin/bash" > /bin/module \
    && chmod +x /bin/module

# install R packages
RUN Rscript -e 'install.packages("vegan")' \
    && Rscript -e 'install.packages("palmerpenguins")' \
    && Rscript -e 'remotes::install_github("G-Thomson/Manu")' \
    && Rscript -e 'install.packages("viridis")' \
    && Rscript -e 'install.packages("ggbeeswarm")' \
    && Rscript -e 'install.packages("ggtext")' \
    && Rscript -e 'install.packages("ggridges")'

# test loading R packages
RUN Rscript -e 'library(vegan)' \
    && Rscript -e 'library(palmerpenguins)' \
    && Rscript -e 'library(Manu)' \
    && Rscript -e 'library(viridis)' \
    && Rscript -e 'library(ggbeeswarm)' \
    && Rscript -e 'library(ggtext)' \
    && Rscript -e 'library(ggridges)'

# copy the introduction to R source into the repo
ARG INTROR_HASH=0729df3d1877356c120612586ca541cf8da17b0e
RUN git clone https://github.com/GenomicsAotearoa/Introduction-to-R /opt/Introduction-to-R \
    && cd /opt/Introduction-to-R \
    && git checkout $INTROR_HASH \
    && rm -rf /opt/Introduction-to-R/.git
