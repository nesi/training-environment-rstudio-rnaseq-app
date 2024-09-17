FROM rocker/geospatial:latest

# install system dependencies
# TODO: move ldap-utils, libnss-ldapd, libpam-ldapd, nscd, nslcd to base image??
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

# install miniconda, create environment and install packages
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh \
    && bash /miniconda.sh -b -p /opt/miniconda3 \
    && rm -f /miniconda.sh \
    && . /opt/miniconda3/etc/profile.d/conda.sh \
    && conda config --set auto_activate_base false \
    && conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda config --set channel_priority strict \
    && conda create --yes --prefix /var/lib/rnaseq \
       python=3.10 \
       fastqc \
       multiqc \
       cutadapt \
       hisat2 \
       samtools \
       subread \
       bamtools \
       salmon \
       star \
       qualimap \
    && echo ". /opt/miniconda3/etc/profile.d/conda.sh" > /etc/profile.d/conda.sh \
    && echo "conda activate /var/lib/rnaseq" >> /etc/profile.d/conda.sh

# install R packages
RUN Rscript -e 'install.packages("BiocManager", repos = "https://cloud.r-project.org")' \
    && Rscript -e 'BiocManager::install("limma")' \
    && Rscript -e 'BiocManager::install("edgeR")' \
    && Rscript -e 'BiocManager::install("DESeq2")' \
    && Rscript -e 'install.packages("beeswarm")' \
    && Rscript -e 'install.packages("knitr")' \
    && Rscript -e 'install.packages("gplots")' \
    && Rscript -e 'BiocManager::install("goseq")' \
    && Rscript -e 'BiocManager::install("GO.db")' \
    && Rscript -e 'BiocManager::install("org.Sc.sgd.db")'

# copy data
COPY docker/RNA_seq/ /var/lib/RNA_seq/
