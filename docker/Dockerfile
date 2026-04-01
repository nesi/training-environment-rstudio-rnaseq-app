FROM rocker/geospatial:4.5.2

# install system dependencies
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
    python3 \
    python3-pip \
    python3-venv \
    python3-dev \
    build-essential \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    default-jre \
    && rm -rf /var/lib/apt/lists/*


# Following libraries are to supplement some of the requirements by more recent R graphics libs
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Make sure RStudio uses the correct R version
RUN echo "rsession-which-r=/usr/local/bin/R" > /etc/rstudio/rserver.conf

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

# Install Python-based tools in a virtual environment
RUN python3 -m venv /opt/venv && \
    /opt/venv/bin/pip install --upgrade pip && \
    /opt/venv/bin/pip install \
    multiqc \
    cutadapt && \
    # Create symlinks to the venv binaries
    ln -s /opt/venv/bin/multiqc /usr/local/bin/ && \
    ln -s /opt/venv/bin/cutadapt /usr/local/bin/

# Install FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc_v0.12.1.zip && \
    mv FastQC /opt/ && \
    chmod +x /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc && \
    rm fastqc_v0.12.1.zip

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
    && Rscript -e 'BiocManager::install("org.Sc.sgd.db")' \
    && Rscript -e 'install.packages("ggiraph")' \
    && Rscript -e 'install.packages("plotly")'

# confirm all R packages load
RUN Rscript -e 'library(BiocManager)' \
    && Rscript -e 'library(limma)' \
    && Rscript -e 'library(edgeR)' \
    && Rscript -e 'library(DESeq2)' \
    && Rscript -e 'library(beeswarm)' \
    && Rscript -e 'library(knitr)' \
    && Rscript -e 'library(gplots)' \
    && Rscript -e 'library(goseq)' \
    && Rscript -e 'library(GO.db)' \
    && Rscript -e 'library(org.Sc.sgd.db)' \
    && Rscript -e 'library(ggiraph)' \
    && Rscript -e 'library(plotly)'

# Set up environment variables for the virtual environment
ENV PATH="/opt/venv/bin:$PATH"

# copy data
COPY docker/RNA_seq/ /var/lib/RNA_seq/
