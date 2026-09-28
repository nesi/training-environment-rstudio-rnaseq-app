# Workshop datasets are prepared in a separate stage so the downloaded zip and
# STAR build tooling never end up in the final image layers.
FROM ubuntu:24.04 AS data

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# STAR is only used here to rebuild the NZ backup genome index
ARG STAR_VERSION=2.7.11b
RUN curl -fsSL -o /tmp/STAR.zip https://github.com/alexdobin/STAR/releases/download/${STAR_VERSION}/STAR_${STAR_VERSION}.zip \
    && unzip -q /tmp/STAR.zip -d /tmp \
    && cp /tmp/STAR_${STAR_VERSION}/Linux_x86_64_static/STAR /usr/local/bin/STAR \
    && rm -rf /tmp/STAR* \
    && STAR --version

# NZ spotty wrasse dataset (GenomicsAotearoa/RNA-seq-data-analysis-workflow-NZ release asset)
ARG NZ_DATA_RELEASE=2026-September
ARG NZ_DATA_SHA256=34a48d817ee23c685db54e3508ce17dfb13a5d36d10a3d4351341f37f0893869
ARG NZ_GENOME=GCF_009762535.1_fNotCel1.pri_genomic
RUN curl -fsSL -o /tmp/RNA_seq.zip https://github.com/GenomicsAotearoa/RNA-seq-data-analysis-workflow-NZ/releases/download/${NZ_DATA_RELEASE}/RNA_seq.zip \
    && echo "${NZ_DATA_SHA256}  /tmp/RNA_seq.zip" | sha256sum -c - \
    && unzip -q /tmp/RNA_seq.zip -x '__MACOSX/*' -d /tmp/nz \
    && rm /tmp/RNA_seq.zip \
    && mkdir -p /var/lib/rnaseq \
    && mv /tmp/nz/RNA_seq /var/lib/rnaseq/nz \
    && find /var/lib/rnaseq/nz -name .DS_Store -delete \
    # the shipped index is built for a full-size genome (SAindex 1.5GB) with sjdbOverhang 100;
    # rebuild it to match the lesson (125bp reads) and the single-chromosome subset genome
    && rm -rf /var/lib/rnaseq/nz/backup-files/genomeindex \
    && mkdir -p /var/lib/rnaseq/nz-genomeindex \
    && cd /var/lib/rnaseq/nz-genomeindex \
    && STAR --runMode genomeGenerate \
        --genomeDir . \
        --genomeFastaFiles /var/lib/rnaseq/nz/Genome/${NZ_GENOME}.fna \
        --sjdbGTFfile /var/lib/rnaseq/nz/Genome/${NZ_GENOME}.gtf \
        --sjdbOverhang 124 \
        --genomeSAindexNbases 11 \
        --runThreadN 4 \
    && rm -rf _STARtmp \
    # the index is large and read-only, so link it rather than copying it into every home dir
    && ln -s /var/lib/rnaseq/nz-genomeindex /var/lib/rnaseq/nz/backup-files/genomeindex \
    # the lesson code refers to both DataFiles/ and data-files/, and writes to results/ and plots/
    && ln -s DataFiles /var/lib/rnaseq/nz/data-files \
    && mkdir -p /var/lib/rnaseq/nz/results /var/lib/rnaseq/nz/plots \
    && echo nz > /var/lib/rnaseq/nz/.dataset

# Yeast dataset (GenomicsAotearoa/RNA-seq-workshop), bundled in this repo
COPY docker/RNA_seq/ /var/lib/rnaseq/yeast/
RUN echo yeast > /var/lib/rnaseq/yeast/.dataset \
    && chmod -R u=rwX,go=rX /var/lib/rnaseq


FROM rocker/geospatial:4.5.2

# install system dependencies
RUN apt update && DEBIAN_FRONTEND=noninteractive apt install -y \
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
RUN apt update && \
    apt install -y --no-install-recommends \
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
    ln -s /opt/venv/bin/cutadapt /usr/local/bin/ && \
    multiqc --version && \
    cutadapt --version

# Install FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc_v0.12.1.zip && \
    mv FastQC /opt/ && \
    chmod +x /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc && \
    rm fastqc_v0.12.1.zip && \
    fastqc --version

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 && \
    tar -xjf samtools-1.18.tar.bz2 && \
    cd samtools-1.18 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm -rf samtools-1.18* && \
    samtools --version

# Install hisat2
RUN wget https://cloud.biohpc.swmed.edu/index.php/s/fE9QCsX3NH4QwBi/download -O hisat2-2.2.1-Linux_x86_64.zip && \
    unzip hisat2-2.2.1-Linux_x86_64.zip && \
    cd hisat2-2.2.1/ && \
    make && \
    cd .. && \
    mv hisat2-2.2.1/hisat2* /usr/local/bin/ && \
    mv hisat2-2.2.1/extract_exons* /usr/local/bin/ && \
    rm -rf hisat2-2.2.1* && \
    hisat2 --version

# Install subread - for featureCounts
RUN wget https://sourceforge.net/projects/subread/files/subread-2.0.3/subread-2.0.3-Linux-x86_64.tar.gz/download -O subread-2.0.3-Linux-x86_64.tar.gz && \
    tar -xzf subread-2.0.3-Linux-x86_64.tar.gz && \
    cp -r subread-2.0.3-Linux-x86_64/bin/* /usr/local/bin/ && \
    rm -rf subread-2.0.3*

# Install bamtools
RUN apt-get update && apt-get install -y cmake && \
    git clone https://github.com/pezmaster31/bamtools.git && \
    cd bamtools && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=/usr/local .. && \
    make && \
    make install && \
    cd ../.. && \
    rm -rf bamtools && \
    apt-get clean && rm -rf /var/lib/apt/lists/* && \
    bamtools --version

# Install STAR with required dependencies
RUN wget https://github.com/alexdobin/STAR/releases/download/2.7.11b/STAR_2.7.11b.zip && \
    unzip STAR_2.7.11b.zip && \
    cp STAR_2.7.11b/Linux_x86_64/STAR* /usr/local/bin/ && \
    rm -rf STAR_2.7.11b* && \
    STAR --version

# Install fastp with required dependencies
RUN wget http://opengene.org/fastp/fastp.1.3.2 && \
    chmod a+x fastp.1.3.2 && \
    mv fastp.1.3.2 /usr/local/bin/fastp && \
    fastp --version

# Install Trimmomatic (NZ workshop), with the NeSI module adapter path the lesson uses
ARG TRIMMOMATIC_VERSION=0.39
RUN wget https://github.com/usadellab/Trimmomatic/releases/download/v${TRIMMOMATIC_VERSION}/Trimmomatic-${TRIMMOMATIC_VERSION}.zip && \
    unzip Trimmomatic-${TRIMMOMATIC_VERSION}.zip && \
    mv Trimmomatic-${TRIMMOMATIC_VERSION} /opt/ && \
    rm Trimmomatic-${TRIMMOMATIC_VERSION}.zip && \
    printf '#!/bin/bash\nexec java -jar /opt/Trimmomatic-%s/trimmomatic-%s.jar "$@"\n' ${TRIMMOMATIC_VERSION} ${TRIMMOMATIC_VERSION} > /usr/local/bin/trimmomatic && \
    chmod +x /usr/local/bin/trimmomatic && \
    mkdir -p /opt/nesi/CS400_centos7_bdw/Trimmomatic/0.39-Java-1.8.0_144 && \
    ln -s /opt/Trimmomatic-${TRIMMOMATIC_VERSION}/adapters /opt/nesi/CS400_centos7_bdw/Trimmomatic/0.39-Java-1.8.0_144/adapters && \
    trimmomatic -version

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

# install R packages for the NZ workshop
RUN Rscript -e 'BiocManager::install(c("clusterProfiler", "enrichplot"))' \
    && Rscript -e 'install.packages("pheatmap")'

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
    && Rscript -e 'library(plotly)' \
    && Rscript -e 'library(tidyverse)' \
    && Rscript -e 'library(clusterProfiler)' \
    && Rscript -e 'library(enrichplot)' \
    && Rscript -e 'library(pheatmap)'

# Set up environment variables for the virtual environment
ENV PATH="/opt/venv/bin:$PATH"

# copy data, one layer per dataset; template/script.sh.erb stages the selected one into ~/RNA_seq
COPY --from=data /var/lib/rnaseq/yeast/ /var/lib/rnaseq/yeast/
COPY --from=data /var/lib/rnaseq/nz-genomeindex/ /var/lib/rnaseq/nz-genomeindex/
COPY --from=data /var/lib/rnaseq/nz/ /var/lib/rnaseq/nz/

# stages the selected dataset into ~/RNA_seq at session start (called from template/script.sh.erb)
COPY --chmod=755 docker/stage-rnaseq-data /usr/local/bin/stage-rnaseq-data
