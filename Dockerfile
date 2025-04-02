FROM ubuntu:20.04

RUN apt-get update && \
    apt-get install -y wget bzip2 ca-certificates git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

ENV PATH="/opt/conda/bin:$PATH"

WORKDIR /app

COPY . /app

RUN conda env create -f GCNPath.yaml && \
    echo "source activate GCNPath" >> ~/.bashrc

CMD ["/bin/bash"]