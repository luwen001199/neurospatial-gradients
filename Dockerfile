FROM python:3.13-slim

WORKDIR /work

RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc g++ libhdf5-dev procps \
    && rm -rf /var/lib/apt/lists/*

COPY environment.yml /tmp/environment.yml
RUN pip install --no-cache-dir \
        scanpy==1.12.4 anndata==0.13.4 numpy==2.5.4 pandas==3.1.0 \
        scipy==1.18.1 scikit-learn==1.9.1 statsmodels==0.15.0 \
        matplotlib==3.11.2 h5py==3.16.0 pillow==13.0.0 \
        pyyaml tqdm snakemake

COPY 00_pipeline /work/00_pipeline
COPY 02_scripts  /work/02_scripts
COPY Snakefile   /work/Snakefile

ENV PYTHONPATH=/work/00_pipeline:/work/02_scripts
CMD ["snakemake", "--cores", "4"]
