FROM nfcore/base
LABEL authors="urmovosa@ut.ee" \
      description="Docker image containing requirements for calculating MDS and visualising population stratification"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/eQTLGenPopAssign/bin:$PATH
ENV TAR="/bin/tar"
RUN ln -s /bin/tar /bin/gtar
RUN apt-get update && apt-get install -y gcc
RUN R -e "install.packages('bigsnpr', dependencies=TRUE, repos = 'http://cran.rstudio.com/')"
RUN wget http://bioconductor.org/packages/3.12/bioc/src/contrib/preprocessCore_1.52.1.tar.gz
RUN R CMD INSTALL --configure-args="--disable-threading" preprocessCore_1.52.1.tar.gz
