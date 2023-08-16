FROM nfcore/base
LABEL authors="urmovosa@ut.ee" \
      description="Docker image containing requirements for calculating MDS and visualising population stratification"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/eQTLGenPopAssign/bin:$PATH
ENV TAR="/bin/tar"
RUN ln -s /bin/tar /bin/gtar
RUN apt-get update && apt-get install -y gcc libcurl4
COPY github_packages/lme4qtl-master.zip /
RUN R -e "remotes::install_local('lme4qtl-master.zip')"
RUN R -e "remotes::install_version('bigsnpr', version = '1.10.8', dependencies = TRUE, repos = 'http://cran.rstudio.com/', upgrade = 'never')"
RUN wget http://bioconductor.org/packages/3.12/bioc/src/contrib/preprocessCore_1.52.1.tar.gz
RUN R CMD INSTALL --configure-args="--disable-threading" preprocessCore_1.52.1.tar.gz
