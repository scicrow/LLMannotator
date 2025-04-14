FROM rocker/r-ver:4.4.3

# Install R packages
RUN /rocker_scripts/install_tidyverse.sh
RUN R -e "install.packages(c('openai', 'foreach', 'do', 'proxy', 'ontologyIndex', 'ontologySimilarity', 'data.table'), repos='https://cloud.r-project.org')"

# Create and set working directory
RUN mkdir -p /home/workdir
COPY . /home/workdir/
WORKDIR /home/workdir/


# Notes on dependencies:
# openai:0.41 # Requires Requires R ≥ 3.5
### https://cran.r-project.org/web/packages/openai/index.html
#
# ontologyIndex:2.12 # Requires R ≥ 3.5
### https://cran.r-project.org/web/packages/ontologyIndex/index.html
#
# ontologySimilarity:2.7 # Requires R ≥ 3.5
### https://cran.r-project.org/web/packages/ontologySimilarity/index.html
#
# data.table:1.17.0 # Requires R ≥ 3.3.0
### https://cran.r-project.org/web/packages/data.table/index.html