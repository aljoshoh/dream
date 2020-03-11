## Start from this Docker image
## FROM ubuntu
FROM rocker/verse:3.6.1

## Install R in Docker image
## RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y r-base

## Install dependency for RCurl
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
	libcurl4-openssl-dev \
	libxml2-dev \
	aptitude

## Install R packages in Docker image
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages(c('tidyr','dplyr','randomForest','glmnet','survival','Hmisc','readr','randomForestSRC','caret'))"
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/h2o/h2o_3.26.0.2.tar.gz', repos=NULL, type='source')"

## Copy run_model.R, models folder, input_data_functions.R & general.R into Docker image
COPY . /home/
RUN mkdir /output
RUN mkdir /input
RUN chmod a+x /home/run_model_sc1.R

## Make Docker container executable
ENTRYPOINT ["Rscript", "/home/run_model_sc1.R"]
