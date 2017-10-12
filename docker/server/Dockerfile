FROM rocker/verse

MAINTAINER Thibaut Jombart <thibautjombart@gmail.com>

RUN apt-get update && apt-get upgrade -y
RUN apt-get install libssl-dev libxml2-dev pandoc pandoc-citeproc libblas-dev liblapack-dev git qpdf -y

## add guest user

RUN adduser --disabled-password --gecos "" guest
RUN usermod -a -G users guest && usermod -a -G staff guest
RUN chmod a+rw /usr/local/lib/R/site-library -R



## install CRAN packages

RUN echo 'options(download.file.method = "libcurl", repos = c(CRAN = "https://cran.ma.imperial.ac.uk"))' > ~/.Rprofile

RUN r -e "install.packages('devtools')" \
 && r -e "install.packages('roxygen2')" \
 && r -e "install.packages('testthat')" \
 && r -e "install.packages('rmarkdown')" \
 && r -e "install.packages('adegenet', dependencies = c('Depends', 'Imports'))" \
 && r -e "install.packages('pegas')" \
 && r -e "install.packages('hierfstat')" \
 && r -e "install.packages('poppr')" \
 && r -e "install.packages('akima')" \
 && r -e "install.packages('maps')" \
 && r -e "install.packages('splancs')" \
 && r -e "install.packages('tripack')"



## install devel packages (github)

RUN r -e "devtools::install_github('thibautjombart/adegenet')"


## clone repos to get sources

RUN su guest
RUN mkdir ~/dev
WORKDIR /home/guest/dev

RUN git clone https://github.com/thibautjombart/adegenet

WORKDIR /home/guest/
RUN ls='ls --color=auto'
