FROM rocker/r-devel-san

MAINTAINER Thibaut Jombart <thibautjombart@gmail.com>

RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y libcurl4-openssl-dev libssl-dev libfontconfig1-dev libxml2-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev gdal-bin proj-bin libgdal-dev libproj-dev libgmp3-dev jags libfftw3-dev

## add guest user

RUN adduser --disabled-password --gecos "" guest
RUN usermod -a -G users guest && usermod -a -G staff guest
RUN chmod a+rw /usr/local/lib/R/site-library -R



## install CRAN packages

# RUN echo 'options(download.file.method = "libcurl", repos = c(CRAN = "https://cran.ma.imperial.ac.uk"))' > ~/.Rprofile

RUN r -e "install.packages('devtools')" \
 && r -e "install.packages('adegenet', dependencies = TRUE)"

RUN r -e 'install.packages("BiocManager")' \
 && r -e 'BiocManager::install("SNPRelate")' \
 && r -e 'install.packages("dartR.base")'


## clone repos to get sources
RUN apt-get install -y git

RUN su guest
RUN mkdir ~/dev
WORKDIR /home/guest/dev

COPY . .

WORKDIR /home/guest/
CMD Rscript -e 'library(dartR.base); gl.pcoa(testset.gl)'
