## Parent Image
FROM ubuntu:20.04


RUN mkdir /opt/genepatt && chmod a+rwx /opt/genepatt

## Install required packages
RUN apt-get update
RUN apt install -y --fix-missing python3.8
RUN apt-get install -y --fix-missing python3-pip
RUN apt-get install -y unzip wget




WORKDIR /opt/genepatt


RUN pip3 install pandas==1.5.3 requests==2.31.0 intervaltree==3.1.0 \
Cython==0.29.28 biopython==1.79 \
reportlab==3.6.8 \ 
pyfaidx==0.6.4 \ 
pysam==0.18.0 \ 
cnvkit==0.9.10 \ 
intervaltree==3.1.0 \ 
Flask==2.2.5 \ 
matplotlib==3.5.1 \ 
numpy==1.22.2 \ 
scipy==1.7.3 \ 
mosek==10.0.38 \ 
future==0.18.3

COPY src* /opt/genepatt

## for testing purposes
RUN mkdir /opt/genepatt/gpunit


## Install Amplicon Classifier
RUN mkdir -p /home/programs
ADD https://github.com/jluebeck/AmpliconClassifier/archive/main.zip /home/programs
RUN cd /home/programs && unzip main.zip
RUN echo export AC_SRC=/home/programs/AmpliconClassifier-main >> ~/.bashrc
RUN mkdir -p /opt/genepatt/.AA_DATA_REPO
RUN echo export AA_DATA_REPO=/opt/genepatt/.AA_DATA_REPO >> ~/.bashrc