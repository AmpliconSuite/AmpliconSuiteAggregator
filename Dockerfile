## Parent Image
FROM ubuntu:20.04


RUN mkdir /opt/genepatt && chmod a+rwx /opt/genepatt
RUN mkdir /opt/genepatt/extracted
RUN mkdir /opt/genepatt/aggregated
RUN mkdir /opt/genepatt/programs
RUN mkdir /opt/genepatt/files
#COPY libs/AmpliconClassifier-main.zip /opt/genepatt/programs




## Install required packages
RUN apt-get update
RUN apt install -y --fix-missing python3.8
RUN apt-get install -y --fix-missing python3-pip
RUN apt-get install -y unzip wget



WORKDIR /opt/genepatt


RUN pip3 install pandas

COPY src* /opt/genepatt

## for testing purposes
RUN mkdir /opt/genepatt/gpunit
COPY gpunit/* /opt/genepatt/gpunit
