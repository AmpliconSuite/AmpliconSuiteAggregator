## Parent Image
FROM ubuntu:20.04


RUN mkdir /opt/genepatt && chmod a+rwx /opt/genepatt

## Install required packages
RUN apt-get update
RUN apt install -y --fix-missing python3.8
RUN apt-get install -y --fix-missing python3-pip
RUN apt-get install -y unzip wget



WORKDIR /opt/genepatt


RUN pip3 install pandas==1.5.3

COPY src* /opt/genepatt

## for testing purposes
RUN mkdir /opt/genepatt/gpunit

