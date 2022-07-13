## Parent Image
FROM ubuntu:20.04


RUN mkdir /opt/genepatt && chmod a+rwx /opt/genepatt
RUN mkdir /opt/genepatt/extracted
RUN mkdir /opt/genepatt/aggregated
WORKDIR /opt/genepatt


## Install required packages
RUN apt-get update
RUN apt install -y --fix-missing python3.8
RUN apt-get install -y --fix-missing python3-pip

RUN pip3 install pandas

COPY src* /opt/genepatt
