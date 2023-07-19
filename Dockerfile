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
matplotlib==3.5.1 \ 
numpy==1.22.2 \ 
scipy==1.7.3


COPY src* /opt/genepatt

## for testing purposes
RUN mkdir /opt/genepatt/gpunit


## Install Amplicon Classifier
RUN mkdir -p /home/programs
ADD https://github.com/AmpliconSuite/AmpliconClassifier/archive/refs/heads/main.zip /home/programs
RUN cd /home/programs && unzip main.zip
#RUN echo export AC_SRC=/home/programs/AmpliconClassifier-main >> ~/.bashrc # handled by .env file
RUN mkdir -p /opt/genepatt/.AA_DATA_REPO
#RUN echo export AA_DATA_REPO=/opt/genepatt/.AA_DATA_REPO >> ~/.bashrc # handled by .env file

ENV AC_SRC=/home/programs/AmpliconClassifier-main
ENV AA_DATA_REPO=/opt/genepatt/.AA_DATA_REPO