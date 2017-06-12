#################################################################
# Source Image
FROM ubuntu:16.04

# Set noninterative mode
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################
LABEL version="1"
LABEL software="resistType"
LABEL software.version="0.1"
LABEL description="Tools for bacterial sequence analysis, focusing on enterobacteriaceae"
LABEL website="https://github.com/hangphan/resistType_docker"
LABEL documentation="https://github.com/hangphan/resistType_docker"
LABEL license="https://github.com/hangphan/resistType_docker"
LABEL tags="Genomics"

# Maintainer
MAINTAINER Hang Phan <hangphan@gmail.com>


# add apt mirror
RUN mv /etc/apt/sources.list /etc/apt/sources.list.bkp && \
    bash -c 'echo -e "deb mirror://mirrors.ubuntu.com/mirrors.txt xenial main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-updates main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-backports main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-security main restricted universe multiverse\n\n" > /etc/apt/sources.list' && \
    cat /etc/apt/sources.list.bkp >> /etc/apt/sources.list && \
        cat /etc/apt/sources.list

# apt update and install global requirements
RUN apt-get clean all
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y  \
    	    git \
	    default-jre \
	    python \
	    python-setuptools \
	    python-dev \
	    build-essential \
	    wget \
	    zip && \
	    apt-get clean && \
	    apt-get purge && \
	    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
	    easy_install pip && \
	    pip install --upgrade pip && \
	    pip install pysam && \
	    pip install BioPython && \
	    pip install numpy


RUN mkdir /data /config

# Install bowtie2

ENV ZIP=bowtie2-2.2.9-linux-x86_64.zip
ENV URL=https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/
ENV FOLDER=bowtie2-2.2.9
ENV DST=/home/bin
RUN mkdir $DST
RUN wget $URL/$ZIP -O $DST/$ZIP && \
    unzip $DST/$ZIP -d $DST && \
    rm $DST/$ZIP && \
    mv $DST/$FOLDER/* $DST && \
    rmdir $DST/$FOLDER

ENV PATH $DST:$PATH

# get resistType scripts
ENV DST=/home/
ENV FOLDER=resistType_docker
ENV URL=https://github.com/hangphan/$FOLDER/
RUN echo "cache-bust" && git clone $URL $DST/$FOLDER
ENV PATH $DST/$FOLDER/bin:$DST/$FOLDER/src/:$PATH

ENV DST=/home/$FOLDER/resources
ENV ZIP=SPAdes-3.10.0-Linux.tar.gz
RUN mkdir -p $DST/ && \
    tar -zxvf $DST/$ZIP  -C /home/$FOLDER/resources/
    ENV PATH $DST/SPAdes-3.10.0-Linux/bin:$PATH

WORKDIR /data/

# CMD ["resistType_v0.1.py -h "]