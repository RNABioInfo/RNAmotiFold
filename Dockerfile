# syntax=docker/dockerfile:1
FROM ubuntu:24.04

RUN apt-get update && apt-get install -y python3 python3-pip python3.12-venv flex bison make libboost-all-dev libgsl-dev build-essential git m4
RUN git clone --branch master https://github.com/RNABioInfo/gapcM.git
RUN git clone --branch main --recurse-submodules https://github.com/RNABioInfo/RNAmotiFold

WORKDIR gapcM
RUN ./configure
RUN make -j 4
RUN make install

WORKDIR ../RNAmotiFold
RUN mkdir Build
RUN python3 -m venv ./Build/ 
ENV PATH="/Build/bin/:$PATH"


RUN . ./Build/bin/activate && pip install .

RUN . ./Build/bin/activate && RNAmotiFold_installer
RUN mv ./Build/bin/RNAmotiFold ../bin/
RUN mv ./Build/bin/RNAmotiCes ../bin/
RUN mv ./Build/bin/RNAmotiAlign ../bin/
RUN mv ./Build/bin/RNAmoSh ../bin/

WORKDIR ..