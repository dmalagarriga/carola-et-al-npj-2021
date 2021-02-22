# Get the base Ubuntu image from Docker Hub
FROM ubuntu:latest

LABEL Name=carola-et-al-njp-2021 Version=0.0.1
# Update apps on the base image
RUN apt-get -y update && apt-get install -y
 
# Install the Clang compiler
RUN apt-get -y install clang

RUN apt-get install libgsl-dev

RUN apt-get install --reinstall make
 
# Copy the current folder which contains C++ source code to the Docker image under /usr/src
COPY . /usr/src/carola-et-al-njp-2021
 
# Specify the working directory
WORKDIR /usr/src/carola-et-al-njp-2021
 