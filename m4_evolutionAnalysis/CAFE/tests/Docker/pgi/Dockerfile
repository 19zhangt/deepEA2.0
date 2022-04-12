# Set the base image
FROM hopobcn/pgi-ce

# Dockerfile author / maintainer 
MAINTAINER Name befulton@iu.edu

RUN apt-get update && apt-get -y upgrade && apt-get -y install git make automake cpputest

COPY cafe_test.sh .

