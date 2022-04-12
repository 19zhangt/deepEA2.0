#!/bin/bash

set -e 

docker run -t registry-snd.docker.iu.edu/befulton/cafe_builds:ubuntu ./cafe_test.sh
docker run -t registry-snd.docker.iu.edu/befulton/cafe_builds:clang ./cafe_test.sh
docker run -t registry-snd.docker.iu.edu/befulton/cafe_builds:gcc54 ./cafe_test.sh
docker run -t registry-snd.docker.iu.edu/befulton/cafe_builds:pgi ./cafe_test.sh

