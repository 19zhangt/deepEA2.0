#!/bin/bash

set -e 

docker run -t registry-snd.docker.iu.edu/befulton/cafe_prebuilds:ubuntu ./cafe_test.sh
docker run -t registry-snd.docker.iu.edu/befulton/cafe_prebuilds:clang ./cafe_test.sh
docker run -t registry-snd.docker.iu.edu/befulton/cafe_prebuilds:gcc54 ./cafe_test.sh
docker run -t registry-snd.docker.iu.edu/befulton/cafe_prebuilds:pgi ./cafe_test.sh

