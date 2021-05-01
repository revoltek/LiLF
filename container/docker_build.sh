#!/bin/bash
docker build . -t revoltek/pill:latest > docker_build.log
#docker run -it revoltek/pill /bin/bash
docker tag 3b2c9c4d8db5 revoltek/pill:latest # or PiLL-0.1
docker push revoltek/pill
singularity build pill-latest.simg docker://revoltek/pill:latest
singularity run --pid --writable-tmpfs --containall --cleanenv -B/home/fdg:/home/fdg,/export/scratch/AG_deGasperin/fdg/:/home/fdg/data pill-latest.simg
