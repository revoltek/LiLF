#!/bin/bash
# docker image prune -a
# docker system prune -a
docker build . -t revoltek/pill:20250113 > docker_build.log 2>&1
docker push revoltek/pill:20250113
singularity build pill.simg docker://revoltek/pill:20250113
singularity run --pid --writable-tmpfs --containall --cleanenv -B/home/fdg:/home/fdg,/export/scratch/AG_deGasperin/fdg/:/home/fdg/data pill.simg
