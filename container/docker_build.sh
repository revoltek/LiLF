#!/bin/bash
docker build . -t revoltek/pill:20240905 > docker_build.log 2>&1
docker push revoltek/pill:20240905
singularity build pill.simg docker://revoltek/pill:20240905
singularity run --pid --writable-tmpfs --containall --cleanenv -B/home/fdg:/home/fdg,/export/scratch/AG_deGasperin/fdg/:/home/fdg/data pill.simg
