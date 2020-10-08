#!/bin/bash

# Add DDF to this dir

docker build . -t lofar_lba
singularity build /home/lofar_lba.sif docker-daemon://lofar_lba:latest
