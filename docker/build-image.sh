#!/bin/sh

# Build image
# --progress=plain: show shell output during build process
# -t feddlib: Give the image a tag by which it can later be referred to (e.g., in "docker images").
# -f dockerfile_feddlib: The configuration/instructions used to build the docker image.
# ./configure-files: Files in this folder are available during the build process of the docker image. Therein, we have configure scripts for Trilinos and FEDDLib.
docker build --progress=plain -t feddlib -f dockerfile_feddlib ./configure-files
