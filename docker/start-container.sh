#!/bin/sh

# Run existing docker container
# -i: Basically have a terminal that you can interact with.
# feddlib: The container name. See "docker ps -a" for a container list.
docker start -i feddlib
# Exit docker container with "exit".
