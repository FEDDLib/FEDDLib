#!/bin/sh

# Create container for docker image
# An image is not changed in the Docker ecosystem. On top of images, containers are created, which use the image as a base. You can then work within the container. Any changes, data etc., are saved inside a container. You can have many containers for one image.
# -it: Basically have a terminal that you can interact with.
# -name: Give the container a name. This also prevents creating the same container twice (in case that is something you want).
# -v: Mount ./shared-files of the host in /ext_docker in the docker container.
#     This will make all files of this folder (where this script etc. reside in) available from within the docker container (read & write).
# feddlib: Docker image that is used to create the container.
# /bin/bash: shell to use for the terminal
docker container create -it --name feddlib -v ./shared-files:/ext_docker feddlib /bin/bash
# To create and run a container in one step, use 
#    docker run -it --name feddlib -v ./shared-files:/ext_docker feddlib /bin/bash
# A container will have been created. Use "docker container ls --all" or (in short) "docker ps -a" to see a list.
# Use the command "docker container rm" or (in short) "docker rm" with the container ID or container name to delete the container, e.g., "docker rm feddlib". Note that this will not delete the image, but it will delete all changes that you as a user have made and data you created.
