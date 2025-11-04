#!/bin/sh

# Delete container created by new-container.sh.
docker container rm feddliba_jk
# The command "docker container rm" or (in short) "docker rm" with the container ID or container name will delete the container (not the image). It will delete all changes that you as a user have made and data you created.
