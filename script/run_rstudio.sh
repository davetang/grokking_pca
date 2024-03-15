#!/usr/bin/env bash

set -euo pipefail

RVER=4.3.3
IMAGE=davetang/rstudio:${RVER}
CONTAINER=rstudio_server_${RVER}
PACKAGEDIR=${HOME}/r_packages_${RVER}
PORT=8889

if [[ ! -d ${PACKAGEDIR} ]]; then
   mkdir ${PACKAGEDIR}
fi

docker run \
   --name ${CONTAINER} \
   -d \
   -p ${PORT}:8787 \
   -v ${PACKAGEDIR}:/packages \
   -v ${HOME}/github/:/home/rstudio/work \
   -e PASSWORD=password \
   -e USERID=$(id -u) \
   -e GROUPID=$(id -g) \
   ${IMAGE}

>&2 echo ${CONTAINER} listening on port ${PORT}

exit 0
