#!/usr/bin/env bash

if [[ ! -e data.tar.gz ]]; then
   wget -c https://www.huber.embl.de/msmb/data.tar.gz
fi

URL=https://www.huber.embl.de/msmb/code
for N in {01..13}; do
   FILE=${N}-chap.R
   if [[ ! -e ${FILE} ]]; then
      wget ${URL}/${FILE}
   fi
done

md5sum -c files.md5
