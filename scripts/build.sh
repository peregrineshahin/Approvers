#!/bin/bash

mkdir -p ../build

cd ../src || exit
make clean
make kaggle
strip cfish

cd ../build || exit

tar --transform 's/.*\///' -cf submission.tar ../src/cfish ../scripts/main.py
zopfli submission.tar
ls -l
