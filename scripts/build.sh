#!/bin/bash

mkdir -p ../build

cd ../src || exit
make clean
make kaggle
strip approvers
chmod +x approvers

cd ../build || exit

tar --transform 's/.*\///' -cf submission.tar ../src/approvers ../scripts/main.py
zopfli submission.tar
ls -l
