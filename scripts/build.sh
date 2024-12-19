#!/bin/bash

mkdir -p ../build

cd ../src || exit
make clean
make -j build CFLAGS="-DKAGGLE -fcommon -Os -s -ffunction-sections -fdata-sections -flto -mavx2" LDFLAGS="-Wl,--gc-sections" ARCH=x86-64-avx2
strip cfish

cd ../build || exit

tar --transform 's/.*\///' -cf submission.tar ../src/cfish ../scripts/main.py
zopfli submission.tar
ls -l
