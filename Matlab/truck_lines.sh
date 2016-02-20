#!/bin/sh

NAME="$1"

convert -size 1920x1080 ../../data/trucks/$NAME.jpg ../../data/trucks/$NAME.pgm

../../src/lsd_1.6/lsd ../../data/trucks/$NAME.pgm $NAME.txt
