#!/bin/sh

../kanatani_codes/projective_reconstruction/projective_reconstruction 2 ../kanatani_codes/projective_reconstruction/truck/truck*.dat \
    --method 2 \
    --output-3D ../kanatani_codes/projective_reconstruction/truck/result/3D.dat \
    --output-proj ../kanatani_codes/projective_reconstruction/truck/result/proj-%03d.dat
