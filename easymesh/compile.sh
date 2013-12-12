#!/bin/bash

gcc -o EasyMesh -O3 easymesh_1_4.c -lm
gcc -o ShowMesh -O3 showmesh_1_0.c -lm -lX11
