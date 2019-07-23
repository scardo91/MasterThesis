VOXEL COLOR CODER
Author: Enrico Polo
Date: Dec. 2017
Oct-tree routine for PCL coding (color)

This repository contains the sources (C and Matlab) to implement a simple lossy color coder for dyniamic pointclouds
based on Octree decomposition and RAHT (Microsoft) and to evaluate Rate Distortion and Computational performances.

Thesis.pdf file contains some theoretical background, the state-of-art in 3D model compression, our implementation in detail and the results obtained.

slides.pdf its a short presentation of the work we have done.

The core of the coder is implemented in C while the preprocessing procedure is performed using Matlab.

PROCEDURE:
1) convert .ply 3D files using the proper Matlab function (see README.txt inside "C code" folder)
2) put the converted files (.bin and .col) inside "C code folder"
3) execute the routine
