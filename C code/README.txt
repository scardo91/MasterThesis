VOXEL COLOR CODER
Author: Enrico Polo
Date: Dec. 2017
Oct-tree routine for PCL coding (color)



1) Compiling

in this folder compile together ac.h, ac.c (the arithmetic coder) and one of the other C files:

- nosortFINAL.c: RAHTv1 coder and decoder for static frames

- sortFINAL.c: RAHTv2 coder and decoder for static frames

- predseq_sort_FINAL.c: dynamic sequence encoding/decoding routine that uses RAHTv2

the files with the suffix 'exec' are the same routines executed 30 times on the same frame to better estimate the execution times.



2) Static routines execution

assuming the sources are compiled in the executable file 'static': 

>>./static <input.bin> <voxel resolution (per axis)> <input.col> <quant. step> <YUV flag> <decomp. levels> 

EXAMPLE: ./static vase.bin 512 vase.col 50 1 9



3) Dynamic routine execution

assuming the sources are compiled in the executable file 'dynamic': 

>>./dynamic <input.bin> <voxel resolution (per axis)> <input.col> <quant. step> <YUV flag> <decomp. levels> <number of frames> <pred. radius> <GOP> 

EXAMPLE: ./static frame0000.bin 512 frame0000.col 50 1 9 150 2 15


All .bin and .col files must be copied in this directory and are obtained by converting .ply 3D models using the functions (convert_ply_bin2(video).m and convert_ply_bin(static).m) in 'Matlab code' folder.
