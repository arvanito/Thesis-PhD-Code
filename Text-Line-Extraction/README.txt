MATLAB code for reproducing the results of our paper

"Binarization-Free Text Line Extraction for Historical Manuscripts"

The folder contains the following files:
- overlay_seams.m
- overlay_text_lines.m
- compute_text_lines.m
- compute_seams.m
- extract_text_lines.m
- run_text_line_extraction.m

and the folder anigaussm/, which contains code for fast 
anisotropic gaussian filtering from the authors of 

"Fast Anisotropic Gauss Filtering", IEEE TIP, vol.12, No.8, 2003.

The code can be also downloaded from the author's website:
http://staff.science.uva.nl/~mark/downloads.html#anigauss


Details of each function are given in the respective files. 
The main files are 'extract_text_lines.m' and 'run_text_line_extraction.m'. 
The file 'extract_text_lines.m' runs the algorithm on one manuscript image. 
The file 'run_text_line_extraction.m' runs the algorithm on 
all manuscript pages of a directory with a specific file extension. 

Parameter selection:
For the Al-Majid and Wadod collections of Saabni et al. dataset we define the following subsets:

- Al-Majid subset 1: majid_(57-64).png
- Al-Majid subset 2: majid_(1-56,65-96).png
- Wadod Spanish: Wadod_(7,9-15,17-19,21-23,25-36,38-42).jpg
- Wadod Arabic subset 1: Wadod_(1-6,8,48,56,59,62,66,68,70).jpg
- Wadod Arabic subset 2: Wadod_(16,20,24,37,43,44,46,47,53,55,58,61,63-65,67,69).png
- Wadod Arabic subset 3: Wadod_(45,49-52,54,57,60).jpg

The selected parameters for the above subsets are shown in Table V of our paper. 

For any questions regarding the code, please contact Nikolaos Arvanitopoulos at 
nick.arvanitopoulos@epfl.ch

Copyright (c) 2014, Nikolaos Arvanitopoulos, EPFL.