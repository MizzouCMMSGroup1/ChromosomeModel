#get our data in, produces a 157x157 matrix

import numpy

if_data = numpy.loadtxt('FragChr7_1MB_Matrix_Lie_nml.txt', delimiter=',')

for i in range(0,157):
    print if_data[i,i]