# prep our data for usage later

import numpy

# get our data in, produces a 157x157 matrix
if_data = numpy.loadtxt('FragChr7_1MB_Matrix_Lie_nml.txt', delimiter=',')

# per paper, drop IF contact values below threshold
IF_CUTOFF = .66
IF_FILENAME = "if_data_stripped.csv"

for i in range(0,157):
    for j in range(0,157):
        if if_data[i,j] < IF_CUTOFF:
            if_data[i,j] = 0

numpy.savetxt(IF_FILENAME, if_data, delimiter=",")

# for the next step
#stripped_if_data = numpy.loadtxt(IF_FILENAME, delimiter=',')
