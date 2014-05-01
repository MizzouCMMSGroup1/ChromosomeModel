'''
Test file
'''

#Numpy/Scipy
import numpy
from scipy import optimize

# matplot lib
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

# Argparser
import argparse

# ChromoUtil
import ChromoUtils as CU

def run(chromo,runtype):
	chromo.init_model()

	args = ()

	opts = {'maxiter' : 100, 'disp' : True }

	results = chromo.optimize(args,runtype,opts)

	print(results)
	print("saving final contact xyz coordinates")

	x = numpy.zeros(chromo.C.NUMBER_CONTACTS)
	y = numpy.zeros(chromo.C.NUMBER_CONTACTS)
	z = numpy.zeros(chromo.C.NUMBER_CONTACTS)

	for i in range(0,chromo.C.NUMBER_CONTACTS):
		x[i] = results.x[i*3]
		y[i] = results.x[i*3+1]
		z[i] = results.x[i*3+2]
		print(results.x[i*3], ',', results.x[i*3+1], ',', results.x[i*3+2])

	mpl.rcParams['legend.fontsize'] = 10

	fig = plt.figure()
	ax = fig.gca(projection='3d')

	ax.plot(x, y, z, label='3d plot of generated contacts')
	ax.legend()

	plt.show()

def main():
	parser = argparse.ArgumentParser(description="Runner for PyChromosomeModeler")
	
	parser.add_argument('-n','--number', help='number of residues to model',type=int,default=5)
	parser.add_argument('-f','--if_file',help='the Interaction Frequency file',type=argparse.FileType('r'))
	
	t_group = parser.add_mutually_exclusive_group()
	t_group.add_argument('-c','--conjugate',action='store_true')
	t_group.add_argument('-a','--anneal',action='store_true')
	
	args = parser.parse_args()

	runtype = 'CG'
	if args.anneal:
		runtype = 'Anneal'

	c = CU.C(args.number,if_file=args.if_file)

	chromo = CU.Chromo(c)

	run(chromo,runtype)

if __name__ == "__main__":
	main()