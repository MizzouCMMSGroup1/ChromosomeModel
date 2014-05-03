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

# Other librs
import copy

# ChromoUtil
import ChromoUtils as CU

def run(regions,runs,runtype,epochs,temp,if_file,outfilename=None):
	args = ()

	opts = {'maxiter' : 100, 'disp' : True }

	x = numpy.zeros(regions)
	y = numpy.zeros(regions)
	z = numpy.zeros(regions)

	chromosomes = []
	for i in range(runs):
		c = CU.C(regions,if_file=open(if_file.name,'r'))
		chromo = CU.Chromo(c)
		chromo.init_model(bounding_box=chromo.C.d_max)
		chromosomes.append(chromo)

	conformations = []
	scores = []

	if runtype == 'SA':
		for i in range(runs):
			temp_outname = "%s-%d" % (outfilename,i)
			linear = CU.Chromo.linear_temperature
			conformation,score = CU.Chromo.simulated_annealing(chromo,epochs,temp,linear,outfilename=temp_outname)
			conformations.append(conformation)
			scores.append(score)
	elif runtype == 'MCMC':
		for i in range(runs):
			temp_outname = "%s-%d" % (outfilename,i)
			conformation,score = CU.Chromo.MCMC(chromo,epochs,outfilename=temp_outname)
			conformations.append(conformation)
			scores.append(score)
	'''
	elif runtype == 'CG':
		results = chromo.optimize(args,runtype,opts)

		print(results)
		print("saving final contact xyz coordinates")

		for i in range(0,chromo.C.NUMBER_CONTACTS):
			x[i] = results.x[i*3]
			y[i] = results.x[i*3+1]
			z[i] = results.x[i*3+2]
			print(x[i], ',', y[i], ',', z[i])
	'''

	'''
	mpl.rcParams['legend.fontsize'] = 10

	fig = plt.figure('Chromo',figsize=(6,6))
	ax = fig.gca(projection='3d')

	ax.plot(x, y, z, label='3d plot of generated contacts')
	ax.legend()

	plt.savefig("%s_%s.png" % (runtype,outfilename), dpi=96, format='png')
	'''

	mpl.rcParams['legend.fontsize'] = 10
	fig = plt.figure(outfilename,figsize=(9,9))
	subplt = fig.add_subplot(111)
	for i in range(runs):
		subplt.plot(range(epochs),scores[i])
	#plt.show()
	plt.savefig("%s_%s-Scores.png" % (runtype,outfilename), dpi=96, format='png')

	for i in range(runs):
		fp = open('%s_%s-%d.pdb' % (runtype,outfilename,i),'w')
		chromo.printPDB(fp)
		fp.close()

def main():
	parser = argparse.ArgumentParser(description="Runner for PyChromosomeModeler")
	
	parser.add_argument('-r','--regions', help='number of regions in the model',type=int,default=20)
	parser.add_argument('-n','--number', help='number of simulations to run',type=int,default=5)
	parser.add_argument('-f','--if_file',help='the Interaction Frequency file',type=argparse.FileType('r'))
	parser.add_argument('-o','--outputname',help='output pdb filename',type=str)

	t_group = parser.add_mutually_exclusive_group()
	t_group.add_argument('-c','--conjugate',action='store_true')
	t_group.add_argument('-a','--anneal',action='store_true')
	t_group.add_argument('-m','--mcmc',action='store_true')

	parser.add_argument('-e','--epochs',help='number of epochs to run for simulated annealing',type=int,default=2000)
	parser.add_argument('-t','--temp',help='maximum temperature for simulated annealing',type=int,default=250)
	
	args = parser.parse_args()

	runtype = 'SA'
	if args.conjugate:
		runtype = 'CG'
	elif args.mcmc:
		runtype = 'MCMC'

	epochs = args.epochs
	temp = args.temp

	run(args.regions,args.number,runtype,epochs,temp,args.if_file,args.outputname)

if __name__ == "__main__":
	main()