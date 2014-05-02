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

def run(chromo,runtype,epochs=1000,temp=2500,outfilename=None):

	chromo.init_model(bounding_box=chromo.C.d_max)

	args = ()

	opts = {'maxiter' : 100, 'disp' : True }

	x = numpy.zeros(chromo.C.NUMBER_CONTACTS)
	y = numpy.zeros(chromo.C.NUMBER_CONTACTS)
	z = numpy.zeros(chromo.C.NUMBER_CONTACTS)

	if runtype == 'SA':
		best_model = CU.Chromo.simulated_annealing(chromo,epochs,temp,CU.Chromo.linear_temperature,outfilename=outfilename)
		
		print(best_model)
		print("saving final contact xyz coordinates")
		
		for i in range(0,chromo.C.NUMBER_CONTACTS):
			x[i] = best_model.coordinate_data[i*3]
			y[i] = best_model.coordinate_data[i*3+1]
			z[i] = best_model.coordinate_data[i*3+2]
			print(x[i], ',', y[i], ',', z[i])
	elif runtype == 'MCMC':
		best_model = CU.Chromo.MCMC(chromo,epochs,outfilename=outfilename)
		
		print(best_model)
		print("saving final contact xyz coordinates")
		
		for i in range(0,chromo.C.NUMBER_CONTACTS):
			x[i] = best_model.coordinate_data[i*3]
			y[i] = best_model.coordinate_data[i*3+1]
			z[i] = best_model.coordinate_data[i*3+2]
			print(x[i], ',', y[i], ',', z[i])
	elif runtype == 'CG':
		results = chromo.optimize(args,runtype,opts)

		print(results)
		print("saving final contact xyz coordinates")

		for i in range(0,chromo.C.NUMBER_CONTACTS):
			x[i] = results.x[i*3]
			y[i] = results.x[i*3+1]
			z[i] = results.x[i*3+2]
			print(x[i], ',', y[i], ',', z[i])

	mpl.rcParams['legend.fontsize'] = 10

	fig = plt.figure('Chromo',figsize=(6,6))
	ax = fig.gca(projection='3d')

	ax.plot(x, y, z, label='3d plot of generated contacts')
	ax.legend()

	plt.show()
	plt.savefig("%s_%s.png" % (runtype,outfilename), dpi=96, format='png')

def main():
	parser = argparse.ArgumentParser(description="Runner for PyChromosomeModeler")
	
	parser.add_argument('-n','--number', help='number of residues to model',type=int,default=5)
	parser.add_argument('-f','--if_file',help='the Interaction Frequency file',type=argparse.FileType('r'))
	parser.add_argument('-o','--outputname',help='output pdb filename',type=str)

	t_group = parser.add_mutually_exclusive_group()
	t_group.add_argument('-c','--conjugate',action='store_true')
	t_group.add_argument('-a','--anneal',action='store_true')
	t_group.add_argument('-m','--mcmc',action='store_true')

	parser.add_argument('-e','--epochs',help='number of epochs to run for simulated annealing',type=int,default=1000)
	parser.add_argument('-t','--temp',help='maximum temperature for simulated annealing',type=int,default=2500)
	
	args = parser.parse_args()

	runtype = 'SA'
	if args.conjugate:
		runtype = 'CG'
	elif args.mcmc:
		runtype = 'MCMC'

	c = CU.C(args.number,if_file=args.if_file)

	epochs = args.epochs
	temp = args.temp

	chromo = CU.Chromo(c)

	run(chromo,runtype,epochs,temp,args.outputname)

	fp = open('%s_%s.pdb' % (runtype,args.outputname),'w')

	chromo.printPDB(fp)

	fp.close()

if __name__ == "__main__":
	main()