'''
ChromoUtils file
Call via:
	import ChromoUtils as CU
Utility files and global variables only
'''

import numpy
import random
import math
from scipy import optimize
import copy

# Generic Chromosome using some Chromosome 7 parameters
class C:
	W1 = 1.0
	W2 = 1.5
	W3 = 1.5
	W4 = 1.5

	d_sq_min = 0.2
	da_sq_max = 1.8
	d_max = 4.5
	d_sq_c = 7.0
	d_sq_max = d_max * d_max

	def __init__(self,num_contacts,if_filename=None,if_file=None):
		self.NUMBER_CONTACTS = num_contacts
		self.NUMBER_CONTACTS_POINTS = self.NUMBER_CONTACTS * 3
		if if_filename is not None:
			self.IF_FILENAME = if_filename
			self.IF_DATA = numpy.loadtxt(self.IF_FILENAME, delimiter=',')
		elif if_file is not None:
			self.IF_FILE = if_file
			self.IF_DATA = numpy.loadtxt(self.IF_FILE, delimiter=',')
		else:
			raise BaseException("Must provide a filename or file")
		self.IF_TOTAL = numpy.sum(self.IF_DATA)

# Chromosome 7 specifc parameters
class C7:
	W1 = 1.0
	W2 = 1.5
	W3 = 1.5
	W4 = 1.5

	d_sq_min = 0.2
	da_sq_max = 1.8
	d_max = 4.5
	d_sq_c = 7.0
	d_sq_max = d_max * d_max

	NUMBER_CONTACTS = 157
	NUMBER_CONTACTS_POINTS = NUMBER_CONTACTS * 3

	def __init__(self,if_filename=None,if_file=None):
		if if_filename is not None:
			self.IF_FILENAME = if_filename
			self.IF_DATA = numpy.loadtxt(self.IF_FILENAME, delimiter=',')
		elif if_file is not None:
			self.IF_FILE = if_file
			self.IF_DATA = numpy.loadtxt(self.IF_FILE, delimiter=',')
		else:
			raise BaseException("Must provide a filename or file")
		self.IF_TOTAL = numpy.sum(self.IF_DATA)

class Chromo:

	def __init__(self,C):
		self.C = C
		self.coordinate_data = numpy.zeros(self.C.NUMBER_CONTACTS_POINTS)

	def init_model(self,bounding_box=0.5):
	    for i in range(0,self.C.NUMBER_CONTACTS_POINTS):
	        self.coordinate_data[i] = bounding_box * (0.5 - random.random())

	def generate_neighborhood(self,bounding_box=0.5):
		random_index = random.randrange(0,self.C.NUMBER_CONTACTS)
		neighborhood = {}
		for dim in ['x','y','z']:
			neighborhood[dim] = []
			plus = copy.copy(self)
			plus.coordinate_data[random_index] += 0.5
			minus = copy.copy(self)
			minus.coordinate_dat[random_index] -= 0.5
		return neighborhood


	def print_model(self):
	    print(self.coordinate_data)


	def max_if(self,i,j):
	    return max(self.C.IF_DATA[i,j], self.C.IF_DATA[j,i])


	def distance_sq(self,i,j):
	    a = [self.coordinate_data[i], self.coordinate_data[i+1], self.coordinate_data[i+2]]
	    b = [self.coordinate_data[j], self.coordinate_data[j+1], self.coordinate_data[j+2]]
	    return (a[0] - b[0])**2 + (a[1] - b[1])**2 + (a[2] - b[2])**2


	def contact_score(self):
	    '''
	    minimize the distance (but keep above min_threshold) between non-sequential pairs that have affinity
	    '''
	    score = 0
	    for i in range(0,self.C.NUMBER_CONTACTS):
	        for j in range(0,self.C.NUMBER_CONTACTS):
	            if i != j or abs(i-j) != 1:
	                d_sq_ij = self.distance_sq(i,j)
	                score += self.C.W1 * math.tanh(self.C.d_sq_c - d_sq_ij) * self.max_if(i,j) + self.C.W2 * math.tanh(d_sq_ij - self.C.d_sq_min) / self.C.IF_TOTAL
	    return score


	def noncontact_score(self):
	    '''
	    maximize the distance (but keep below max_threshold) between non-sequential pairs that don't have affinity
	    '''
	    score = 0
	    for i in range(0,self.C.NUMBER_CONTACTS):
	        for j in range(0,self.C.NUMBER_CONTACTS):
	            if i != j or abs(i-j) != 1:
	                d_sq_ij = self.distance_sq(i,j)
	                score += self.C.W3 * math.tanh(self.C.d_sq_max - d_sq_ij) / self.C.IF_TOTAL + self.C.W4 * math.tanh(d_sq_ij - self.C.d_sq_c) / self.C.IF_TOTAL
	    return score


	def pair_smoothing(self):
	    '''
	    keep adjacent contacts (eg |i-j|==1) with slightly lower score than above so they are prioritized for optimization
	    '''
	    score = 0
	    for i in range(0,self.C.NUMBER_CONTACTS):
	        for j in range(0,self.C.NUMBER_CONTACTS):
	            if abs(i-j) == 1:
	                d_sq_ij = self.distance_sq(i,j)
	                score += self.C.W1 * self.max_if(i,j) * math.tanh(self.C.da_sq_max - d_sq_ij) + self.C.W2 * math.tanh(d_sq_ij - self.C.d_sq_min) / self.C.IF_TOTAL
	    return score


	def model_score(self):
	    return self.contact_score() + self.noncontact_score() + self.pair_smoothing()


	# shim between skeleton and cg code
	iter_tracker = 0
	old_score = 0

	def f(self,x, *args):
	    #print x
	    self.iter_tracker += 1
	    
	    for i in range(0,self.C.NUMBER_CONTACTS_POINTS):
	        self.coordinate_data[i] = x[i]
	    
	    current_score = self.model_score()

	    print("iter:", self.iter_tracker, "score:", current_score, "change:", current_score - self.old_score)

	    self.old_score = current_score
	    
	    return current_score

	def optimize(self,args,runtype,opts):
		random_start = copy.copy(self.coordinate_data)
		return optimize.minimize(self.f, random_start, args=args, method=runtype, options=opts)
