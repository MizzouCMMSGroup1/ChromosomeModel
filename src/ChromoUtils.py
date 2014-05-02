'''
ChromoUtils file
Call via:
	import ChromoUtils as CU
Utility files and global variables only
'''

import numpy
from numpy import random
import random as rand
import math
from scipy import optimize
import copy
import math

# Generic Chromosome using some Chromosome 7 parameters
class C:
	W1 = 1.0
	W2 = 1.5
	W3 = 1.5
	W4 = 1.5

	d_sq_min = 0.2
	d_min = math.sqrt(d_sq_min)
	da_sq_max = 1.8
	d_max = 4.5
	d_sq_c = 7.0
	d_c = math.sqrt(d_sq_c)
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
		self.IF_MAX = self.IF_DATA.max()

# Chromosome 7 specifc parameters
class C7:
	W1 = 1.0
	W2 = 1.5
	W3 = 1.5
	W4 = 1.5

	d_sq_min = 0.2
	d_min = math.sqrt(d_sq_min)
	da_sq_max = 1.8
	d_max = 4.5
	d_sq_c = 7.0
	d_c = math.sqrt(d_sq_c)
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
		self.IF_MAX = self.IF_DATA.max()

class Chromo:

	def __init__(self,C):
		self.C = C
		self.coordinate_data = numpy.zeros(self.C.NUMBER_CONTACTS_POINTS)

	def init_model(self,bounding_box=0.5):
	    for i in range(0,self.C.NUMBER_CONTACTS_POINTS):
	        self.coordinate_data[i] = bounding_box * (0.5 - random.random())

	def generate_neighborhood(self,scale=0.5):
		index = random.randint(0,self.C.NUMBER_CONTACTS)
		neighborhood = {}
		movement = scale
		for i,dim in enumerate(['x','y','z']):
			neighborhood[dim] = []
			plus = copy.copy(self)
			plus.coordinate_data = copy.copy(self.coordinate_data)
			plus.coordinate_data[index*3+i] += movement
			neighborhood[dim].append(plus)
			minus = copy.copy(self)
			minus.coordinate_data = copy.copy(self.coordinate_data)
			minus.coordinate_data[index*3+i] -= movement
			neighborhood[dim].append(minus)
		return (neighborhood,index)

	def random_neighbor(self,scale=0.5):
		(neighborhood,index) = self.generate_neighborhood(scale)
		dim = random.choice(list(neighborhood.keys()))
		coin = random.randint(2)
		return (neighborhood[dim][coin],index)


	def print_model(self):
	    print(self.coordinate_data)


	def max_if(self,i,j):
	    return max(self.C.IF_DATA[i,j], self.C.IF_DATA[j,i])/self.C.IF_TOTAL


	def distance_sq(self,i,j):
	    a = [self.coordinate_data[i], self.coordinate_data[i+1], self.coordinate_data[i+2]]
	    b = [self.coordinate_data[j], self.coordinate_data[j+1], self.coordinate_data[j+2]]
	    return (a[0] - b[0])**2 + (a[1] - b[1])**2 + (a[2] - b[2])**2


	def single_contact_score(self,i):
		score = 0
		for j in range(0,self.C.NUMBER_CONTACTS):
			if i != j and abs(i-j) != 1:
				d_sq_ij = self.distance_sq(i,j)
				if d_sq_ij < self.C.d_sq_c:
					score += self.C.W1 * math.tanh(self.C.d_sq_c - d_sq_ij) * self.max_if(i,j) + self.C.W2 * math.tanh(d_sq_ij - self.C.d_sq_min) / self.C.IF_TOTAL
		return score

	def contact_score(self):
	    '''
	    minimize the distance (but keep above min_threshold) between non-sequential pairs that have affinity
	    '''
	    score = 0
	    for i in range(0,self.C.NUMBER_CONTACTS):
	    	'''
	        for j in range(0,self.C.NUMBER_CONTACTS):
	            if i != j and abs(i-j) != 1:
	                d_sq_ij = self.distance_sq(i,j)
	                if d_sq_ij < self.C.d_sq_c:
	                	score += self.C.W1 * math.tanh(self.C.d_sq_c - d_sq_ij) * self.max_if(i,j) + self.C.W2 * math.tanh(d_sq_ij - self.C.d_sq_min) / self.C.IF_TOTAL
	    	'''
	    	score += self.single_contact_score(i)
	    return score

	def single_noncontact_score(self,i):
		score = 0
		for j in range(0,self.C.NUMBER_CONTACTS):
			if i != j and abs(i-j) != 1:
				d_sq_ij = self.distance_sq(i,j)
				if d_sq_ij < self.C.d_sq_c:
					score += self.C.W3 * math.tanh(self.C.d_sq_max - d_sq_ij) / self.C.IF_TOTAL + self.C.W4 * math.tanh(d_sq_ij - self.C.d_sq_c) / self.C.IF_TOTAL
		return score
	
	def noncontact_score(self):
	    '''
	    maximize the distance (but keep below max_threshold) between non-sequential pairs that don't have affinity
	    '''
	    score = 0
	    for i in range(0,self.C.NUMBER_CONTACTS):
	    	'''
	        for j in range(0,self.C.NUMBER_CONTACTS):
	            if i != j and abs(i-j) != 1:
	                d_sq_ij = self.distance_sq(i,j)
	                if d_sq_ij > self.C.d_sq_c:
	                	score += self.C.W3 * math.tanh(self.C.d_sq_max - d_sq_ij) / self.C.IF_TOTAL + self.C.W4 * math.tanh(d_sq_ij - self.C.d_sq_c) / self.C.IF_TOTAL
	    	'''
	    	score += self.single_noncontact_score(i)
	    return score

	def single_pair_smoothing(self,i):
		score = 0
		for j in range(0,self.C.NUMBER_CONTACTS):
			if abs(i-j) == 1:
				d_sq_ij = self.distance_sq(i,j)
				score += self.C.W1 * (self.C.IF_MAX/self.C.IF_TOTAL) * math.tanh(self.C.da_sq_max - d_sq_ij) + self.C.W2 * math.tanh(d_sq_ij - self.C.d_sq_min) / self.C.IF_TOTAL
		return score

	def pair_smoothing(self):
	    '''
	    keep adjacent contacts (eg |i-j|==1) with slightly lower score than above so they are prioritized for optimization
	    '''
	    score = 0
	    for i in range(0,self.C.NUMBER_CONTACTS):
	    	'''
	        for j in range(0,self.C.NUMBER_CONTACTS):
	            if abs(i-j) == 1:
	                d_sq_ij = self.distance_sq(i,j)
	                score += self.C.W1 * (self.C.IF_MAX/self.C.IF_TOTAL) * math.tanh(self.C.da_sq_max - d_sq_ij) + self.C.W2 * math.tanh(d_sq_ij - self.C.d_sq_min) / self.C.IF_TOTAL
	    	'''
	    	score += self.single_pair_smoothing(i)
	    return score


	def model_score(self):
	    return self.contact_score() + self.noncontact_score() + self.pair_smoothing()

	'''
	Print to PDB
	'''
	def printPDB(self,fp):
		lines = []
		for i in range(self.C.NUMBER_CONTACTS):
			x = self.coordinate_data[i*3]
			y = self.coordinate_data[i*3+1]
			z = self.coordinate_data[i*3+2]
			lines.append('ATOM  %5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n' % (
				i,"C",' ',"XXX",'A',i,' ',
				x,y,z,0.0,0.0,'C ','  '
			))
		fp.writelines(lines)

	'''
	Optimizations
	'''
	def update_score(old_conformation,old_score,new_conformation,index):
		minus_score = old_conformation.single_contact_score(index)+old_conformation.single_noncontact_score(index)+old_conformation.single_pair_smoothing(index)
		plus_score = new_conformation.single_contact_score(index)+new_conformation.single_noncontact_score(index)+new_conformation.single_pair_smoothing(index)
		return old_score - minus_score + plus_score

	#temperature calculator. non-linear decrease
	def sigmoid_temperature(e,t,k):
		return -t*2/(1 + math.exp(-5*k/t)) + t*2

	def linear_temperature(e,t,k):
		return (-t/e)*k + t

	def simulated_annealing(seed,epochs,temp,temp_func):
		current_conformation = seed
		current_score = current_conformation.model_score()

		for i in range(1,epochs+1):
			T = temp_func(epochs,temp,i)
			T = T if T > 1e-6 else 1e6 # Prevent divide-by-zero
			# Generate neighbor and score/diff
			(new_conformation,index) = current_conformation.random_neighbor(current_conformation.C.d_min/2)#*T/temp)
			new_score = Chromo.update_score(current_conformation,current_score,new_conformation,index)
			score_diff = (new_score - current_score) * 1e4
			# New conformation is better, accept
			if score_diff > 0:
				current_conformation = new_conformation
				current_score = new_score
				print("Iter",i," - Accepting higher score")
			else:
				# Limit score diff. Small negative moves may mean
				# drastic degredation of quality due to hyper-tangent
				# score_diff = 0 if score_diff > -0.005 else score_diff
				score_diff = 0 if score_diff > -0.5 else score_diff
				# Scale diff, as changes are normally in the 1e-2 range or below
				prob_to_accept = math.exp(score_diff/T)
				# Don't accept 0 score changes due to hyper-tangent above
				if score_diff != 0 and random.random() < prob_to_accept:
					current_conformation = new_conformation
					current_score = new_score
					print("Iter",i," - Accepting lower score")
				else:
					print("Iter",i," - Ignoring lower score")
			print("Iter",i," - Current score:",current_score)
			#print(current_conformation.coordinate_data)
		return current_conformation

	# MCMC
	def MCMC(seed,epochs):
		current_conformation = seed
		current_score = current_conformation.model_score()

		for i in range(1,epochs+1):
			(neighborhood,index) = current_conformation.generate_neighborhood(current_conformation.C.d_min*i/epochs)
			candidates = [[current_conformation,current_score,0]]
			
			min_score = current_score
			max_score = current_score
			
			for k,neighbors in neighborhood.items():
				for neighbor in neighbors:
					neighbor_score = Chromo.update_score(current_conformation,current_score,neighbor,index)
					candidates.append([neighbor,neighbor_score,0])
					if neighbor_score < min_score:
						min_score = neighbor_score
					if neighbor_score > max_score:
						max_score = neighbor_score
			score_range = max_score - min_score
			
			if score_range == 0:
				for j in range(len(candidates)):
					candidates[j][2] = 1/len(candidates)
			else:
				for j in range(len(candidates)):
					candidates[j][2] = abs((candidates[j][1] - min_score)/score_range)
			
			new_conformation = None
			while new_conformation is None:
				index = random.randint(len(candidates))
				if random.random() < candidates[index][2]:
					new_conformation = candidates[index]
			
			current_conformation = new_conformation[0]
			current_score = new_conformation[1]
			
			print("Iter",i," - Current score:",current_score)
		return current_conformation
			

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
